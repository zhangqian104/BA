#include <fstream>
#include <sys/stat.h>
#include <map>
#include<set>
#include<random>
#include "imu.h"
#include "utilities.h"
#include "triangulation.h"

using namespace std;

std::string houseModelPath = "C:\\Users\\Administrator\\Desktop\\simulation\\bin\\house_model\\house.txt";
std::string outputDataPath = "C:\\Users\\Administrator\\Desktop\\simulation\\simulation\\OutputData";
std::string camPoseAndObsPath = "CamPoseAndObs.txt";

extern IdType g_Id = 0;
extern IdType NewId()
{
	return g_Id++;
}
extern void ResetId()
{
	g_Id = 0;
}

struct MapPoint
{
	IdType id;

	double x;
	double y;
	double z;

	double u;
	double v;
};

struct CamPose
{
	IdType id;
	double t;

	double qw;
	double qx;
	double qy;
	double qz;
	double tx;
	double ty;
	double tz;

	Eigen::Quaterniond q;
	Eigen::Vector3d tran;
	Eigen::Matrix3d R;

	std::vector<MapPoint> mpInCam;
	std::vector<RowPose> rsRowPoses;
};

struct KeyFrame
{
	double timestamp;
	std::vector<Eigen::Vector3d> imuAccs;
	std::vector<Eigen::Vector3d> imuGyros;
	std::vector<RowPose> rsPoses;
	PointsType pointsCam;    // 3D points in current cam
	ObsersType featuresCam;  // corresponding 2D observations

	Eigen::Matrix3d Rwc;
	Eigen::Vector3d twc;
	Eigen::Vector3d twcNoise;
};

bool IsNewPoint(PointsType &points, Eigen::Vector4d &point)
{
	for (unsigned int idx = 0; idx < points.size();idx++) {
		if (points[idx].second == point) {
			return false;
		}
	}
	return true;
}

void InsertNewPoint(PointsType &points, Eigen::Vector4d &point)
{
	if (IsNewPoint(points, point)) {
		std::pair<IdType, Eigen::Vector4d> newPoint(NewId(), point);
		points.push_back(newPoint);
	}
}

void ReadCamPoseAndObs(std::vector<CamPose> &allCamPose)
{
	std::ifstream f;
	f.open(camPoseAndObsPath);

	ResetId();

	allCamPose.clear();
	std::vector<MapPoint> currObsInOneCam;
	std::vector<MapPoint> prevObsInOneCam;
	bool isNewPose = true;
	unsigned int numCam = 0;

	while (!f.eof()) {
		std::string s;
		std::getline(f, s);
		if (!s.empty()) {
			std::stringstream ss;
			ss << s;
			int firstNum;
			ss >> firstNum;
			if (firstNum < 0) {
				CamPose camPose;
				camPose.id = NewId();
				ss >> camPose.t;
				ss >> camPose.qw;
				ss >> camPose.qx;
				ss >> camPose.qy;
				ss >> camPose.qz;
				ss >> camPose.tx;
				ss >> camPose.ty;
				ss >> camPose.tz;

				camPose.q = Eigen::Quaterniond (camPose.qw, camPose.qx, camPose.qy, camPose.qz);
				camPose.R = camPose.q.toRotationMatrix();
				camPose.tran = Eigen::Vector3d (camPose.tx, camPose.ty, camPose.tz);

				prevObsInOneCam = currObsInOneCam;
				currObsInOneCam.clear();
				if (allCamPose.size() > 0 && numCam>=1) {
					allCamPose[numCam-1].mpInCam = prevObsInOneCam;
				}
				allCamPose.push_back(camPose);
				isNewPose = true;
				numCam++;
			} else {
				MapPoint mp;
				double tmp;
				mp.id = static_cast<IdType>(firstNum);
				ss >> mp.x;
				ss >> mp.y;
				ss >> mp.z;
				ss >> tmp;
				ss >> mp.u;
				ss >> mp.v;
				if (isNewPose) {
					currObsInOneCam.clear();
					currObsInOneCam.push_back(mp);
					isNewPose = false;
				} else {
					currObsInOneCam.push_back(mp);
				}
			}
		}
	}
}

void CreateOrReadPoints(PointsType& points)
{
    std::ifstream f;
    f.open(houseModelPath);

    while(!f.eof()){
        std::string s;
        std::getline(f,s);
        if(!s.empty()){
            std::stringstream ss;
            ss << s;
            double x,y,z;
            ss >> x;
            ss >> y;
            ss >> z;
            Eigen::Vector4d pt0( x, y, z, 1 );
            ss >> x;
            ss >> y;
            ss >> z;
            Eigen::Vector4d pt1( x, y, z, 1 );

			InsertNewPoint(points, pt0);
			InsertNewPoint(points, pt1);
        }
    }

    // save points
    SavePoints("all_points.txt", points);
}

void CreatePointObs(std::vector<MotionData> &camData, PointsType &points)
{
	Param param;
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> noise(0.0, param.pixelNoiseSigma);

	std::ofstream saveFile;
	saveFile.open("CamPoseAndObs.txt");

	for (int n = 0; n < camData.size(); ++n){
		MotionData data = camData[n];
		Eigen::Matrix4d Twc = Eigen::Matrix4d::Identity();
		Twc.block(0, 0, 3, 3) = data.Rwb;
		Twc.block(0, 3, 3, 1) = data.twb;

		// 遍历所有的特征点，看哪些特征点在视野里
		PointsType pointsCam;    // ３维点在当前cam视野里
		ObsersType featuresCam;  // 对应的２维图像坐标
		for (int i = 0; i < points.size(); ++i) {
			IdType pointId = points[i].first;
			Eigen::Vector4d pw = points[i].second;          // 最后一位存着feature id
			pw[3] = 1;                               //改成齐次坐标最后一位
			Eigen::Vector4d pc1 = Twc.inverse() * pw; // T_wc.inverse() * Pw  -- > point in cam frame

			if (pc1(2) < 0) continue; // z必须大于０,在摄像机坐标系前方

			Eigen::Vector2d obsNorm(pc1(0) / pc1(2), pc1(1) / pc1(2));
			double u = obsNorm[0] * param.fx + param.cx;
			double v = obsNorm[1] * param.fy + param.cy;

			//Create rolling shutter observation
			if (param.rsT > 0.0 && v > 0. && v < param.image_h) {
				unsigned int row = static_cast<unsigned int>(std::floor(v));
				Eigen::Matrix4d iTwc = Eigen::Matrix4d::Identity();
				iTwc.block(0, 0, 3, 3) = data.rsPoses[row].Rwb;
				iTwc.block(0, 3, 3, 1) = data.rsPoses[row].twb;
				Eigen::Vector4d ipc1 = iTwc.inverse()*pw;
				if (ipc1(2) < 0) { 
					continue; 
				}
				Eigen::Vector2d iObsNorm(ipc1(0) / ipc1(2), ipc1(1) / ipc1(2));
				u = iObsNorm[0] * param.fx + param.cx;
				v = iObsNorm[1] * param.fy + param.cy;
			}

			if (param.addPixedNoise) {
				u += noise(generator);
				v += noise(generator);
			}
			if (u>0. && u<param.image_w && v>0. && v<param.image_h ){
				Eigen::Vector2d pt(u, v);
				ObserType obs(pointId, pt);
				pointsCam.push_back(points[i]);
				featuresCam.push_back(obs);
			}
		}

		// save points
		camData[n].pointsCam = pointsCam;
		camData[n].featuresCam = featuresCam;

		std::stringstream filename1;
		filename1 << outputDataPath <<"\\"<< "AllObsInImg_time_"<< n << ".txt";
		SaveFeatures(filename1.str(), data, pointsCam, featuresCam);
		SavaPoseAndFeatures(saveFile, data, pointsCam, featuresCam);
	}
}

void CreateIMUPoseAndObs(std::vector< MotionData > &imuData, std::vector< MotionData > &imuDataNoise)
{
	// IMU model
	Param params; //Initial parameters in Param
	IMU imuGen(params); //Create imu pose and measures

	// create imu data
	// imu pose gyro acc
	/*std::vector< MotionData > imuData;
	std::vector< MotionData > imuDataNoise;*/

	for (float t = params.tStart; t < params.tEnd;) {
		MotionData imu = imuGen.MotionModel(t);
		imuData.push_back(imu);

		// add imu noise
		MotionData imuNoise = imu;
		imuGen.AddIMUnoise(imuNoise);
		imuDataNoise.push_back(imuNoise);

		t += 1.0 / params.imuFrequency;
	}
	imuGen.initVelocity = imuData[0].imuVelocity;
	imuGen.init_twb = imuData.at(0).twb;
	imuGen.initRwb = imuData.at(0).Rwb;
	savePose("imu_pose.txt", imuData);
	savePose("imu_pose_noise.txt", imuDataNoise);

	imuGen.testImu("imu_pose.txt", "imu_integ_pose.txt");     // test the imu data, integrate the imu data to generate the imu trajecotry
	imuGen.testImu("imu_pose_noise.txt", "imu_integ_pose_noise.txt");
}

void CreateCamPoseAndObs(PointsType &points, std::vector< MotionData > &camData)
{
	// IMU model
	Param params; //Initial parameters in Param
	IMU imuGen(params); //Create imu pose and measures

	// cam pose
	//std::vector< MotionData > camData;
	for (float t = params.tStart; t < params.tEnd;) {
		MotionData imu = imuGen.MotionModel(t);   // imu body frame to world frame motion
		MotionData cam;

		cam.timestamp = imu.timestamp;
		cam.Rwb = imu.Rwb * params.Rbc;    // cam frame in world frame: Rwb actually is Rwc
		cam.twb = imu.twb + imu.Rwb * params.tbc; //  Twc = Twb * Tbc ,  t = Rwb * tbc + twb

		//Create rolling shutter row poses
		if (params.rsT > 0.0 && params.image_h>0.0) {
			double readoutTime = params.rsT / params.image_h;
			for (unsigned int row = 0; row < static_cast<unsigned int>(params.image_h); row++) {
				double iTime = t + row * readoutTime;
				RowPose iPose;
				MotionData iImu = imuGen.MotionModel(iTime);
				iPose.row = row;
				iPose.timestamp = iTime;
				iPose.Rwb = iImu.Rwb*params.Rbc;
				iPose.twb = iImu.twb + iImu.Rwb*params.tbc;
				iPose.qwb = Eigen::Quaterniond (iPose.Rwb);
				cam.rsPoses.push_back(iPose);
			}
		}

		camData.push_back(cam);
		t += 1.0 / params.camFrequency;
	}
	savePose("cam_pose.txt", camData);

	// points obs in image
	CreatePointObs(camData, points);
}

int main(){
    // Read Model and Get 3d points and lines
	PointsType points;
	CreateOrReadPoints(points);
	std::vector< MotionData > imuData;
	std::vector< MotionData > imuDataNoise;
	CreateIMUPoseAndObs(imuData,imuDataNoise);
	std::vector< MotionData > camData;
	CreateCamPoseAndObs(points,camData);

	//Key frames for lba unit test
	Param params;
	std::random_device rd;
	std::default_random_engine generator(rd());
	std::normal_distribution<double> noise(0.0, params.camPoseNoise);
	std::vector<KeyFrame> keyFrames;
	for (unsigned int idx1 = 1; idx1 < camData.size(); idx1++) {
		KeyFrame keyFrame;
		MotionData prevCam = camData[idx1 - 1];
		MotionData currCam = camData[idx1];
		double prevCamTime = prevCam.timestamp;
		double currCamTime = currCam.timestamp;
		std::vector<Eigen::Vector3d> imuAccs;
		std::vector<Eigen::Vector3d> imuGyros;
		for (unsigned int idx2 = 0; idx2 < imuDataNoise.size(); idx2++) {
			MotionData imu = imuDataNoise[idx2];
			double imuTime = imu.timestamp;
			if (imuTime > prevCamTime && imuTime <= currCamTime) {
				imuAccs.push_back(imu.imuAcc);
				imuGyros.push_back(imu.imuGyro);
			}
		}

		keyFrame.timestamp = currCamTime;
		keyFrame.imuAccs = imuAccs;
		keyFrame.imuGyros = imuGyros;
		keyFrame.rsPoses = currCam.rsPoses;
		keyFrame.pointsCam = currCam.pointsCam;
		keyFrame.featuresCam = currCam.featuresCam;
		keyFrame.Rwc = currCam.Rwb;
		keyFrame.twc = currCam.twb;
		Eigen::Vector3d poseNoise(noise(generator), noise(generator), noise(generator));
		keyFrame.twcNoise = Eigen::Vector3d(keyFrame.twc+poseNoise);
		keyFrames.push_back(keyFrame);
		std::cout << noise(generator) << std::endl;
	}

	std::vector<CamPose> allCamPose;
	ReadCamPoseAndObs(allCamPose);

	Param param;
	Triangulation triang;
	CamPose cam1 = allCamPose[0];
	Eigen::Matrix<double, 3, 4> Twc1;
	Twc1 << cam1.R, cam1.tran;
	std::vector<MapPoint> mpInCam1 = cam1.mpInCam;
	for (unsigned int camIdx = 1; camIdx < allCamPose.size(); camIdx++) {
		CamPose cam2 = allCamPose[camIdx];
		Eigen::Matrix<double, 3, 4> Twc2;
		Twc2 << cam2.R, cam2.tran;
		std::vector<MapPoint> mpInCam2 = cam2.mpInCam;

		//Get matchers between cam1 and cam2
		std::vector<std::pair<IdType, IdType>> matchers;
		for (unsigned int idx1 = 0; idx1 < mpInCam1.size(); idx1++) {
			for (unsigned int idx2 = 0; idx2 < mpInCam2.size(); idx2++) {
				if (mpInCam1[idx1].id == mpInCam2[idx2].id) {
					std::pair<IdType, IdType> match(idx1, idx2);
					matchers.push_back(match);
				}
			}
		}

		//Calculate parallax angle
		double avgParaAngle = 0.;
		for (unsigned int mIdx = 0; mIdx < matchers.size(); mIdx++) {
			std::pair<IdType, IdType> match = matchers[mIdx];
			double x, y, z;
			x = cam1.tx;
			y = cam1.ty;
			z = cam1.tz;
			Eigen::Vector3d leftPt(x, y, z);
			x = cam2.tx;
			y = cam2.ty;
			z = cam2.tz;
			Eigen::Vector3d rightPt(x, y, z);
			x = mpInCam1[match.first].x;
			y = mpInCam1[match.first].y;
			z = mpInCam1[match.first].z;
			Eigen::Vector3d midPt(x, y, z);
			double parallaxAngle = triang.CalParallaxAngle(leftPt, midPt, rightPt);
			avgParaAngle += parallaxAngle;
			cout <<"cam: "<<cam1.id<<"_"<<cam2.id<<"; "<<"mappoint: "<< mpInCam1 [match.first].id<<"; "
				<<"parallaxAngle: "<< parallaxAngle << std::endl;
		}
		if (matchers.size() > 0) {
			avgParaAngle /= matchers.size();
			cout << "avgParaAngle: " << avgParaAngle << std::endl;
		}

		//Triangulation Evaluation
		double avgSvdError3D = 0.0;
		int totalSvdNum = 0;
		double avgImprMidpointError3D = 0.0;
		int totalImprMidpointNum = 0;

		for (unsigned int mIdx = 0; mIdx < matchers.size(); mIdx++) {
			std::pair<IdType, IdType> match = matchers[mIdx];
			Eigen::Vector3d mpTrue(mpInCam1[match.first].x, mpInCam1[match.first].y, mpInCam1[match.first].z);
			Eigen::Vector3d obs1(mpInCam1[match.first].u, mpInCam1[match.first].v, 1.0);
			Eigen::Vector3d obs2(mpInCam2[match.second].u, mpInCam2[match.second].v, 1.0);
			Eigen::Vector4d triMP;
			Eigen::Vector3d triMp;
			bool isSVDValid = false;

			//Linear SVD method
			isSVDValid = triang.LinearSVDMethod(param.K, Twc1, Twc2, obs1, obs2, triMP);
			if (isSVDValid) {
				triMp = triMP.segment(0, 3);
				double error = CalculateError3D(mpTrue, triMp);
				if (error < 20.0) {
					avgSvdError3D += error;
					totalSvdNum++;
				}
				cout << "cam: " << cam1.id << "_" << cam2.id << "; " << "mappoint: " << mpInCam1[match.first].id << "; "
					<< "svdTria: " <<triMP[0]<<", "<<triMP[1]<<", "<<triMP[2]<< std::endl;
			}

			//Midpoint method（根据论文公式，三角化出的点不正确）
			/*isSVDValid = triang.MidPointMethod(param.K, Twc1, Twc2, obs1, obs2, triMP);
			if (isSVDValid) {
				cout << "cam: " << cam1.id << "_" << cam2.id << "; " << "mappoint: " << mpInCam1[match.first].id << "; "
					<< "midpointTria: " << triMP[0] << ", " << triMP[1] << ", " << triMP[2] << std::endl;
			}*/

			//Improved midpoint method
			isSVDValid = triang.ImprovedMidPointMethod(param.K, Twc1, Twc2, obs1, obs2, triMP);
			if (isSVDValid) {
				triMp = triMP.segment(0, 3);
				double error = CalculateError3D(mpTrue, triMp);
				if (error < 20.0) {
					avgImprMidpointError3D += error;
						totalImprMidpointNum++;
				}
				cout << "cam: " << cam1.id << "_" << cam2.id << "; " << "mappoint: " << mpInCam1[match.first].id << "; "
					<< "improvedMidpointTria: " << triMP[0] << ", " << triMP[1] << ", " << triMP[2] << std::endl;
			}
		}
		if (totalSvdNum > 0) {
			avgSvdError3D /= totalSvdNum;
			cout << "avgSvdError3D: " << avgSvdError3D << std::endl;
		}
		if (totalImprMidpointNum > 0) {
			avgImprMidpointError3D /= totalImprMidpointNum;
			cout << "avgImprMidpointError3D: " << avgImprMidpointError3D << std::endl;
		}
	}

    return 0;
}
