#include <random>
#include "gene_simulate_data.h"
#include"utility.h"

namespace SLAM_SIMULATION {
	GeneSimulateData::GeneSimulateData(Param p_, std::string inputPath_, std::string outputPath_)
		:params(p_), inputPath(inputPath_), outputPath(outputPath_)
	{
		gyroBias = Eigen::Vector3d::Zero();
		accBias = Eigen::Vector3d::Zero();
	}

	GeneSimulateData::~GeneSimulateData()
	{

	}

	void GeneSimulateData::AddIMUnoise(MotionData& data)
	{
		std::random_device rd;
		std::default_random_engine generator_(rd());
		std::normal_distribution<double> noise(0.0, 1.0);

		Eigen::Vector3d noise_gyro(noise(generator_), noise(generator_), noise(generator_));
		Eigen::Matrix3d gyro_sqrt_cov = params.gyroNoiseSigma * Eigen::Matrix3d::Identity();
		data.imuGyro = data.imuGyro + gyro_sqrt_cov * noise_gyro / sqrt(params.imuTimestep) + gyroBias;

		Eigen::Vector3d noise_acc(noise(generator_), noise(generator_), noise(generator_));
		Eigen::Matrix3d acc_sqrt_cov = params.accNoiseSigma * Eigen::Matrix3d::Identity();
		data.imuAcc = data.imuAcc + acc_sqrt_cov * noise_acc / sqrt(params.imuTimestep) + accBias;

		// gyro_bias update
		Eigen::Vector3d noise_gyro_bias(noise(generator_), noise(generator_), noise(generator_));
		gyroBias += params.gyroBiasSigma * sqrt(params.imuTimestep) * noise_gyro_bias;
		data.imuGyroBias = gyroBias;

		// acc_bias update
		Eigen::Vector3d noise_acc_bias(noise(generator_), noise(generator_), noise(generator_));
		accBias += params.accBiasSigma * sqrt(params.imuTimestep) * noise_acc_bias;
		data.imuAccBias = accBias;
	}

	MotionData GeneSimulateData::MotionModel(double t)
	{
		MotionData data;
		// param
		float ellipseX = 20;
		float ellipseY = 20;
		float z = 0;           // z轴做sin运动
		float K1 = 10;          // z轴的正弦频率是x，y的k1倍
		float K = M_PI / 5.f;    // t * K = 2pi 由于仿真时间是t s, 系数K控制yaw正好旋转一圈，运动一周

		// translation
		// twb:  body frame in world frame
		Eigen::Vector3d position(ellipseX * cos(K * t) + 5, ellipseY * sin(K * t) + 5, z * sin(K1 * K * t) + 5);
		Eigen::Vector3d dp(-K * ellipseX * sin(K*t), K * ellipseY * cos(K*t), z*K1*K * cos(K1 * K * t));              // position导数　in world frame
		double K2 = K * K;
		Eigen::Vector3d ddp(-K2 * ellipseX * cos(K*t), -K2 * ellipseY * sin(K*t), -z * K1*K1*K2 * sin(K1 * K * t));     // position二阶导数

		// Rotation
		double kRoll = 0.1;
		double kPitch = 0.2;
		Eigen::Vector3d eulerAngles(kRoll * cos(t), kPitch * sin(t), K*t);   // roll ~ [-0.2, 0.2], pitch ~ [-0.3, 0.3], yaw ~ [0,2pi]
		Eigen::Vector3d eulerAnglesRates(-kRoll * sin(t), kPitch * cos(t), K);      // euler angles 的导数

	//    Eigen::Vector3d eulerAngles(0.0,0.0, K*t );   // roll ~ 0, pitch ~ 0, yaw ~ [0,2pi]
	//    Eigen::Vector3d eulerAnglesRates(0.,0. , K);      // euler angles 的导数

		Eigen::Matrix3d Rwb = euler2Rotation(eulerAngles);         // body frame to world frame
		Eigen::Vector3d imuGyro = eulerRates2bodyRates(eulerAngles) * eulerAnglesRates;   //  euler rates trans to body gyro

		Eigen::Vector3d gn(0, 0, -9.8); //gravity in navigation frame(ENU)   ENU (0,0,-9.81)  NED(0,0,9,81)
		Eigen::Vector3d imuAcc = Rwb.transpose() * (ddp - gn);  //  Rbw * Rwn * gn = gs

		data.imuGyro = imuGyro;
		data.imuAcc = imuAcc;
		data.Rwb = Rwb;
		data.twb = position;
		data.imuVelocity = dp;
		data.timestamp = t;
		return data;
	}

	void GeneSimulateData::TestImu(std::string src, std::string dist)
	{
		std::vector<MotionData>imudata;
		LoadPose(src, imudata);

		std::ofstream save_points;
		save_points.open(dist);

		double dt = params.imuTimestep;
		Eigen::Vector3d Pwb = init_twb;              // position :    from  imu measurements
		Eigen::Quaterniond Qwb(initRwb);            // quaterniond:  from imu measurements
		Eigen::Vector3d Vw = initVelocity;          // velocity  :   from imu measurements
		Eigen::Vector3d gw(0, 0, -9.81);    // ENU frame
		Eigen::Vector3d temp_a;
		Eigen::Vector3d theta;
		for (int i = 1; i < imudata.size(); ++i) {

			MotionData imupose = imudata[i];

			//delta_q = [1 , 1/2 * thetax , 1/2 * theta_y, 1/2 * theta_z]
			Eigen::Quaterniond dq;
			Eigen::Vector3d dtheta_half = imupose.imuGyro * dt / 2.0;
			dq.w() = 1;
			dq.x() = dtheta_half.x();
			dq.y() = dtheta_half.y();
			dq.z() = dtheta_half.z();
			dq.normalize();

			/// imu 动力学模型 欧拉积分
			Eigen::Vector3d acc_w = Qwb * (imupose.imuAcc) + gw;  // aw = Rwb * ( acc_body - acc_bias ) + gw
			Qwb = Qwb * dq;
			Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
			Vw = Vw + acc_w * dt;

			/// 中值积分

			//　按着imu postion, imu quaternion , cam postion, cam quaternion 的格式存储，由于没有cam，所以imu存了两次
			save_points << imupose.timestamp << " "
				<< Qwb.w() << " "
				<< Qwb.x() << " "
				<< Qwb.y() << " "
				<< Qwb.z() << " "
				<< Pwb(0) << " "
				<< Pwb(1) << " "
				<< Pwb(2) << " "
				<< Qwb.w() << " "
				<< Qwb.x() << " "
				<< Qwb.y() << " "
				<< Qwb.z() << " "
				<< Pwb(0) << " "
				<< Pwb(1) << " "
				<< Pwb(2) << " "
				<< std::endl;
		}

		std::cout << "IMU integration test end" << std::endl;

	}

	bool GeneSimulateData::IsNewPoint(PointsType &points, Eigen::Vector4d &point)
	{
		for (unsigned int idx = 0; idx < points.size(); idx++) {
			if (points[idx].second == point) {
				return false;
			}
		}
		return true;
	}

	void GeneSimulateData::InsertNewPoint(PointsType &points, Eigen::Vector4d &point)
	{
		if (IsNewPoint(points, point)) {
			std::pair<IdType, Eigen::Vector4d> newPoint(NewId(), point);
			points.push_back(newPoint);
		}
	}

	void GeneSimulateData::ReadCamPoseAndObs(std::vector<CamPose> &allCamPose)
	{
		std::ifstream f;
		f.open(inputPath);

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

					camPose.q = Eigen::Quaterniond(camPose.qw, camPose.qx, camPose.qy, camPose.qz);
					camPose.R = camPose.q.toRotationMatrix();
					camPose.tran = Eigen::Vector3d(camPose.tx, camPose.ty, camPose.tz);

					prevObsInOneCam = currObsInOneCam;
					currObsInOneCam.clear();
					if (allCamPose.size() > 0 && numCam >= 1) {
						allCamPose[numCam - 1].mpInCam = prevObsInOneCam;
					}
					allCamPose.push_back(camPose);
					isNewPose = true;
					numCam++;
				}
				else {
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
					}
					else {
						currObsInOneCam.push_back(mp);
					}
				}
			}
		}
	}

	void GeneSimulateData::CreateOrReadPoints(PointsType& points)
	{
		std::ifstream f;
		f.open(inputPath);

		while (!f.eof()) {
			std::string s;
			std::getline(f, s);
			if (!s.empty()) {
				std::stringstream ss;
				ss << s;
				double x, y, z;
				ss >> x;
				ss >> y;
				ss >> z;
				Eigen::Vector4d pt0(x, y, z, 1);
				ss >> x;
				ss >> y;
				ss >> z;
				Eigen::Vector4d pt1(x, y, z, 1);

				InsertNewPoint(points, pt0);
				InsertNewPoint(points, pt1);
			}
		}

		// save points
		SavePoints("all_points.txt", points);
	}

	void GeneSimulateData::CreatePointObs(std::vector<MotionData> &camData, PointsType &points)
	{
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::normal_distribution<double> noise(0.0, params.pixelNoiseSigma);

		std::ofstream saveFile;
		saveFile.open("CamPoseAndObs.txt");

		for (int n = 0; n < camData.size(); ++n) {
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
				double u = obsNorm[0] * params.fx + params.cx;
				double v = obsNorm[1] * params.fy + params.cy;

				//Create rolling shutter observation
				if (params.rsT > 0.0 && v > 0. && v < params.image_h) {
					unsigned int row = static_cast<unsigned int>(std::floor(v));
					Eigen::Matrix4d iTwc = Eigen::Matrix4d::Identity();
					iTwc.block(0, 0, 3, 3) = data.rsPoses[row].Rwb;
					iTwc.block(0, 3, 3, 1) = data.rsPoses[row].twb;
					Eigen::Vector4d ipc1 = iTwc.inverse()*pw;
					if (ipc1(2) < 0) {
						continue;
					}
					Eigen::Vector2d iObsNorm(ipc1(0) / ipc1(2), ipc1(1) / ipc1(2));
					u = iObsNorm[0] * params.fx + params.cx;
					v = iObsNorm[1] * params.fy + params.cy;
				}

				if (params.addPixedNoise) {
					u += noise(generator);
					v += noise(generator);
				}
				if (u > 0. && u<params.image_w && v>0. && v < params.image_h) {
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
			filename1 << outputPath << "\\" << "AllObsInImg_time_" << n << ".txt";
			SaveFeatures(filename1.str(), data, pointsCam, featuresCam);
			SavaPoseAndFeatures(saveFile, data, pointsCam, featuresCam);
		}
	}

	void GeneSimulateData::CreateIMUPoseAndObs(std::vector< MotionData > &imuData, std::vector< MotionData > &imuDataNoise)
	{
		for (float t = params.tStart; t < params.tEnd;) {
			MotionData imu = MotionModel(t);
			imuData.push_back(imu);

			// add imu noise
			MotionData imuNoise = imu;
			AddIMUnoise(imuNoise);
			imuDataNoise.push_back(imuNoise);

			t += 1.0 / params.imuFrequency;
		}

		initVelocity = imuData[0].imuVelocity;
		init_twb = imuData[0].twb;
		initRwb = imuData[0].Rwb;

		savePose("imu_pose.txt", imuData);
		savePose("imu_pose_noise.txt", imuDataNoise);

		TestImu("imu_pose.txt", "imu_integ_pose.txt");     // test the imu data, integrate the imu data to generate the imu trajecotry
		TestImu("imu_pose_noise.txt", "imu_integ_pose_noise.txt");
	}

	void GeneSimulateData::CreateCamPoseAndObs(PointsType &points, std::vector< MotionData > &camData)
	{
		// cam pose
		//std::vector< MotionData > camData;
		for (float t = params.tStart; t < params.tEnd;) {
			MotionData imu = MotionModel(t);   // imu body frame to world frame motion
			MotionData cam;

			cam.timestamp = imu.timestamp;
			cam.Rwb = imu.Rwb * params.Rbc;    // cam frame in world frame: Rwb actually is Rwc
			cam.twb = imu.twb + imu.Rwb * params.tbc; //  Twc = Twb * Tbc ,  t = Rwb * tbc + twb
			cam.imuVelocity = imu.imuVelocity;

			//Create rolling shutter row poses
			if (params.rsT > 0.0 && params.image_h > 0.0) {
				double readoutTime = params.rsT / params.image_h;
				for (unsigned int row = 0; row < static_cast<unsigned int>(params.image_h); row++) {
					double iTime = t + row * readoutTime;
					RowPose iPose;
					MotionData iImu = MotionModel(iTime);
					iPose.row = row;
					iPose.timestamp = iTime;
					iPose.Rwb = iImu.Rwb*params.Rbc;
					iPose.twb = iImu.twb + iImu.Rwb*params.tbc;
					iPose.qwb = Eigen::Quaterniond(iPose.Rwb);
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
}