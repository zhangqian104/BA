#include <fstream>
#include <sys/stat.h>
#include <map>
#include<set>
#include<random>
#include "common.h"
#include "gene_simulate_data.h"
#include "utility.h"
#include "triangulation.h"
#include "nonlinear_optimization_solver.h"

using namespace SLAM_SIMULATION;

std::string houseModelPath = "C:\\Users\\Administrator\\Desktop\\simulation\\bin\\house_model\\house.txt";
std::string outputDataPath = "C:\\Users\\Administrator\\Desktop\\simulation\\simulation\\OutputData";
std::string camPoseAndObsPath = "CamPoseAndObs.txt";

bool NonlinearOptimizationSolver(std::vector<KeyFrame> &keyFrames)
{
	int numKFs = keyFrames.size();
	if (numKFs <= 0) {
		return false;
	}

	NonlinearSolver::LossFunction *lossFunction;
	lossFunction = new NonlinearSolver::CauchyLoss(1.0);

	// step1. 构建 problem
	NonlinearSolver::Problem problem(NonlinearSolver::Problem::ProblemType::SLAM_PROBLEM);
	std::vector<std::shared_ptr<NonlinearSolver::VertexPose>> vertexCams_vec;
	std::vector<std::shared_ptr<NonlinearSolver::VertexSpeedBias>> vertexVB_vec;
	int pose_dim = 0;

	for (int i = 0; i < numKFs + 1; i++){
		KeyFrame kf = keyFrames[i];
		std::shared_ptr<NonlinearSolver::VertexPose> vertexCam(new NonlinearSolver::VertexPose());
		Eigen::VectorXd pose(7);
		pose << kf.twcNoise[0], kf.twcNoise[1], kf.twcNoise[2], kf.qwc.x, kf.qwc.y, kf.qwc.z, kf.qwc.w;
		vertexCam->SetParameters(pose);
		vertexCams_vec.push_back(vertexCam);
		problem.AddVertex(vertexCam);
		pose_dim += vertexCam->LocalDimension();

		std::shared_ptr<NonlinearSolver::VertexSpeedBias> vertexVB(new NonlinearSolver::VertexSpeedBias());
		Eigen::VectorXd vb(9);
		vb << kf.vel[0], kf.vel[1], kf.vel[2], kf.accBiass.back()[0], kf.accBiass.back()[1], kf.accBiass.back()[2],
			kf.gyroBiass.back()[0], kf.gyroBiass.back()[1], kf.gyroBiass.back()[2];
		vertexVB->SetParameters(vb);
		vertexVB_vec.push_back(vertexVB);
		problem.AddVertex(vertexVB);
		pose_dim += vertexVB->LocalDimension();
	}

	// IMU
	for (int i = 1; i < numKFs; i++)
	{
		KeyFrame prevKF = keyFrames[i-1];
		KeyFrame currKF = keyFrames[i];
		int j = i + 1;
		if (currKF.timestamp-prevKF.timestamp > 10.0)
			continue;

		std::shared_ptr<NonlinearSolver::EdgeImu> imuEdge(new NonlinearSolver::EdgeImu(currKF.imuPreintgB));
		std::vector<std::shared_ptr<NonlinearSolver::Vertex>> edge_vertex;
		edge_vertex.push_back(vertexCams_vec[i-1]);
		edge_vertex.push_back(vertexVB_vec[i-1]);
		edge_vertex.push_back(vertexCams_vec[i]);
		edge_vertex.push_back(vertexVB_vec[i]);
		imuEdge->SetVertex(edge_vertex);
		problem.AddEdge(imuEdge);
	}

	// Visual Factor
	std::vector<std::shared_ptr<NonlinearSolver::VertexInverseDepth>> vertexPt_vec;
	{
		int feature_index = -1;
		// 遍历每一个特征
		for (auto &it_per_id : f_manager.feature)
		{
			it_per_id.used_num = it_per_id.feature_per_frame.size();
			if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < numKFs - 2))
				continue;

			++feature_index;

			int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
			Vector3d pts_i = it_per_id.feature_per_frame[0].point;

			std::shared_ptr<NonlinearSolver::VertexInverseDepth> verterxPoint(new NonlinearSolver::VertexInverseDepth());
			Eigen::Matrix<double, Eigen::Dynamic, 1> inv_d(1);
			inv_d << para_Feature[feature_index][0];
			verterxPoint->SetParameters(inv_d);
			problem.AddVertex(verterxPoint);
			vertexPt_vec.push_back(verterxPoint);

			// 遍历所有的观测
			for (auto &it_per_frame : it_per_id.feature_per_frame)
			{
				imu_j++;
				if (imu_i == imu_j)
					continue;

				Vector3d pts_j = it_per_frame.point;

				std::shared_ptr<NonlinearSolver::EdgeReprojection> edge(new NonlinearSolver::EdgeReprojection(pts_i, pts_j));
				std::vector<std::shared_ptr<NonlinearSolver::Vertex>> edge_vertex;
				edge_vertex.push_back(verterxPoint);
				edge_vertex.push_back(vertexCams_vec[imu_i]);
				edge_vertex.push_back(vertexCams_vec[imu_j]);

				edge->SetVertex(edge_vertex);
				edge->SetInformation(project_sqrt_info_.transpose() * project_sqrt_info_);

				edge->SetLossFunction(lossFunction);
				problem.AddEdge(edge);
			}
		}
	}


	problem.Solve(10);

	// update parameter
	for (int i = 0; i < numKFs + 1; i++)
	{
		Eigen::Matrix<double, Eigen::Dynamic, 1> p = vertexCams_vec[i]->Parameters();
		for (int j = 0; j < 7; ++j)
		{
			para_Pose[i][j] = p[j];
		}

		Eigen::Matrix<double, Eigen::Dynamic, 1> vb = vertexVB_vec[i]->Parameters();
		for (int j = 0; j < 9; ++j)
		{
			para_SpeedBias[i][j] = vb[j];
		}
	}

	// 遍历每一个特征
	for (int i = 0; i < vertexPt_vec.size(); ++i)
	{
		Eigen::Matrix<double, Eigen::Dynamic, 1> f = vertexPt_vec[i]->Parameters();
		para_Feature[i][0] = f[0];
	}

	return true;
}

int main(){
    // Read Model and Get 3d points and lines
	Param params;
	GeneSimulateData geneData(params,houseModelPath,outputDataPath);
	PointsType points;
	geneData.CreateOrReadPoints(points);
	std::vector< MotionData > imuData;
	std::vector< MotionData > imuDataNoise;
	geneData.CreateIMUPoseAndObs(imuData,imuDataNoise);
	std::vector< MotionData > camData;
	geneData.CreateCamPoseAndObs(points,camData);

	//Key frames for lba unit test
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
		std::vector<Eigen::Vector3d> imuAccBiass;
		std::vector<Eigen::Vector3d> imuGyroBiass;
		for (unsigned int idx2 = 0; idx2 < imuDataNoise.size(); idx2++) {
			MotionData imu = imuDataNoise[idx2];
			double imuTime = imu.timestamp;
			if (imuTime > prevCamTime && imuTime <= currCamTime) {
				imuAccs.push_back(imu.imuAcc);
				imuGyros.push_back(imu.imuGyro);
				imuAccBiass.push_back(imu.imuAccBias);
				imuGyroBiass.push_back(imu.imuGyroBias);
			}
		}

		keyFrame.timestamp = currCamTime;
		keyFrame.imuAccs = imuAccs;
		keyFrame.imuGyros = imuGyros;
		keyFrame.accBiass = imuAccBiass;
		keyFrame.gyroBiass = imuGyroBiass;
		keyFrame.rsPoses = currCam.rsPoses;
		keyFrame.pointsCam = currCam.pointsCam;
		keyFrame.featuresCam = currCam.featuresCam;
		keyFrame.Rwc = currCam.Rwb;
		keyFrame.qwc = Eigen::Quaterniond (currCam.Rwb);
		keyFrame.twc = currCam.twb;
		keyFrame.vel = currCam.imuVelocity;
		Eigen::Vector3d poseNoise(noise(generator), noise(generator), noise(generator));
		keyFrame.twcNoise = Eigen::Vector3d(keyFrame.twc+poseNoise);
		keyFrames.push_back(keyFrame);
	}

	NonlinearOptimizationSolver(keyFrames);

//#define TRIANGULATION_EVALUATION
#ifdef TRIANGULATION_EVALUATION
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
#endif

    return 0;
}
