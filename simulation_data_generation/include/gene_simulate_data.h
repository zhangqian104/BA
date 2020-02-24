#ifndef GENE_SIMULATE_DATA_H
#define GENE_SIMULATE_DATA_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <vector>
#include"common.h"
#include "param.h"

namespace SLAM_SIMULATION {
	class GeneSimulateData
	{
	public:
		GeneSimulateData(Param p, std::string inputPath, std::string outputPath);
		~GeneSimulateData();

		MotionData MotionModel(double t);

		void AddIMUnoise(MotionData& data);
		void TestImu(std::string src, std::string dist);        // imu数据进行积分，用来看imu轨迹

		bool IsNewPoint(PointsType &points, Eigen::Vector4d &point);
		void InsertNewPoint(PointsType &points, Eigen::Vector4d &point);
		void ReadCamPoseAndObs(std::vector<CamPose> &allCamPose);
		void CreateOrReadPoints(PointsType& points);
		void CreatePointObs(std::vector<MotionData> &camData, PointsType &points);
		void CreateIMUPoseAndObs(std::vector< MotionData > &imuData, std::vector< MotionData > &imuDataNoise);
		void CreateCamPoseAndObs(PointsType &points, std::vector< MotionData > &camData);

	private:
		Param params;
		std::string inputPath;
		std::string outputPath;

		Eigen::Vector3d gyroBias;
		Eigen::Vector3d accBias;
		Eigen::Vector3d initVelocity;
		Eigen::Vector3d init_twb;
		Eigen::Matrix3d initRwb;
	};
}
#endif //IMUSIMWITHPOINTLINE_IMU_H
