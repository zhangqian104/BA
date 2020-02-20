#ifndef IMUSIMWITHPOINTLINE_IMU_H
#define IMUSIMWITHPOINTLINE_IMU_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <vector>

#include "param.h"

using IdType = unsigned long;
using PointType = std::pair<IdType, Eigen::Vector4d>;
using PointsType = std::vector<PointType, Eigen::aligned_allocator<PointType> >;
using ObserType = std::pair<IdType, Eigen::Vector2d>;
using ObsersType = std::vector<ObserType, Eigen::aligned_allocator<ObserType>>;
using LineType = std::pair<Eigen::Vector4d, Eigen::Vector4d>;
using LinesType = std::vector<LineType, Eigen::aligned_allocator<LineType> >;

struct RowPose
{
	unsigned int row;
	double timestamp;
	Eigen::Quaterniond qwb;
	Eigen::Vector3d twb;
	Eigen::Matrix3d Rwb;
};

struct MotionData
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    double timestamp;

    Eigen::Matrix3d Rwb;
    Eigen::Vector3d twb;

    Eigen::Vector3d imuAcc;
    Eigen::Vector3d imuGyro;

    Eigen::Vector3d imuGyroBias;
    Eigen::Vector3d imuAccBias;

    Eigen::Vector3d imuVelocity;

	std::vector<RowPose> rsPoses;

	PointsType pointsCam;    // ３维点在当前cam视野里
	ObsersType featuresCam;  // 对应的２维图像坐标
};

// euler2Rotation:   body frame to interitail frame
Eigen::Matrix3d euler2Rotation( Eigen::Vector3d  eulerAngles);
Eigen::Matrix3d eulerRates2bodyRates(Eigen::Vector3d eulerAngles);


class IMU
{
public:
    IMU(Param p);
    Param param;
    Eigen::Vector3d gyroBias;
    Eigen::Vector3d accBias;

    Eigen::Vector3d initVelocity;
    Eigen::Vector3d init_twb;
    Eigen::Matrix3d initRwb;

    MotionData MotionModel(double t);

    void AddIMUnoise(MotionData& data);
    void testImu(std::string src, std::string dist);        // imu数据进行积分，用来看imu轨迹

};

#endif //IMUSIMWITHPOINTLINE_IMU_H
