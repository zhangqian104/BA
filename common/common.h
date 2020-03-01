#ifndef COMMON_H
#define COMMON_H

#include<vector>
#include <map>
#include<set>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Geometry>

const float M_PI = 3.14159f;
const Eigen::Vector3d G = Eigen::Vector3d (0.0, 0.0, 9.8);

enum StateOrder
{
	O_P = 0,
	O_R = 3,
	O_V = 6,
	O_BA = 9,
	O_BG = 12
};

typedef unsigned long IdType;
typedef std::pair<IdType, Eigen::Vector4d> PointType;
typedef std::vector<PointType, Eigen::aligned_allocator<PointType> > PointsType;
typedef std::pair<IdType, Eigen::Vector2d> ObserType;
typedef std::vector<ObserType, Eigen::aligned_allocator<ObserType>> ObsersType;
typedef std::pair<Eigen::Vector4d, Eigen::Vector4d> LineType;
typedef std::vector<LineType, Eigen::aligned_allocator<LineType> > LinesType;

struct MapPoint
{
	IdType id;

	double x = -1.0;
	double y = -1.0;
	double z = -1.0;

	double u=-1.0;
	double v=-1.0;
};

struct RowPose
{
	unsigned int row;
	double timestamp;
	Eigen::Quaterniond qwb;
	Eigen::Vector3d twb;
	Eigen::Matrix3d Rwb;
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



struct IMUData
{
	double gyroNoiseSigma = 1.7e-4;
	double accNoiseSigma = 2.0e-3;

	// covariance of measurement
	Eigen::Matrix3d gyroNoiseCov = Eigen::Matrix3d::Identity()*gyroNoiseSigma*gyroNoiseSigma / 0.005/**100*/;       // sigma_g * sigma_g / dt, ~6e-6*10
	Eigen::Matrix3d accNoiseCov = Eigen::Matrix3d::Identity()*accNoiseSigma*accNoiseSigma / 0.005 * 100;       // sigma_a * sigma_a / dt, ~8e-4*10

	double gyroBiasRw = 2.0e-5/**10*/;  //2e-12*1e3
	double accoBiasRw = 5.0e-3/**10*/;  //4.5e-8*1e2

	// covariance of bias random walk
	Eigen::Matrix3d gyrBiasRWCov = Eigen::Matrix3d::Identity()*gyroBiasRw*gyroBiasRw;     // sigma_gw * sigma_gw * dt, ~2e-12
	Eigen::Matrix3d accBiasRWCov = Eigen::Matrix3d::Identity()*accoBiasRw*accoBiasRw;     // sigma_aw * sigma_aw * dt, ~4.5e-8

	// Raw data of imu's
	Eigen::Vector3d gyro;    //gyr data
	Eigen::Vector3d acc;    //acc data
	double t;      //timestamp
};

#endif