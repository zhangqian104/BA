#ifndef COMMON_H
#define COMMON_H

#include<vector>
#include <map>
#include<set>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Geometry>
#include "imu_integration.h"

constexpr float M_PI = 3.14159;
#define G Eigen::Vector3d (0.0, 0.0, 9.8);

typedef unsigned long IdType;
typedef std::pair<IdType, Eigen::Vector4d> PointType;
typedef std::vector<PointType, Eigen::aligned_allocator<PointType> > PointsType;
typedef std::pair<IdType, Eigen::Vector2d> ObserType;
typedef std::vector<ObserType, Eigen::aligned_allocator<ObserType>> ObsersType;
typedef std::pair<Eigen::Vector4d, Eigen::Vector4d> LineType;
typedef std::vector<LineType, Eigen::aligned_allocator<LineType> > LinesType;

extern IdType g_Id = 0;
extern void ResetId()
{
	g_Id = 0;
}
extern IdType NewId()
{
	return g_Id++;
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

struct KeyFrame
{
	double timestamp;
	std::vector<Eigen::Vector3d> imuAccs;
	std::vector<Eigen::Vector3d> imuGyros;
	std::vector<Eigen::Vector3d> accBiass;
	std::vector<Eigen::Vector3d> gyroBiass;
	std::vector<RowPose> rsPoses;
	PointsType pointsCam;    // 3D points in current cam
	ObsersType featuresCam;  // corresponding 2D observations

	Eigen::Matrix3d Rwc;
	Eigen::Quaterniond qwc;
	Eigen::Vector3d twc;
	Eigen::Vector3d twcNoise;
	Eigen::Vector3d vel;

	/*IMUPreintegRecursive imuPreintgR;
	IMUPreintegBatch imuPreintgB;*/
};

struct IMUData
{
	// covariance of measurement
	Eigen::Matrix3d gyrMeasCov = Eigen::Matrix3d::Identity()*1.7e-4*1.7e-4 / 0.005/**100*/;       // sigma_g * sigma_g / dt, ~6e-6*10
	Eigen::Matrix3d accMeasCov = Eigen::Matrix3d::Identity()*2.0e-3*2.0e-3 / 0.005 * 100;       // sigma_a * sigma_a / dt, ~8e-4*10

	double gyrBiasRw2 = 2.0e-5*2.0e-5/**10*/;  //2e-12*1e3
	double accBiasRw2 = 5.0e-3*5.0e-3/**10*/;  //4.5e-8*1e2

	// covariance of bias random walk
	Eigen::Matrix3d gyrBiasRWCov = Eigen::Matrix3d::Identity()*gyrBiasRw2;     // sigma_gw * sigma_gw * dt, ~2e-12
	Eigen::Matrix3d accBiasRWCov = Eigen::Matrix3d::Identity()*accBiasRw2;     // sigma_aw * sigma_aw * dt, ~4.5e-8

	// Raw data of imu's
	Eigen::Vector3d gyro;    //gyr data
	Eigen::Vector3d acc;    //acc data
	double t;      //timestamp
};

#endif