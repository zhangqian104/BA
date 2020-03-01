#ifndef MYSLAM_BACKEND_IMUEDGE_H
#define MYSLAM_BACKEND_IMUEDGE_H

#include <memory>
#include <string>
#include "imu_integration.h"
#include "edge.h"

namespace NonlinearSolver {

/**
 * 此边是IMU误差，此边为4元边，与之相连的顶点有：Pi Mi Pj Mj
 */
class EdgeImu : public Edge {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    explicit EdgeImu(SLAM_SIMULATION::IMUPreintegRecursive* _pre_integration):pre_integration_(_pre_integration),
          Edge(15, 4, std::vector<std::string>{"VertexPose", "VertexSpeedBias", "VertexPose", "VertexSpeedBias"}) {
//        if (pre_integration_) {
//            pre_integration_->GetJacobians(dr_dbg_, dv_dbg_, dv_dba_, dp_dbg_, dp_dba_);
//            Mat99 cov_meas = pre_integration_->GetCovarianceMeasurement();
//            Mat66 cov_rand_walk = pre_integration_->GetCovarianceRandomWalk();
//            Mat1515 cov = Mat1515::Zero();
//            cov.block<9, 9>(0, 0) = cov_meas;
//            cov.block<6, 6>(9, 9) = cov_rand_walk;
//            SetInformation(cov.inverse());
//        }
    }

    /// 返回边的类型信息
    virtual std::string TypeInfo() const override { return "EdgeImu"; }

    /// 计算残差
    virtual void ComputeResidual() override;

    /// 计算雅可比
    virtual void ComputeJacobians() override;

//    static void SetGravity(const Vec3 &g) {
//        gravity_ = g;
//    }

private:
	SLAM_SIMULATION::IMUPreintegRecursive* pre_integration_;
    static Eigen::Vector3d gravity_;

	Eigen::Matrix3d dp_dba_ = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dp_dbg_ = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dr_dbg_ = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dv_dba_ = Eigen::Matrix3d::Zero();
	Eigen::Matrix3d dv_dbg_ = Eigen::Matrix3d::Zero();
};

}

#endif
