#include "vertex_pose.h"

namespace NonlinearSolver {

void VertexPose::Plus(const Eigen::Matrix<double, Eigen::Dynamic, 1> &delta) {
	Eigen::Matrix<double, Eigen::Dynamic, 1> &parameters = Parameters();
    parameters.head<3>() += delta.head<3>();
	Eigen::Quaterniond q(parameters[6], parameters[3], parameters[4], parameters[5]);
    q = q * Sophus::SO3d::exp(Vec3(delta[3], delta[4], delta[5])).unit_quaternion();  // right multiplication with so3
    q.normalized();
    parameters[3] = q.x();
    parameters[4] = q.y();
    parameters[5] = q.z();
    parameters[6] = q.w();
//    Qd test = Sophus::SO3d::exp(Vec3(0.2, 0.1, 0.1)).unit_quaternion() * Sophus::SO3d::exp(-Vec3(0.2, 0.1, 0.1)).unit_quaternion();
//    std::cout << test.x()<<" "<< test.y()<<" "<<test.z()<<" "<<test.w() <<std::endl;
}

}
