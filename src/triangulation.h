#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#define M_PI 3.14159

class Triangulation {
public:
	Triangulation();
	~Triangulation();

	bool MidPointMethod(const Eigen::Matrix3d &K, const Eigen::Matrix<double, 3, 4> &T1, const Eigen::Matrix<double, 3, 4> &T2,
		const Eigen::Vector3d &obs1, const Eigen::Vector3d &obs2, Eigen::Vector4d &x);
	bool ImprovedMidPointMethod(const Eigen::Matrix3d &K, const Eigen::Matrix<double, 3, 4> &T1, const Eigen::Matrix<double, 3, 4> &T2,
		const Eigen::Vector3d &obs1, const Eigen::Vector3d &obs2, Eigen::Vector4d &x);
	bool LinearSVDMethod(const Eigen::Matrix3d &K, const Eigen::Matrix<double, 3, 4> &T1, const Eigen::Matrix<double,3,4> &T2,
		const Eigen::Vector3d &obs1, const Eigen::Vector3d &obs2,Eigen::Vector4d &x);

	double CalParallaxAngle(const Eigen::Vector3d &leftPt, const Eigen::Vector3d &midPt, const Eigen::Vector3d &rightPt);

private:
	
};

#endif