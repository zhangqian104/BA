#include "triangulation.h"
#include "utilities.h"
#include <eigen3/Eigen/SVD>
#include <eigen3/Eigen/Dense>

Triangulation::Triangulation() {}

Triangulation::~Triangulation() {}

bool Triangulation::MidPointMethod(const Eigen::Matrix3d &K, const Eigen::Matrix<double, 3, 4> &T1, const Eigen::Matrix<double, 3, 4> &T2,
	const Eigen::Vector3d &obs1, const Eigen::Vector3d &obs2, Eigen::Vector4d &x)
{
	Eigen::Matrix3d R1 = T1.block<3, 3>(0, 0);
	Eigen::Vector3d t1 = T1.block<3, 1>(0, 3);
	Eigen::Matrix3d R2 = T2.block<3, 3>(0, 0);
	Eigen::Vector3d t2 = T2.block<3, 1>(0, 3);

	Eigen::Matrix3d R21 = R2.transpose()*R1;
	Eigen::Vector3d t21 = R2.transpose()*(t1 - t2);

	Eigen::Vector3d f1 = K.inverse()*obs1;
	Eigen::Vector3d f2 = K.inverse()*obs2;

	Eigen::Vector3d Rxf1 = R21 * f1;
	Eigen::Vector3d p = Rxf1.cross(f2);
	Eigen::Vector3d q = Rxf1.cross(t21);
	Eigen::Vector3d r = f2.cross(t21);

	if (p.norm() < 10e-6) {
		return false;
	}

	double d1 = abs(p.dot(r)) / p.norm();
	double d2 = abs(p.dot(q)) / p.norm();

	Eigen::Vector3d x0 = 0.5*(t21 + d1 * R21*f1 + d2 * f2);

	x0 = R2 * x0 + t2;

	x.block<3, 1>(0, 0) = x0;
	x[3] = 1.0;
	return true;
}

bool Triangulation::ImprovedMidPointMethod(const Eigen::Matrix3d &K, const Eigen::Matrix<double, 3, 4> &T1, const Eigen::Matrix<double, 3, 4> &T2,
	const Eigen::Vector3d &obs1, const Eigen::Vector3d &obs2, Eigen::Vector4d &x)
{
	Eigen::Matrix3d R1 = T1.block<3, 3>(0, 0);
	Eigen::Vector3d t1 = T1.block<3, 1>(0, 3);
	Eigen::Matrix3d R2 = T2.block<3, 3>(0, 0);
	Eigen::Vector3d t2 = T2.block<3, 1>(0, 3);

	Eigen::Matrix3d R21 = R2.transpose()*R1;
	Eigen::Vector3d t21 = R2.transpose()*(t1 - t2);

	Eigen::Vector3d f1 = K.inverse()*obs1;
	Eigen::Vector3d f2 = K.inverse()*obs2;

	Eigen::Vector3d Rxf1 = R21 * f1;
	Eigen::Vector3d p = Rxf1.cross(f2);
	Eigen::Vector3d q = Rxf1.cross(t21);
	Eigen::Vector3d r = f2.cross(t21);

	if (p.norm() < 10e-6) {
		return false;
	}

	double d1 = r.norm() / p.norm();
	double d2 = q.norm() / p.norm();

	Eigen::Vector3d x0 = 0.5*(t21 + d1 * R21*f1 + d2 * f2);

	x0 = R2 *x0+t2;

	x.block<3, 1>(0, 0) = x0;
	x[3] = 1.0;
	return true;
}

bool Triangulation::LinearSVDMethod(const Eigen::Matrix3d &K, const Eigen::Matrix<double, 3, 4> &T1,
	const Eigen::Matrix<double, 3, 4> &T2, const Eigen::Vector3d &obs1, const Eigen::Vector3d &obs2, Eigen::Vector4d &x)
{
	Eigen::Matrix3d R1 = T1.block<3, 3>(0, 0);
	Eigen::Vector3d t1 = T1.block<3, 1>(0, 3);
	Eigen::Matrix3d R2 = T2.block<3, 3>(0, 0);
	Eigen::Vector3d t2 = T2.block<3, 1>(0, 3);
	Eigen::Matrix<double, 3, 4> T1tran;
	T1tran << R1.transpose(), -R1.transpose()*t1;
	Eigen::Matrix<double, 3, 4> T2tran;
	T2tran << R2.transpose(), -R2.transpose()*t2;

	Eigen::Matrix3d obsSkewSymmMat1 = GetSkewSymmMatrix(obs1);
	Eigen::Matrix3d obsSkewSymmMat2 = GetSkewSymmMatrix(obs2);
	Eigen::Matrix<double, 3, 4> obsKxT1 = obsSkewSymmMat1 * K * T1tran;
	Eigen::Matrix<double, 3, 4> obsKxT2 = obsSkewSymmMat2 * K * T2tran;
	Eigen::Matrix4d A;
	A.block<2, 4>(0, 0) = obsKxT1.block<2, 4>(0, 0);
	A.block<2, 4>(2, 0) = obsKxT2.block<2, 4>(0, 0);
	Eigen::JacobiSVD<Eigen::Matrix4d> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	auto V = svd.matrixV();
	//auto U = svd.matrixU();
	//auto S = svd.singularValues();

	/*Eigen::Vector4d triangulatedPoint;
	triangulatedPoint = A.jacobiSvd(Eigen::ComputeFullV).matrixV().rightCols<1>();*/

	x = V.col(3);
	if (std::fabs(x[3]) < 10e-6) {
		return false;
	}
	x /= x[3];

	return true;
}

double Triangulation::CalParallaxAngle(const Eigen::Vector3d &leftPt, const Eigen::Vector3d &midPt, const Eigen::Vector3d &rightPt)
{
	Eigen::Vector3d v1 = leftPt - midPt;
	Eigen::Vector3d v2 = rightPt - midPt;

	double parallaxAngle = atan2(v1.cross(v2).norm(), v1.transpose() * v2);
	if (v1.cross(v2).z() < 0) {
		parallaxAngle = 2 * M_PI - parallaxAngle;
	}
	return parallaxAngle * 180 / M_PI;
}