#ifndef IMU_INTEGRATION_H
#define IMU_INTEGRATION_H

#include <vector>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include "so3.h"
#include "common.h"
#include "utility.h"
#include "param.h"

namespace SLAM_SIMULATION {

class IMUPreintegRecursive
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	IMUPreintegRecursive(IMUData imuNoise_, const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
		const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg);

	void AddNewImu(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr);

	void Repropagate(const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg);

	void MidpointIntegration(double _dt,
		const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
		const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
		const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
		const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg,
		Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
		Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg, bool update_jacobian);

	void EulerIntegration(double _dt, const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
		const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
		const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
		const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg,
		Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
		Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg, bool update_jacobian);

	void Propagate(double _dt, const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1);

	Eigen::Matrix<double, 15, 1> Evaluate(const Eigen::Vector3d &Pi, const Eigen::Quaterniond &Qi, const Eigen::Vector3d &Vi, const Eigen::Vector3d &Bai, const Eigen::Vector3d &Bgi,
		const Eigen::Vector3d &Pj, const Eigen::Quaterniond &Qj, const Eigen::Vector3d &Vj, const Eigen::Vector3d &Baj, const Eigen::Vector3d &Bgj)
	{
		Eigen::Matrix<double, 15, 1> residuals;

		Eigen::Matrix3d dp_dba = jacobian.block<3, 3>(O_P, O_BA);
		Eigen::Matrix3d dp_dbg = jacobian.block<3, 3>(O_P, O_BG);

		Eigen::Matrix3d dq_dbg = jacobian.block<3, 3>(O_R, O_BG);

		Eigen::Matrix3d dv_dba = jacobian.block<3, 3>(O_V, O_BA);
		Eigen::Matrix3d dv_dbg = jacobian.block<3, 3>(O_V, O_BG);

		Eigen::Vector3d dba = Bai - linearized_ba;
		Eigen::Vector3d dbg = Bgi - linearized_bg;

		Eigen::Quaterniond corrected_delta_q = delta_q * Utility::deltaQ(dq_dbg * dbg);
		Eigen::Vector3d corrected_delta_v = delta_v + dv_dba * dba + dv_dbg * dbg;
		Eigen::Vector3d corrected_delta_p = delta_p + dp_dba * dba + dp_dbg * dbg;

		residuals.block<3, 1>(O_P, 0) = Qi.inverse() * (0.5 * G * sum_dt * sum_dt + Pj - Pi - Vi * sum_dt) - corrected_delta_p;
		residuals.block<3, 1>(O_R, 0) = 2 * (corrected_delta_q.inverse() * (Qi.inverse() * Qj)).vec();
		residuals.block<3, 1>(O_V, 0) = Qi.inverse() * (G * sum_dt + Vj - Vi) - corrected_delta_v;
		residuals.block<3, 1>(O_BA, 0) = Baj - Bai;
		residuals.block<3, 1>(O_BG, 0) = Bgj - Bgi;
		return residuals;
	}

public:
	double dt;
	Eigen::Vector3d acc_0, gyr_0;
	Eigen::Vector3d acc_1, gyr_1;

	Eigen::Vector3d linearized_acc, linearized_gyr;
	Eigen::Vector3d linearized_ba, linearized_bg;

	Eigen::Matrix<double, 15, 15> jacobian, covariance;
	Eigen::Matrix<double, 15, 15> step_jacobian;
	Eigen::Matrix<double, 15, 18> step_V;
	Eigen::Matrix<double, 18, 18> noise;

	double sum_dt;
	Eigen::Vector3d delta_p;
	Eigen::Quaterniond delta_q;
	Eigen::Vector3d delta_v;

	std::vector<double> dt_buf;
	std::vector<Eigen::Vector3d> acc_buf;
	std::vector<Eigen::Vector3d> gyr_buf;

private:
	IMUData imuNoise;
};

class IMUPreintegBatch
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	IMUPreintegBatch(IMUData imuNoise_);
	IMUPreintegBatch(const IMUPreintegBatch& pre);

	//IMUPreintegrator& operator = (const IMUPreintegrator& pre);

	// reset to initial state
	void reset();

	// incrementally update 1)delta measurements, 2)jacobians, 3)covariance matrix
	void update(const Eigen::Vector3d& omega, const Eigen::Vector3d& acc, const double& dt);

	// delta measurements, position/velocity/rotation(matrix)
	inline Eigen::Vector3d getDeltaP() const    // P_k+1 = P_k + V_k*dt + R_k*a_k*dt*dt/2
	{
		return _delta_P;
	}
	inline Eigen::Vector3d getDeltaV() const    // V_k+1 = V_k + R_k*a_k*dt
	{
		return _delta_V;
	}
	inline Eigen::Matrix3d getDeltaR() const   // R_k+1 = R_k*exp(w_k*dt).     NOTE: Rwc, Rwc'=Rwc*[w_body]x
	{
		return _delta_R;
	}

	// jacobian of delta measurements w.r.t bias of gyro/acc
	inline Eigen::Matrix3d getJPBiasg() const     // position / gyro
	{
		return _J_P_Biasg;
	}
	inline Eigen::Matrix3d getJPBiasa() const     // position / acc
	{
		return _J_P_Biasa;
	}
	inline Eigen::Matrix3d getJVBiasg() const     // velocity / gyro
	{
		return _J_V_Biasg;
	}
	inline Eigen::Matrix3d getJVBiasa() const     // velocity / acc
	{
		return _J_V_Biasa;
	}
	inline Eigen::Matrix3d getJRBiasg() const  // rotation / gyro
	{
		return _J_R_Biasg;
	}

	// noise covariance propagation of delta measurements
	// note: the order is rotation-velocity-position here
	inline Eigen::Matrix<double, 9, 9> getCovPVPhi() const
	{
		return _cov_P_V_Phi;
	}

	inline double getDeltaTime() const {
		return deltaTime;
	}

	// skew-symmetric matrix
	static Eigen::Matrix3d skew(const Eigen::Vector3d& v)
	{
		return Sophus::SO3::hat(v);
	}

	// exponential map from vec3 to mat3x3 (Rodrigues formula)
	static Eigen::Matrix3d Expmap(const Eigen::Vector3d& v)
	{
		return Sophus::SO3::exp(v).matrix();
	}

	// right jacobian of SO(3)
	static Eigen::Matrix3d JacobianR(const Eigen::Vector3d& w)
	{
		Eigen::Matrix3d Jr = Eigen::Matrix3d::Identity();
		double theta = w.norm();
		if (theta < 0.00001)
		{
			return Jr;// = Matrix3d::Identity();
		}
		else
		{
			Eigen::Vector3d k = w.normalized();  // k - unit direction vector of w
			Eigen::Matrix3d K = skew(k);
			Jr = Eigen::Matrix3d::Identity()
				- (1 - cos(theta)) / theta * K
				+ (1 - sin(theta) / theta)*K*K;
		}
		return Jr;
	}
	static Eigen::Matrix3d JacobianRInv(const Eigen::Vector3d& w)
	{
		Eigen::Matrix3d Jrinv = Eigen::Matrix3d::Identity();
		double theta = w.norm();

		// very small angle
		if (theta < 0.00001)
		{
			return Jrinv;
		}
		else
		{
			Eigen::Vector3d k = w.normalized();  // k - unit direction vector of w
			Eigen::Matrix3d K = Sophus::SO3::hat(k);
			Jrinv = Eigen::Matrix3d::Identity()
				+ 0.5*Sophus::SO3::hat(w)
				+ (1.0 - (1.0 + cos(theta))*theta / (2.0*sin(theta))) *K*K;
		}

		return Jrinv;
		/*
	 * in gtsam:
	 *
	 *   double theta2 = omega.dot(omega);
	 *  if (theta2 <= std::numeric_limits<double>::epsilon()) return I_3x3;
	 *  double theta = std::sqrt(theta2);  // rotation angle
	 *  * Right Jacobian for Log map in SO(3) - equation (10.86) and following equations in
	 *   * G.S. Chirikjian, "Stochastic Models, Information Theory, and Lie Groups", Volume 2, 2008.
	 *   * logmap( Rhat * expmap(omega) ) \approx logmap( Rhat ) + Jrinv * omega
	 *   * where Jrinv = LogmapDerivative(omega);
	 *   * This maps a perturbation on the manifold (expmap(omega))
	 *   * to a perturbation in the tangent space (Jrinv * omega)
	 *
	 *  const Matrix3 W = skewSymmetric(omega); // element of Lie algebra so(3): W = omega^
	 *  return I_3x3 + 0.5 * W +
	 *         (1 / (theta * theta) - (1 + cos(theta)) / (2 * theta * sin(theta))) *
	 *             W * W;
	 *
	 * */
	}

	// left jacobian of SO(3), Jl(x) = Jr(-x)
	static Eigen::Matrix3d JacobianL(const Eigen::Vector3d& w)
	{
		return JacobianR(-w);
	}
	// left jacobian inverse
	static Eigen::Matrix3d JacobianLInv(const Eigen::Vector3d& w)
	{
		return JacobianRInv(-w);
	}


	inline Eigen::Quaterniond normalizeRotationQ(const Eigen::Quaterniond& r)
	{
		Eigen::Quaterniond _r(r);
		if (_r.w() < 0)
		{
			_r.coeffs() *= -1;
		}
		return _r.normalized();
	}

	inline Eigen::Matrix3d normalizeRotationM(const Eigen::Matrix3d& R)
	{
		Eigen::Quaterniond qr(R);
		return normalizeRotationQ(qr).toRotationMatrix();
	}
private:
	/*
	 * NOTE:
	 * don't add pointer as member variable.
	 * operator = is used in g2o
	*/

	// delta measurements, position/velocity/rotation(matrix)
	Eigen::Vector3d _delta_P;    // P_k+1 = P_k + V_k*dt + R_k*a_k*dt*dt/2
	Eigen::Vector3d _delta_V;    // V_k+1 = V_k + R_k*a_k*dt
	Eigen::Matrix3d _delta_R;    // R_k+1 = R_k*exp(w_k*dt).     note: Rwc, Rwc'=Rwc*[w_body]x

	// jacobian of delta measurements w.r.t bias of gyro/acc
	Eigen::Matrix3d _J_P_Biasg;     // position / gyro
	Eigen::Matrix3d _J_P_Biasa;     // position / acc
	Eigen::Matrix3d _J_V_Biasg;     // velocity / gyro
	Eigen::Matrix3d _J_V_Biasa;     // velocity / acc
	Eigen::Matrix3d _J_R_Biasg;   // rotation / gyro

	// noise covariance propagation of delta measurements
	Eigen::Matrix<double, 9, 9> _cov_P_V_Phi;

	double deltaTime;
	IMUData imuNoise;
};
}

#endif