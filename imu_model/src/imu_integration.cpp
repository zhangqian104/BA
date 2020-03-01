#include "imu_integration.h"

namespace SLAM_SIMULATION {
IMUPreintegRecursive::IMUPreintegRecursive(IMUData imuNoise_, const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
	const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg)
{
	imuNoise = imuNoise_;
	acc_0 = _acc_0;
	gyr_0 = _gyr_0;
	linearized_acc = _acc_0;
	linearized_gyr = _gyr_0;
	linearized_ba = _linearized_ba;
	linearized_bg = _linearized_bg;
	jacobian = Eigen::Matrix<double, 15, 15>::Identity();
	covariance = Eigen::Matrix<double, 15, 15>::Zero();
	sum_dt = 0.0;
	delta_p = Eigen::Vector3d::Zero();
	delta_q = Eigen::Quaterniond::Identity();
	delta_v = Eigen::Vector3d::Zero();

	noise = Eigen::Matrix<double, 18, 18>::Zero();
	noise.block<3, 3>(0, 0) = imuNoise.accNoiseCov;//(ACC_N * ACC_N) * Eigen::Matrix3d::Identity();/
	noise.block<3, 3>(3, 3) = imuNoise.gyroNoiseCov;//(GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(6, 6) = imuNoise.accNoiseCov;//(ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(9, 9) = imuNoise.gyroNoiseCov;//(GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(12, 12) = imuNoise.accBiasRWCov;//(ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(15, 15) = imuNoise.gyrBiasRWCov;//(GYR_W * GYR_W) * Eigen::Matrix3d::Identity();
}

void IMUPreintegRecursive::AddNewImu(double dt, const Eigen::Vector3d &acc, const Eigen::Vector3d &gyr)
{
	dt_buf.push_back(dt);
	acc_buf.push_back(acc);
	gyr_buf.push_back(gyr);
	Propagate(dt, acc, gyr);
}

void IMUPreintegRecursive::Repropagate(const Eigen::Vector3d &_linearized_ba, const Eigen::Vector3d &_linearized_bg)
{
sum_dt = 0.0;
acc_0 = linearized_acc;
gyr_0 = linearized_gyr;
delta_p.setZero();
delta_q.setIdentity();
delta_v.setZero();
linearized_ba = _linearized_ba;
linearized_bg = _linearized_bg;
jacobian.setIdentity();
covariance.setZero();
for (int i = 0; i < static_cast<int>(dt_buf.size()); i++)
	Propagate(dt_buf[i], acc_buf[i], gyr_buf[i]);
}

void IMUPreintegRecursive::MidpointIntegration(double _dt,
const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg,
Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg, bool update_jacobian)
{
//ROS_INFO("midpoint integration");
	Eigen::Vector3d un_acc_0 = delta_q * (_acc_0 - linearized_ba);
	Eigen::Vector3d un_gyr = 0.5 * (_gyr_0 + _gyr_1) - linearized_bg;
result_delta_q = delta_q * Eigen::Quaterniond(1, un_gyr(0) * _dt / 2, un_gyr(1) * _dt / 2, un_gyr(2) * _dt / 2);
Eigen::Vector3d un_acc_1 = result_delta_q * (_acc_1 - linearized_ba);
Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
result_delta_p = delta_p + delta_v * _dt + 0.5 * un_acc * _dt * _dt;
result_delta_v = delta_v + un_acc * _dt;
result_linearized_ba = linearized_ba;
result_linearized_bg = linearized_bg;

if (update_jacobian)
{
	Eigen::Vector3d w_x = 0.5 * (_gyr_0 + _gyr_1) - linearized_bg;
	Eigen::Vector3d a_0_x = _acc_0 - linearized_ba;
	Eigen::Vector3d a_1_x = _acc_1 - linearized_ba;
	Eigen::Matrix3d R_w_x, R_a_0_x, R_a_1_x;

	R_w_x << 0, -w_x(2), w_x(1),
		w_x(2), 0, -w_x(0),
		-w_x(1), w_x(0), 0;
	R_a_0_x << 0, -a_0_x(2), a_0_x(1),
		a_0_x(2), 0, -a_0_x(0),
		-a_0_x(1), a_0_x(0), 0;
	R_a_1_x << 0, -a_1_x(2), a_1_x(1),
		a_1_x(2), 0, -a_1_x(0),
		-a_1_x(1), a_1_x(0), 0;

	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(15, 15);
	F.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
	F.block<3, 3>(0, 3) = -0.25 * delta_q.toRotationMatrix() * R_a_0_x * _dt * _dt +
		-0.25 * result_delta_q.toRotationMatrix() * R_a_1_x * (Eigen::Matrix3d::Identity() - R_w_x * _dt) * _dt * _dt;
	F.block<3, 3>(0, 6) = Eigen::MatrixXd::Identity(3, 3) * _dt;
	F.block<3, 3>(0, 9) = -0.25 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * _dt * _dt;
	F.block<3, 3>(0, 12) = -0.25 * result_delta_q.toRotationMatrix() * R_a_1_x * _dt * _dt * -_dt;
	F.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity() - R_w_x * _dt;
	F.block<3, 3>(3, 12) = -1.0 * Eigen::MatrixXd::Identity(3, 3) * _dt;
	F.block<3, 3>(6, 3) = -0.5 * delta_q.toRotationMatrix() * R_a_0_x * _dt +
		-0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * (Eigen::Matrix3d::Identity() - R_w_x * _dt) * _dt;
	F.block<3, 3>(6, 6) = Eigen::Matrix3d::Identity();
	F.block<3, 3>(6, 9) = -0.5 * (delta_q.toRotationMatrix() + result_delta_q.toRotationMatrix()) * _dt;
	F.block<3, 3>(6, 12) = -0.5 * result_delta_q.toRotationMatrix() * R_a_1_x * _dt * -_dt;
	F.block<3, 3>(9, 9) = Eigen::Matrix3d::Identity();
	F.block<3, 3>(12, 12) = Eigen::Matrix3d::Identity();
	//cout<<"A"<<endl<<A<<endl;

	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(15, 18);
	V.block<3, 3>(0, 0) = 0.25 * delta_q.toRotationMatrix() * _dt * _dt;
	V.block<3, 3>(0, 3) = 0.25 * -result_delta_q.toRotationMatrix() * R_a_1_x  * _dt * _dt * 0.5 * _dt;
	V.block<3, 3>(0, 6) = 0.25 * result_delta_q.toRotationMatrix() * _dt * _dt;
	V.block<3, 3>(0, 9) = V.block<3, 3>(0, 3);
	V.block<3, 3>(3, 3) = 0.5 * Eigen::MatrixXd::Identity(3, 3) * _dt;
	V.block<3, 3>(3, 9) = 0.5 * Eigen::MatrixXd::Identity(3, 3) * _dt;
	V.block<3, 3>(6, 0) = 0.5 * delta_q.toRotationMatrix() * _dt;
	V.block<3, 3>(6, 3) = 0.5 * -result_delta_q.toRotationMatrix() * R_a_1_x  * _dt * 0.5 * _dt;
	V.block<3, 3>(6, 6) = 0.5 * result_delta_q.toRotationMatrix() * _dt;
	V.block<3, 3>(6, 9) = V.block<3, 3>(6, 3);
	V.block<3, 3>(9, 12) = Eigen::MatrixXd::Identity(3, 3) * _dt;
	V.block<3, 3>(12, 15) = Eigen::MatrixXd::Identity(3, 3) * _dt;

	//step_jacobian = F;
	//step_V = V;
	jacobian = F * jacobian;
	covariance = F * covariance * F.transpose() + V * noise * V.transpose();
}

}

void IMUPreintegRecursive::EulerIntegration(double _dt, const Eigen::Vector3d &_acc_0, const Eigen::Vector3d &_gyr_0,
const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1,
const Eigen::Vector3d &delta_p, const Eigen::Quaterniond &delta_q, const Eigen::Vector3d &delta_v,
const Eigen::Vector3d &linearized_ba, const Eigen::Vector3d &linearized_bg,
Eigen::Vector3d &result_delta_p, Eigen::Quaterniond &result_delta_q, Eigen::Vector3d &result_delta_v,
Eigen::Vector3d &result_linearized_ba, Eigen::Vector3d &result_linearized_bg, bool update_jacobian)
{
result_delta_p = delta_p + delta_v * _dt + 0.5 * (delta_q * (_acc_1 - linearized_ba)) * _dt * _dt;
result_delta_v = delta_v + delta_q * (_acc_1 - linearized_ba) * _dt;
Eigen::Vector3d omg = _gyr_1 - linearized_bg;
omg = omg * _dt / 2;
Eigen::Quaterniond dR(1, omg(0), omg(1), omg(2));
result_delta_q = (delta_q * dR);
result_linearized_ba = linearized_ba;
result_linearized_bg = linearized_bg;

if (update_jacobian)
{
	Eigen::Vector3d w_x = _gyr_1 - linearized_bg;
	Eigen::Vector3d a_x = _acc_1 - linearized_ba;
	Eigen::Matrix3d R_w_x, R_a_x;

	R_w_x << 0, -w_x(2), w_x(1),
		w_x(2), 0, -w_x(0),
		-w_x(1), w_x(0), 0;
	R_a_x << 0, -a_x(2), a_x(1),
		a_x(2), 0, -a_x(0),
		-a_x(1), a_x(0), 0;

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(15, 15);
	// one step euler 0.5
	A.block<3, 3>(0, 3) = 0.5 * (-1 * delta_q.toRotationMatrix()) * R_a_x * _dt;
	A.block<3, 3>(0, 6) = Eigen::MatrixXd::Identity(3, 3);
	A.block<3, 3>(0, 9) = 0.5 * (-1 * delta_q.toRotationMatrix()) * _dt;
	A.block<3, 3>(3, 3) = -R_w_x;
	A.block<3, 3>(3, 12) = -1 * Eigen::MatrixXd::Identity(3, 3);
	A.block<3, 3>(6, 3) = (-1 * delta_q.toRotationMatrix()) * R_a_x;
	A.block<3, 3>(6, 9) = (-1 * delta_q.toRotationMatrix());
	//cout<<"A"<<endl<<A<<endl;

	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(15, 12);
	U.block<3, 3>(0, 0) = 0.5 * delta_q.toRotationMatrix() * _dt;
	U.block<3, 3>(3, 3) = Eigen::MatrixXd::Identity(3, 3);
	U.block<3, 3>(6, 0) = delta_q.toRotationMatrix();
	U.block<3, 3>(9, 6) = Eigen::MatrixXd::Identity(3, 3);
	U.block<3, 3>(12, 9) = Eigen::MatrixXd::Identity(3, 3);

	// put outside
	Eigen::Matrix<double, 12, 12> noise = Eigen::Matrix<double, 12, 12>::Zero();
	noise.block<3, 3>(0, 0) = imuNoise.accNoiseCov;//(ACC_N * ACC_N) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(3, 3) = imuNoise.gyroNoiseCov;//(GYR_N * GYR_N) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(6, 6) = imuNoise.accBiasRWCov;//(ACC_W * ACC_W) * Eigen::Matrix3d::Identity();
	noise.block<3, 3>(9, 9) = imuNoise.gyrBiasRWCov;//(GYR_W * GYR_W) * Eigen::Matrix3d::Identity();

	//write F directly
	Eigen::MatrixXd F, V;
	F = (Eigen::MatrixXd::Identity(15, 15) + _dt * A);
	V = _dt * U;
	step_jacobian = F;
	step_V = V;
	jacobian = F * jacobian;
	covariance = F * covariance * F.transpose() + V * noise * V.transpose();
}

}

void IMUPreintegRecursive::Propagate(double _dt, const Eigen::Vector3d &_acc_1, const Eigen::Vector3d &_gyr_1)
{
	dt = _dt;
	acc_1 = _acc_1;
	gyr_1 = _gyr_1;
	Eigen::Vector3d result_delta_p;
	Eigen::Quaterniond result_delta_q;
	Eigen::Vector3d result_delta_v;
	Eigen::Vector3d result_linearized_ba;
	Eigen::Vector3d result_linearized_bg;

	MidpointIntegration(_dt, acc_0, gyr_0, _acc_1, _gyr_1, delta_p, delta_q, delta_v,
		linearized_ba, linearized_bg,
		result_delta_p, result_delta_q, result_delta_v,
		result_linearized_ba, result_linearized_bg, 1);

	//checkJacobian(_dt, acc_0, gyr_0, acc_1, gyr_1, delta_p, delta_q, delta_v,
	//                    linearized_ba, linearized_bg);
	delta_p = result_delta_p;
	delta_q = result_delta_q;
	delta_v = result_delta_v;
	linearized_ba = result_linearized_ba;
	linearized_bg = result_linearized_bg;
	delta_q.normalize();
	sum_dt += dt;
	acc_0 = acc_1;
	gyr_0 = gyr_1;
}

IMUPreintegBatch::IMUPreintegBatch(const IMUPreintegBatch& pre) :
	_delta_P(pre._delta_P),
	_delta_V(pre._delta_V),
	_delta_R(pre._delta_R),
	_J_P_Biasg(pre._J_P_Biasg),
	_J_P_Biasa(pre._J_P_Biasa),
	_J_V_Biasg(pre._J_V_Biasg),
	_J_V_Biasa(pre._J_V_Biasa),
	_J_R_Biasg(pre._J_R_Biasg),
	_cov_P_V_Phi(pre._cov_P_V_Phi),
	deltaTime(pre.deltaTime),
	imuNoise(pre.imuNoise)
{

}

IMUPreintegBatch::IMUPreintegBatch(IMUData imuNoise_)
{
	imuNoise = imuNoise_;
	// delta measurements, position/velocity/rotation(matrix)
	_delta_P.setZero();    // P_k+1 = P_k + V_k*dt + R_k*a_k*dt*dt/2
	_delta_V.setZero();    // V_k+1 = V_k + R_k*a_k*dt
	_delta_R.setIdentity();    // R_k+1 = R_k*exp(w_k*dt).     note: Rwc, Rwc'=Rwc*[w_body]x

	// jacobian of delta measurements w.r.t bias of gyro/acc
	_J_P_Biasg.setZero();     // position / gyro
	_J_P_Biasa.setZero();     // position / acc
	_J_V_Biasg.setZero();     // velocity / gyro
	_J_V_Biasa.setZero();     // velocity / acc
	_J_R_Biasg.setZero();   // rotation / gyro

	// noise covariance propagation of delta measurements
	_cov_P_V_Phi.setZero();

	deltaTime = 0;
}

void IMUPreintegBatch::reset()
{
	// delta measurements, position/velocity/rotation(matrix)
	_delta_P.setZero();    // P_k+1 = P_k + V_k*dt + R_k*a_k*dt*dt/2
	_delta_V.setZero();    // V_k+1 = V_k + R_k*a_k*dt
	_delta_R.setIdentity();    // R_k+1 = R_k*exp(w_k*dt).     note: Rwc, Rwc'=Rwc*[w_body]x

	// jacobian of delta measurements w.r.t bias of gyro/acc
	_J_P_Biasg.setZero();     // position / gyro
	_J_P_Biasa.setZero();     // position / acc
	_J_V_Biasg.setZero();     // velocity / gyro
	_J_V_Biasa.setZero();     // velocity / acc
	_J_R_Biasg.setZero();   // rotation / gyro

	// noise covariance propagation of delta measurements
	_cov_P_V_Phi.setZero();

	deltaTime = 0;

}

// incrementally update 1)delta measurements, 2)jacobians, 3)covariance matrix
// acc: acc_measurement - bias_a, last measurement!! not current measurement
// omega: gyro_measurement - bias_g, last measurement!! not current measurement
void IMUPreintegBatch::update(const Eigen::Vector3d& omega, const Eigen::Vector3d& acc, const double& dt)
{
	double dt2 = dt * dt;

	Eigen::Matrix3d dR = Expmap(omega*dt);
	Eigen::Matrix3d Jr = JacobianR(omega*dt);

	// noise covariance propagation of delta measurements
	// err_k+1 = A*err_k + B*err_gyro + C*err_acc
	Eigen::Matrix3d I3x3 = Eigen::Matrix3d::Identity();
	Eigen::Matrix<double, 9, 9> A = Eigen::Matrix<double, 9, 9>::Identity();
	A.block<3, 3>(6, 6) = dR.transpose();
	A.block<3, 3>(3, 6) = -_delta_R * skew(acc)*dt;
	A.block<3, 3>(0, 6) = -0.5*_delta_R*skew(acc)*dt2;
	A.block<3, 3>(0, 3) = I3x3 * dt;
	Eigen::Matrix<double, 9, 3> Bg = Eigen::Matrix<double, 9, 3>::Zero();
	Bg.block<3, 3>(6, 0) = Jr * dt;
	Eigen::Matrix<double, 9, 3> Ca = Eigen::Matrix<double, 9, 3>::Zero();
	Ca.block<3, 3>(3, 0) = _delta_R * dt;
	Ca.block<3, 3>(0, 0) = 0.5*_delta_R*dt2;
	_cov_P_V_Phi = A * _cov_P_V_Phi*A.transpose() +
		Bg * imuNoise.gyroNoiseCov*Bg.transpose() +//IMUData::getGyrMeasCov()
		Ca * imuNoise.accNoiseCov*Ca.transpose();//IMUData::getAccMeasCov()


	// jacobian of delta measurements w.r.t bias of gyro/acc
	// update P first, then V, then R
	_J_P_Biasa += _J_V_Biasa * dt - 0.5*_delta_R*dt2;
	_J_P_Biasg += _J_V_Biasg * dt - 0.5*_delta_R*skew(acc)*_J_R_Biasg*dt2;
	_J_V_Biasa += -_delta_R * dt;
	_J_V_Biasg += -_delta_R * skew(acc)*_J_R_Biasg*dt;
	_J_R_Biasg = dR.transpose()*_J_R_Biasg - Jr * dt;

	// delta measurements, position/velocity/rotation(matrix)
	// update P first, then V, then R. because P's update need V&R's previous state
	_delta_P += _delta_V * dt + 0.5*_delta_R*acc*dt2;    // P_k+1 = P_k + V_k*dt + R_k*a_k*dt*dt/2
	_delta_V += _delta_R * acc*dt;
	_delta_R = normalizeRotationM(_delta_R*dR);  // normalize rotation, in case of numerical error accumulation


//    // noise covariance propagation of delta measurements
//    // err_k+1 = A*err_k + B*err_gyro + C*err_acc
//    Matrix3d I3x3 = Matrix3d::Identity();
//    MatrixXd A = MatrixXd::Identity(9,9);
//    A.block<3,3>(6,6) = dR.transpose();
//    A.block<3,3>(3,6) = -_delta_R*skew(acc)*dt;
//    A.block<3,3>(0,6) = -0.5*_delta_R*skew(acc)*dt2;
//    A.block<3,3>(0,3) = I3x3*dt;
//    MatrixXd Bg = MatrixXd::Zero(9,3);
//    Bg.block<3,3>(6,0) = Jr*dt;
//    MatrixXd Ca = MatrixXd::Zero(9,3);
//    Ca.block<3,3>(3,0) = _delta_R*dt;
//    Ca.block<3,3>(0,0) = 0.5*_delta_R*dt2;
//    _cov_P_V_Phi = A*_cov_P_V_Phi*A.transpose() +
//        Bg*IMUData::getGyrMeasCov*Bg.transpose() +
//        Ca*IMUData::getAccMeasCov()*Ca.transpose();

	// delta time
	deltaTime += dt;

}

}