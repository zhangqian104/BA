#include "utility.h"

Eigen::Matrix3d Utility::g2R(const Eigen::Vector3d &g)
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d ng1 = g.normalized();
    Eigen::Vector3d ng2{0, 0, 1.0};
    R0 = Eigen::Quaterniond::FromTwoVectors(ng1, ng2).toRotationMatrix();
    double yaw = Utility::R2ypr(R0).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    // R0 = Utility::ypr2R(Eigen::Vector3d{-90, 0, 0}) * R0;
    return R0;
}

void SavePoints(std::string filename, PointsType points)
{
	std::ofstream savePoints;
	savePoints.open(filename.c_str());

	for (int i = 0; i < points.size(); ++i) {
		IdType id = points[i].first;
		Eigen::Vector4d p = points[i].second;

		savePoints << id << " "
			<< p(0) << " "
			<< p(1) << " "
			<< p(2) << " "
			<< p(3) << std::endl;
	}
}

void SavaPoseAndFeatures(std::ofstream &saveFile,
	MotionData data, PointsType points, ObsersType features)
{
	double time = data.timestamp;
	Eigen::Quaterniond q(data.Rwb);
	Eigen::Vector3d t = data.twb;

	saveFile << -1 << " "
		<< time << " "
		<< q.w() << " "
		<< q.x() << " "
		<< q.y() << " "
		<< q.z() << " "
		<< t(0) << " "
		<< t(1) << " "
		<< t(2) << " "
		<< std::endl;

	for (int i = 0; i < points.size(); ++i) {
		IdType id = points[i].first;
		Eigen::Vector4d p = points[i].second;
		Eigen::Vector2d f = features[i].second;
		saveFile << id << " "
			<< p(0) << " "
			<< p(1) << " "
			<< p(2) << " "
			<< p(3) << " "
			<< f(0) << " "
			<< f(1) << " "
			<< std::endl;
	}
}

void SaveFeatures(std::string filename, MotionData data, PointsType points, ObsersType features)
{
	std::ofstream savePoints;
	savePoints.open(filename.c_str());
	double time = data.timestamp;
	Eigen::Quaterniond q(data.Rwb);
	Eigen::Vector3d t = data.twb;

	savePoints << -1 << " "
		<< time << " "
		<< q.w() << " "
		<< q.x() << " "
		<< q.y() << " "
		<< q.z() << " "
		<< t(0) << " "
		<< t(1) << " "
		<< t(2) << " "
		<< std::endl;

	for (int i = 0; i < points.size(); ++i) {
		IdType id = points[i].first;
		Eigen::Vector4d p = points[i].second;
		Eigen::Vector2d f = features[i].second;
		savePoints << id << " "
			<< p(0) << " "
			<< p(1) << " "
			<< p(2) << " "
			<< p(3) << " "
			<< f(0) << " "
			<< f(1) << " "
			<< std::endl;
	}
}
void SaveLines(std::string filename,
	std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > features)
{
	std::ofstream savePoints;
	savePoints.open(filename.c_str());

	for (int i = 0; i < features.size(); ++i) {
		Eigen::Vector4d f = features[i];
		savePoints << f(0) << " "
			<< f(1) << " "
			<< f(2) << " "
			<< f(3) << " "
			<< std::endl;
	}
}

void LoadPose(std::string filename, std::vector<MotionData>& pose)
{
	std::ifstream f;
	f.open(filename.c_str());

	if (!f.is_open()) {
		std::cerr << " can't open LoadFeatures file " << std::endl;
		return;
	}

	while (!f.eof()) {
		std::string s;
		std::getline(f, s);

		if (!s.empty()) {
			std::stringstream ss;
			ss << s;

			MotionData data;
			double time;
			Eigen::Quaterniond q;
			Eigen::Vector3d t;
			Eigen::Vector3d gyro;
			Eigen::Vector3d acc;

			ss >> time;
			ss >> q.w();
			ss >> q.x();
			ss >> q.y();
			ss >> q.z();
			ss >> t(0);
			ss >> t(1);
			ss >> t(2);
			ss >> gyro(0);
			ss >> gyro(1);
			ss >> gyro(2);
			ss >> acc(0);
			ss >> acc(1);
			ss >> acc(2);

			data.timestamp = time;
			data.imuGyro = gyro;
			data.imuAcc = acc;
			data.twb = t;
			data.Rwb = Eigen::Matrix3d(q);
			pose.push_back(data);
		}
	}
}

void savePose(std::string filename, std::vector<MotionData> pose)
{
	std::ofstream f;
	f.open(filename.c_str());

	for (int i = 0; i < pose.size(); ++i) {
		MotionData data = pose[i];
		double time = data.timestamp;
		Eigen::Quaterniond q(data.Rwb);
		Eigen::Vector3d t = data.twb;
		Eigen::Vector3d gyro = data.imuGyro;
		Eigen::Vector3d acc = data.imuAcc;

		f << time << " " << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << " "
			<< t(0) << " " << t(1) << " " << t(2) << " "
			<< gyro(0) << " " << gyro(1) << " " << gyro(2) << " "
			<< acc(0) << " " << acc(1) << " " << acc(2) << " " << std::endl;
		if (data.rsPoses.size() > 0) {
			for (unsigned int row = 0; row < data.rsPoses.size(); row++) {
				RowPose iPose = data.rsPoses[row];
				f << iPose.row << " " << iPose.timestamp << " " << iPose.qwb.w() << " " << iPose.qwb.x() << " " << iPose.qwb.y() << " " << iPose.qwb.z() << " "
					<< iPose.twb(0) << " " << iPose.twb(1) << " " << iPose.twb(2) << " " << std::endl;
			}
		}
	}
}

Eigen::Matrix3d GetSkewSymmMatrix(Eigen::Vector3d v)
{
	Eigen::Matrix3d mat;
	mat << 0.0, -v[2], v[1],
		v[2], 0.0, -v[0],
		-v[1], v[0], 0.0;
	return mat;
}

double CalculateError3D(const Eigen::Vector3d pt1, const Eigen::Vector3d pt2)
{
	return (pt1 - pt2).norm();
}

// euler2Rotation:   body frame to interitail frame
Eigen::Matrix3d euler2Rotation(Eigen::Vector3d  eulerAngles)
{
	double roll = eulerAngles(0);
	double pitch = eulerAngles(1);
	double yaw = eulerAngles(2);

	double cr = cos(roll); double sr = sin(roll);
	double cp = cos(pitch); double sp = sin(pitch);
	double cy = cos(yaw); double sy = sin(yaw);

	Eigen::Matrix3d RIb;
	RIb << cy * cp, cy*sp*sr - sy * cr, sy*sr + cy * cr*sp,
		sy*cp, cy *cr + sy * sr*sp, sp*sy*cr - cy * sr,
		-sp, cp*sr, cp*cr;
	return RIb;
}

Eigen::Matrix3d eulerRates2bodyRates(Eigen::Vector3d eulerAngles)
{
	double roll = eulerAngles(0);
	double pitch = eulerAngles(1);

	double cr = cos(roll); double sr = sin(roll);
	double cp = cos(pitch); double sp = sin(pitch);

	Eigen::Matrix3d R;
	R << 1, 0, -sp,
		0, cr, sr*cp,
		0, -sr, cr*cp;

	return R;
}
