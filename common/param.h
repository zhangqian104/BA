#ifndef IMUSIM_PARAM_H
#define IMUSIM_PARAM_H

#include <eigen3/Eigen/Core>

class Param{
public:
    Param();//里面初始化外参
	~Param();
    // time
    int imuFrequency = 200;
    int camFrequency = 30;
    double imuTimestep = 1./imuFrequency;
    double camTimestep = 1./camFrequency;
    double tStart = 0.;
    double tEnd = 10.;  //  20 s
	double rsT = 0.01;

    // noise
    double gyroBiasSigma = 1.0e-5;
    double accBiasSigma = 0.0001;

    double gyroNoiseSigma = 0.015;    // rad/s * 1/sqrt(hz)
    double accNoiseSigma = 0.019;      //　m/(s^2) * 1/sqrt(hz)

	const bool addPixedNoise = false;
    double pixelNoiseSigma = 3.;              // 1 pixel noise

	double camPoseNoise = 0.05;

    // cam f
    double fx = 460.;
    double fy = 460.;
    double cx = 255.;
    double cy = 255.;
    double image_w = 640.;
    double image_h = 480.;
	Eigen::Matrix3d K;

    // 外参数
    Eigen::Matrix3d Rbc;   // cam to body
    Eigen::Vector3d tbc;     // cam to body
};


#endif //IMUSIM_PARAM_H
