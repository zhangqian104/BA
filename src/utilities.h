#ifndef IMUSIMWITHPOINTLINE_UTILITIES_H
#define IMUSIMWITHPOINTLINE_UTILITIES_H

#include "imu.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <iostream>
#include <vector>
#include <fstream>

// save 3d points to file
void SavePoints(std::string filename, PointsType points);

// save 3d points and it's obs in image
void SaveFeatures(std::string filename, MotionData data, PointsType points, ObsersType features);

void SavaPoseAndFeatures(std::ofstream &saveFile, MotionData data,
	PointsType points, ObsersType features);

// save line obs
void SaveLines(std::string filename,
                std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > features);


void LoadPose(std::string filename, std::vector<MotionData>& pose);

// save imu body data
void savePose(std::string filename, std::vector<MotionData> pose);

Eigen::Matrix3d GetSkewSymmMatrix(Eigen::Vector3d v);

double CalculateError3D(const Eigen::Vector3d pt1, const Eigen::Vector3d pt2);

#endif //IMUSIMWITHPOINTLINE_UTILITIES_H


