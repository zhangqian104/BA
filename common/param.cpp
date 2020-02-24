#include "param.h"

Param::Param()
{
    Eigen::Matrix3d R;   // 把body坐标系朝向旋转一下,得到相机坐标系，好让它看到landmark,  相机坐标系的轴在body坐标系中的表示
    // 相机朝着轨迹里面看， 特征点在轨迹外部， 这里我们采用这个
    R << 0, 0, -1,
            -1, 0, 0,
            0, 1, 0;
	//R.setIdentity();
    Rbc = R;
    tbc = Eigen::Vector3d(0.05,0.04,0.03);

	Eigen::Matrix3d mK;
	mK << fx, 0.0, cx,
				0.0, fy, cy,
				0.0, 0.0, 1.0;
	K = mK;
}

Param::~Param()
{

}