#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */

	VectorXd y = z - H_ * x_;
	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd h = VectorXd(3);
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
	float rho = sqrt(pow(px,2)+pow(py,2));
  float theta = atan2(py,px);
  float rho_dot = (px * vx + py * vy)/rho;
  h << rho, theta, rho_dot;
	VectorXd y = z - h;

  float Pi = 3.14159265358979323846;
  if( y(1) > Pi )
      y(1) -= 2*Pi;
  if( y(1) < -Pi )
      y(1) += 2*Pi;

	MatrixXd S = H_ * P_ * H_.transpose() + R_;
	MatrixXd K = P_ * H_.transpose() * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
