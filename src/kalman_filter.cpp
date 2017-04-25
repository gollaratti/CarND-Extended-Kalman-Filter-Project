#include "kalman_filter.h"

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
  P_ = (F_ * P_) * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = (H_ * P_) * Ht + R_;
  MatrixXd K = P_*Ht*S.inverse();
  x_= x_ + K*y;
  unsigned int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  #define PI 3.1415926535897
  VectorXd h_x = VectorXd(3);
  double ro = atan2(x_(1,0),x_(0,0));
  double sqt_px_2__py_2 = sqrt(x_(0,0)*x_(0,0) + x_(1,0)*x_(1,0));
  
  //set-up h_x
  h_x(0,0) = sqt_px_2__py_2;
  h_x(1,0) = ro;
  if(sqt_px_2__py_2 > 0.000001) {
    h_x(2,0) = (x_(0,0)*x_(2,0) + x_(1,0)*x_(3,0))/sqt_px_2__py_2;
    }
  else {
	if((x_(0,0)*x_(2,0) + x_(1,0)*x_(3,0)) > 0.0) {
		h_x(2,0) = LONG_MAX;
	}
	else {
		h_x(2,0) = -LONG_MAX;
	}
    }
         
  VectorXd y = z - h_x;
	  //limit ro to [-PI,+PI]
	  while (y(1,0) > PI) {
		  y(1,0) -= 2*PI;
	  }
	  while (y(1,0) < -PI) {
		  y(1,0) += 2*PI;
	  }
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = (H_ * P_) * Ht + R_;
  MatrixXd K = P_*Ht*S.inverse();
  x_= x_ + K*y;
  unsigned int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_)*P_;
}
