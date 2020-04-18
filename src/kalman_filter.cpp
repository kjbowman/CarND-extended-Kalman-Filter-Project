#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using std::cout;
using std::endl;
using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::UpdateQ(double dt) {
  Q_ = MatrixXd(4,4);
  double dt2 = dt * dt;
  double dt3_2 = (dt * dt2)/2.0;
  double dt4_4 = (dt2 * dt2)/4.0;

  Q_ << dt4_4*noise_ax,       0,        dt3_2*noise_ax,       0,
              0,        dt4_4*noise_ay,       0,        dt3_2*noise_ay,
        dt3_2*noise_ax,       0,        dt2*noise_ax,         0,
              0,        dt3_2*noise_ay,       0,        dt2*noise_ay;
}

void KalmanFilter::UpdateF(double dt) {
  F_(0,2) = dt;
  F_(1,3) = dt;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

// for Lidar measurements
void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - (H_ * x_);

  UpdateCommon(y, H_);
}

// for Radar measurements
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  Tools tools;

  VectorXd hx = VectorXd(3);
  hx = tools.Cartesian_to_Polar(x_);

  VectorXd y = z - hx;
  y(1) = atan2(sin(y(1)), cos(y(1))); // normalize angle to [-pi, pi]

  MatrixXd Hj = tools.CalculateJacobian(x_);

  UpdateCommon(y, Hj);
}

// common to both Lidar and Radar, given the appropriate H matrix
void KalmanFilter::UpdateCommon(const VectorXd &y, const MatrixXd &H) {
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R_;
  MatrixXd K = (P_ * Ht) * S.inverse();

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}
