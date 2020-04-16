#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

// for Lidar measurements
void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

// for Radar measurements
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  Tools tools;

  VectorXd Hx = tools.Cartesian_to_Polar(x_);
  VectorXd y = z - Hx;
  // normalize angle in y vector (limit to +/- pi)
  double phi = y(1);
  y(1) = atan2(sin(phi), cos(phi));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht * R_;
  MatrixXd K = P_ * Ht * S.inverse();

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
