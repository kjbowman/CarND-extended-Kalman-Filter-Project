#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::Cartesian_to_Polar(const VectorXd& cartesian) {
  const double epsilon = 0.0001;
  VectorXd polar(3);

  double px = cartesian[0];
  double py = cartesian[1];
  double vx = cartesian[2];
  double vy = cartesian[3];

  double rho = sqrt(px*px + py*py);
  double phi = atan2(py, px);
  double rho_dot = (rho > epsilon) ? (px*vx + py*vy) / rho : 0.0;

  polar << rho, phi, rho_dot;
  return polar;
}

VectorXd Tools::Polar_to_Cartesian(const VectorXd& polar) {
  VectorXd cartesian(4);

  double rho = polar[0];
  double phi = polar[1];
  double rho_dot = polar[2];

  double cos_phi = cos(phi);
  double sin_phi = sin(phi);

  double px = rho * cos_phi;
  double py = rho * sin_phi;
  double vx = rho_dot * cos_phi;
  double vy = rho_dot * sin_phi;

  cartesian << px, py, vx, vy;
  return cartesian;
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if((estimations.size() == 0) || (estimations.size() != ground_truth.size())) {
      cout << "vector size error" << endl;
      return rmse;
  }

  // accumulate squared residuals
  VectorXd accum(4);
  accum << 0,0,0,0;
  VectorXd residual(4);
  residual << 0,0,0,0;

  for (int i=0; i < estimations.size(); ++i) {
    residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();   // square it
    accum += residual;
  }

  // calculate the mean
  VectorXd mean(4);
  mean = accum.array() / estimations.size();
  // calculate the squared root
  rmse = mean.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // compute the Jacobian matrix
  float px2_py2 = px*px + py*py;
  float sqrt_px2_py2 = pow(px2_py2, 0.5);
  float root3_2_px2_py2 = pow(px2_py2, 1.5);

  // check division by zero
  if(fabs(px2_py2) < 0.0001) {
      cout << "CalculateJacobian() - Error - Division by zero" << endl;
      return Hj;
  }

  Hj(0,0) = px / sqrt_px2_py2;
  Hj(0,1) = py / sqrt_px2_py2;

  Hj(1,0) = -py / px2_py2;
  Hj(1,1) = px / px2_py2;

  Hj(2,0) = py * (vx*py - vy*px) / root3_2_px2_py2;
  Hj(2,1) = px * (vy*px - vx*py) / root3_2_px2_py2;
  Hj(2,2) = px / sqrt_px2_py2;
  Hj(2,3) = py / sqrt_px2_py2;

  return Hj;
}
