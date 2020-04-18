#include "tools.h"
#include <iostream>
#include <algorithm>  // std::max

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;
using std::max;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::Cartesian_to_Polar(const VectorXd& cartesian) {
  VectorXd polar(3);

  double px = cartesian(0);
  double py = cartesian(1);
  double vx = cartesian(2);
  double vy = cartesian(3);

  double rho = sqrt(px*px + py*py);
  double phi = atan2(py, px);
  double rho_dot = (rho > epsilon_) ? (px*vx + py*vy) / rho : 0.0;

  polar << rho, phi, rho_dot;
  return polar;
}

VectorXd Tools::Polar_to_Cartesian(const VectorXd& polar) {
  VectorXd cartesian(4);

  double rho = polar(0);
  double phi = polar(1);
  double rho_dot = polar(2);

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
  VectorXd residual(4);
  residual << 0,0,0,0;

  // check for valid input vectors
  if((estimations.size() == 0) || (estimations.size() != ground_truth.size())) {
      cout << "vector size error" << endl;
      return rmse;
  }

  // RMS Error calculation
  for (int i=0; i < estimations.size(); ++i) {
    residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();   // square it
    rmse += residual;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj = MatrixXd::Zero(3, 4);
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // compute the Jacobian matrix
  double px2_py2 = max(px*px + py*py, epsilon_);
  double sqrt_px2_py2 = max(sqrt(px2_py2), epsilon_);
  double root3_2_px2_py2 = max(px2_py2 * sqrt_px2_py2, epsilon_);

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
