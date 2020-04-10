#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

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
  if(fabs(px2_px2) < 0.0001) {
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
