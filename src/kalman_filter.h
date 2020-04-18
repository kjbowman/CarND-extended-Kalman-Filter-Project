#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"
#include "tools.h"

class KalmanFilter {
public:
  // constructor
  KalmanFilter();

  // destructor
  virtual ~KalmanFilter();

  // update the process noise covariance matrix
  void UpdateQ(double dt);

  // update the transition matrix
  void UpdateF(double dt);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

  /**
   * Part of the update function common to both the standard standard
   * and the extended filters
   * @param y the measurement update
   * @param H the measurement matrix
   */
  void UpdateCommon(const Eigen::VectorXd &y, const Eigen::MatrixXd &H);

  // public attributes ...
  Eigen::VectorXd x_;   // state vector
  Eigen::MatrixXd P_;   // state covariance matrix
  Eigen::MatrixXd F_;   // state transition matrix
  Eigen::MatrixXd Q_;   // process covariance matrix
  Eigen::MatrixXd H_;   // measurement matrix (lidar)
  // Eigen::MatrixXd Hj_;  // measurement matrix (radar)
  Eigen::MatrixXd R_;   // measurement covariance matrix

private:
  float noise_ax {9.0};
  float noise_ay {9.0};

};

#endif // KALMAN_FILTER_H_
