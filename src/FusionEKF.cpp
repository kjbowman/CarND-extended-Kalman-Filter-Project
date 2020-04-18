#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>   // cos

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

// Constructor
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process (Q) and measurement noises (R)
   */
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ <<  1,    0,    0,    0,
              0,    1,    0,    0,
              0,    0, 1000,    0,
              0,    0,    0, 1000;

  ekf_.F_ = MatrixXd::Identity(4, 4);
  ekf_.UpdateF(1);

  ekf_.Q_ = MatrixXd::Zero(4, 4);

  // ekf H will always be for laser - radar will compute and use Hj instead of H
  ekf_.H_ = H_laser_;
}

// Destructor
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  Tools tools;
  // ----- INITIALIZATION -----
  if (!is_initialized_) {
    previous_timestamp_ = measurement_pack.timestamp_;

    // first measurement - initialize state vector
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
      ekf_.x_ = tools.Polar_to_Cartesian(measurement_pack.raw_measurements_);
      // don't use vx, vy from initial radar measurement
      ekf_.x_[2] = 0;
      ekf_.x_[3] = 0;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
    }

    // done initializing
    is_initialized_ = true;
    return;
  }

  // ----- PREDICTION -----
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) * 1.0e-6;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.UpdateF(dt);
  ekf_.UpdateQ(dt);

  ekf_.Predict();

  // ----- UPDATE -----
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {  // MeasurementPackage::LIDAR
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
