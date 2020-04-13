#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  // constructor
  Tools();

  // destructor
  virtual ~Tools();

  // Cartesian <-> Polar conversion helpers
  Eigen::VectorXd Cartesian_to_Polar(const Eigen::VectorXd& v);
  Eigen::VectorXd Polar_to_Cartesian(const Eigen::VectorXd& v);

  // helper method to calculate RMSE
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                const std::vector<Eigen::VectorXd> &ground_truth);

  // helper method to calculate Jacobians
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

};

#endif  // TOOLS_H_
