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

  // helper method to convert Cartesian to polor
  // @param v vector(x, y, vx, vy)
  // @return vector(rho, phi, rho_dot)
  Eigen::VectorXd Cartesian_to_Polar(const Eigen::VectorXd& cartesian);

  // helper method to convert polar to Cartesian
  // @param v vector(rho, phi, rho_dot)
  // @return vector(x, y, vx, vy)
  Eigen::VectorXd Polar_to_Cartesian(const Eigen::VectorXd& polar);

  // helper method to calculate RMSE
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                const std::vector<Eigen::VectorXd> &ground_truth);

  // helper method to calculate Jacobians
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

};

#endif  // TOOLS_H_
