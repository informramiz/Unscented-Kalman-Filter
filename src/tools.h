#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to normalize angle between [-pi, pi]
   */
  static float NormalizeAngle(float angle_rad);

  /**
   * A helper method to calculate NIS value
   * which can be used to verify consistency of filter
   */
  static double CalculateNIS(const Eigen::VectorXd& z_predicted, const Eigen::VectorXd& z, const Eigen::MatrixXd& S);
};

#endif /* TOOLS_H_ */
