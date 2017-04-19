#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
    std::cout << "calculateRMSE() - Error - Invalid Estimations vector size" << std::endl;
    return rmse;
  }

  //accumulate squared residuals
  VectorXd sum(4);
  sum << 0, 0, 0, 0;

  for(int i = 0; i < estimations.size(); ++i) {
    VectorXd estimated = estimations[i];
    VectorXd actual = ground_truth[i];

    VectorXd diff = (estimated - actual);
    //coefficient-wise multiplication
    diff = diff.array() * diff.array();
    sum = sum + diff;
  }

  //calculate the mean
  Eigen::VectorXd mean = sum / estimations.size();

  //calculate the squared root
  rmse = mean.array().sqrt();

  return rmse;
}

float Tools::NormalizeAngle(float angle_rad) {
  angle_rad = angle_rad - 2*M_PI*floor((angle_rad+ M_PI)/(2 * M_PI ));
  return angle_rad;
}

double Tools::CalculateNIS(const VectorXd& z_predicted, const VectorXd& z, const MatrixXd& S) {
  VectorXd z_diff = z - z_predicted;
  double NIS = z_diff.transpose() * S.inverse() * z_diff;

  return NIS;
}

//float Tools::NormalizeAngle(float angle_rad) {
//  //angle normalization
//    while (angle_rad > M_PI)
//      angle_rad -= 2.*M_PI;
//
//    while (angle_rad < -M_PI)
//      angle_rad += 2. * M_PI;
//
//    return angle_rad;
//}
