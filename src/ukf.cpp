/*
 * ukf.cpp
 *
 *  Created on: Apr 6, 2017
 *  Author: ramiz
 */

#include <iostream>
#include <cmath>
#include "ukf.h"

UKF::UKF() {
  //TODO
}

UKF::~UKF() {
  //TODO
}

void UKF::Init() {

}

void UKF::GenerateSigmaPoints(MatrixXd * Xsig_out) {
  //set state dimension
  int n_x = 5;
  //set augmented state dimension
  int n_aug = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;

  //Process noise variance longitudinal acceleration
  double variance_a = std::pow(std_a, 2);

  //Process noise variance yaw acceleration
  double variance_yawdd = std::pow(std_yawdd, 2);

  //set noise covariance matrix
  MatrixXd Q = MatrixXd::Zero(2, 2);
  Q <<    variance_a, 0,
          0, variance_yawdd;

  //define spreading parameter lamda (design parameter of UKF)
  double lambda = 3 - n_aug;

  //set example state vector
  VectorXd x = VectorXd(n_x);
  x << 5.7441,
      1.3800,
      2.2049,
      0.5015,
      0.3528;

  //set example covariance matrix
  MatrixXd P = MatrixXd(n_x, n_x);
  P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug);
  //set augmented mean vector, considering noise mean is zero
  x_aug.head(5) = x;

  //create augmented covariance matrix with Process Covariance matrix P
  //and process noise covariance matrix Q
  MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);
  P_aug.topLeftCorner(n_x, n_x) = P;
  P_aug.bottomRightCorner(2, 2) = Q;

  //calculate total sigma points to generate
  int total_sigma_points = 2 * n_aug + 1;

  //create a sigma point matrix
  MatrixXd Xsigma = MatrixXd(n_aug, total_sigma_points);

  //calculate square root of P_aug by using Cholesky decomposition
  MatrixXd A = P_aug.llt().matrixL();

  //calculate sigma points
  //rule part1: first point is mean x so
  Xsigma.col(0) = x_aug;

  //rule part2: calculate next 1 to n_aug points
  Xsigma.block(0, 1, n_aug, n_aug) = (sqrt(lambda + n_aug) * A).colwise() + x_aug;

  //rule part3: calculate next n_aug to 2n_aug points
  Xsigma.block(0, n_aug + 1, n_aug, n_aug) = ( -1 * sqrt(lambda + n_aug) * A).colwise() + x_aug;

  *Xsig_out = Xsigma;
  /* expected result:
     Xsig_aug =
    5.7441  5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441
      1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38
    2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049
    0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015
    0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528
         0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0
         0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641
  */
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {
  //set state dimension
  int n_x = 5;
  //set augmented state dimension
  int n_aug = 7;

  //create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug <<
      5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
        1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
      2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
      0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
      0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
           0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
           0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //time diff in secs
  double delta_t = 0.1;

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for(int i = 0; i < Xsig_aug.cols(); ++i) {
//    if (i == 7) {
//      std::cout << "i == " << i << ", x = " << std::endl << Xsig_aug.col(i) << std::endl;
//    }
    Xsig_pred.col(i) = PredictSingleSigmaPoint(Xsig_aug.col(i), delta_t);
  }

  *Xsig_out = Xsig_pred;

// Expected result
//  Xsig_pred =
//   5.93553  6.06251  5.92217   5.9415  5.92361  5.93516  5.93705  5.93553  5.80832  5.94481  5.92935  5.94553  5.93589  5.93401  5.93553
//   1.48939  1.44673  1.66484  1.49719    1.508  1.49001  1.49022  1.48939   1.5308  1.31287  1.48182  1.46967  1.48876  1.48855  1.48939
//    2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.23954   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049  2.17026   2.2049
//   0.53678 0.473387 0.678098 0.554557 0.643644 0.543372  0.53678 0.538512 0.600173 0.395462 0.519003 0.429916 0.530188  0.53678 0.535048
//    0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528 0.387441 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528 0.318159
}

VectorXd UKF::PredictSingleSigmaPoint(const VectorXd & x_aug, double delta_t) {
  double v = x_aug(2);
  double yaw_angle = x_aug(3);
  double yaw_rate = x_aug(4);
  double noise_longitudinal_a = x_aug(5);
  double noise_yaw_a = x_aug(6);

  //do some common calculations
  double sin_yaw_angle = std::sin(yaw_angle);
  double cos_yaw_angle = std::cos(yaw_angle);
  double delta_t_square = delta_t * delta_t;

  double px_change = 0;
  double py_change = 0;

  if (std::fabs(yaw_rate) > 0.001) {
    //calculations specific to this case
    double c1 = v / yaw_rate;
    double c2 = yaw_angle + yaw_rate * delta_t;
    double sin_c2 = std::sin(c2);
    double cos_c2 = std::cos(c2);

    px_change = (c1 * (sin_c2 - sin_yaw_angle));
    py_change = (c1 * (-cos_c2 + cos_yaw_angle));
  }
  else {
    px_change = (v * cos_yaw_angle * delta_t);
    py_change = (v * sin_yaw_angle * delta_t);
  }

  VectorXd x_change = VectorXd(5);
  x_change << px_change,
              py_change,
              0,
              (yaw_rate * delta_t),
              0;

  VectorXd x_noise = VectorXd(5);
  x_noise <<  (0.5 * delta_t_square * cos_yaw_angle * noise_longitudinal_a),
              (0.5 * delta_t_square * sin_yaw_angle * noise_longitudinal_a),
              (delta_t * noise_longitudinal_a),
              (0.5 * delta_t_square * noise_yaw_a),
              (delta_t * noise_yaw_a);

  //add noise to change in x
  x_change += x_noise;

  //extract state vector x from augmented state vector (contains last 2 elements as noise values after augmentation)
  VectorXd x_old = x_aug.head(5);

  //calculate total change in x state
  VectorXd x_prediction = x_old + x_change;

  return x_prediction;
}



