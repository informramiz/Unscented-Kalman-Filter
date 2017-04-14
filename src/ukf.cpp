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



