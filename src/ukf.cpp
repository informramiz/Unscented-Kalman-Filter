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

  //define spreading parameter lamda (design parameter of UKF)
  double lambda = n_x - 3;

  //set example state
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

  //calculate total sigma points to generate
  int total_sigma_points = 2 * n_x + 1;

  //create a sigma point matrix
  MatrixXd Xsigma = MatrixXd(n_x, total_sigma_points);

  //calculate square root of P by using Cholesky decomposition
  MatrixXd A = P.llt().matrixL();

  //calculate sigma points
  //rule part1: first point is mean x so
  Xsigma.col(0) = x;

  //rule part2: calculate next 1 to n_x points
  Xsigma.block(0, 1, n_x, n_x) = (sqrt(lambda + n_x) * A).colwise() + x;

  //rule part3: calculate next n_x to 2n_x points
  Xsigma.block(0, n_x + 1, n_x, n_x) = ( -1 * sqrt(lambda + n_x) * A).colwise() + x;

  *Xsig_out = Xsigma;
}



