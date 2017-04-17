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

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
  //set state dimension
    int n_x = 5;

    //set augmented dimension
    int n_aug = 7;

    //define spreading parameter
    double lambda = 3 - n_aug;

    //create example matrix with predicted sigma points
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
    Xsig_pred <<
             5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
               1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
              2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
             0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
              0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

    //create vector for weights
    VectorXd weights = VectorXd(2 * n_aug + 1);

    //create vector for predicted state mean
    VectorXd x = VectorXd(n_x);

    //create matrix for predicted state covariance matrix
    MatrixXd P = MatrixXd(n_x, n_x);

    //set weights
    weights.fill(1 / (2 * (lambda + n_aug)));
    weights(0) = lambda / (lambda + n_aug);

    //for ease of calculation let's make a matrix of weights
    MatrixXd W = weights.transpose().replicate(Xsig_pred.rows(), 1);

    //-----predict state mean-------
    //multiply weights with sigma points element wise
    MatrixXd Xsig_pred_weighted = (Xsig_pred.array() * W.array()).matrix();
    //take element wise (rowwise) sum of all sigma points
    x = Xsig_pred_weighted.rowwise().sum();

    //-----predict state covariance matrix-----
    MatrixXd X_mean_distance = (Xsig_pred.colwise() - x);
    //angle normalization
    for (int i = 0; i < X_mean_distance.cols(); ++i) {
      X_mean_distance.col(i)(3) = NormalizeAngle(X_mean_distance.col(i)(3));
    }
    //multiply each sigma point mean distance with it's weight
    MatrixXd X_weighted_mean_distance = (X_mean_distance.array() * W.array()).matrix();
    P = X_weighted_mean_distance * X_mean_distance.transpose();

    *x_out = x;
    *P_out = P;

    /*
     expected result x:
     x =
        5.93637

        1.49035

        2.20528

        0.536853

        0.353577

      expected result p:
      P =
          0.00543425 -0.0024053 0.00341576 -0.00348196 -0.00299378

          -0.0024053 0.010845 0.0014923 0.00980182 0.00791091

          0.00341576 0.0014923 0.00580129 0.000778632 0.000792973

          -0.00348196 0.00980182 0.000778632 0.0119238 0.0112491

          -0.00299378 0.00791091 0.000792973 0.0112491 0.0126972
     */
}

double UKF::NormalizeAngle(double angle) {
  //angle normalization
  while (angle > M_PI)
    angle -= 2.*M_PI;

  while (angle < -M_PI)
    angle += 2. * M_PI;

  return angle;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //calculate total sigma points
  int total_sigma_points = 2 * n_aug + 1;

  //create vector for weights
  VectorXd weights = VectorXd(total_sigma_points);
  //set weights
  weights.fill(1 / (2 * (lambda + n_aug)));
  weights(0) = lambda / (lambda + n_aug);

  //radar measurement noise standard deviation radius/range in m
  double std_radar_range = 0.3;

  //radar measurement noise standard deviation angle in radians
  double std_radar_phi = 0.0175;

  //radar measurement noise standard deviation radius change (range rate) in m/s
  double std_radar_range_rate = 0.1;

  //calculate Measurement noise matrix
  MatrixXd R = MatrixXd::Zero(n_z, n_z);
  R <<  std_radar_range * std_radar_range, 0, 0,
        0, std_radar_phi * std_radar_phi, 0,
        0, 0, std_radar_range_rate * std_radar_range_rate;

  //A sample matrix for predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, total_sigma_points);
  Xsig_pred <<
           5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
             1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
            2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
           0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
            0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //Matrix for sigma points in measurement space
  MatrixXd Zsig_pred = MatrixXd(n_z, total_sigma_points);

  //Create a vector for predicted measurement mean
  VectorXd z_predicted = VectorXd(n_z);

  //Create a Matrix for predicted measurement covariance
  MatrixXd S = MatrixXd::Zero(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred.cols(); ++i) {
    Zsig_pred.col(i) = MapToPolar(Xsig_pred.col(i));
  }

  //for ease of calculation let's make a matrix of weights
  MatrixXd W = weights.transpose().replicate(Zsig_pred.rows(), 1);

  //-----predict state mean-------
  //multiply weights with sigma points element wise
  MatrixXd Zsig_pred_weighted = (Zsig_pred.array() * W.array()).matrix();
  //take element wise (rowwise) sum of all sigma points
  z_predicted = Zsig_pred_weighted.rowwise().sum();

  //-----predict state covariance matrix-----
  MatrixXd Z_mean_distance = (Zsig_pred.colwise() - z_predicted);
  //angle normalization
  for (int i = 0; i < Z_mean_distance.cols(); ++i) {
    Z_mean_distance.col(i)(1) = NormalizeAngle(Z_mean_distance.col(i)(1));
  }
  //multiply each sigma point mean distance with it's weight
  MatrixXd Z_weighted_mean_distance = (Z_mean_distance.array() * W.array()).matrix();
  S = Z_weighted_mean_distance * Z_mean_distance.transpose() + R;

  *z_out = z_predicted;
  *S_out = S;
}

Eigen::VectorXd UKF::MapToPolar(const Eigen::VectorXd& x) {

  float px = x(0);
  float py = x(1);
  float v = x(2);
  float yaw_angle = x(3);
//  float yaw_rate = x(4);

  //calculate speed vx and vy components
  float vx = v * cos(yaw_angle);
  float vy = v * sin(yaw_angle);

  float px2_py2_sum = px * px + py * py;
  float px2_py2_sum_sqrt = sqrt(px2_py2_sum);

  Eigen::VectorXd z_predicted(3);
  //apply non-linear transformation h(x)
  z_predicted <<  px2_py2_sum_sqrt,
                  atan2(py, px),
                  ((px * vx + py * vy) / px2_py2_sum_sqrt);

  return z_predicted;
}



