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
  //if set false means filter is not yet initialized
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //TODO: Fine this noise parameter value
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  //TODO: Fine this noise parameter value
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  /**********************
   * ---Start---
   * Below noise values provided by manufacturer and should not be changed
   ***********************/

  // Laser measurement noise standard deviation position1 in m
  std_laser_px_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laser_py_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radar_r_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radar_phi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radar_rd_ = 0.3;

  /**********************
   * ---End---
   ***********************/

  NIS_laser_ = 0;
  NIS_radar_ = 0;

  //augmented state mean dimension
  n_aug_ = 7;

  //state mean dimension
  n_x_ = 5;

  time_us_ = 0;

  //calculate total sigma points count
  total_sigma_points_ = 2 * n_aug_ + 1;

  //set spreading parameter
  lambda_ = 3 - n_aug_;

  //initialize weights
  weights_ = VectorXd(total_sigma_points_);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  //initialize predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, total_sigma_points_);

  //initialize state mean and convariance
  x_ = VectorXd::Zero(n_x_);
  P_ = MatrixXd::Zero(n_x_, n_x_);
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
      X_mean_distance.col(i)(3) = Tools::NormalizeAngle(X_mean_distance.col(i)(3));
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
    Z_mean_distance.col(i)(1) = Tools::NormalizeAngle(Z_mean_distance.col(i)(1));
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

void UKF::UpdateRadarState(VectorXd* x_out, MatrixXd* P_out) {
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

  VectorXd weights = VectorXd(total_sigma_points);
  //set weights
  weights.fill(1 / (2 * (lambda + n_aug)));
  weights(0) = lambda / (lambda + n_aug);

  //create example matrix with predicted sigma points in state space
  MatrixXd Xsig_pred = MatrixXd(n_x, total_sigma_points);
  Xsig_pred <<
      5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
      1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
      2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
      0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
      0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
      5.93637,
      1.49035,
      2.20528,
      0.536853,
      0.353577;

  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
      0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
      -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
      0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
      -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
      -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
      0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred <<
      6.12155,
      0.245993,
      2.10313;

  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S <<
      0.0946171, -0.000139448,   0.00407016,
      -0.000139448,  0.000617548, -0.000770652,
      0.00407016, -0.000770652,    0.0180917;

  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
      5.9214,   //rho in m
      0.2187,   //phi in rad
      2.0062;   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //for ease of calculation let's make a matrix of weights
  MatrixXd W = weights.transpose().replicate(Xsig_pred.rows(), 1);

  //calculate Cross correlation
  MatrixXd X_mean_distance = (Xsig_pred.colwise() - x);
  //angle normalization
  for (int i = 0; i < X_mean_distance.cols(); ++i) {
    X_mean_distance.col(i)(3) = Tools::NormalizeAngle(X_mean_distance.col(i)(3));
  }
  MatrixXd X_mean_distance_weighted = (X_mean_distance.array() * W.array()).matrix();

  MatrixXd Z_mean_distance = Zsig.colwise() - z_pred;
  //angle normalization
  for (int i = 0; i < Z_mean_distance.cols(); ++i) {
    Z_mean_distance.col(i)(1) = Tools::NormalizeAngle(Z_mean_distance.col(i)(1));
  }
  Tc = X_mean_distance_weighted * Z_mean_distance.transpose();

  //calculate Kalman Gain
  MatrixXd K = Tc * S.inverse();

  //update state mean
  x = x + K * (z - z_pred);

  //update state covariance matrix
  P = P - K * S * K.transpose();

  *x_out = x;
  *P_out = P;

  /*
   * expected result x:
     x =
      5.92276

      1.41823

      2.15593

      0.489274

      0.321338

    expected result P:
    P =
      0.00361579 -0.000357881 0.00208316 -0.000937196 -0.00071727

      -0.000357881 0.00539867 0.00156846 0.00455342 0.00358885

      0.00208316 0.00156846 0.00410651 0.00160333 0.00171811

      -0.000937196 0.00455342 0.00160333 0.00652634 0.00669436

      -0.00071719 0.00358884 0.00171811 0.00669426 0.00881797
   * */
}

/**
 * @param {MeasurementPackage} measurement_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
  /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
      /**
       * Initialize the state ekf_.x_ with the first measurement.
       * Create the covariance matrix.
       * Remember: you'll need to convert radar from polar to cartesian coordinates.
       */
      //first measurement
      std::cout << "EKF: " << std::endl;
      x_ = VectorXd(5);
      x_ << 1, 1, 1, 1, 1;

      if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
        //Convert radar from polar to cartesian coordinates and initialize state.
        double rho = measurement_package.raw_measurements_[0];
        double phi = measurement_package.raw_measurements_[1];
        double px = rho * cos(phi);
        double py = rho * sin(phi);

        x_ <<  px,
              py,
              0,
              0,
              0;
      }
      else if (measurement_package.sensor_type_ == MeasurementPackage::LASER) {
        x_ <<  measurement_package.raw_measurements_[0],
                    measurement_package.raw_measurements_[1],
                    0,
                    0,
                    0;
      }

      // done initializing, no need to predict or update
      timestamp_ = measurement_package.timestamp_;
      is_initialized_ = true;
      return;
    }
}

Eigen::VectorXd UKF::GetMeanState() const {
  return x_;
}

double UKF::GetLaserNIS() const {
  return NIS_laser_;
}

double UKF::GetRadarNIS() const {
  return NIS_radar_;
}



