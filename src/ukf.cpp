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
  std_a_ = 0.5;

  //TODO: Fine this noise parameter value
  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 5;

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

  //initialize measurement covariance matrix for Laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laser_px_ * std_laser_px_, 0,
              0, std_laser_py_ * std_laser_py_;

  //initialize measurement covariance matrix for Radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radar_r_ * std_radar_r_, 0, 0,
              0, std_radar_phi_ * std_radar_phi_, 0,
              0, 0, std_radar_rd_ * std_radar_rd_;

  //initialize noise covariance matrix
  Q_ = MatrixXd(2, 2);
  Q_ << std_a_ * std_a_, 0,
        0,  std_yawdd_ * std_yawdd_;

  //initialize the state transition matrix H for laser
  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  NIS_laser_ = 0;
  NIS_radar_ = 0;

  //augmented state mean dimension
  n_aug_ = 7;

  //state mean dimension
  n_x_ = 5;

  //radar measurement state dimension
  n_z_radar_ = 3;

  //laser measurement state dimension
  n_z_laser_ = 2;

  timestamp_ = 0;

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
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_.diagonal().fill(0.5);
}

UKF::~UKF() {
  //TODO
}

void UKF::Init(const MeasurementPackage& measurement_package) {

  /**
   * Initialize the state ekf_.x_ with the first measurement.
   * Create the covariance matrix.
   * Remember: you'll need to convert radar from polar to cartesian coordinates.
   */
  //first measurement
  std::cout << "EKF: " << std::endl;
  x_ = VectorXd(5);

  if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
    //Convert radar from polar to cartesian coordinates and initialize state.
    double rho = measurement_package.raw_measurements_[0];
    double phi = measurement_package.raw_measurements_[1];
    double px = rho * cos(phi);
    double py = rho * sin(phi);

    x_ <<   px,
            py,
            0,
            0,
            0;
  }
  else if (measurement_package.sensor_type_ == MeasurementPackage::LASER) {
    x_ <<   measurement_package.raw_measurements_[0],
            measurement_package.raw_measurements_[1],
            0,
            0,
            0;
  }

  // done initializing, no need to predict or update
  timestamp_ = measurement_package.timestamp_;
  is_initialized_ = true;
}

void UKF::GenerateSigmaPoints(MatrixXd * Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  //set augmented mean vector, considering noise mean is zero
  x_aug.head(5) = x_;

  //create augmented covariance matrix with Process Covariance matrix P
  //and process noise covariance matrix Q
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;

  //create a sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, total_sigma_points_);

  //calculate square root of P_aug by using Cholesky decomposition
  MatrixXd A = P_aug.llt().matrixL();

  //calculate sigma points
  //rule part1: first point is mean x so
  Xsig_aug.col(0) = x_aug;

  //rule part2: calculate next 1 to n_aug points
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = (std::sqrt(lambda_ + n_aug_) * A).colwise() + x_aug;

  //rule part3: calculate next n_aug to 2n_aug points
  Xsig_aug.block(0, n_aug_ + 1, n_aug_, n_aug_) = ( -1 * sqrt(lambda_ + n_aug_) * A).colwise() + x_aug;

  *Xsig_out = Xsig_aug;
}

void UKF::PredictSigmaPoints(const MatrixXd & Xsig_aug, double delta_t, MatrixXd* Xsig_out) {
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, total_sigma_points_);

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  for(int i = 0; i < Xsig_aug.cols(); ++i) {
    Xsig_pred.col(i) = PredictSingleSigmaPoint(Xsig_aug.col(i), delta_t);
  }

  *Xsig_out = Xsig_pred;
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

void UKF::PredictMeanAndCovariance() {
    //for ease of calculation let's make a matrix of weights
    MatrixXd W = weights_.transpose().replicate(Xsig_pred_.rows(), 1);

    //-----predict state mean-------
    //multiply weights with sigma points element wise
    MatrixXd Xsig_pred_weighted = (Xsig_pred_.array() * W.array()).matrix();
    //take element wise (rowwise) sum of all sigma points
    x_ = Xsig_pred_weighted.rowwise().sum();

    //-----predict state covariance matrix-----
    MatrixXd X_mean_distance = (Xsig_pred_.colwise() - x_);
    //angle normalization
    for (int i = 0; i < X_mean_distance.cols(); ++i) {
      X_mean_distance.col(i)(3) = Tools::NormalizeAngle(X_mean_distance.col(i)(3));
    }
    //multiply each sigma point mean distance with it's weight
    MatrixXd X_weighted_mean_distance = (X_mean_distance.array() * W.array()).matrix();
    P_ = X_weighted_mean_distance * X_mean_distance.transpose();
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //Matrix for sigma points in measurement space
  MatrixXd Zsig_pred = MatrixXd(n_z, total_sigma_points_);

  //Create a vector for predicted measurement mean
  VectorXd z_predicted = VectorXd(n_z);

  //Create a Matrix for predicted measurement covariance
  MatrixXd S = MatrixXd::Zero(n_z, n_z);

  //transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); ++i) {
    Zsig_pred.col(i) = MapToPolar(Xsig_pred_.col(i));
  }

  //for ease of calculation let's make a matrix of weights
  MatrixXd W = weights_.transpose().replicate(Zsig_pred.rows(), 1);

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
  S = Z_weighted_mean_distance * Z_mean_distance.transpose() + R_radar_;

  *z_out = z_predicted;
  *S_out = S;
  *Zsig_out = Zsig_pred;
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

void UKF::UpdateStateWithRadar(const VectorXd & z) {

  VectorXd z_pred = VectorXd(n_z_radar_);
  MatrixXd Zsig = MatrixXd(n_z_radar_, total_sigma_points_);
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);

  //predict radar measurement
  PredictRadarMeasurement(&z_pred, &S, &Zsig);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

  //for ease of calculation let's make a matrix of weights
  MatrixXd W = weights_.transpose().replicate(Xsig_pred_.rows(), 1);

  //calculate Cross correlation
  MatrixXd X_mean_distance = (Xsig_pred_.colwise() - x_);
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
  x_ = x_ + K * (z - z_pred);

  //update state covariance matrix
  P_ = P_ - K * S * K.transpose();

  //update Radar NIS
  NIS_radar_ = Tools::CalculateNIS(z_pred, z, S);
}

void UKF::UpdateStateWithLaser(const VectorXd& z) {
  //calculate the difference between
  //predicted and actual measurement
  VectorXd z_predicted = H_laser_ * x_;
  VectorXd y = z - z_predicted;

  //calculate Kalman Gain
  MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_;
  MatrixXd K = P_ * H_laser_.transpose() * S.inverse();

  //update the state based on Kalman gain and difference between our belief (x)
  //and measurement received.
  x_ = x_ + K * y;

  // update the state covariance/uncertainty based on measurement
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_laser_) * P_;

  //update Laser NIS
  NIS_radar_ = Tools::CalculateNIS(z_predicted, z, S);
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
     Init(measurement_package);
     is_initialized_ = true;
     return;
   }

//   std::cout << "measurement received: " << std::endl;

   //compute the time elapsed between the current and previous measurements in seconds
   double detla_t = (measurement_package.timestamp_ - timestamp_) / 1000000.0;

   //update timestamp to new measurement received timestamp
   timestamp_ = measurement_package.timestamp_;

   Predict(detla_t);

   if(measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
     double rho = measurement_package.raw_measurements_[0];
     double phi = measurement_package.raw_measurements_[1];
     double rho_dot = measurement_package.raw_measurements_[2];

     //check to avoid division by zero
     double px = x_(0);
     double py = x_(1);
     if(px == 0 || py == 0) {
       return;
     } else if(fabs(sqrt(px * px + py * py)) < 0.0001) {
       return;
     }

     VectorXd z = VectorXd(3);
     z <<  rho,
           phi,
           rho_dot;

     Update(z, measurement_package.sensor_type_);

   } else if(measurement_package.sensor_type_ == MeasurementPackage::LASER) {

     VectorXd z = VectorXd(2);
     z << measurement_package.raw_measurements_[0],
          measurement_package.raw_measurements_[1];
     Update(z, measurement_package.sensor_type_);
   }
}

void UKF::Predict(double delta_t) {
  //Generate augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, total_sigma_points_);
  GenerateSigmaPoints(&Xsig_aug);

  //Predict sigma points
  PredictSigmaPoints(Xsig_aug, delta_t, &Xsig_pred_);

  //Predict mean and covariance from sigma points
  PredictMeanAndCovariance();
}

void UKF::Update(const VectorXd& z, MeasurementPackage::SensorType sensor_type) {
  //check laser type
  //call PredictMeasurement____(Radar|Lidar) method
    //call UpdateStateWith__(Radar|Lidar) method
  if(sensor_type == MeasurementPackage::RADAR) {
    UpdateStateWithRadar(z);
  }
  else if(sensor_type == MeasurementPackage::LASER) {
    UpdateStateWithLaser(z);
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



