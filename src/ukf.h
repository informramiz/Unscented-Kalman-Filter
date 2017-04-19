/*
 * ukf.h
 *
 *  Created on: Apr 6, 2017
 *  Author: ramiz
 *  Description: Header file which contains all the functions for UKF
 */

#ifndef UKF_H_
#define UKF_H_

#include <vector>
#include "Eigen/Dense"
#include "tools.h"
#include "measurement_package.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class UKF {
public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Init Initializes Unscented Kalman filter
   */
  void Init();

  /**
   * Following functions define structure of UKF
   * and I will add code to each of them incrementally
   */
  void GenerateSigmaPoints(MatrixXd* Xsig_out);
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);
  void SigmaPointPrediction(const MatrixXd & Xsig_aug, double delta_t, MatrixXd* Xsig_out);
  void PredictMeanAndCovariance();
  void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out);
  void UpdateRadarState(VectorXd* x_out, MatrixXd* P_out);

  /**
   * ProcessMeasurement
   * @param measurement_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage measurement_package);

  /**
   * Returns mean state x
   */
  Eigen::VectorXd GetMeanState() const;

  /**
   * Returns Laser NIS
   */
  double GetLaserNIS() const;

  /**
   * Returns Radar NIS
   */
  double GetRadarNIS() const;
private:
  /**
   * Predicts a single sigma point based on augmented sigma point
   * provided
   * @param x, augmented 7D sigma point
   */
  VectorXd PredictSingleSigmaPoint(const VectorXd & x_aug, double delta_t);

  /**
   * Non-linear function h(x) that maps cartesian coordinates @param x =(px, py, vx, vy)
   * to polar coordinates (range=rho, angle=phi, range_rate=rho_dot)
   */
  Eigen::VectorXd MapToPolar(const Eigen::VectorXd& x);

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laser_px_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laser_py_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radar_r_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radar_phi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radar_rd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Total sigma points count
  int total_sigma_points_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* Measurement Covariance Matrix for Radar
  MatrixXd R_radar_;

  ///* Measurement covariance matrix for Laser
  MatrixXd R_laser_;

  ///* Noise covariance matrix
  MatrixXd Q_;
};



#endif /* UKF_H_ */
