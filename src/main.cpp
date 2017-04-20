//============================================================================
// Name        : UnscentedKalmanFilter.cpp
// Author      : Ramiz Raja
// Version     :
// Copyright   : MIT Licensed
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void CheckArguments(int argc, char* argv[]);
void ReadMeasurements(string in_file_name, vector<MeasurementPackage> &measurement_pack_list,
                       vector<GroundTruthPackage>& gt_pack_list);
void CheckFiles(const string& in_name, const string& out_name);
void TestUkf(const string& in_file_name, const string& out_file_name);
void RunUkf(const vector<MeasurementPackage>& measurement_pack_list,
            const vector<GroundTruthPackage>& gt_pack_list,
            const string& out_file_name);

void RunUkfWithFileOutputOptimizedForPythonVisualization(const vector<MeasurementPackage>& measurement_pack_list,
            const vector<GroundTruthPackage>& gt_pack_list,
            const string& out_file_name);

using namespace std;

int main(int argc, char* argv[]) {
  //check validity of arguments
  CheckArguments(argc, argv);

  string in_file_name = argv[1];
  string out_file_name = argv[2];
  //validate files existence
  CheckFiles(in_file_name, out_file_name);

  TestUkf(in_file_name, out_file_name);

	return 0;
}

void TestUkf(const string& in_file_name, const string& out_file_name) {
  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;
  ReadMeasurements(in_file_name, measurement_pack_list, gt_pack_list);

  RunUkf(measurement_pack_list, gt_pack_list, out_file_name);
}

void RunUkf(const vector<MeasurementPackage>& measurement_pack_list,
            const vector<GroundTruthPackage>& gt_pack_list,
            const string& out_file_name) {

  ofstream out_file(out_file_name, ofstream::out);
  if (!out_file.is_open()) {
    std::cerr << "Error opening file: " << out_file_name << std::endl;
    exit(EXIT_FAILURE);
  }

  // column names for output file
  out_file << "time_stamp" << "\t";
  out_file << "px_state" << "\t";
  out_file << "py_state" << "\t";
  out_file << "v_state" << "\t";
  out_file << "yaw_angle_state" << "\t";
  out_file << "yaw_rate_state" << "\t";
  out_file << "sensor_type" << "\t";
  out_file << "NIS" << "\t";
  out_file << "px_measured" << "\t";
  out_file << "py_measured" << "\t";
  out_file << "px_ground_truth" << "\t";
  out_file << "py_ground_truth" << "\t";
  out_file << "vx_ground_truth" << "\t";
  out_file << "vy_ground_truth" << "\n";

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  //Call the UKF-based fusion
  size_t N = measurement_pack_list.size();
  for (size_t k = 0; k < N; ++k) {
    // start filtering from the second frame (the speed is unknown in the first
    // frame)
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // timestamp
    out_file << measurement_pack_list[k].timestamp_ << "\t"; // pos1 - est

    // output the estimation
    VectorXd x = ukf.GetMeanState();
    out_file << x(0) << "\t"; //px
    out_file << x(1) << "\t"; //py
    out_file << x(2) << "\t"; //speed
    out_file << x(3) << "\t"; //yaw
    out_file << x(4) << "\t"; //yaw_rate

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // sensor type
      out_file << "lidar" << "\t";
      // NIS value
      out_file << ukf.GetLaserNIS() << "\t";

      // output the estimation
      out_file << measurement_pack_list[k].raw_measurements_(0) << "\t"; //p1_meas
      out_file << measurement_pack_list[k].raw_measurements_(1) << "\t"; //p2_meas

    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // sensor type
      out_file << "radar" << "\t";

      // NIS value
      out_file << ukf.GetRadarNIS() << "\t";

      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file << ro * cos(phi) << "\t"; // p1_meas
      out_file << ro * sin(phi) << "\t"; // ps_meas
    }

    // output the ground truth packages
    out_file << gt_pack_list[k].gt_values_(0) << "\t"; //p1
    out_file << gt_pack_list[k].gt_values_(1) << "\t"; //p2
    out_file << gt_pack_list[k].gt_values_(2) << "\t"; //vx
    out_file << gt_pack_list[k].gt_values_(3) << "\n"; //vy

    // convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    x = ukf.GetMeanState();
    float x_estimate_ = x(0);
    float y_estimate_ = x(1);
    float vx_estimate_ = x(2) * cos(x(3));
    float vy_estimate_ = x(2) * sin(x(3));

    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }

  // compute the accuracy (RMSE)
  cout << "Accuracy - RMSE:" << endl << Tools::CalculateRMSE(estimations, ground_truth) << endl;

  if (out_file.is_open()) {
    out_file.close();
  }
}


void RunUkfWithFileOutputOptimizedForPythonVisualization(const vector<MeasurementPackage>& measurement_pack_list,
            const vector<GroundTruthPackage>& gt_pack_list,
            const string& out_file_name) {

  ofstream out_file(out_file_name, ofstream::out);
  if (!out_file.is_open()) {
    std::cerr << "Error opening file: " << out_file_name << std::endl;
    exit(EXIT_FAILURE);
  }

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  //Call the UKF-based fusion
  size_t N = measurement_pack_list.size();
  for (size_t k = 0; k < N; ++k) {
    // start filtering from the second frame (the speed is unknown in the first
    // frame)
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // timestamp
    out_file << measurement_pack_list[k].timestamp_ << "\t"; // pos1 - est

    // output the estimation
    VectorXd x = ukf.GetMeanState();
    out_file << x(0) << "\t"; //px
    out_file << x(1) << "\t"; //py
    out_file << x(2) << "\t"; //speed
    out_file << x(3) << "\t"; //yaw
    out_file << x(4) << "\t"; //yaw_rate

    double NIS_laser = 0;
    double NIS_radar = 0;

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // NIS value
      NIS_laser = ukf.GetLaserNIS();

      // output the estimation
      out_file << measurement_pack_list[k].raw_measurements_(0) << "\t"; //p1_meas
      out_file << measurement_pack_list[k].raw_measurements_(1) << "\t"; //p2_meas

    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // NIS value
      NIS_radar = ukf.GetRadarNIS();

      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file << ro * cos(phi) << "\t"; // p1_meas
      out_file << ro * sin(phi) << "\t"; // ps_meas
    }

    // output the ground truth packages
    out_file << gt_pack_list[k].gt_values_(0) << "\t"; //p1
    out_file << gt_pack_list[k].gt_values_(1) << "\t"; //p2
    out_file << gt_pack_list[k].gt_values_(2) << "\t"; //vx
    out_file << gt_pack_list[k].gt_values_(3) << "\t"; //vy

    out_file << NIS_laser << "\t";
    out_file << NIS_radar << "\n";

    // convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    x = ukf.GetMeanState();
    float x_estimate_ = x(0);
    float y_estimate_ = x(1);
    float vx_estimate_ = x(2) * cos(x(3));
    float vy_estimate_ = x(2) * sin(x(3));

    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }

  // compute the accuracy (RMSE)
  cout << "Accuracy - RMSE:" << endl << Tools::CalculateRMSE(estimations, ground_truth) << endl;

  if (out_file.is_open()) {
    out_file.close();
  }
}

void ReadMeasurements(string in_file_name, vector<MeasurementPackage> &measurement_pack_list,
                       vector<GroundTruthPackage>& gt_pack_list) {
  ifstream in_file(in_file_name.c_str(), ifstream::in);

  if (!in_file.is_open()) {
    std::cerr << "Can not open file: " << in_file_name << std::endl;
    exit(EXIT_FAILURE);
  }

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  string line;
  while (getline(in_file, line) && !in_file.eof()) {

    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long timestamp;

    // reads first element from the current line
    iss >> sensor_type;
    if (sensor_type.compare("L") == 0) {
      // LASER MEASUREMENT

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float x;
      float y;
      iss >> x;
      iss >> y;
      meas_package.raw_measurements_ << x, y;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // RADAR MEASUREMENT
      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  if (in_file.is_open()) {
    in_file.close();
  }
}

void CheckArguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void CheckFiles(const string& in_name, const string& out_name) {
  ifstream in_file(in_name.c_str(), ifstream::in);

  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }
  in_file.close();

  ofstream out_file(out_name.c_str(), ofstream::out);
  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
  out_file.close();
}
