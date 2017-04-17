//============================================================================
// Name        : UnscentedKalmanFilter.cpp
// Author      : Ramiz Raja
// Version     :
// Copyright   : MIT Licensed
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "ukf.h"

using namespace std;

int main() {
	UKF ukf;

	VectorXd z = VectorXd::Zero(3);
	MatrixXd S = MatrixXd::Zero(3, 3);
	ukf.PredictRadarMeasurement(&z, &S);

	//print result
  std::cout << "Predicted measurement mean z" << std::endl;
  std::cout << z << std::endl;
  std::cout << "Predicted measurement covariance matrix S" << std::endl;
  std::cout << S << std::endl;

	return 0;
}
