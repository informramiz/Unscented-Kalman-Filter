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

	VectorXd x = VectorXd(5);
	MatrixXd P = MatrixXd(5, 5);
	ukf.PredictMeanAndCovariance(&x, &P);

	//print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

	return 0;
}
