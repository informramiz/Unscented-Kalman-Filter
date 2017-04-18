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

	VectorXd x = VectorXd::Zero(5);
	MatrixXd P = MatrixXd::Zero(5, 5);
	ukf.UpdateState(&x, &P);

	//print result
  std::cout << "Updated mean x" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Updated covariance matrix P" << std::endl;
  std::cout << P << std::endl;

	return 0;
}
