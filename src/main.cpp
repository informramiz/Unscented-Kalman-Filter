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
	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!

	MatrixXd Xsig_pred = MatrixXd(15, 5);
	UKF ukf;
	ukf.SigmaPointPrediction(&Xsig_pred);

	//print result
  std::cout << "Xsig = " << std::endl << Xsig_pred << std::endl;

	return 0;
}
