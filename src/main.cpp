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

	MatrixXd Xsig = MatrixXd(11, 5);
	UKF ukf;
	ukf.GenerateSigmaPoints(&Xsig);

	//print result
  std::cout << "Xsig = " << std::endl << Xsig << std::endl;

	return 0;
}
