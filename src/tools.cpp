#include <iostream>
#include "tools.h"

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
		VectorXd rmse(4);
	vector<VectorXd> estm;
	rmse << 0,0,0,0;

	// check for empty estimations
	if (estimations.size() == 0) {
		cout << "CalculateRMSE ERROR: estimation vector size should not be zero" << endl;
		return rmse;
	}

	// check for estimations and ground_truth vectors being equal size.
	if (estimations.size() != ground_truth.size()) {
		cout << "CalculateRMSE ERROR: estimation vector size should equal ground_truth size" << endl;
		return rmse;
	}

	for (int i=0; i < estimations.size(); ++i) {
		VectorXd converted(4);
		converted << estimations[i][0],						// px
				 estimations[i][1],							// py
				 cos(estimations[i][3])*estimations[i][2],	// vx
				 sin(estimations[i][3])*estimations[i][2];	// vy
		estm.push_back(converted);
	}

	//accumulate squared residuals errors
	for (int i=0; i < estm.size(); ++i) {
		VectorXd errors = estm[i] - ground_truth[i];
		errors = errors.array() * errors.array();
		rmse += errors;
	}

	//calculate the mean from the rolling sum
	rmse /= estm.size();

	//calculate the squared root from the mean
	rmse = rmse.array().sqrt();

	//return the RMSE for this estimation.
	return rmse;
}
}