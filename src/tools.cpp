#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// Accumulate squared results
	VectorXd residual(4);
	for (unsigned int i = 0; i < estimations.size(); ++i) {
		residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();

		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse / estimations.size();

	// Square root the results
	rmse = rmse.array().sqrt();

	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3, 4);

	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	// pre-computed terms
	float term1 = (px*px) + (py*py);
	float term2 = sqrt(term1);
	float term3 = pow(term1, 1.5);

	if (fabs(term1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by zero" << endl;
		return Hj;
	}

		//compute the Jacobian matrix
	Hj << px/term2, py/term2, 0, 0,
	     -py/term1, px/term1, 0, 0,
	     py*(vx*py - vy*px)/term3, px*(vy*px - vx*py)/term3, px/term2, py/term2;

	return Hj;
}
