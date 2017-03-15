#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 30;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 30;

	// Laser measurement noise standard deviation position1 in m
	std_laspx_ = 0.15;

	// Laser measurement noise standard deviation position2 in m
	std_laspy_ = 0.15;

	// Radar measurement noise standard deviation radius in m
	std_radr_ = 0.3;

	// Radar measurement noise standard deviation angle in rad
	std_radphi_ = 0.03;

	// Radar measurement noise standard deviation radius change in m/s
	std_radrd_ = 0.3;

	// Set state dimension
	int n_x = 5;

	// Set augmented dimension
	int n_aug_ = 7; 

	// Sigma point spreading parameter
	double lambda = 3 - n_aug;

	// Weights vector
  	VectorXd weights = VectorXd(2*n_aug+1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

	}

	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {

	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	/**
	TODO:

	Complete this function! Estimate the object's location. Modify the state
	vector, x_. Predict sigma points, the state, and the state covariance matrix.
	*/
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
	/**
	TODO:

	Complete this function! Use lidar data to update the belief about the object's
	position. Modify the state vector, x_, and covariance, P_.

	You'll also need to calculate the lidar NIS.
	*/
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	// TODO: Calulate the NIS

	// Transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug + 1; i++) {    //2n+1 simga points

		// extract values for better readibility
		double p_x = Xsig_pred(0,i);
		double p_y = Xsig_pred(1,i);
		double v    = Xsig_pred(2,i);
		double yaw = Xsig_pred(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);						//r
		Zsig(1,i) = atan2(p_y,p_x);									//phi
		Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);	//r_dot
	}

	// Mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug+1; i++) {
			z_pred = z_pred + weights(i) * Zsig.col(i);
	}

	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug + 1; i++) {    // 2n+1 simga points
		// Residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// Angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights(i) * z_diff * z_diff.transpose();
	}

	// Add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z,n_z);
	R << std_radr*std_radr, 0, 0,
		 0, std_radphi*std_radphi, 0,
		 0, 0,std_radrd*std_radrd;
	S = S + R;
}
