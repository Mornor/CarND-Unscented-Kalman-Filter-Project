#include <iostream>
#include "ukf.h"

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

	is_initialized_ = false;

	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.003;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.5;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	std_a_ = 3;

	// Process noise standard deviation yaw acceleration in rad/s^2
	std_yawdd_ = 0.09;

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

	// Set state dimension (px, py, speed v (magnitude of the velocty), Si (angle of the orientation toward which the tracking object move
	// and yaw rate Si dot)
	n_x = 5;

	// Set measurement dimension, radar can measure r, phi, and r_dot
  	n_z = 3;

	// Set augmented dimension
	n_aug_ = 7; 

	// Sigma point spreading parameter
	lambda = 3 - n_aug_;

	// initial state vector
	x_ = VectorXd(n_x);

	// Augmented mean vector
  	VectorXd x_aug = VectorXd(n_aug_);

  	// Sigma point matrix
    Xsig_pred = MatrixXd(n_x, 2 * n_aug_ + 1); 		// predicted 
    Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); 	// augmented

  	// Initial covariance matrix
	P_ = MatrixXd(n_x, n_x);

	// Sensor matrices
	R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;

    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
	/*****************************************************************************
	*  Initialization
	****************************************************************************/

	if(!is_initialized_){
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			x_ = InitRadar(measurement_pack);
		}
    
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			x_ = InitLaser(measurement_pack);
		}

		// Init covariance matrix
		P_ <<	1, 0, 0, 0, 0,
				0, 1, 0, 0, 0,
				0, 0, 1, 0, 0,
				0, 0, 0, 1, 0,
				0, 0, 0, 0, 1;

		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;
		return; 
	}

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	Prediction(dt);
	previous_timestamp_ = measurement_pack.timestamp_;


	/*****************************************************************************
	*  Update
	****************************************************************************/
	/*if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(measurement_pack);
	}

	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(measurement_pack);
	}*/
}

/**
* Intial Radar measurement
*/
VectorXd UKF::InitRadar(const MeasurementPackage& measurement_pack){
	VectorXd x = VectorXd(n_x);

	float rho = measurement_pack.raw_measurements_[0]; 
	float phi = measurement_pack.raw_measurements_[1]; 
	float range_rate = measurement_pack.raw_measurements_[2]; 

	float px = rho * cos(phi);  
	float py = rho * sin(phi);  
	float v = range_rate;  
	float yaw = phi;  
	float yawd = 0.0;

	x << px, py, v, yaw, yawd;

	return x;
}

/**
* Intial Laser (Lidar) measurement
*/
VectorXd UKF::InitLaser(const MeasurementPackage &measurement_pack){
	VectorXd x = VectorXd(n_x);
	float rho = measurement_pack.raw_measurements_[0]; 
	float phi = measurement_pack.raw_measurements_[1]; 
	x << rho, phi, 0, 0, 0;
	return x;
}

MatrixXd UKF::AugmentedSigmaPoints() {
	// Augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	// Augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	// Sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	// Augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	// Augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_ * std_a_;
	P_aug(6, 6) = std_yawdd_ * std_yawdd_;

	// Square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	// Augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i < n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug_) * L.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda + n_aug_) * L.col(i);
	}

  return Xsig_aug;
}


MatrixXd UKF::PredictAugmentedSigmaPoints(double dt){
	//create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug_ + 1);

	//predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    	//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double px_p, py_p;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd * dt) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * dt));
    	} else {
      		px_p = p_x + v * dt * cos(yaw);
      		py_p = p_y + v * dt * sin(yaw);
    	}

		double v_p = v;
		double yaw_p = yaw + yawd * dt;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5 * nu_a * dt * dt * cos(yaw);
		py_p = py_p + 0.5 * nu_a * dt * dt * sin(yaw);
		v_p = v_p + nu_a * dt;

		yaw_p = yaw_p + 0.5 * nu_yawdd * dt * dt;
		yawd_p = yawd_p + nu_yawdd * dt;

		//write predicted sigma point into right column
		Xsig_pred(0, i) = px_p;
		Xsig_pred(1, i) = py_p;
		Xsig_pred(2, i) = v_p;
		Xsig_pred(3, i) = yaw_p;
		Xsig_pred(4, i) = yawd_p;
	}

	return Xsig_pred;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
	// Generate augmented sigma points matrix 
	MatrixXd Xsig_aug = AugmentedSigmaPoints();

	// Predict the future one usingd delta_t (dt)
	Xsig_pred = PredictAugmentedSigmaPoints(dt);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
	
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
	// TODO: Calulate the NIS

}
