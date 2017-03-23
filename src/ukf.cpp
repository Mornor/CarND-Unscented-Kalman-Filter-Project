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

	// Sensor matrices
	R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_, 0,
                0, std_laspy_;

    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

	// initial state vector
	x_ = VectorXd(5);

	// initial covariance matrix
	P_ = MatrixXd(5, 5);
	P_ << 	1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1000, 0, 0,
            0, 0, 0, 1000, 0,
            0, 0, 0, 0, 1000;

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

	// Augmented state covariance
  	MatrixXd P_aug = MatrixXd(7, 7);
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_ * std_a_;
	P_aug(6,6) = std_yawdd_ * std_yawdd_;

	// Set state dimension (px, py, speed v (magnitude of the velocty), Si (angle of the orientation toward which the tracking object move
	// and yaw rate Si dot)
	int n_x = 5;

	// Set measurement dimension, radar can measure r, phi, and r_dot
  	int n_z = 3;

	// Set augmented dimension
	int n_aug_ = 7; 

	// Sigma point spreading parameter
	double lambda = 3 - n_aug_;

	// Augmented mean vector
  	VectorXd x_aug = VectorXd(n_aug_);

	// Sigma point matrix
	Xsig = MatrixXd(n_x, 2 * n_x + 1); 				// initial
    Xsig_pred = MatrixXd(n_x, 2 * n_aug_ + 1); 		// predicted 
    Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); 	// augmented

	// Weights vector
  	VectorXd weights = VectorXd(2 * n_aug_ +  1);
  	double weight_0 = lambda / (lambda + n_aug_);
  	weights(0) = weight_0;
  	for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
    	double weight = 0.5 / (n_aug_ + lambda);
    	weights(i) = weight;
  	}

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

		// Generatte sigma points
		GenerateSigmaPoints(&Xsig);

		// Generate augmented sigma points matrix
		AugmentedSigmaPoints(&Xsig_aug);

		// first measurement, we do not now px and py (position x and y), neither vx and vy
        float p_x = 0; 
        float p_y = 0;
        float v_x = 0; 
        float v_y = 0; 
        float v = 0;  // Speed, magnitude of the velocity
        float yaw_rate = 0; 
        float yaw_rate_dot = 0; 

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            float rho = measurement_pack.raw_measurements_[0];
            float phi = measurement_pack.raw_measurements_[1];
            float range_rate = measurement_pack.raw_measurements_[2];

            p_x = rho * cos(phi);
            p_y = rho * sin(phi);

            // Compute v (magnitude of the speed), by using v_x and v_y first
            v_x = range_rate * cos(phi);
            v_y = range_rate * sin(phi);
            v = sqrt(v_x * v_x + v_y * v_y);

            // Compute yaw_rate, make sure it is possible by checking v_x
            yaw_rate = fabs(v_x) > 0.0001 ? atan(v_y / v_x) : 0;

            if(fabs(p_x) < 0.0001){
                p_x = 1;
                P_(0,0) = 1000;
            }
            
            if(fabs(p_y) < 0.0001){
                p_y = 1;
                P_(1,1) = 1000;
            }

        }
    
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            p_x = measurement_pack.raw_measurements_[0];
            p_y = measurement_pack.raw_measurements_[1];
        }

    	x_ << p_x, p_y, v, yaw_rate, yaw_rate_dot;

	}

	previous_timestamp_ = measurement_pack.timestamp_;
	is_initialized_ = true;
	return; 

	/*****************************************************************************
	*  Prediction
	****************************************************************************/

	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
	Prediction(dt);
	previous_timestamp_ = measurement_pack.timestamp_;


	/*****************************************************************************
	*  Update
	****************************************************************************/
	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(measurement_pack);
	}

	else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(measurement_pack);
	}
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
	//set first column of sigma point matrix
	Xsig.col(0) = x_;

	// Calculate square root of P
  	MatrixXd A = P_.llt().matrixL();

  	// set remaining sigma points
  	for (int i = 0; i < n_x; i++){
    	Xsig.col(i+1)     = x_ + sqrt(lambda+n_x) * A.col(i);
    	Xsig.col(i+1+n_x) = x_ - sqrt(lambda+n_x) * A.col(i);
  	}
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i < n_aug_; i++){
		Xsig_aug.col(i+1)        = x_aug + sqrt(lambda+n_aug_) * L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
	}
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
	// Predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++){

		// Extract values for better readability
		double p_x = Xsig_aug(0,i);
		double p_y = Xsig_aug(1,i);
		double v = Xsig_aug(2,i);
		double yaw_rate = Xsig_aug(3,i);
		double yaw_rate_dot = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		// Predicted state values
		double px_p, py_p;

		// Avoid division by zero
		if (fabs(yaw_rate_dot) > 0.001) {
			px_p = p_x + v/yaw_rate_dot * (sin(yaw_rate + yaw_rate_dot * dt) - sin(yaw_rate));
			py_p = p_y + v/yaw_rate_dot * (cos(yaw_rate) - cos(yaw_rate + yaw_rate_dot * dt));
		} else {
			px_p = p_x + v * dt * cos(yaw_rate);
			py_p = p_y + v * dt * sin(yaw_rate);
		}

		double v_p = v;
		double yaw_p = yaw_rate + yaw_rate_dot * dt ;
		double yawd_p = yaw_rate_dot;

		// Add noise
		px_p = px_p + 0.5 * nu_a * dt * dt  * cos(yaw_rate);
		py_p = py_p + 0.5 * nu_a * dt * dt  * sin(yaw_rate);
		v_p = v_p + nu_a * dt ;

		yaw_p = yaw_p + 0.5 * nu_yawdd * dt * dt ;
		yawd_p = yawd_p + nu_yawdd * dt ;

		// Write predicted sigma point into right column
		Xsig_pred(0,i) = px_p;
		Xsig_pred(1,i) = py_p;
		Xsig_pred(2,i) = v_p;
		Xsig_pred(3,i) = yaw_p;
		Xsig_pred(4,i) = yawd_p;
	}
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
	// TODO: Calulate the NIS
	VectorXd z = measurement_pack.raw_measurements_;

    VectorXd z_pred = H_laser_ * x_;
    VectorXd z_diff = z - z_pred;
    MatrixXd Ht = H_laser_.transpose();
    MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    //new estimate
    x_ = x_ + (K * z_diff);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_laser_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
	// TODO: Calulate the NIS

	// Transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {    //2n+1 simga points

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
	for (int i=0; i < 2 * n_aug_ + 1; i++) {
			z_pred = z_pred + weights(i) * Zsig.col(i);
	}

	// Measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {    // 2n+1 simga points
		// Residual
		VectorXd z_diff = Zsig.col(i) - z_pred;

		// Angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights(i) * z_diff * z_diff.transpose();
	}

	// Add measurement noise covariance matrix
	MatrixXd R = MatrixXd(n_z,n_z);
	R << std_radr_ * std_radr_, 0, 0,
		 0, std_radphi_ * std_radphi_, 0,
		 0, 0, std_radrd_ * std_radrd_;
	S = S + R;

	// Calculate cross correlation matrix
  	MatrixXd Tc = MatrixXd(n_x, n_z); // Cross correlation Tc
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
		// Residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// Angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2. * M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;

		// State difference
		VectorXd x_diff = Xsig_pred.col(i) - x_;
		
		// Angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2. * M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;

		Tc = Tc + weights(i) * x_diff * z_diff.transpose();
	}

	// Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// Residual
	VectorXd z = measurement_pack.raw_measurements_;
	VectorXd z_diff = z - z_pred;

	// Angle normalization
	while (z_diff(1)> M_PI) z_diff(1) -= 2. * M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;

	// Update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
}
