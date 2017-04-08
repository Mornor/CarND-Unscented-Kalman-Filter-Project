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
	std_yawdd_ = 0.9;

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

	// Radar measurement dimension can measure r, phi, and r_dot
	n_z = 3;

	// Lidar can only measure px and py
  	n_z_lidar = 2;

	// Set augmented dimension
	n_aug_ = 7; 

	// Sigma point spreading parameter
	lambda = 3 - n_aug_;

	// initial state vector
	x_ = VectorXd(n_x);

	// Augmented mean vector
  	x_aug = VectorXd(n_aug_);

  	// Sigma point matrix
  	Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1); 	// augmented 
    Xsig_pred = MatrixXd(n_x, 2 * n_aug_ + 1); 		// predicted (augmented one)

  	// Initial covariance matrix
	P_ = MatrixXd(n_x, n_x);

	// Sensor matrices
	R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_ * std_laspx_, 0,
                0, std_laspy_ * std_laspy_;

    R_radar_ = MatrixXd(3, 3);
    R_radar_ <<	std_radr_ * std_radr_, 0, 0,
				0, std_radphi_ * std_radphi_, 0,
				0, 0, std_radrd_ * std_radrd_;

    // Vector for weights
  	weights = VectorXd(2 * n_aug_ + 1);
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

	// Dividing large time steps into smaller prediction intervals helps to maintain numerical stability.
	while (dt > 0.1){
		Prediction(0.05);
		dt -= 0.05;
	}
	
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

/**
* Intial Radar measurement
*/
VectorXd UKF::InitRadar(const MeasurementPackage& measurement_pack){
	VectorXd x = VectorXd(n_x);

	double rho = measurement_pack.raw_measurements_[0]; 
	double phi = measurement_pack.raw_measurements_[1]; 
	double range_rate = measurement_pack.raw_measurements_[2]; 

	double px = rho * cos(phi);  
	double py = rho * sin(phi);  
	double v = range_rate;  
	double yaw = phi;  
	double yawd = 0.0;

	x << px, py, v, yaw, yawd;

	return x;
}

void UKF::PredictMeanAndCovariance(){
	// Set weights
	double weight_0 = lambda / (lambda + n_aug_);
	weights(0) = weight_0;
	
	for (int i=1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
		double weight = 0.5 / (n_aug_ + lambda);
		weights(i) = weight;
	}

	// Predicted state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_ = x_+ weights(i) * Xsig_pred.col(i);
	}

	// Predicted state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		VectorXd x_diff = Xsig_pred.col(i) - x_; // state difference VectorXd x_diff = Xsig_pred.col(i) - x_
		//angle normalization
		while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
		while (x_diff(3) <- M_PI) x_diff(3) += 2.* M_PI;
		P_ = P_ + weights(i) * x_diff * x_diff.transpose();
	}
}

/**
* Intial Laser (Lidar) measurement
*/
VectorXd UKF::InitLaser(const MeasurementPackage &measurement_pack){
	VectorXd x = VectorXd(n_x);
	double px = measurement_pack.raw_measurements_[0]; 
	double py = measurement_pack.raw_measurements_[1]; 
	x << px, py, 0, 0, 0;
	return x;
}

void UKF::AugmentedSigmaPoints() {
	// Augmented mean vector
	x_aug = VectorXd(n_aug_);

	// Augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

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
}


void UKF::PredictAugmentedSigmaPoints(double dt){
	// Predict sigma points
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    	//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		// Predicted state values
		double px_p, py_p;

		// Avoid division by zero
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

		// Add noise
		px_p = px_p + 0.5 * nu_a * dt * dt * cos(yaw);
		py_p = py_p + 0.5 * nu_a * dt * dt * sin(yaw);
		v_p = v_p + nu_a * dt;

		yaw_p = yaw_p + 0.5 * nu_yawdd * dt * dt;
		yawd_p = yawd_p + nu_yawdd * dt;

		// Write predicted sigma point into right column
		Xsig_pred(0, i) = px_p;
		Xsig_pred(1, i) = py_p;
		Xsig_pred(2, i) = v_p;
		Xsig_pred(3, i) = yaw_p;
		Xsig_pred(4, i) = yawd_p;
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} dt the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
	// Generate augmented sigma points matrix
	AugmentedSigmaPoints(); // Xsig_aug 

	// Predict the future one usingd delta_t (dt)
	PredictAugmentedSigmaPoints(dt); // Xsig_pred

	// Predict mean and covariance
	PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
	double px = measurement_pack.raw_measurements_[0]; 
	double py = measurement_pack.raw_measurements_[1]; 

	VectorXd z = VectorXd(n_z_lidar);
	z << px, py; 

	// Predict Lidar measurement
	MatrixXd Zsig;
	VectorXd z_pred;
	MatrixXd S;
	PredictLidarMeasurement(&Zsig, &z_pred, &S);

	// Update the state matrix x_ and the covariance matrix P_
	UpdateState(z, Zsig, z_pred, S, n_z_lidar);

	// update NIS
  	NIS_laser_ = (measurement_pack.raw_measurements_-z_pred).transpose()*S.inverse()*(measurement_pack.raw_measurements_-z_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {

	double rho = measurement_pack.raw_measurements_[0]; 
	double phi = measurement_pack.raw_measurements_[1]; 
	double range_rate = measurement_pack.raw_measurements_[2]; 

	VectorXd z = VectorXd(n_z);
	z << rho, phi, range_rate;

	// Predict Radar measurement
	MatrixXd Zsig;
	VectorXd z_pred;
	MatrixXd S; 
	PredictRadarMeasurement(&Zsig, &z_pred, &S);

	// Update the state matrix x_ and the covariance matrix P_
	UpdateState(z, Zsig, z_pred, S, n_z);

	// Update NIS
	NIS_radar_ = (measurement_pack.raw_measurements_-z_pred).transpose()*S.inverse()*(measurement_pack.raw_measurements_-z_pred);
}

void UKF::PredictRadarMeasurement(MatrixXd *Zsig_out, VectorXd *z_pred_out, MatrixXd *S_out) {
	// Matrix for sigma points in measurement space
  	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  	VectorXd z_pred = VectorXd(n_z);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
		// Extract values for better readibility
		double p_x = Xsig_pred(0,i);
		double p_y = Xsig_pred(1,i);
		double v  = Xsig_pred(2,i);
		double yaw = Xsig_pred(3,i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		// measurement model
		Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);	//r
		Zsig(1,i) = atan2(p_y,p_x);				//phi 			
		// Avoid division by zero                                
		if (Zsig(0, i) < 0.001) {
      		Zsig(2, i) = (p_x * v1 + p_y * v2) / 0.001;  //r_dot
    	} else {
      		Zsig(2, i) = (p_x * v1 + p_y * v2) / Zsig(0, i);  //r_dot;
    	}
	}

	// Mean predicted measurement
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred = z_pred + weights(i) * Zsig.col(i);
	}

	// Measurement covariance matrix 
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  
		VectorXd z_diff = Zsig.col(i) - z_pred; //residual
		// Angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights(i) * z_diff * z_diff.transpose();
	}

	// Add measurement noise covariance matrix
	S = S + R_radar_;

	*z_pred_out = z_pred;
	*Zsig_out = Zsig; 
	*S_out = S; 
}

void UKF::PredictLidarMeasurement(MatrixXd *Zsig_out, VectorXd *z_pred_out, MatrixXd *S_out) {
	MatrixXd Zsig = MatrixXd(n_z_lidar, 2 * n_aug_ + 1);
	VectorXd z_pred = VectorXd(n_z_lidar);

	for (int i = 0; i < 2 * n_aug_ + 1; i++) { 
		// Measurement model
		Zsig(0,i) = Xsig_pred(0,i); // px
		Zsig(1,i) = Xsig_pred(1,i); // py
	}

	// Mean predicted measurement
	z_pred.fill(0.0);
	for (int i=0; i < 2 * n_aug_+1; i++) {
		z_pred = z_pred + weights(i) * Zsig.col(i);
	}

	// Measurement covariance matrix
	MatrixXd S = MatrixXd(n_z_lidar, n_z_lidar);
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig.col(i) - z_pred;
		// Angle normalization
		while (z_diff(1) > M_PI) z_diff(1)-=2. * M_PI;
		while (z_diff(1) < -M_PI) z_diff(1)+=2. * M_PI;

		S = S + weights(i) * z_diff * z_diff.transpose();
	}

	// Add measurement noise covariance matrix
	S = S + R_laser_;

	*z_pred_out = z_pred;
	*Zsig_out = Zsig;
	*S_out = S;
}

void UKF::UpdateState(const VectorXd &z, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, int src_dimension) {
	// Cross correlation matrix Tc
	MatrixXd Tc = MatrixXd(n_x, src_dimension);

	// Calculate cross correlation matrix
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd z_diff = Zsig.col(i) - z_pred; //residual
		
		// Angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		// State difference
		VectorXd x_diff = Xsig_pred.col(i) - x_;
		
		// Angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		Tc = Tc + weights(i) * x_diff * z_diff.transpose(); // Buggy line in UpdateStae for Lidar data
	}

	// Kalman gain K;
	MatrixXd K = Tc * S.inverse();

	// Residual
	VectorXd z_diff = z - z_pred;

	// Angle normalization
	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	// Update state mean and covariance matrix
	x_ = x_ + K * z_diff;
	P_ = P_ - K * S * K.transpose();
}
