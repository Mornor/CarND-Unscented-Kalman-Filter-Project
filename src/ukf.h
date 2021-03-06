#ifndef UKF_H
#define UKF_H
#include "Eigen/Dense"
#include "measurement_package.h"
#include "ground_truth_package.h"
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* Predicted sigma points matrix
  MatrixXd Xsig_pred;

  ///* predicted sigma points matrix
  MatrixXd Xsig_aug;

  // Augmented mean vector
  VectorXd x_aug;

  // Sensor Noise 
  MatrixXd R_laser_;
  MatrixXd R_radar_;
  
  ///* time when the state is true, in us
  long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  double previous_timestamp_;

  ///* State dimension
  int n_x;

  // Set measurement dimension, radar can measure r, phi, and r_dot
  int n_z;

  // Set measurement dimension, lidar
  int n_z_lidar;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  VectorXd weights;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   * @param gt_package The ground truth of the state x at measurement time
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double dt);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void PredictMeanAndCovariance();
  void PredictRadarMeasurement(MatrixXd *Zsig, VectorXd *z_pred, MatrixXd *S);
  void PredictLidarMeasurement(MatrixXd *Zsig, VectorXd *z_pred, MatrixXd *S);
  void UpdateState(const VectorXd &z, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, int src_dimension);

  void AugmentedSigmaPoints();
  void PredictAugmentedSigmaPoints(double delta_t);

  VectorXd InitRadar(const MeasurementPackage &measurement_pack);
  VectorXd InitLaser(const MeasurementPackage &measurement_pack);


};

#endif /* UKF_H */
