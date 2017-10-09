#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    previous_timestamp_ = 0;

    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    R_laser_ << 0.0225, 0,
            0, 0.0225;

    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    /**
    TODO:
      * Finish initializing the FusionEKF.
      * Set the process and measurement noises
    */
    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;
    ekf_.Q_ = MatrixXd(4, 4);
    noise_px_ = 9;
    noise_py_ = 9;
    noise_vx_ = 100;
    noise_vy_ = 100;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        cout << "EKF:" << endl;
        // initialize timestamp
        previous_timestamp_ = measurement_pack.timestamp_;
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            const float px = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1));
            const float py = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1));
            const float vx = measurement_pack.raw_measurements_(2) * cos(measurement_pack.raw_measurements_(1));
            const float vy = measurement_pack.raw_measurements_(2) * sin(measurement_pack.raw_measurements_(1));

            ekf_.x_ << px, py, vx, vy;
        } else {
            ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
        }
        ekf_.P_ << 1000, 0, 0, 0,
                0, 1000, 0, 0,
                0, 0, 1000, 0,
                0, 0, 0, 1000;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/


    /**
     TODO:
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
     */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt in seconds

    if (dt > 0.0001) {
        previous_timestamp_ = measurement_pack.timestamp_;
        ekf_.F_(0, 2) = dt;
        ekf_.F_(1, 3) = dt;
        float c1 = dt * dt;
        float c2 = c1 * dt / 2;
        float c3 = c1 * c1 / 4;
        ekf_.Q_ << c3 * noise_px_, 0, c2 * noise_vx_, 0,
                0, c3 * noise_py_, 0, c2 * noise_vy_,
                c2 * noise_px_, 0, c1 * noise_vx_, 0,
                0, c2 * noise_py_, 0, c1 * noise_vy_;
        ekf_.Predict();
    }

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.R_ = R_laser_;
        ekf_.H_ = H_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
