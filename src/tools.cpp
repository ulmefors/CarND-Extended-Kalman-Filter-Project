#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

/**
 * @param estimations Estimated cartesian position, velocity
 * @param ground_truth True cartesian position, velocity
 * @return Root Mean Squared Error of estimations compared to ground truth
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    VectorXd residual(4);
    VectorXd sum_sq_error(4);
    sum_sq_error << 0,0,0,0;

    long nb_meas = estimations.size();
    long nb_truth = ground_truth.size();

    // Verify that number of measurements matches number of ground truth entries
    if (nb_meas == 0 || nb_meas != nb_truth) {
        std::cout << "Error! Tools::CalculateRMSE() - Invalid Measurement/Ground Truth length" << std::endl;
        return rmse;
    }

    // Sum residuals and calculate mean squared error
    for (int i = 0; i < nb_meas; i++) {
        residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        sum_sq_error = sum_sq_error + residual;
    }
    rmse = (sum_sq_error/nb_meas).array().sqrt();

    return rmse;
}

/**
 * @param x_state Cartesian position, velocity
 * @return Jacobian matrix. Element (0,0) set to NaN if division by zero.
 */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state){

    float RAD_THRESHOLD = 0.00001;

    // Jacobian Matrix
    MatrixXd Hj(3,4);

    // Cartesian position and velocities
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // Radius^2, radius, and radius^3
    float rho2 = std::pow(px, 2) + std::pow(py, 2);
    float rho1 = std::pow(rho2, 0.5);
    float rho3 = std::pow(rho1, 3);

    // Verify that radius is not (close to) zero
    if (rho1 < RAD_THRESHOLD) {
			std::cout << "Error! Tools::CalculateJacobian() - Division by zero" << std::endl;
			Hj.setZero();
			Hj(0,0) = nan(0);
    } else {
      Hj<<px/rho1,    py/rho1,    0,  0,
          -py/rho2,   px/rho2,    0,  0,
          py*(vx*py-vy*px)/rho3, px*(vy*px-vx*py)/rho3, px/rho1,  py/rho1;
    }

    return Hj;
}
