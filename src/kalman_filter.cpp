#include "kalman_filter.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Vector4d;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // Noise average is zero
  VectorXd u = VectorXd::Zero(x_.size());

  x_ = F_*x_ + u;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::UpdateKF(const VectorXd &z) {
  VectorXd y = z - H_*x_;
	Update(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	VectorXd y = z - CartesianToPolar(x_);
	Update(y);
}

void KalmanFilter::Update(const Eigen::VectorXd &y) {
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_*P_*Ht + R_;
	MatrixXd K = P_*Ht*S.inverse();
	MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

	// Update the state and covariance matrices.
	x_ = x_ + K*y;
	P_ = (I - K*H_)*P_;
}

Vector3d KalmanFilter::CartesianToPolar(const Vector4d &x){
	double px = x(0);
	double py = x(1);
	double vx = x(2);
	double vy = x(3);

	double rho = sqrt(px*px + py*py);
	double phi = atan2(py, px);
	double rho_dot = (px*vx + py*vy)/rho;

	Vector3d h_x {rho, phi, rho_dot};
	return h_x;
}
