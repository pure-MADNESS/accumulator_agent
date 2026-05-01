#include "acc_ekf.hpp"

AccEKF::AccEKF(double max_capacity_wh, double efficiency) 
    : EKF(2, 1), _cap_wh(max_capacity_wh), _eta(efficiency) {
    
    x << 0.5, 0.0;

    Q.resize(2, 2);
    Q << 1e-3, 0,    // slow soc
         0,    0.1;   // varying bias
    R.resize(1, 1);
    R << 0.01; 
}

void AccEKF::set_input(double p_net) {
    _p_net = p_net;
}

VectorXd AccEKF::f(const VectorXd& x, double dt) {
    VectorXd x_new(2);
    double soc = x(0);
    double p_bias = x(1);


    double p_eff = (_p_net + p_bias);
    double power_flow = (p_eff > 0) ? (p_eff * _eta) : (p_eff / _eta);
    
    x_new(0) = soc + (power_flow * (dt / 3600.0)) / _cap_wh;
    x_new(1) = p_bias; // bias modeled as a constant

    x_new(0) = std::clamp(x_new(0), 0.0, 1.0);
    
    return x_new;
}

MatrixXd AccEKF::F(const VectorXd& x, double dt) {
    MatrixXd Fj(2, 2);
    
    // d(f1)/d(soc) = 1
    // d(f1)/d(pbias) = (dt / 3600) * eta_term / _cap_wh
    double eta_term = ((_p_net + x(1)) > 0) ? _eta : (1.0 / _eta);
    double df1_dpb = (dt / 3600.0) * eta_term / _cap_wh;

    Fj << 1.0, df1_dpb,
          0.0, 1.0;
    
    return Fj;
}

VectorXd AccEKF::h(const VectorXd& x) {
    VectorXd z_pred(1);
    z_pred(0) = x(0); // the measure is directly the SOC
    return z_pred;
}

MatrixXd AccEKF::H(const VectorXd& x) {
    MatrixXd Hj(1, 2);
    Hj << 1.0, 0.0;
    return Hj;
}