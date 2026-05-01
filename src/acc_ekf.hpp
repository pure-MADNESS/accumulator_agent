#ifndef __ACC_EKF_HPP__
#define __ACC_EKF_HPP__

#include "ekf.hpp"
#include <algorithm>

using namespace Eigen;

class AccEKF : public EKF {
public:
    /**
     * state x: [SoC (0-1), P_bias (W)]^T
     * measure z: [SoC_fmu]^T
     */
    AccEKF(double max_capacity_wh, double efficiency);

    void set_input(double p_net);

    // SoC_dot = (P_net + P_bias) * eta / Capacity
    VectorXd f(const VectorXd& x, double dt) override;
    MatrixXd F(const VectorXd& x, double dt) override;

    VectorXd h(const VectorXd& x) override;
    MatrixXd H(const VectorXd& x) override;

private:
    double _cap_wh;
    double _eta;
    double _p_net;
};

#endif