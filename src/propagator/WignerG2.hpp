#ifndef _WIGNER_G2_HPP
#define _WIGNER_G2_HPP

#include <armadillo>

#include "../utilities.hpp"
#include "propagator.hpp"

class WignerG2 : public Propagator {
    double    gamma;
    unsigned  nx;
    arma::vec mass;

public:
    WignerG2() {
        num_total_beads = 2;
        mass            = 1;
    }
    double get_tau(unsigned atom) override { return (2. * mass(atom) + tau * gamma) / gamma; }
    void   set_params(std::shared_ptr<Potential> pot, double Tau, pt::ptree::value_type p, arma::vec mass,
                      double beta) override;
    double operator()(const arma::cube &conf) override;
};

REGISTER_TYPE_GENERAL(WignerG2, Propagator)

#endif
