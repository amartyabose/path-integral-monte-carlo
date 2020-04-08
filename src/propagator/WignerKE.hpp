#ifndef _WIGNER_KE_HPP
#define _WIGNER_KE_HPP

#include <armadillo>

#include "../utilities.hpp"
#include "propagator.hpp"

class WignerKE : public Propagator {
    double    beta;
    unsigned  np;
    arma::vec mass;

public:
    WignerKE() {
        num_total_beads = 2;
        mass            = 1;
    }
    void   set_params(std::shared_ptr<Potential> pot, double Tau, pt::ptree::value_type p, arma::vec mass,
                      double beta) override;
    double operator()(const arma::cube &momentum) override;
};

REGISTER_TYPE_GENERAL(WignerKE, Propagator)

#endif
