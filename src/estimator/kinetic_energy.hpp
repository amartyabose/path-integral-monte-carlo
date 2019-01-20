#ifndef _KE_H_
#define _KE_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class KineticEnergy : public Estimator {
public:
    double beta;
    int num_prop;
    void setup(pt::ptree params, unsigned nblocks=10) {
        beta = params.get<double>("beta");
        num_prop = params.get<int>("num_propagators");
        values_real_plus = values_real_minus = values_imag_plus = values_imag_minus = arma::zeros<arma::vec>(nblocks);
        i_real_plus = i_real_minus = i_imag_plus = i_imag_minus = arma::zeros<arma::vec>(nblocks);
    }
    double eval(boost::shared_ptr<Configuration> x) {
        double val = 0;
        double tau = beta/x->num_beads();
        for (unsigned i=0; i<x->num_beads(); i++)
            val += arma::accu((x->time_slice(i) - x->time_slice((i+1)%x->num_beads())) % (x->time_slice(i) - x->time_slice((i+1)%x->num_beads())));
        return x->num_dims()*x->num_atoms()/(2.*tau) - val / (2 * x->num_beads() * tau * tau);
    }
};

REGISTER_TYPE_GENERAL(KineticEnergy, Estimator)

#endif
