#ifndef _WIGNER_KE_HPP
#define _WIGNER_KE_HPP

#include <boost/shared_ptr.hpp>
#include "../utilities.hpp"
#include "propagator.hpp"

class WignerKE : public Propagator {
    double mass, beta;
    unsigned np;
public:
    WignerKE() { num_total_beads = 2; mass = 1;}
    virtual void set_params(boost::shared_ptr<Potential> pot, double Tau, pt::ptree::value_type p, pt::ptree params) {
        V   = pot;
        np = p.second.get<unsigned>("<xmlattr>.np", 0);
        tau = p.second.get<double>("<xmlattr>.dt") / (num_total_beads - 1.);
        //beta = params.get<double>("beta");
        beta = p.second.get<double>("<xmlattr>.alpha");
        mass = params.get<double>("mass", 1);
    }
    double operator()(const arma::cube &momentum) {
        double ke = arma::accu(momentum%momentum)/(2.*mass);
        double var = (beta-tau) * ke;
        return exp_series(var, np);
    }
};

REGISTER_TYPE_GENERAL(WignerKE, Propagator)

#endif
