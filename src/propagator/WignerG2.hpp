#ifndef _WIGNER_G2_HPP
#define _WIGNER_G2_HPP

#include <boost/shared_ptr.hpp>
#include "../utilities.hpp"
#include "propagator.hpp"

class WignerG2 : public Propagator {
    double mass, gamma;
    unsigned nx;
public:
    WignerG2() { num_total_beads = 2; mass = 1;}
    virtual void set_params(boost::shared_ptr<Potential> pot, double Tau, pt::ptree::value_type p, pt::ptree params) {
        V   = pot;
        tau = Tau / (num_total_beads - 1.);
        nx = p.second.get<unsigned>("<xmlattr>.nx", 0);
        gamma = p.second.get<unsigned>("<xmlattr>.gamma", 1);
        mass = params.get<double>("mass", 1);
    }
    virtual double get_tau() {
        return (2.*mass + tau*gamma)/gamma;
    }
    double operator()(const arma::cube &conf) {
        double total_amplitude = - tau/2. * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(1), arma::span::all)));
        arma::mat x0 = conf(arma::span::all, arma::span(0), arma::span::all);
        arma::mat xf = conf(arma::span::all, arma::span(1), arma::span::all);
        double dist2 = arma::accu((x0-xf)%(x0-xf));
        double var = dist2 * mass / (2. * get_tau());
        double poly_expansion = 1;
        double term = var;
        for(unsigned long i=1; i<nx; i++) {
            poly_expansion += term;
            term *= var/(i+1);
        }

        return std::exp(total_amplitude) * poly_expansion;
    }
};

REGISTER_TYPE_GENERAL(WignerG2, Propagator)

#endif
