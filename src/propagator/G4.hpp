#ifndef _G4_HPP_
#define _G4_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G4 : public Propagator {
public:
    G4() { num_total_beads = 3; }
    double operator()(const arma::cube &conf) {
        double total_amplitude = 1;
        total_amplitude *=
            4./3. * std::exp(-tau/2. * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + 2. * (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(2), arma::span::all))))
            - 1./3. * std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(2), arma::span::all))));
        return total_amplitude;
    }
};

REGISTER_TYPE_GENERAL(G4, Propagator)

#endif
