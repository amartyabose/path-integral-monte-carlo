#ifndef _G2_HPP_
#define _G2_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G2 : public Propagator {
public:
    G2() { num_total_beads = 2; }
    double operator()(const arma::cube &conf) {
        double total_amplitude = - tau/2. * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(1), arma::span::all)));
        return std::exp(total_amplitude);
    }
};

REGISTER_TYPE_GENERAL(G2, Propagator)

#endif
