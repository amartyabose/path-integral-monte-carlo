#ifndef _G6_HPP_
#define _G6_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G6 : public Propagator {
public:
    G6() { num_total_beads = 5; }
    double operator()(const arma::cube &conf) {
        double total_amplitude = 1;
        total_amplitude *=
            64./45. * std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all))/2. + (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(2), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(3), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(4), arma::span::all))/2.))
            - 4./9. * std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + 2.*(*V)(conf(arma::span::all, arma::span(2), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(4), arma::span::all))))
            + 1./45. * std::exp(-2*tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(4), arma::span::all))));
        return total_amplitude;
    }
};

REGISTER_TYPE_GENERAL(G6, Propagator)

#endif
