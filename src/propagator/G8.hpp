#ifndef _G8_HPP_
#define _G8_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G8 : public Propagator {
public:
    G8() { num_total_beads = 7; }
    double operator()(const arma::cube &conf) {
        double total_amplitude = 1;
        total_amplitude *=
            54./35. * std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all))/2. + (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(2), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(3), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(4), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(5), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(6), arma::span::all))/2.))
            - 27./40. * std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + 2.*(*V)(conf(arma::span::all, arma::span(2), arma::span::all)) + 2.*(*V)(conf(arma::span::all, arma::span(4), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(6), arma::span::all))))
            + 2./15. * std::exp(-tau * (1.5*(*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + 3*(*V)(conf(arma::span::all, arma::span(3), arma::span::all)) + 1.5*(*V)(conf(arma::span::all, arma::span(6), arma::span::all))))
            - 1./840. * std::exp(-3 * tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) + (*V)(conf(arma::span::all, arma::span(6), arma::span::all))));
        return total_amplitude;
    }
};

REGISTER_TYPE_GENERAL(G8, Propagator)

#endif
