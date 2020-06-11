#ifndef _G2_HPP_
#define _G2_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G2 : public Propagator {
public:
    G2() { num_total_beads = 2; }
    double     operator()(const arma::cube &conf) override;
    double     operator()(const arma::cube &conf, unsigned index) override;
    arma::cube derivative(const arma::cube &conf) override;
};

REGISTER_TYPE_GENERAL(G2, Propagator)

#endif
