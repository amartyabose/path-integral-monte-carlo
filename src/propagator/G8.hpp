#ifndef _G8_HPP_
#define _G8_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G8 : public Propagator {
public:
    G8() { num_total_beads = 7; }
    double     operator()(const arma::cube &conf) override;
    arma::cube derivative(const arma::cube &conf) override;
};

REGISTER_TYPE_GENERAL(G8, Propagator)

#endif
