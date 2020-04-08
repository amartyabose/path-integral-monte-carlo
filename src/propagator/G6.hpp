#ifndef _G6_HPP_
#define _G6_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G6 : public Propagator {
public:
    G6() { num_total_beads = 5; }
    double     operator()(const arma::cube &conf) override;
    arma::cube derivative(const arma::cube &conf) override;
};

REGISTER_TYPE_GENERAL(G6, Propagator)

#endif
