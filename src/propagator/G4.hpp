#ifndef _G4_HPP_
#define _G4_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G4 : public Propagator {
public:
    G4() { num_total_beads = 3; }
    double     operator()(const arma::cube &conf) override;
    arma::cube derivative(const arma::cube &conf) override;
};

REGISTER_TYPE_GENERAL(G4, Propagator)

#endif
