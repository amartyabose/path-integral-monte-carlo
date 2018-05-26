#ifndef _HARMONIC_HPP_
#define _HARMONIC_HPP_

#include "potential.hpp"
#include "../utilities.hpp"

class HarmonicOscillator : public Potential {
    double omega;
public:
    void setup(pt::ptree node) {
        omega = node.get<double>("frequency");
    }

    double operator() (arma::mat x) {
        return 0.5 * omega * omega * arma::accu(x%x);
    }
};

REGISTER_TYPE_GENERAL(HarmonicOscillator, Potential)

#endif
