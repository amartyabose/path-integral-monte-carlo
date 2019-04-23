#ifndef _MORSE_HPP_
#define _MORSE_HPP_

#include "potential.hpp"
#include "../utilities.hpp"

class Morse : public Potential {
    double rmin, alpha, cutoff_radius, D;
public:
    void setup(pt::ptree node) {
        rmin = node.get<double>("rmin");
        alpha = node.get<double>("alpha");
        D = node.get<double>("D", 1);
        cutoff_radius = node.get<double>("cutoff_radius", -1.);
    }

    double operator()(arma::mat x) {
        double pe = 0;
        double exp_dist = (1 - std::exp(-alpha * (x(0,0)-rmin)));
        pe += D * exp_dist * exp_dist;

        return pe;
    }
};

REGISTER_TYPE_GENERAL(Morse, Potential)

#endif
