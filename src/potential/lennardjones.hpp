#ifndef _LENNARDJONES_HPP_
#define _LENNARDJONES_HPP_

#include "potential.hpp"
#include "../utilities.hpp"

class LennardJones : public Potential {
    double rmin, epsilon;
public:
    void setup(pt::ptree node) {
        rmin = node.get<double>("rmin");
        epsilon = node.get<double>("epsilon");
    }

    double operator()(arma::mat x) {
        double pe = 0;
        for(unsigned i=0; i<x.n_rows; i++)
            for(unsigned j=0; j<i; j++) {
                double dist = arma::norm(x.row(i) - x.row(j), 2);
                double rmin_over_dist_6 = std::pow(rmin/dist, 6.);
                pe += epsilon * rmin_over_dist_6 * (rmin_over_dist_6 - 2.);
            }
        return pe;
    }
};

REGISTER_TYPE_GENERAL(LennardJones, Potential)

#endif
