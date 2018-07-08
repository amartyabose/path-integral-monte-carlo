#ifndef _LENNARDJONES_HPP_
#define _LENNARDJONES_HPP_

#include "potential.hpp"
#include "../utilities.hpp"

class LennardJones : public Potential {
    double rmin, epsilon, cutoff_radius;
public:
    void setup(pt::ptree node) {
        rmin = node.get<double>("rmin");
        epsilon = node.get<double>("epsilon");
        cutoff_radius = node.get<double>("cutoff_radius", -1.);
    }

    double operator()(arma::mat x) {
        double pe = 0;
        for(unsigned i=0; i<x.n_rows; i++)
            for(unsigned j=0; j<i; j++) {
                double dist = arma::norm(x.row(i) - x.row(j), 2);
                double rmin_over_dist_6 = std::pow(rmin/dist, 6.);
                pe += epsilon * rmin_over_dist_6 * (rmin_over_dist_6 - 2.);
            }

        if(cutoff_radius>0) {
            arma::mat com = arma::mean(x, 0);
            for(unsigned i=0; i<x.n_rows; i++) {
                double dist = arma::norm(x.row(i) - com.row(0), 2);
                pe += epsilon * std::pow(dist/cutoff_radius, 20.);
            }
        }
        return pe;
    }
};

REGISTER_TYPE_GENERAL(LennardJones, Potential)

#endif
