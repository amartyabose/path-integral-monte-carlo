#ifndef _LENNARDJONES_HPP_
#define _LENNARDJONES_HPP_

#include "../utilities.hpp"
#include "potential.hpp"

class LennardJones : public Potential {
    double rmin, epsilon, cutoff_radius;

public:
    void setup(pt::ptree node) override;

    double               operator()(arma::mat const &x, unsigned index) const override;
    std::complex<double> operator()(arma::cx_mat const &x) const override { return 0; }

    arma::mat derivative(arma::mat const &x) const override;
};

REGISTER_TYPE_GENERAL(LennardJones, Potential)

#endif
