#ifndef _MORSE_HPP_
#define _MORSE_HPP_

#include "../utilities.hpp"
#include "potential.hpp"

class Morse : public Potential {
    double rmin, alpha, cutoff_radius, D;

public:
    void setup(pt::ptree node) override;

    double               operator()(arma::mat const &x, unsigned index) const override;
    std::complex<double> operator()(arma::cx_mat const &x) const override;
};

REGISTER_TYPE_GENERAL(Morse, Potential)

#endif
