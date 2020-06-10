#ifndef _HARMONIC_HPP_
#define _HARMONIC_HPP_

#include "../utilities.hpp"
#include "potential.hpp"

class HarmonicOscillator : public Potential {
    double omega;

public:
    void setup(pt::ptree node) override;

    double               operator()(arma::mat const &x, unsigned index) override;
    std::complex<double> operator()(arma::cx_mat const &x) override;
};

REGISTER_TYPE_GENERAL(HarmonicOscillator, Potential)

#endif
