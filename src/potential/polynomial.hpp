#ifndef _POLYNOMIAL_HPP_
#define _POLYNOMIAL_HPP_

#include <vector>

#include "../utilities.hpp"
#include "potential.hpp"

class PolynomialPotential : public Potential {
    std::vector<double> coeffs;

public:
    void setup(pt::ptree node) override;

    double               operator()(arma::mat const &x) override;
    std::complex<double> operator()(arma::cx_mat const &x) override;
};

REGISTER_TYPE_GENERAL(PolynomialPotential, Potential)

#endif
