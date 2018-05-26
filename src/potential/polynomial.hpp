#ifndef _POLYNOMIAL_HPP_
#define _POLYNOMIAL_HPP_

#include <string>
#include <vector>
using std::string;
using std::vector;

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "potential.hpp"
#include "../utilities.hpp"

class PolynomialPotential : public Potential {
    vector<double> coeffs;
public:
    void setup(pt::ptree node) {
        string coeffs_string = node.get<string>("coeffs");
        vector<string> coeff_vec;
        boost::split(coeff_vec, coeffs_string, boost::is_any_of("\t "));
        for(unsigned i=0; i<coeff_vec.size(); i++)
            coeffs.push_back(boost::lexical_cast<double>(coeff_vec[i]));
    }

    double operator()(arma::mat x) {
        double pe = 0;
        arma::mat temp = arma::ones<arma::mat>(arma::size(x));
        for(unsigned i=0; i<coeffs.size(); i++) {
            pe += coeffs[i] * arma::accu(temp);
            temp = temp % x;
        }
        return pe;
    }
};

REGISTER_TYPE_GENERAL(PolynomialPotential, Potential)

#endif
