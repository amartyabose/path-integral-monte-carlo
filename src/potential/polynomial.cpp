#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "polynomial.hpp"

void PolynomialPotential::setup(pt::ptree node) {
    std::string              coeffs_string = node.get<std::string>("coeffs");
    std::vector<std::string> coeff_vec;
    boost::split(coeff_vec, coeffs_string, boost::is_any_of("\t "));
    for (unsigned i = 0; i < coeff_vec.size(); i++)
        coeffs.push_back(boost::lexical_cast<double>(coeff_vec[i]));
}

double PolynomialPotential::operator()(arma::mat const &x, unsigned index) const {
    double pe = coeffs[0];
    for (unsigned d = 0; d < x.n_cols; d++) {
        double dist = x(index, d);
        double temp = dist;
        for (unsigned i = 1; i < coeffs.size(); i++) {
            pe += coeffs[i] * temp;
            temp = temp * dist;
        }
    }
    return pe;
}

std::complex<double> PolynomialPotential::operator()(arma::cx_mat const &x) const {
    std::complex<double> pe = coeffs[0];
    for (unsigned atom = 0; atom < x.n_rows; atom++) {
        for (unsigned d = 0; d < x.n_cols; d++) {
            std::complex<double> dist = x(atom, d);
            std::complex<double> temp = dist;
            for (unsigned i = 1; i < coeffs.size(); i++) {
                pe += coeffs[i] * temp;
                temp = temp * dist;
            }
        }
    }
    return pe;
}
