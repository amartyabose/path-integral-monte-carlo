#include "morse.hpp"

void Morse::setup(pt::ptree node) {
    rmin          = node.get<double>("rmin");
    alpha         = node.get<double>("alpha");
    D             = node.get<double>("D", 1);
    cutoff_radius = node.get<double>("cutoff_radius", -1.);
}

double Morse::operator()(arma::mat const &x) {
    double pe = 0;
    for (unsigned atom = 0; atom < x.n_rows; atom++) {
        double exp_dist = (1. - std::exp(-alpha * (x(0, 0) - rmin)));
        pe += D * exp_dist * exp_dist;
    }
    return pe;
}

std::complex<double> Morse::operator()(arma::cx_mat const &x) {
    std::complex<double> pe = 0;
    for (unsigned atom = 0; atom < x.n_rows; atom++) {
        std::complex<double> exp_dist = (1. - std::exp(-alpha * (x(0, 0) - rmin)));
        pe += D * exp_dist * exp_dist;
    }
    return pe;
}