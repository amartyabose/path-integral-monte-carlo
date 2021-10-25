#include "morse.hpp"
#include "../boundary_conditions/boundary_conditions.hpp"

void Morse::setup(pt::ptree node) {
    rmin          = node.get<double>("rmin");
    alpha         = node.get<double>("alpha");
    D             = node.get<double>("D", 1);
    cutoff_radius = node.get<double>("cutoff_radius", -1.);
}

double Morse::operator()(arma::mat const &x, unsigned index) const {
    arma::vec disp     = bc->wrap_vector(x.row(index));
    double    dist     = arma::norm(disp);
    double    exp_dist = (1. - std::exp(-alpha * (dist - rmin)));
    double    pe       = D * exp_dist * exp_dist;
    return pe;
}

std::complex<double> Morse::operator()(arma::cx_mat const &x) const {
    std::complex<double> pe = 0;
    for (unsigned atom = 0; atom < x.n_rows; atom++) {
        std::complex<double> exp_dist = (1. - std::exp(-alpha * (x(0, 0) - rmin)));
        pe += D * exp_dist * exp_dist;
    }
    return pe;
}
