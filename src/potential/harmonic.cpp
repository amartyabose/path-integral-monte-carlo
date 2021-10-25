#include "harmonic.hpp"

void HarmonicOscillator::setup(pt::ptree node) { omega = node.get<double>("frequency"); }

double HarmonicOscillator::operator()(arma::mat const &x, unsigned atom) const {
    double pe = 0.5 * omega * omega * arma::accu(x.row(atom) % x.row(atom));
    return pe;
}

std::complex<double> HarmonicOscillator::operator()(arma::cx_mat const &x) const {
    std::complex<double> pe = 0.0;
    for (unsigned atom = 0; atom < x.n_rows; atom++)
        pe += 0.5 * omega * omega * arma::accu(x.row(atom) % x.row(atom));
    return pe;
}
