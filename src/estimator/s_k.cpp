#include <cmath>
#include <complex>
#include <string>

#include "../boundary_conditions/boundary_conditions.hpp"
#include "s_k.hpp"

std::complex<double> I(0, 1);

void Sk::setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) {
    type        = type_;
    double Rmax = node.get<double>("<xmlattr>.rmax");
    kmax        = node.get<double>("<xmlattr>.kmax");
    delta_k     = 2 * arma::datum::pi / Rmax;
    ks_1d       = arma::regspace(delta_k, delta_k, kmax);
    n_cols      = 2;
    n_rows      = ks_1d.n_rows;

    isotropic = true;

    boost::optional<std::string> orientation = node.get_optional<std::string>("<xmlattr>.axis");
    if (orientation) {
        isotropic = false;
        axis      = arma::normalise(arma::vec(orientation.get()));
    }
}

double Sk::calc_sk(arma::mat const &pos, arma::vec const &k) {
    std::complex<double> amplitude = 0.0;
    for (unsigned atom = 0; atom < pos.n_rows; atom++) {
        arma::vec atomic_pos = pos.row(atom);
        amplitude += std::exp(-I * arma::dot(k, atomic_pos));
    }
    double ampl_norm = std::abs(amplitude);
    return ampl_norm * ampl_norm / pos.n_rows;
}

double Sk::calc_sk(std::shared_ptr<Configuration> const &x, arma::vec const &ks) {
    double val = 0.0;
    if (type == "pimc" || type == "pimd") {
        for (unsigned i = 0; i < x->num_beads(); ++i)
            val += calc_sk(x->time_slice(i), ks);
        val /= x->num_beads();
    } else if (type == "pimc" || type == "pimd")
        val = calc_sk(x->pos(), ks);

    return val;
}

arma::mat Sk::eval(std::shared_ptr<Configuration> const &x) {
    arma::mat ans = arma::zeros<arma::mat>(n_rows, n_cols);
    ans.col(0)    = ks_1d;
    if (x->num_dims() == 1) {
        for (unsigned ki = 0; ki < ks_1d.n_rows; ++ki)
            ans(ki, 1) = calc_sk(x, arma::ones<arma::vec>(1) * ks_1d(ki));
    } else if (x->num_dims() == 2) {
        arma::vec nums = arma::zeros<arma::vec>(n_rows);
        for (unsigned kx = 0; kx < n_rows; ++kx)
            for (unsigned ky = 0; ky < n_rows; ++ky) {
                unsigned len = std::sqrt(kx * kx + ky * ky);
                if (len * len == kx * kx + ky * ky && len > 0 && len <= n_rows) {
                    nums[len - 1]++;
                    arma::vec k = arma::zeros<arma::vec>(2);
                    k(0)        = kx;
                    k(1)        = ky;
                    ans(len - 1, 1) += calc_sk(x, k * delta_k);
                }
            }
        ans.col(1) /= nums;
    } else if (x->num_dims() == 3) {
        arma::vec nums = arma::zeros<arma::vec>(n_rows);
        for (unsigned kx = 0; kx < n_rows; ++kx)
            for (unsigned ky = 0; ky < n_rows; ++ky)
                for (unsigned kz = 0; kz < n_rows; ++kz) {
                    unsigned len = std::sqrt(kx * kx + ky * ky + kz * kz);
                    if (len * len == kx * kx + ky * ky + kz * kz && len > 0 && len <= n_rows) {
                        nums[len - 1]++;
                        arma::vec k = arma::zeros<arma::vec>(3);
                        k(0)        = kx;
                        k(1)        = ky;
                        k(2)        = kz;
                        ans(len - 1, 1) += calc_sk(x, k * delta_k);
                    }
                }
        ans.col(1) /= nums;
    }
    return ans;
}
