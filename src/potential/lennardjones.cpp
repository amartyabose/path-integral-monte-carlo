#include "lennardjones.hpp"
#include "../units.hpp"

void LennardJones::setup(pt::ptree node) {
    epsilon       = node.get<double>("epsilon") / units.energy_to_non_base_units;
    cutoff_radius = node.get<double>("cutoff_radius", -1.);

    auto rmin_  = node.get_optional<double>("rmin");
    auto sigma_ = node.get_optional<double>("sigma");
    if (!rmin_ && !sigma_)
        throw std::runtime_error("Neither Rmin nor Sigma provided.");
    else if (rmin_ && sigma_)
        throw std::runtime_error("Both Rmin or Sigma provided.");
    else if (rmin_)
        rmin = rmin_.get();
    else if (sigma_)
        rmin = sigma_.get() * 1.12246204831; // the silly value is 2^(1/6)
}

double LennardJones::operator()(arma::mat const &x, unsigned i) {
    double pe = 0;
    for (unsigned j = 0; j < x.n_rows; j++) {
        if (j == i)
            continue;

        arma::rowvec disp = bc->wrap_vector(x.row(i) - x.row(j));
        double       dist = std::sqrt(arma::accu(disp % disp));

        double rmin_over_dist_6 = std::pow(rmin / dist, 6.);
        pe += epsilon * rmin_over_dist_6 * (rmin_over_dist_6 - 2.);
    }

    if (cutoff_radius > 0) {
        arma::mat    com  = arma::mean(x, 0);
        arma::rowvec disp = bc->wrap_vector(x.row(i) - com.row(0));
        double       dist = std::sqrt(arma::accu(disp % disp));
        pe += epsilon * std::pow(dist / cutoff_radius, 20.);
    }
    return pe;
}

arma::mat LennardJones::derivative(arma::mat const &x) {
    arma::mat grad = arma::zeros<arma::mat>(arma::size(x));
    for (unsigned i = 0; i < x.n_rows; i++)
        for (unsigned j = 0; j < i; j++) {
            arma::rowvec disp = bc->wrap_vector(x.row(i) - x.row(j));
            double       dist = std::sqrt(arma::accu(disp % disp));

            double rmin_over_dist_6 = std::pow(rmin / dist, 6.);

            double deriv_r = -12 * epsilon * (rmin_over_dist_6 - 1.0) * rmin_over_dist_6 / dist;
            grad.row(i) += deriv_r * disp / dist;
            grad.row(j) -= deriv_r * disp / dist;
        }

    if (cutoff_radius > 0) {
        arma::mat com = arma::mean(x, 0);
        for (unsigned i = 0; i < x.n_rows; i++) {
            arma::rowvec disp    = bc->wrap_vector(x.row(i) - com.row(0));
            double       dist    = std::sqrt(arma::accu(disp % disp));
            double       deriv_r = 20.0 * epsilon * std::pow(dist / cutoff_radius, 20) / dist;
            grad.row(i) += deriv_r * disp / dist;
        }
    }
    return grad / 2;
}