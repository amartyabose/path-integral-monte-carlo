#include "lennardjones.hpp"

void LennardJones::setup(pt::ptree node) {
    rmin          = node.get<double>("rmin");
    epsilon       = node.get<double>("epsilon");
    cutoff_radius = node.get<double>("cutoff_radius", -1.);
}

double LennardJones::operator()(arma::mat const &x) {
    double pe = 0;
    for (unsigned i = 0; i < x.n_rows; i++)
        for (unsigned j = 0; j < i; j++) {
            arma::vec disp = x.row(i) - x.row(j);

            double dist             = arma::accu(disp % disp);
            double rmin_over_dist_6 = std::pow(rmin / dist, 6.);
            pe += epsilon * rmin_over_dist_6 * (rmin_over_dist_6 - 2.);
        }

    if (cutoff_radius > 0) {
        arma::mat com = arma::mean(x, 0);
        for (unsigned i = 0; i < x.n_rows; i++) {
            arma::vec disp = x.row(i) - com.row(0);

            double dist = arma::accu(disp % disp);
            pe += epsilon * std::pow(dist / cutoff_radius, 20.);
        }
    }
    return pe;
}

std::complex<double> LennardJones::operator()(arma::cx_mat const &x) {
    std::complex<double> pe = 0;
    for (unsigned i = 0; i < x.n_rows; i++)
        for (unsigned j = 0; j < i; j++) {
            arma::cx_vec disp = x.row(i) - x.row(j);

            std::complex<double> dist             = arma::accu(disp % disp);
            std::complex<double> rmin_over_dist_6 = std::pow(rmin / dist, 6.);
            pe += epsilon * rmin_over_dist_6 * (rmin_over_dist_6 - 2.);
        }

    if (cutoff_radius > 0) {
        arma::cx_mat com = arma::mean(x, 0);
        for (unsigned i = 0; i < x.n_rows; i++) {
            arma::cx_vec disp = x.row(i) - com.row(0);

            std::complex<double> dist = arma::accu(disp % disp);
            pe += epsilon * std::pow(dist / cutoff_radius, 20.);
        }
    }
    return pe;
}
