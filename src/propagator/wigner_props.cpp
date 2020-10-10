#include "WignerG2.hpp"
#include "WignerKE.hpp"

void WignerG2::set_params(double Tau, pt::ptree::value_type p, arma::vec mass_, double beta_) {
    tau   = Tau / (num_total_beads - 1.);
    nx    = p.second.get<unsigned>("<xmlattr>.nx", 0);
    gamma = p.second.get<unsigned>("<xmlattr>.gamma", 1);
    mass  = mass_;
}

double WignerG2::operator()(const arma::cube &conf) {
    double total_amplitude = -tau / 2. *
                             ((*pot)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                              (*pot)(conf(arma::span::all, arma::span(1), arma::span::all)));
    arma::mat x0  = conf(arma::span::all, arma::span(0), arma::span::all);
    arma::mat xf  = conf(arma::span::all, arma::span(1), arma::span::all);
    double    var = 1;
    for (unsigned atom = 0; atom < x0.n_rows; atom++) {
        arma::vec pos0  = x0.row(atom);
        arma::vec posf  = xf.row(atom);
        auto      disp  = bc->wrap_vector(pos0 - posf);
        double    dist2 = arma::accu(disp % disp);
        var *= utilities::exp_series(dist2 * mass(atom) / (2. * get_tau(atom)), nx);
    }
    return total_amplitude * var;
}

void WignerKE::set_params(double Tau, pt::ptree::value_type p, arma::vec mass_, double beta_) {
    np   = p.second.get<unsigned>("<xmlattr>.np", 0);
    tau  = p.second.get<double>("<xmlattr>.dt_frac") * beta_ / (num_total_beads - 1.);
    beta = p.second.get<double>("<xmlattr>.alpha");
    mass = mass_;
}

double WignerKE::operator()(const arma::cube &momentum) {
    double    ans = 1;
    arma::mat p   = momentum.slice(0);
    for (unsigned atom = 0; atom < p.n_rows; atom++) {
        double ke  = 0.5 * arma::dot(p.row(atom), p.row(atom)) / mass(atom);
        double var = (beta - tau) * ke;
        ans *= utilities::exp_series(var, np);
    }
    return ans;
}
