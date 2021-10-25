#include "wigner_momentum.hpp"

void WignerMomentum::setup(pt::ptree::value_type node, double beta_, arma::vec mass_) {
    mass  = mass_;
    beta  = node.second.get<double>("<xmlattr>.alpha"); // * arma::ones<arma::vec>(mass.n_rows);
    sigma = arma::sqrt(mass / beta);
}

void WignerMomentum::operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) const {
    unsigned             dimensions     = conf->num_dims();
    WignerConfiguration *wigner_newconf = dynamic_cast<WignerConfiguration *>(conf->duplicate());
    for (unsigned atom = 0; atom < atom_nums.n_rows; atom++) {
        for (unsigned dim = 0; dim < dimensions; dim++)
            wigner_newconf->set_momentum(atom_nums(atom), dim, random_normal(0, sigma(atom)));
    }
    check_amplitude(conf, std::shared_ptr<Configuration>(wigner_newconf));
}

void WignerMomentum::check_amplitude(std::shared_ptr<Configuration> &conf_old,
                                     std::shared_ptr<Configuration>  conf_new) const {
    moves_tried++;

    WignerConfiguration *wigner_newconf = dynamic_cast<WignerConfiguration *>(conf_new.get());
    WignerConfiguration *wigner_oldconf = dynamic_cast<WignerConfiguration *>(conf_old.get());
    arma::mat            temp_pold      = wigner_oldconf->get_momentum();
    arma::mat            temp_pnew      = wigner_newconf->get_momentum();
    arma::cube           pold           = arma::cube(temp_pold.n_rows, temp_pold.n_cols, 1);
    pold.slice(0)                       = temp_pold;
    arma::cube pnew                     = arma::cube(temp_pnew.n_rows, temp_pnew.n_cols, 1);
    pnew.slice(0)                       = temp_pnew;

    double new_weight = 1, old_weight = 1;
    for (unsigned p = 0; p < propagator.size(); p++) {
        new_weight *= (*propagator[p])(pnew);
        old_weight *= (*propagator[p])(pold);
    }

    if (new_weight / old_weight > random_float(0, 1)) {
        conf_old.reset(conf_new->duplicate());
        moves_accepted++;
    }
}
