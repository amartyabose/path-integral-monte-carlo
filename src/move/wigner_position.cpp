#include "wigner_position.hpp"

void WignerPosition::setup(pt::ptree::value_type node, double beta_, arma::vec mass_) {
    delta_beta = node.second.get<double>("<xmlattr>.dt_frac") * beta_;
    mass       = mass_;
    sigma      = arma::sqrt(delta_beta / (4. * mass));
}

void WignerPosition::operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) const {
    arma::mat            xmean       = (conf->time_slice(0) + conf->time_slice(conf->num_beads(0) - 1)) / 2.;
    unsigned             dimensions  = conf->num_dims();
    WignerConfiguration *wigner_conf = dynamic_cast<WignerConfiguration *>(conf.get());
    for (unsigned atom = 0; atom < atom_nums.n_rows; atom++) {
        for (unsigned dim = 0; dim < dimensions; dim++)
            wigner_conf->set_position(atom_nums(atom), dim, random_normal(xmean(atom_nums(atom), dim), sigma(atom)));
        moves_tried++;
        moves_accepted++;
    }
}
