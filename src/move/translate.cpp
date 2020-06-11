#include "translate.hpp"

void Translate::setup(pt::ptree::value_type node, double beta_, arma::vec mass_) {
    max_step = node.second.get<double>("length");
    // arma::ones<arma::mat>(params.get<unsigned>("num_atoms"), params.get<unsigned>("num_dimensions"));
}

void Translate::operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) {
    unsigned  dims = conf->num_dims();
    arma::vec step(dims);
    for (unsigned atom = 0; atom < atom_nums.n_rows; atom++) {
        for (unsigned i = 0; i < dims; i++)
            step(i) = random_float(-max_step, max_step);

        std::shared_ptr<Configuration> conf_new(conf->duplicate());
        conf_new->shift(atom_nums(atom), step);

        check_amplitude(conf, conf_new, atom);
    }
}
