#include "brownian_bridge.hpp"
#include "../units.hpp"

double BrownianBridge::get_temp(unsigned bead_num, unsigned atom) {
    double   temp          = 0;
    unsigned beads_counted = 0;
    for (unsigned p = 0; beads_counted < bead_num && p < propagator.size(); p++)
        for (unsigned b = 0; beads_counted < bead_num && b < propagator[p]->total_beads() - 1; b++) {
            beads_counted++;
            temp += propagator[p]->get_tau(atom);
        }
    return temp;
}

void BrownianBridge::setup(pt::ptree::value_type node, double beta_, arma::vec mass_) {
    num_beads_moved = node.second.get<unsigned>("length");
    num_attempts    = node.second.get<unsigned>("attempts");
    lambda          = units.hbar * units.hbar / (2. * mass_);
}

void BrownianBridge::set_beta(arma::vec beta_) { beta = beta_; }

void BrownianBridge::operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) {
    for (unsigned atom_ind = 0; atom_ind < atom_nums.n_rows; atom_ind++) {
        unsigned nbeads = conf->num_augmented_beads();
        unsigned atom   = atom_nums(atom_ind);
        for (unsigned attempt = 0; attempt < num_attempts; attempt++) {
            unsigned start = random_integer(0, nbeads - 1);
            unsigned end   = (start + num_beads_moved + 1) % nbeads;

            std::shared_ptr<Configuration> conf_new(conf->duplicate());
            for (unsigned b = (start + 1) % nbeads, nmoved = 0; nmoved < num_beads_moved;
                 b = (b + 1) % nbeads, nmoved++) {
                double tau0         = get_temp((b + nbeads - 1) % nbeads, atom);
                double tau1         = get_temp(b, atom);
                double taue         = get_temp(end, atom);
                double start_to_new = tau1 - tau0;
                double new_to_end   = taue - tau1;
                if (new_to_end < 0)
                    new_to_end = taue + (beta(atom) - tau1);
                if (start_to_new < 0)
                    start_to_new = tau1 + (beta(atom) - tau0);
                start_to_new *= 2. * lambda(atom);
                new_to_end *= 2. * lambda(atom);

                double start_to_end = start_to_new + new_to_end;
                auto   sigma_bb     = std::sqrt(start_to_new * new_to_end / start_to_end);

                auto start_bead = conf_new->augmented_bead_position(atom, (nbeads + b - 1) % nbeads);
                auto end_bead   = conf_new->augmented_bead_position(atom, end);
                auto disp       = bc->wrap_vector(end_bead - start_bead);

                arma::vec mean_loc = start_bead + start_to_new * disp / start_to_end;

                arma::vec new_pos = arma::zeros<arma::vec>(conf->num_dims());
                for (unsigned d = 0; d < conf->num_dims(); d++)
                    new_pos(d) = random_normal(mean_loc(d), sigma_bb);

                conf_new->augmented_set(atom, b, bc->wrap_vector(new_pos));
            }
            check_amplitude(conf, conf_new, atom);
        }
    }
}
