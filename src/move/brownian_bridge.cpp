#include "brownian_bridge.hpp"
#include "../units.hpp"

double BrownianBridge::get_temp(unsigned bead_num, unsigned atom) const {
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

void BrownianBridge::operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) const {
    for (unsigned atom_ind = 0; atom_ind < atom_nums.n_rows; atom_ind++) {
        unsigned atom   = atom_nums(atom_ind);
        unsigned nbeads = conf->num_augmented_beads(atom);
        if (conf->type_of_polymers[atom] == 'o')
            nbeads += 1;
        for (unsigned attempt = 0; attempt < num_attempts; attempt++) {
            unsigned start = random_integer(0, nbeads - 1);
            unsigned end   = (start + num_beads_moved + 1) % nbeads;
            if (conf->type_of_polymers[atom] == 'o') {
                start = random_integer(0, nbeads - 2 - num_beads_moved);
                end   = start + num_beads_moved + 1;
            }

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

            if (conf->type_of_polymers[atom] == 'o') {
                std::shared_ptr<Configuration> conf_new_zero(conf->duplicate());

                double tau0     = get_temp(0, atom);
                double tau1     = get_temp(1, atom);
                double tau      = (tau1 - tau0) * 2 * lambda(atom);
                double sigma_bb = std::sqrt(tau);

                arma::vec new_pos = arma::zeros<arma::vec>(conf->num_dims());
                for (unsigned d = 0; d < conf->num_dims(); d++)
                    new_pos(d) = random_normal(conf_new->augmented_bead_position(atom, 1, d), sigma_bb);

                conf_new_zero->augmented_set(atom, 0, bc->wrap_vector(new_pos));
                check_amplitude(conf, conf_new_zero, atom);

                std::shared_ptr<Configuration> conf_new_last(conf->duplicate());

                tau0     = get_temp(nbeads - 2, atom);
                tau1     = get_temp(nbeads - 1, atom);
                tau      = (tau1 - tau0) * 2 * lambda(atom);
                sigma_bb = std::sqrt(tau);

                new_pos = arma::zeros<arma::vec>(conf->num_dims());
                for (unsigned d = 0; d < conf->num_dims(); d++)
                    new_pos(d) = random_normal(conf_new->augmented_bead_position(atom, nbeads - 2, d), sigma_bb);

                conf_new_last->augmented_set(atom, nbeads - 1, bc->wrap_vector(new_pos));
                check_amplitude(conf, conf_new_last, atom);
            }
        }
    }
}
