#include <algorithm>
#include <string>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "spdlog/spdlog.h"

#include "boundary_conditions/boundary_conditions.hpp"
#include "configuration.hpp"
#include "random.hpp"

Configuration::Configuration(unsigned natoms_, unsigned ndimensions_, std::vector<unsigned> bead_nums, arma::vec m,
                             arma::vec open_chains) {
    bead_num  = bead_nums;
    positions = arma::zeros<arma::cube>(natoms_, bead_nums.back() + 1, ndimensions_);
    type_of_polymers.resize(natoms_, 'c');
    for (unsigned i = 0; i < open_chains.n_rows; ++i)
        type_of_polymers[open_chains[i]] = 'o';

    atom_names.resize(natoms_);
    mass   = m;
    natoms = natoms_;
    ndims  = ndimensions_;
}

void Configuration::load_config(std::string filename) {
    std::ifstream ifs(filename);

    unsigned num_atoms_from_config;
    ifs >> num_atoms_from_config;
    if (num_atoms_from_config != natoms)
        throw std::runtime_error("Number of atoms in the initial configuration file, " +
                                 std::to_string(num_atoms_from_config) + " does not match the input file, " +
                                 std::to_string(natoms) + ".");

    arma::mat pos(natoms, ndims);
    for (unsigned i = 0; i < natoms; i++) {
        std::string name;
        ifs >> name;
        atom_names[i] = name;
        for (unsigned d = 0; d < ndims; d++) {
            double val;
            ifs >> val;
            pos(i, d) = val;
        }
    }

    for (unsigned b = 0; b < positions.n_cols; b++)
        positions(arma::span::all, arma::span(b), arma::span::all) = pos;
}

void Configuration::random_config() {
    for (unsigned d = 0; d < num_dims(); d++)
        for (unsigned atom = 0; atom < num_atoms(); atom++) {
            if (type_of_polymers[atom] == 'c')
                for (unsigned b = 0; b < num_augmented_beads(atom); b++)
                    augmented_set(atom, b, d, random_normal(0, 1));
            else
                for (unsigned b = 0; b <= num_augmented_beads(atom); b++)
                    augmented_set(atom, b, d, random_normal(0, 1));
        }
}

void Configuration::augmented_set(unsigned atom_num, unsigned time_ind, unsigned dim, double value) {
    positions(atom_num, time_ind, dim) = value;
    if (type_of_polymers[atom_num] == 'c' && time_ind == 0)
        positions(atom_num, num_augmented_beads(atom_num), dim) = value;
}

void Configuration::augmented_set(unsigned atom_num, unsigned time_ind, arma::vec value) {
    for (unsigned d = 0; d < num_dims(); d++) {
        positions(atom_num, time_ind, d) = value(d);
        if (type_of_polymers[atom_num] == 'c' && time_ind == 0)
            positions(atom_num, num_augmented_beads(atom_num), d) = value(d);
    }
}

void Configuration::shift(unsigned atom_num, arma::vec shift_amt) {
    for (unsigned time_ind = 0; time_ind < positions.n_cols; time_ind++) {
        arma::vec pos     = positions.tube(atom_num, time_ind);
        arma::vec new_pos = bc->wrap_vector(pos + shift_amt);

        for (unsigned d = 0; d < num_dims(); d++)
            positions(atom_num, time_ind, d) = new_pos(d);
    }
}

unsigned Configuration::num_dims() const { return ndims; }

unsigned Configuration::num_beads() const { return bead_num.size(); }

unsigned Configuration::num_augmented_beads(unsigned atom_num) const { return positions.n_cols - 1; }

unsigned Configuration::num_atoms() const { return natoms; }

double Configuration::augmented_bead_position(unsigned atom_num, unsigned time_ind, unsigned dim) const {
    return arma::as_scalar(positions(atom_num, time_ind, dim));
}

arma::vec Configuration::augmented_bead_position(unsigned atom_num, unsigned time_ind) const {
    return positions.tube(atom_num, time_ind);
}

arma::vec Configuration::bead_position(unsigned atom_num, unsigned time_ind) const {
    return positions.tube(atom_num, bead_num[time_ind]);
}

arma::mat Configuration::necklace(unsigned atom_num) const {
    return positions(arma::span(atom_num), arma::span::all, arma::span::all);
}

arma::mat Configuration::time_slice(unsigned time_ind) const {
    return positions(arma::span::all, arma::span(bead_num[time_ind]), arma::span::all);
}

arma::mat Configuration::augmented_time_slice(unsigned time_ind) const {
    return positions(arma::span::all, arma::span(time_ind), arma::span::all);
}

arma::cube Configuration::get_augmented_segment(unsigned t1, unsigned t2) const {
    return positions(arma::span::all, arma::span(bead_num[t1], bead_num[t2]), arma::span::all);
}

arma::mat Configuration::CoM() const {
    arma::mat ans = arma::zeros<arma::mat>(num_atoms(), num_dims());
    for (unsigned b = 0; b < num_beads() - 1; b++)
        ans += time_slice(b);

    return ans / (num_beads() - 1);
}

arma::mat Configuration::pos() const { return time_slice(0); }

arma::mat Configuration::get_momentum() const { return arma::zeros<arma::mat>(1, num_dims()); }

std::string Configuration::repr(int frame_cnt) const {
    std::string output = std::to_string(num_atoms() * (num_beads() - 1)) + "\n" + std::to_string(frame_cnt) + "\n";
    for (unsigned slices = 0; slices < num_beads() - 1; slices++) {
        arma::mat slice = time_slice(slices);
        for (unsigned r = 0; r < slice.n_rows; r++) {
            output += atom_names[r] + '\t';
            for (unsigned c = 0; c < slice.n_cols; c++)
                output += (boost::format("%.5e\t") % slice(r, c)).str();
            for (int k = 0; k < 3 - slice.n_cols; k++)
                output += (boost::format("%.5e\t") % 0).str();
            output += "slice" + std::to_string(slices + 1) + '\n';
        }
    }
    return output;
}
