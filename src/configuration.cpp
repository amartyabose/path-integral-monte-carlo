#include <algorithm>
#include <string>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "spdlog/spdlog.h"

#include "boundary_conditions/boundary_conditions.hpp"
#include "configuration.hpp"
#include "random.hpp"

Configuration::Configuration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums, arma::vec m,
                             arma::mat bead_specific_mass) {
    bead_num  = bead_nums;
    positions = arma::zeros<arma::cube>(natoms, bead_nums.back() + 1, ndimensions);
    type_of_polymers.resize(natoms, 'c');
    mass = m;
}

void Configuration::load_config(std::string filename) {
    arma::mat pos;
    pos.load(filename);
    arma::mat wrapped_pos = bc->wrap_coordinates(pos);
    for (unsigned b = 0; b < positions.n_cols; b++)
        positions(arma::span::all, arma::span(b), arma::span::all) = wrapped_pos;
}

void Configuration::random_config() {
    for (unsigned d = 0; d < num_dims(); d++)
        for (unsigned b = 0; b < num_augmented_beads(); b++)
            for (unsigned atom = 0; atom < num_atoms(); atom++)
                augmented_set(atom, b, d, random_normal(0, 1));
}

void Configuration::augmented_set(unsigned atom_num, unsigned time_ind, unsigned dim, double value) {
    positions(atom_num, time_ind, dim) = value;
    if (type_of_polymers[atom_num] == 'c' && time_ind == 0)
        positions(atom_num, num_augmented_beads(), dim) = value;
}

void Configuration::augmented_set(unsigned atom_num, unsigned time_ind, arma::vec value) {
    for (unsigned d = 0; d < num_dims(); d++) {
        positions(atom_num, time_ind, d) = value(d);
        if (type_of_polymers[atom_num] == 'c' && time_ind == 0)
            positions(atom_num, num_augmented_beads(), d) = value(d);
    }
}

// Not used. Is this correct?????
void Configuration::shift_time_bead(unsigned time_ind, arma::mat vals) {
    positions(arma::span::all, arma::span(time_ind), arma::span::all) += vals;
    if (!time_ind)
        positions(arma::span::all, arma::span(num_augmented_beads()), arma::span::all) += vals;
}

void Configuration::shift(unsigned atom_num, arma::vec shift_amt) {
    for (unsigned time_ind = 0; time_ind < positions.n_cols; time_ind++) {
        arma::vec pos     = positions.tube(atom_num, time_ind);
        arma::vec new_pos = bc->wrap_vector(pos + shift_amt);

        positions.tube(atom_num, time_ind) = new_pos;
    }
}

unsigned Configuration::num_dims() const { return positions.n_slices; }

unsigned Configuration::num_beads() const { return bead_num.size(); }

unsigned Configuration::num_augmented_beads() const {
    // return (type_of_polymers[atom_num]=='o') ? positions.n_cols :
    // positions.n_cols-1;
    return positions.n_cols - 1;
}

unsigned Configuration::num_atoms() const { return positions.n_rows; }

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
        arma::mat slice = bc->center_box(time_slice(slices));
        for (unsigned r = 0; r < slice.n_rows; r++) {
            output += "AT\t";
            for (unsigned c = 0; c < slice.n_cols; c++)
                output += (boost::format("%.5e\t") % slice(r, c)).str();
            for (int k = 0; k < 3 - slice.n_cols; k++)
                output += (boost::format("%.5e\t") % 0).str();
            output += "\n";
        }
    }
    return output;
}
