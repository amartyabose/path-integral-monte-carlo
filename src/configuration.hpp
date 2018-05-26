#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

#include <armadillo>

struct Configuration {
    arma::cube positions;
    Configuration(unsigned natoms, unsigned nbeads, unsigned ndimensions) {
        positions = arma::zeros<arma::cube>(natoms, nbeads, ndimensions);
    }
    unsigned num_dims() const {
        return positions.n_slices;
    }
    unsigned num_beads() const {
        return positions.n_cols;
    }
    unsigned num_atoms() const {
        return positions.n_rows;
    }
    arma::vec bead_position(unsigned atom_num, unsigned time_ind) const {
        return positions.tube(atom_num, time_ind);
    }
    arma::mat necklace(unsigned atom_num) const {
        return positions(arma::span(atom_num), arma::span::all, arma::span::all);
    }
    arma::mat time_slice(unsigned time_ind) const {
        return positions(arma::span::all, arma::span(time_ind), arma::span::all);
    }
    arma::mat CoM() const {
        return arma::mean(positions, 1);
    }
};

#endif
