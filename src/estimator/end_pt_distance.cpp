#include <string>

#include "end_pt_distance.hpp"

arma::mat End2EndDist::eval(std::shared_ptr<Configuration> const &x) {
    arma::vec zeroth_pos = x->bead_position(atom_num, 0);
    arma::vec final_pos  = x->augmented_bead_position(atom_num, x->num_augmented_beads(atom_num));
    arma::mat diff       = arma::zeros<arma::mat>(n_rows, n_cols);
    arma::vec vecdiff    = bc->wrap_vector(zeroth_pos - final_pos);
    for (unsigned i=0; i<vecdiff.n_rows; i++)
        diff(0, i) = vecdiff(i);
    return diff;
}
