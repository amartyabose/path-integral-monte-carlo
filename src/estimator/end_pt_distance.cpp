#include <string>

#include "end_pt_distance.hpp"

arma::mat End2EndDist::eval(std::shared_ptr<Configuration> const &x) {
    arma::vec zeroth_pos = x->bead_position(atom_num, 0);
    arma::vec final_pos  = x->augmented_bead_position(atom_num, x->num_augmented_beads(atom_num));
    return zeroth_pos - final_pos;
}
