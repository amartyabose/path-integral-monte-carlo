#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

#include <vector>

#include <armadillo>

class Configuration {
protected:
    arma::cube positions;
    std::vector<char> type_of_polymers;   // 'c' for closed; 'o' for open
    // if type_of_polymers is 'c' for a particular atom,
    // then the first and last beads would be identical
    // else the first and last beads can differ
    std::vector<unsigned> bead_num;
public:
    Configuration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums);
    void augmented_set(unsigned atom_num, unsigned time_ind, unsigned dim, double value);
    void shift(unsigned atom_num, arma::vec shift_amt);
    unsigned num_dims() const;
    unsigned num_beads(unsigned atom_num) const;
    unsigned num_atoms() const;
    double augmented_bead_position(unsigned atom_num, unsigned time_ind, unsigned dim) const;
    arma::vec bead_position(unsigned atom_num, unsigned time_ind) const;
    arma::mat necklace(unsigned atom_num) const;
    arma::mat time_slice(unsigned time_ind) const;
    arma::cube get_augmented_segment(unsigned t1, unsigned t2) const;
    arma::mat CoM() const;
};

#endif
