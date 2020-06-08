#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

#include <vector>

#include <armadillo>

#include "potential/potential.hpp"

class Configuration {
protected:
    arma::vec  mass;
    arma::cube positions;

    unsigned natoms, ndims;

    std::vector<std::string> atom_names;

    std::vector<char> type_of_polymers; // 'c' for closed; 'o' for open
    // if type_of_polymers is 'c' for a particular atom,
    // then the first and last beads would be identical
    // else the first and last beads can differ
    std::vector<unsigned> bead_num;

public:
    Configuration() = default;
    Configuration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums, arma::vec mass,
                  arma::mat bead_specific_mass = arma::zeros<arma::mat>(0, 0));
    virtual void           load_config(std::string filename);
    virtual void           random_config();
    virtual Configuration *duplicate() { return new Configuration(*this); }
    void                   augmented_set(unsigned atom_num, unsigned time_ind, unsigned dim, double value);
    void                   augmented_set(unsigned atom_num, unsigned time_ind, arma::vec value);
    virtual void           shift(unsigned atom_num, arma::vec shift_amt);
    virtual void           shift_time_bead(unsigned time_ind, arma::mat vals);
    unsigned               num_dims() const;
    unsigned               num_beads() const;
    unsigned               num_augmented_beads() const;
    unsigned               num_atoms() const;
    double                 augmented_bead_position(unsigned atom_num, unsigned time_ind, unsigned dim) const;
    arma::vec              augmented_bead_position(unsigned atom_num, unsigned time_ind) const;
    arma::vec              bead_position(unsigned atom_num, unsigned time_ind) const;
    arma::mat              necklace(unsigned atom_num) const;
    arma::mat              time_slice(unsigned time_ind) const;
    arma::mat              augmented_time_slice(unsigned time_ind) const;
    arma::cube             get_augmented_segment(unsigned t1, unsigned t2) const;
    arma::mat              CoM() const;
    virtual arma::vec      get_mass() const { return mass; }
    virtual arma::mat      pos() const;
    virtual arma::mat      get_momentum() const;
    virtual std::string    repr(int frame_cnt) const;

    virtual std::complex<double> weight() const { return 1.; }
};

#endif
