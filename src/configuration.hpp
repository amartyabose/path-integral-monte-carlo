#ifndef _CONFIGURATION_HPP_
#define _CONFIGURATION_HPP_

#include <vector>

#include <armadillo>

#include "potentials.hpp"

class Configuration {
protected:
    arma::cube positions;
    std::vector<char> type_of_polymers;   // 'c' for closed; 'o' for open
    // if type_of_polymers is 'c' for a particular atom,
    // then the first and last beads would be identical
    // else the first and last beads can differ
    std::vector<unsigned> bead_num;
public:
    Configuration() {}
    Configuration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums);
    virtual Configuration* duplicate() {return new Configuration(*this);}
    void augmented_set(unsigned atom_num, unsigned time_ind, unsigned dim, double value);
    virtual void shift(unsigned atom_num, arma::vec shift_amt);
    unsigned num_dims() const;
    unsigned num_beads() const;
    unsigned num_augmented_beads() const;
    unsigned num_atoms() const;
    double augmented_bead_position(unsigned atom_num, unsigned time_ind, unsigned dim) const;
    arma::vec bead_position(unsigned atom_num, unsigned time_ind) const;
    arma::mat necklace(unsigned atom_num) const;
    arma::mat time_slice(unsigned time_ind) const;
    arma::cube get_augmented_segment(unsigned t1, unsigned t2) const;
    arma::mat CoM() const;
    virtual std::complex<double> weight() const { return 1.; }
    virtual arma::mat pos() const;
    virtual std::string header() const;
    virtual std::string repr(const boost::shared_ptr<Potential> &V) const;
};

#endif
