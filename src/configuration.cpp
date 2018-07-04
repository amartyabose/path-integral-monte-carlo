#include <algorithm>
#include <string>

#include <boost/lexical_cast.hpp>

#include "configuration.hpp"

Configuration::Configuration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums) {
    bead_num = bead_nums;
    positions = arma::zeros<arma::cube>(natoms, bead_nums.back()+1, ndimensions);
    type_of_polymers.resize(natoms, 'c');
}

void Configuration::augmented_set(unsigned atom_num, unsigned time_ind, unsigned dim, double value) {
    positions(atom_num, time_ind, dim) = value;
    if (type_of_polymers[atom_num] == 'c' && time_ind == 0)
        //positions(atom_num, num_beads(atom_num), dim) = value;
        positions(atom_num, num_augmented_beads(), dim) = value;
}

void Configuration::shift(unsigned atom_num, arma::vec shift_amt) {
    for(unsigned time_ind=0; time_ind<positions.n_cols; time_ind++)
        positions.tube(atom_num, time_ind) += shift_amt;
}

unsigned Configuration::num_dims() const {
    return positions.n_slices;
}

unsigned Configuration::num_beads() const {
    return bead_num.size();
}

unsigned Configuration::num_augmented_beads() const {
    //return (type_of_polymers[atom_num]=='o') ? positions.n_cols : positions.n_cols-1;
    return positions.n_cols-1;
}

unsigned Configuration::num_atoms() const {
    return positions.n_rows;
}

double Configuration::augmented_bead_position(unsigned atom_num, unsigned time_ind, unsigned dim) const {
    return positions(atom_num, time_ind, dim);
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

arma::cube Configuration::get_augmented_segment(unsigned t1, unsigned t2) const {
    return positions(arma::span::all, arma::span(bead_num[t1], bead_num[t2]), arma::span::all);
}

arma::mat Configuration::CoM() const {
    return arma::mean(positions, 1);
}

arma::mat Configuration::pos() const {
    return time_slice(0);
}

std::string Configuration::header() const {
    std::string head = "weight\t";
    if(positions.n_slices==1)
        for(unsigned atom=0; atom<positions.n_rows; atom++)
            head += "pos_" + boost::lexical_cast<std::string>(atom) + "\t";
    head += "potential\n";
    return head;
}

std::string Configuration::repr(const boost::shared_ptr<Potential> &V) const {
    std::string data = boost::lexical_cast<string>(weight().real()) + "\t";
    arma::mat t0 = time_slice(0);
    if(t0.n_cols==1)
        for(unsigned r=0; r<t0.n_rows; r++)
            for(unsigned c=0; c<t0.n_cols; c++)
                data += boost::lexical_cast<std::string>(t0(r, c)) + "\t";

    data += boost::lexical_cast<std::string>((*V)(pos()));

    return data + "\n";
}
