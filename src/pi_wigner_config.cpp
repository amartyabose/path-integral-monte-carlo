#include <cmath>
#include <complex>
#include <string>

#include <boost/lexical_cast.hpp>

#include "pi_wigner_config.hpp"

WignerConfiguration::WignerConfiguration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums) : Configuration(natoms, ndimensions, bead_nums) {
    momentum = arma::zeros<arma::mat>(natoms, ndimensions);
    wigner_pos = arma::zeros<arma::mat>(natoms, ndimensions);
}

void WignerConfiguration::set_position(unsigned atom_num, unsigned dim, double value) {
    wigner_pos(atom_num, dim) = value;
}

void WignerConfiguration::set_momentum(unsigned atom_num, unsigned dim, double value) {
    momentum(atom_num, dim) = value;
}

arma::mat WignerConfiguration::pos() const {
    return wigner_pos;
}

void WignerConfiguration::shift(unsigned atom_num, arma::vec shift_amt) {
    Configuration::shift(atom_num, shift_amt);
    wigner_pos.row(atom_num) += shift_amt;
}

std::complex<double> WignerConfiguration::weight() const {
    std::complex<double> I(0,1);
    return std::exp(-I*arma::accu(momentum % (time_slice(0)-time_slice(num_beads()-2))));
}

std::string WignerConfiguration::header() const {
    std::string head = "weight\t";
    if(positions.n_slices==1)
        for(unsigned atom=0; atom<positions.n_rows; atom++)
            head += "pos_" + boost::lexical_cast<std::string>(atom) + "\tmom_" + boost::lexical_cast<std::string>(atom) + "\t";
    head += "potential\n";
    return head;
}

std::string WignerConfiguration::repr(const boost::shared_ptr<Potential> &V) const {
    std::string data = boost::lexical_cast<std::string>(weight().real()) + "\t";

    if(wigner_pos.n_cols==1)
        for(unsigned r=0; r<wigner_pos.n_rows; r++)
            for(unsigned c=0; c<wigner_pos.n_cols; c++)
                data += boost::lexical_cast<std::string>(wigner_pos(r, c)) + "\t" + boost::lexical_cast<std::string>(momentum(r,c)) + "\t";

    data += boost::lexical_cast<std::string>((*V)(pos()));

    return data + "\n";
}
