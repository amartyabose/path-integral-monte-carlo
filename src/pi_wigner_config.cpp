#include <cmath>
#include <complex>
#include <string>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "pi_wigner_config.hpp"

WignerConfiguration::WignerConfiguration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums,
                                         arma::vec mass, arma::mat bead_specific_mass)
    : Configuration(natoms, ndimensions, bead_nums, mass, bead_specific_mass) {
    momentum   = arma::zeros<arma::mat>(natoms, ndimensions);
    wigner_pos = arma::zeros<arma::mat>(natoms, ndimensions);
}

void WignerConfiguration::load_config(std::string filename) {
    Configuration::load_config(filename);
    wigner_pos.load(filename);
}

void WignerConfiguration::set_position(unsigned atom_num, unsigned dim, double value) {
    wigner_pos(atom_num, dim) = value;
}

void WignerConfiguration::set_momentum(unsigned atom_num, unsigned dim, double value) {
    momentum(atom_num, dim) = value;
}

arma::mat WignerConfiguration::pos() const { return wigner_pos; }

void WignerConfiguration::shift(unsigned atom_num, arma::vec shift_amt) {
    Configuration::shift(atom_num, shift_amt);
    wigner_pos.row(atom_num) += shift_amt.st();
}

std::complex<double> WignerConfiguration::weight() const {
    std::complex<double> I(0, 1);
    return std::exp(-I * arma::accu(momentum % (time_slice(0) - time_slice(num_beads() - 2))));
}

std::string WignerConfiguration::repr(int frame_cnt) const {
    std::string output = std::to_string(num_atoms()) + "\n" + std::to_string(frame_cnt) + "\n";
    for (unsigned r = 0; r < wigner_pos.n_rows; r++) {
        output += "AT\t";
        for (unsigned c = 0; c < wigner_pos.n_cols; c++)
            output += (boost::format("%.5e\t") % wigner_pos(r, c)).str();
        for (int k = 0; k < 3 - wigner_pos.n_cols; k++)
            output += (boost::format("%.5e\t") % 0).str();
        for (unsigned c = 0; c < momentum.n_cols; c++)
            output += (boost::format("%.5e\t") % momentum(r, c)).str();
        for (int k = 0; k < 3 - wigner_pos.n_cols; k++)
            output += (boost::format("%.5e\t") % 0).str();
        output += "\n";
    }
    return output;
}
