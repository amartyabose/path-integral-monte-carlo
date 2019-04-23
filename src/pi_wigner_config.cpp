#include <cmath>
#include <complex>
#include <string>

#include <boost/lexical_cast.hpp>

#include "pi_wigner_config.hpp"

WignerConfiguration::WignerConfiguration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums, double mass) : Configuration(natoms, ndimensions, bead_nums, mass) {
    momentum = arma::zeros<arma::mat>(natoms, ndimensions);
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

arma::mat WignerConfiguration::pos() const {
    return wigner_pos;
}

void WignerConfiguration::shift(unsigned atom_num, arma::vec shift_amt) {
    Configuration::shift(atom_num, shift_amt);
    wigner_pos.row(atom_num) += shift_amt.st();
}

std::complex<double> WignerConfiguration::weight() const {
    std::complex<double> I(0,1);
    return std::exp(-I*arma::accu(momentum % (time_slice(0)-time_slice(num_beads()-2))));
}

std::string WignerConfiguration::header() const {
    std::string names = "Re(weight), Im(weight)";
    for(unsigned a=0; a<num_atoms(); a++)
        for(unsigned d=0; d<num_dims(); d++) {
            names += ", position atom"+boost::lexical_cast<std::string>(a)+" dim"+boost::lexical_cast<std::string>(d);
            names += ", momentum atom"+boost::lexical_cast<std::string>(a)+" dim"+boost::lexical_cast<std::string>(d);
        }
    names += ", potential energy, kinetic energy";
    return names;
}

std::string WignerConfiguration::repr(const boost::shared_ptr<Potential> &V) const {
    std::string data = boost::lexical_cast<std::string>(weight().real()) + "\t";

    if(wigner_pos.n_cols==1)
        for(unsigned r=0; r<wigner_pos.n_rows; r++)
            for(unsigned c=0; c<wigner_pos.n_cols; c++)
                data += boost::lexical_cast<std::string>(wigner_pos(r, c)) + "\t" + boost::lexical_cast<std::string>(momentum(r,c)) + "\t";

    data += boost::lexical_cast<std::string>((*V)(pos())) + "\t" + boost::lexical_cast<std::string>(arma::accu(momentum%momentum)/(2.*mass));

    return data + "\n";
}

std::vector<double> WignerConfiguration::to_vec() const {
    std::vector<double> vec;
    arma::mat t0 = time_slice(0);
    vec.push_back(weight().real());
    vec.push_back(weight().imag());
    if(t0.n_cols==1)
        for(unsigned r=0; r<t0.n_rows; r++)
            for(unsigned c=0; c<t0.n_cols; c++) {
                vec.push_back(wigner_pos(r, c));
                vec.push_back(momentum(r, c));
            }
    return vec;
}
