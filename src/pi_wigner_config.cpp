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

std::string WignerConfiguration::repr() const {
    std::string data;
    for(unsigned r=0; r<wigner_pos.n_rows; r++)
        for(unsigned c=0; c<wigner_pos.n_cols; c++)
            data += boost::lexical_cast<std::string>(wigner_pos(r, c)) + "\t" + boost::lexical_cast<std::string>(momentum(r,c)) + "\t";

    //arma::mat t0 = time_slice(0);
    //arma::mat tm1 = time_slice(num_beads()-2);
    //for(unsigned r=0; r<t0.n_rows; r++)
    //    for(unsigned c=0; c<t0.n_cols; c++)
    //        data += boost::lexical_cast<std::string>(t0(r, c)) + "\t";
    //for(unsigned r=0; r<tm1.n_rows; r++)
    //    for(unsigned c=0; c<tm1.n_cols; c++)
    //        data += boost::lexical_cast<std::string>(tm1(r, c)) + "\t";

    std::complex<double> I(0,1);
    std::complex<double> phase = std::exp(-I*arma::accu(momentum % (time_slice(0)-time_slice(num_beads()-2))));
    data += boost::lexical_cast<std::string>(phase.real());

    return data + "\n";
}
