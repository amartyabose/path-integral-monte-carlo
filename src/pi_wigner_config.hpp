#ifndef _PIWIGNER_CONFIGURATION_HPP_
#define _PIWIGNER_CONFIGURATION_HPP_

#include "configuration.hpp"

class WignerConfiguration : public Configuration {
protected:
    arma::mat wigner_pos, momentum;
public:
    WignerConfiguration() {}
    WignerConfiguration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums);
    virtual Configuration* duplicate() {return new WignerConfiguration(*this);}
    arma::mat get_momentum() const { return momentum; }
    void set_position(unsigned atom_num, unsigned dim, double value);
    void set_momentum(unsigned atom_num, unsigned dim, double value);
    virtual std::string repr() const;
};

#endif
