#ifndef _PIWIGNER_CONFIGURATION_HPP_
#define _PIWIGNER_CONFIGURATION_HPP_

#include "configuration.hpp"

class WignerConfiguration : public Configuration {
protected:
    arma::mat wigner_pos, momentum;
public:
    WignerConfiguration() {}
    WignerConfiguration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums, double mass);
    virtual void load_config(std::string filename);
    virtual Configuration* duplicate() {return new WignerConfiguration(*this);}
    virtual void shift(unsigned atom_num, arma::vec shift_amt);
    arma::mat get_momentum() const { return momentum; }
    void set_position(unsigned atom_num, unsigned dim, double value);
    void set_momentum(unsigned atom_num, unsigned dim, double value);
    virtual std::complex<double> weight() const;
    virtual arma::mat pos() const;
    virtual std::string header() const;
    virtual std::string repr(const boost::shared_ptr<Potential> &V) const;
    virtual std::vector<double> to_vec() const;
};

#endif
