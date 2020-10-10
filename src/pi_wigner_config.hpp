#ifndef _PIWIGNER_CONFIGURATION_HPP_
#define _PIWIGNER_CONFIGURATION_HPP_

#include "configuration.hpp"

class WignerConfiguration : public Configuration {
protected:
    arma::mat wigner_pos, momentum;

public:
    WignerConfiguration() = default;
    WignerConfiguration(unsigned natoms, unsigned ndimensions, std::vector<unsigned> bead_nums, arma::vec mass);
    void           load_config(std::string filename) override;
    Configuration *duplicate() override { return new WignerConfiguration(*this); }
    void           shift(unsigned atom_num, arma::vec shift_amt) override;
    arma::mat      get_momentum() const override { return momentum; }
    void           set_position(unsigned atom_num, unsigned dim, double value);
    void           set_momentum(unsigned atom_num, unsigned dim, double value);
    arma::mat      pos() const override;
    std::string    repr(int frame_cnt, int slice_num) const override;

    std::complex<double> weight() const override;
};

#endif
