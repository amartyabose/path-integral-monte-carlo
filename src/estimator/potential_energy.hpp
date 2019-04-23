#ifndef _PE_H_
#define _PE_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class PotentialEnergy : public Estimator {
public:
    shared_ptr<Potential> pot;
    void setup(std::string type_, pt::ptree params, unsigned nblocks=10) {
        pt::ptree pot_tree = params.get_child("potential");
        pot = Potential::create(pot_tree.get<string>("<xmlattr>.name"));
        pot->setup(pot_tree);
        values_real_plus = values_real_minus = values_imag_plus = values_imag_minus = arma::zeros<arma::vec>(nblocks);
        i_real_plus = i_real_minus = i_imag_plus = i_imag_minus = arma::zeros<arma::vec>(nblocks);
    }
    double eval(boost::shared_ptr<Configuration> x) {
        return (*pot)(x->pos());
    }
};

REGISTER_TYPE_GENERAL(PotentialEnergy, Estimator)

#endif
