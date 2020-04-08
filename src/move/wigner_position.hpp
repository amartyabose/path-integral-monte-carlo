#ifndef _WIGNER_POSITION_HPP_
#define _WIGNER_POSITION_HPP_

#include <armadillo>

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class WignerPosition : public Move {
    double    delta_beta;
    arma::vec mass, sigma;

public:
    void setup(pt::ptree::value_type node, double beta_, arma::vec mass_) override;
    void operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) override;
};

REGISTER_TYPE_GENERAL(WignerPosition, Move)

#endif
