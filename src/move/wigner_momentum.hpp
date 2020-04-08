#ifndef _WIGNER_MOMENTUM_HPP_
#define _WIGNER_MOMENTUM_HPP_

#include <armadillo>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class WignerMomentum : public Move {
    double    beta;
    arma::vec mass, sigma;

public:
    void setup(pt::ptree::value_type node, double beta_, arma::vec mass_) override;
    void operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) override;
    void check_amplitude(std::shared_ptr<Configuration> &conf_old, std::shared_ptr<Configuration> conf_new) override;
};

REGISTER_TYPE_GENERAL(WignerMomentum, Move)

#endif
