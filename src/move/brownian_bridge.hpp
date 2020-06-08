#ifndef _BROWNIAN_BRIDGE_HPP_
#define _BROWNIAN_BRIDGE_HPP_

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class BrownianBridge : public Move {
    unsigned  num_beads_moved;
    unsigned  num_attempts;
    arma::vec beta;
    arma::vec lambda;
    double    get_temp(unsigned bead_num, unsigned atom);

public:
    void setup(pt::ptree::value_type node, double beta_, arma::vec mass_) override;
    void set_beta(arma::vec beta) override;
    void operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) override;
};

REGISTER_TYPE_GENERAL(BrownianBridge, Move)

#endif
