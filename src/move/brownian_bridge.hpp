#ifndef _BROWNIAN_BRIDGE_HPP_
#define _BROWNIAN_BRIDGE_HPP_

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class BrownianBridge : public Move {
    unsigned num_beads_moved;
    unsigned num_attempts;
    double tau;
public:
    void setup(pt::ptree::value_type node, pt::ptree params) {
        num_beads_moved = node.second.get<unsigned>("length");
        num_attempts = node.second.get<unsigned>("attempts");
        tau = params.get<double>("beta")/params.get<unsigned>("num_beads");
    }

    Configuration operator()(Configuration conf, arma::uvec atom_nums) {
        unsigned nbeads = conf.num_beads();
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++)
            for(unsigned attempt=0; attempt<num_attempts; attempt++) {
                unsigned start = random_integer(0, nbeads-1);
                unsigned end = (start + num_beads_moved + 1) % nbeads;
                Configuration conf_new = conf;
                for(unsigned b = (start+1)%nbeads, nmoved=0; nmoved < num_beads_moved; b = (b+1)%nbeads, nmoved++) {
                    unsigned L = num_beads_moved - nmoved + 1;
                    for(unsigned d=0; d<conf.num_dims(); d++)
                        conf_new.positions(atom_nums(atom), b, d) = random_normal(((L-1.)*conf_new.positions(atom_nums(atom), (nbeads+b-1)%nbeads, d) + conf_new.positions(atom_nums(atom), end, d))/L, std::sqrt((L-1.)/L*tau));
                }
                check_amplitude(conf, conf_new);
            }
        return conf;
    }
};

REGISTER_TYPE_GENERAL(BrownianBridge, Move)

#endif
