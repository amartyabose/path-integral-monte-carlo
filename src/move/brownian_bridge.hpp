#ifndef _BROWNIAN_BRIDGE_HPP_
#define _BROWNIAN_BRIDGE_HPP_

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class BrownianBridge : public Move {
    unsigned num_beads_moved;
    unsigned num_attempts;
    double beta, tau;
    double get_temp(unsigned bead_num) {
        double temp = 0;
        unsigned beads_counted = 0;
        for(unsigned p=0; beads_counted<bead_num && p<propagator.size(); p++)
            for(unsigned b=0; beads_counted<bead_num && b<propagator[p]->total_beads()-1; b++) {
                beads_counted++;
                temp += propagator[p]->get_tau();
            }
        return temp;
    }
public:
    void setup(pt::ptree::value_type node, pt::ptree params) {
        num_beads_moved = node.second.get<unsigned>("length");
        num_attempts = node.second.get<unsigned>("attempts");
        beta = params.get<double>("beta");//params.get<unsigned>("num_propagators");
        tau = params.get<double>("beta")/params.get<unsigned>("num_propagators");
    }

    Configuration operator()(Configuration conf, arma::uvec atom_nums) {
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++) {
            unsigned nbeads = conf.num_beads(atom_nums(atom));
            for(unsigned attempt=0; attempt<num_attempts; attempt++) {
                unsigned start = random_integer(0, nbeads-1);
                unsigned end = (start + num_beads_moved + 1) % nbeads;
                Configuration conf_new = conf;
                for(unsigned b = (start+1)%nbeads, nmoved=0; nmoved < num_beads_moved; b = (b+1)%nbeads, nmoved++) {
                    double tau0 = get_temp((b+nbeads-1)%nbeads);
                    double tau1 = get_temp(b);
                    double taue = get_temp(end);
                    double start_to_new = tau1 - tau0;
                    double new_to_end = taue - tau1;
                    if(new_to_end<0)
                        new_to_end = taue + (beta - tau1);
                    if(start_to_new<0)
                        start_to_new = tau1 + (beta - tau0);
                    double start_to_end = start_to_new + new_to_end;

                    for(unsigned d=0; d<conf.num_dims(); d++)
                        conf_new.augmented_set(atom_nums(atom), b, d, random_normal(((new_to_end)*conf_new.augmented_bead_position(atom_nums(atom), (nbeads+b-1)%nbeads, d) + start_to_new*conf_new.augmented_bead_position(atom_nums(atom), end, d))/start_to_end, std::sqrt(start_to_new * new_to_end / start_to_end)));
                }
                check_amplitude(conf, conf_new);
            }
        }
        return conf;
    }
};

REGISTER_TYPE_GENERAL(BrownianBridge, Move)

#endif
