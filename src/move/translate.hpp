#ifndef _TRANSLATE_HPP_
#define _TRANSLATE_HPP_

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class Translate : public Move {
    arma::mat max_step;
public:
    void setup(pt::ptree::value_type node, pt::ptree params) {
        max_step = node.second.get<double>("length") * arma::ones<arma::mat>(params.get<unsigned>("num_atoms"), params.get<unsigned>("num_dimensions"));
    }

    Configuration operator()(Configuration conf, arma::uvec atom_nums) {
        arma::vec step(max_step.n_cols);
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++) {
            for(unsigned i=0; i<max_step.n_cols; i++)
                step(i) = random_float(-1,1) * max_step(atom_nums(atom), i);

            Configuration conf_new = conf;
            for(unsigned time_ind=0; time_ind < conf.num_beads(atom); time_ind++)
                for(unsigned d=0; d<max_step.n_cols; d++)
                    //conf_new.positions.tube(atom_nums(atom), time_ind) += step;
                    conf_new.augmented_set(atom_nums(atom), time_ind, d,
                                           conf_new.augmented_bead_position(atom_nums(atom), time_ind, d) + step(d));

            check_amplitude(conf, conf_new);
        }
        return conf;
    }
};

REGISTER_TYPE_GENERAL(Translate, Move)

#endif
