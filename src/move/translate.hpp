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
            for(unsigned time_ind=0; time_ind < conf.positions.n_cols; time_ind++)
                conf_new.positions.tube(atom_nums(atom), time_ind) += step;

            check_amplitude(conf, conf_new);
        }
        return conf;
    }
};

REGISTER_TYPE_GENERAL(Translate, Move)

#endif
