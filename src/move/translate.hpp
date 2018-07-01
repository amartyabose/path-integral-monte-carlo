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

    void operator()(boost::shared_ptr<Configuration> &conf, arma::uvec atom_nums) {
        arma::vec step(max_step.n_cols);
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++) {
            for(unsigned i=0; i<max_step.n_cols; i++)
                step(i) = random_float(-1,1) * max_step(atom_nums(atom), i);

            boost::shared_ptr<Configuration> conf_new(conf->duplicate());
            conf_new->shift(atom_nums(atom), step);

            check_amplitude(conf, conf_new);
        }
    }
};

REGISTER_TYPE_GENERAL(Translate, Move)

#endif
