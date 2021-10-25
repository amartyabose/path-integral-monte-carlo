#ifndef _TRANSLATE_HPP_
#define _TRANSLATE_HPP_

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class Translate : public Move {
    double max_step;

public:
    void setup(pt::ptree::value_type node, double beta_, arma::vec mass_) override;
    void operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) const override;
};

REGISTER_TYPE_GENERAL(Translate, Move)

#endif
