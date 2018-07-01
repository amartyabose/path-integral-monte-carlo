#ifndef _WIGNER_POSITION_HPP_
#define _WIGNER_POSITION_HPP_

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class WignerPosition : public Move {
    double delta_beta, mass, sigma;
public:
    void setup(pt::ptree::value_type node, pt::ptree params) {
        delta_beta = node.second.get<double>("<xmlattr>.dt");
        mass = params.get<double>("mass", 1);
        sigma = std::sqrt(delta_beta/(4.*mass));
    }

    void operator()(boost::shared_ptr<Configuration> &conf, arma::uvec atom_nums) {
        arma::mat xmean = (conf->time_slice(0) + conf->time_slice(conf->num_beads()-2))/2.;
        unsigned dimensions = conf->num_dims();
        WignerConfiguration *wigner_conf = dynamic_cast<WignerConfiguration*>(conf.get());
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++) {
            for(unsigned dim=0; dim<dimensions; dim++)
                wigner_conf->set_position(atom_nums(atom), dim, random_normal(xmean(atom_nums(atom), dim), sigma));
            moves_tried++;
            moves_accepted++;
        }
    }
};

REGISTER_TYPE_GENERAL(WignerPosition, Move)

#endif
