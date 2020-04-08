#ifndef _MOVE_HPP_
#define _MOVE_HPP_

#include <map>
#include <memory>

#include <armadillo>

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include "../configuration.hpp"
#include "../pi_wigner_config.hpp"
#include "../propagator/propagator.hpp"
#include "../random.hpp"

class Move;

struct MoveFactory {
    virtual std::shared_ptr<Move> create() = 0;
};

class Move {
    static std::map<std::string, MoveFactory *> &get_factory();

public:
    std::vector<std::shared_ptr<Propagator>> propagator;

    unsigned moves_accepted, moves_tried;

    Move() { moves_accepted = moves_tried = 0; }

    virtual void setup(pt::ptree::value_type node, double beta, arma::vec mass) = 0;

    virtual void set_beta(arma::vec beta) {}

    virtual void operator()(std::shared_ptr<Configuration> &conf, arma::uvec atom_nums) = 0;

    virtual void check_amplitude(std::shared_ptr<Configuration> &conf_old, std::shared_ptr<Configuration> conf_new) {
        moves_tried++;
        double new_weight = 1, old_weight = 1;
        for (unsigned p = 0; p < propagator.size(); p++) {
            new_weight *= (*propagator[p])(conf_new->get_augmented_segment(p, p + 1));
            if (new_weight < 0)
                return;
            old_weight *= (*propagator[p])(conf_old->get_augmented_segment(p, p + 1));
        }
        if (new_weight / old_weight > random_float(0, 1)) {
            conf_old.reset(conf_new->duplicate());
            moves_accepted++;
        }
    }

    static void                  registerType(const std::string &name, MoveFactory *factory);
    static std::shared_ptr<Move> create(const std::string &name);
};

#endif
