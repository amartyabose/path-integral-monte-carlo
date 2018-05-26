#ifndef _MOVE_HPP_
#define _MOVE_HPP_

#include <map>

#include <armadillo>

#include <boost/shared_ptr.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
namespace pt = boost::property_tree;

#include "../random.hpp"
#include "../configuration.hpp"
#include "../propagators.hpp"

class Move;

struct MoveFactory {
    virtual Move* create() = 0;
};

class Move {
    static std::map<std::string, MoveFactory*> factories;
public:
    boost::shared_ptr<Propagator> propagator;
    unsigned moves_accepted, moves_tried;

    Move() {
        moves_accepted = moves_tried = 0;
    }

    virtual void setup(pt::ptree::value_type node, pt::ptree params) = 0;

    virtual Configuration operator() (Configuration conf, arma::uvec atom_nums) = 0;

    void check_amplitude(Configuration &conf_old, Configuration &conf_new) {
        moves_tried++;
        if((*propagator)(conf_new)/(*propagator)(conf_old) > random_float(0,1)) {
            conf_old = conf_new;
            moves_accepted++;
        }
    }

    static void registerType(const std::string &name, MoveFactory *factory) {
        factories[name] = factory;
    }

    static boost::shared_ptr<Move> create(const std::string &name) {
        if(factories.find(name)==factories.end()) {
            std::cerr<<"Not a valid move"<<std::endl;
            exit(1);
        }
        return boost::shared_ptr<Move>(factories[name]->create());
    }
};
std::map<std::string, MoveFactory*> Move::factories;

#endif
