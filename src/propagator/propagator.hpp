#ifndef _PROPAGATOR_HPP_
#define _PROPAGATOR_HPP_

#include <map>

#include <boost/shared_ptr.hpp>

#include "../configuration.hpp"
#include "../potential/potential.hpp"

class Propagator;

struct PropagatorFactory {
    virtual Propagator* create() = 0;
};

class Propagator {
    static std::map<std::string, PropagatorFactory*> factories;
public:
    boost::shared_ptr<Potential> V;
    double tau;
    virtual double operator()(Configuration conf) = 0;
    static void registerType(const std::string &name, PropagatorFactory *factory) {
        factories[name] = factory;
    }
    static boost::shared_ptr<Propagator> create(const std::string &name) {
        if(factories.find(name)==factories.end()) {
            std::cerr<<"Not a valid propagator"<<std::endl;
            exit(1);
        }
        return boost::shared_ptr<Propagator>(factories[name]->create());
    }
};
std::map<std::string, PropagatorFactory*> Propagator::factories;

#endif
