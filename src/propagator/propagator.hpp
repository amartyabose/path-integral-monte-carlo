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
    static std::map<std::string, PropagatorFactory*>& get_factory();
protected:
    boost::shared_ptr<Potential> V;
    double tau;
    unsigned num_total_beads;
    // This is the total number of beads required to describe a single
    // high temperature propagator. There would be two real beads (the
    // extreme ones) and then whatever number of virtual beads...
public:
    virtual void set_params(boost::shared_ptr<Potential> pot, double Tau, pt::ptree::value_type p, pt::ptree params) {
        V = pot;
        tau = Tau / (num_total_beads - 1.);
    }
    virtual double get_tau() {
        return tau;
    }
    virtual double operator()(const arma::cube &conf) = 0;
    unsigned total_beads() {
        return num_total_beads;
    }
    static void registerType(const std::string &name, PropagatorFactory *factory);
    static boost::shared_ptr<Propagator> create(const std::string &name);
};

#endif
