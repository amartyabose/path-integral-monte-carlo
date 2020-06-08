#ifndef _PROPAGATOR_HPP_
#define _PROPAGATOR_HPP_

#include <map>
#include <memory>

#include <armadillo>

#include "../configuration.hpp"
#include "../potential/potential.hpp"

class Propagator;

struct PropagatorFactory {
    virtual std::shared_ptr<Propagator> create() = 0;
};

class Propagator {
    static std::map<std::string, PropagatorFactory *> &get_factory();

protected:
    double   tau;
    unsigned num_total_beads;
    // This is the total number of beads required to describe a single
    // high temperature propagator. There would be two real beads (the
    // extreme ones) and then whatever number of virtual beads...
public:
    virtual void set_params(double Tau, pt::ptree::value_type p, arma::vec mass_, double beta) {
        tau = Tau / (num_total_beads - 1.);
    }
    virtual double     get_tau(unsigned atom) { return tau; }
    virtual double     operator()(const arma::cube &conf) = 0;
    virtual arma::cube derivative(const arma::cube &conf) { return arma::zeros<arma::cube>(arma::size(conf)); }
    unsigned           total_beads() { return num_total_beads; }

    static void                        registerType(const std::string &name, PropagatorFactory *factory);
    static std::shared_ptr<Propagator> create(const std::string &name);
};

#endif
