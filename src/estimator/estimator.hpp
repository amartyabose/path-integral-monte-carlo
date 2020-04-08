#ifndef _ESTIMATOR_H_
#define _ESTIMATOR_H_

#include <memory>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <armadillo>

#include "../configuration.hpp"
#include "../potential/potential.hpp"
#include "../propagator/propagator.hpp"
#include "../units.hpp"

class Estimator;

struct EstimatorFactory {
    virtual std::shared_ptr<Estimator> create() = 0;
};

class Estimator {
    static std::map<std::string, EstimatorFactory *> &get_factory();

public:
    unsigned  n_rows, n_cols;
    arma::vec mass;
    arma::mat bead_specific_mass;

    std::shared_ptr<Potential> pot;

    virtual void      setup(std::string type, arma::vec masses, double beta, unsigned num_beads, pt::ptree node) = 0;
    virtual arma::mat eval(std::shared_ptr<Configuration> const &x)                                              = 0;

    static void                       registerType(const std::string &name, EstimatorFactory *factory);
    static std::shared_ptr<Estimator> create(const std::string &name);
};

#endif
