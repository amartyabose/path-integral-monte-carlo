#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <map>
#include <memory>

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <armadillo>

#include "../boundary_conditions/boundary_conditions.hpp"

class Potential;

struct PotentialFactory {
    virtual std::shared_ptr<Potential> create() = 0;
};

class Potential {
    static std::map<std::string, PotentialFactory *> &get_factory();

public:
    virtual void                 setup(pt::ptree node)             = 0;
    virtual double               operator()(arma::mat const &x)    = 0;
    virtual std::complex<double> operator()(arma::cx_mat const &x) = 0;

    virtual arma::mat derivative(arma::mat const &x);

    static void registerType(const std::string &name, PotentialFactory *factory);

    static std::shared_ptr<Potential> create(const std::string &name);
};

extern std::shared_ptr<Potential> pot;

#endif
