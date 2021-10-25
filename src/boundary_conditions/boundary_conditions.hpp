#ifndef _BOUNDARY_CONDS_H_
#define _BOUNDARY_CONDS_H_

#include <map>
#include <memory>

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <armadillo>

class BoundaryConditions;

struct BoundaryConditionsFactory {
    virtual std::shared_ptr<BoundaryConditions> create() = 0;
};

class BoundaryConditions {
    static std::map<std::string, BoundaryConditionsFactory *> &get_factory();

public:
    virtual void         setup(pt::ptree p)                          = 0;
    virtual arma::mat    center_box(arma::mat const &x) const        = 0;
    virtual arma::mat    wrap_coordinates(arma::mat const &x) const  = 0;
    virtual arma::rowvec wrap_rowvector(arma::rowvec const &x) const = 0;
    virtual arma::vec    wrap_vector(arma::vec const &x) const       = 0;

    static void registerType(const std::string &name, BoundaryConditionsFactory *factory);
    static std::shared_ptr<BoundaryConditions> create(const std::string &name);
};

extern std::shared_ptr<BoundaryConditions> bc;

#endif
