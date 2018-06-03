#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_

#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
namespace pt = boost::property_tree;

#include <armadillo>

class Potential;

struct PotentialFactory {
    virtual Potential* create() = 0;
};

class Potential {
    static std::map<std::string, PotentialFactory*>& get_factory();
public:
    virtual void setup(pt::ptree node) = 0;
    virtual double operator()(arma::mat x) = 0;
    static void registerType(const std::string &name, PotentialFactory *factory);
    static boost::shared_ptr<Potential> create(const std::string &name);
};

#endif
