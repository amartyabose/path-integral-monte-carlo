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
    static std::map<std::string, PotentialFactory*> factories;
public:
    virtual void setup(pt::ptree node) = 0;
    virtual double operator()(arma::mat x) = 0;
    static void registerType(const std::string &name, PotentialFactory *factory) {
        factories[name] = factory;
    }
    static boost::shared_ptr<Potential> create(const std::string &name) {
        if(factories.find(name) == factories.end()) {
            std::cerr<<"Not a valid potential"<<std::endl;
            exit(1);
        }
        return boost::shared_ptr<Potential>(factories[name]->create());
    }
};
std::map<std::string, PotentialFactory*> Potential::factories;

#endif
