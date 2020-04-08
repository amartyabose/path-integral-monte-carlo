#ifndef UNITS_H_
#define UNITS_H_

#include <map>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <armadillo>

struct Dimensions {
    std::map<std::string, int> dim;
    Dimensions(std::map<std::string, int> _dim) : dim(_dim) {}
    friend bool operator<(Dimensions const &, Dimensions const &);
};

// This is an arbitrary ordering function to allow us to use Dimensions as
// keys in a map. Change this with impunity if required --- doesn't affect
// the code anywhere else.
bool operator<(Dimensions const &d1, Dimensions const &d2);

extern Dimensions const Mass, Length, Time, Force, Energy, Temperature;

struct Units {
    std::map<Dimensions, double> scaling; // scaling with respect to SI units.
                                          // Temperature in Kelvin.

    double hbar, kB;
    double force_to_base_units, energy_to_non_base_units;

    void setup(pt::ptree node);

    double scale(Dimensions const &dimension, bool constrain_to_base_units = false) const {
        double val = 1;

        if ((!constrain_to_base_units) && scaling.find(dimension) != scaling.end())
            val = scaling.at(dimension);
        else
            val = std::pow(scaling.at(Mass), dimension.dim.at("mass")) *
                  std::pow(scaling.at(Length), dimension.dim.at("length")) *
                  std::pow(scaling.at(Time), dimension.dim.at("time"));

        return val;
    }

    double specialify(Dimensions const &dimension, double const val) const {
        return val * scale(dimension, true) / scale(dimension);
    }

    double non_specialify(Dimensions const &dimension, double const val) const {
        return val * scale(dimension) / scale(dimension, true);
    }

    double convert_to(Units const u2, Dimensions const &dimension, double const val) const {
        return val * scale(dimension) / u2.scale(dimension);
    }

    arma::vec specialify(Dimensions const &dimension, arma::vec const val) const {
        return val * scale(dimension, true) / scale(dimension);
    }

    arma::vec non_specialify(Dimensions const &dimension, arma::vec const val) const {
        return val * scale(dimension) / scale(dimension, true);
    }

    arma::vec convert_to(Units const u2, Dimensions const &dimension, arma::vec const val) const {
        return val * scale(dimension) / u2.scale(dimension);
    }

    arma::mat specialify(Dimensions const &dimension, arma::mat const val) const {
        return val * scale(dimension, true) / scale(dimension);
    }

    arma::mat non_specialify(Dimensions const &dimension, arma::mat const val) const {
        return val * scale(dimension) / scale(dimension, true);
    }

    arma::mat convert_to(Units const u2, Dimensions const &dimension, arma::mat const val) const {
        return val * scale(dimension) / u2.scale(dimension);
    }
};

extern Units units;

#endif
