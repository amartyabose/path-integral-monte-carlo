#ifndef _RANDOM_HPP_
#define _RANDOM_HPP_

#include <boost/random.hpp>
namespace bRand = boost::random;

extern bRand::mt19937 generator;

int    random_integer(int low, int high);
double random_float(double low, double high);
double random_normal(double mean, double sigma);

#endif
