#include "random.hpp"
bRand::mt19937 generator;

int random_integer(int low, int high) {
    bRand::uniform_int_distribution<> dist(low, high);
    return dist(generator);
}

double random_float(double low, double high) {
    bRand::uniform_real_distribution<double> dist(low, high);
    return dist(generator);
}

double random_normal(double mean, double sigma) {
    bRand::normal_distribution<double> dist(mean, sigma);
    return dist(generator);
}
