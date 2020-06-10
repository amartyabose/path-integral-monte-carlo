#ifndef _GR_H_
#define _GR_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class Gr : public Estimator {
protected:
    std::string type;

    double    Rmax, bin_size;
    unsigned  nbins;
    arma::vec Rs, normalization;

public:
    void      setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override;
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(Gr, Estimator)

#endif
