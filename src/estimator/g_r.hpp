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
    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        type   = type_;
        Rmax   = node.get<double>("<xmlattr>.rmax");
        n_cols = 2;
        n_rows = nbins = node.get<unsigned>("<xmlattr>.nbins", 20);

        Rs            = arma::linspace(0, Rmax, nbins);
        bin_size      = Rs(1) - Rs(0);
        normalization = arma::datum::nan *
                        arma::ones<arma::vec>(arma::size(Rs)); // To be constructed the first time eval() is called.
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(Gr, Estimator)

#endif
