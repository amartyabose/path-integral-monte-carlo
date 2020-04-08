#ifndef _POS_CORR_H_
#define _POS_CORR_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class PositionCorrelation : public Estimator {
    double   beta;
    unsigned num_beads;

public:
    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        if (type_ == "wigner")
            throw std::runtime_error("Imaginary time correlations not available with Wigner calculations.");
        beta      = beta_;
        num_beads = num_beads_;
        n_rows    = num_beads;
        n_cols    = 2;
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(PositionCorrelation, Estimator)

#endif