#ifndef _CV_H_
#define _CV_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class CV : public Estimator {
public:
    std::string type;
    double      beta;

    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        type   = type_;
        beta   = beta_;
        n_cols = 1;
        n_rows = 2;
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(CV, Estimator)

#endif
