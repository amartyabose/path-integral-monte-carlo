#ifndef _KE_H_
#define _KE_H_

#include "../configuration.hpp"
#include "../utilities.hpp"
#include "estimator.hpp"

class KineticEnergy : public Estimator {
public:
    double      beta;
    std::string type;

    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        type   = type_;
        beta   = beta_;
        mass   = mass_;
        n_rows = n_cols = 1;
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(KineticEnergy, Estimator)

#endif
