#ifndef _PE_H_
#define _PE_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class PotentialEnergy : public Estimator {
public:
    std::string type;

    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        type   = type_;
        n_rows = n_cols = 1;
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(PotentialEnergy, Estimator)

#endif
