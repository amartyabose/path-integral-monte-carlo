#ifndef _VKE_H_
#define _VKE_H_

#include "../configuration.hpp"
#include "../utilities.hpp"
#include "estimator.hpp"

class VirialKineticEnergy : public Estimator {
public:
    double      beta;
    std::string type;

    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        if (type_ == "wigner")
            throw std::runtime_error("Virial observables are not available with Wigner calculations.");
        type   = type_;
        beta   = beta_;
        n_rows = n_cols = 1;
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(VirialKineticEnergy, Estimator)

#endif
