#ifndef _SK_H_
#define _SK_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class Sk : public Estimator {
protected:
    std::string type;

    bool      isotropic;
    arma::vec axis;

    double    kmax, delta_k;
    arma::vec ks_1d;

    double calc_sk(arma::mat const &x, arma::vec const &ks);
    double calc_sk(std::shared_ptr<Configuration> const &x, arma::vec const &ks);

public:
    void      setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override;
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(Sk, Estimator)

#endif
