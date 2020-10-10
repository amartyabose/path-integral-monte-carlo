#ifndef _END_DIST_H_
#define _END_DIST_H_

#include "../utilities.hpp"
#include "estimator.hpp"

class End2EndDist : public Estimator {
public:
    std::string type;
    unsigned    atom_num;

    void setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) override {
        type = type_;
        if (type == "wigner")
            throw std::runtime_error("End2EndDist not defined for Wigner calculations.");
        n_rows   = 1;
        n_cols   = node.get<unsigned>("<xmlattr>.dims", 3);
        atom_num = node.get<unsigned>("<xmlattr>.atom_num");
    }
    arma::mat eval(std::shared_ptr<Configuration> const &x) override;
};

REGISTER_TYPE_GENERAL(End2EndDist, Estimator)

#endif
