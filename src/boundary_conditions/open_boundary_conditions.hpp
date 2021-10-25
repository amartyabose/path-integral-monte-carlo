#ifndef _OBC_HPP_
#define _OBC_HPP_

#include "../utilities.hpp"
#include "boundary_conditions.hpp"

class OpenBoundaryConditions : public BoundaryConditions {
public:
    void         setup(pt::ptree node) override;
    arma::mat    center_box(arma::mat const &x) const override;
    arma::mat    wrap_coordinates(arma::mat const &x) const override;
    arma::rowvec wrap_rowvector(arma::rowvec const &x) const override;
    arma::vec    wrap_vector(arma::vec const &x) const override;
};

REGISTER_TYPE_GENERAL(OpenBoundaryConditions, BoundaryConditions)

#endif
