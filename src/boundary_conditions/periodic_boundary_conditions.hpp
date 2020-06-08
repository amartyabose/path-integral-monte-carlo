#ifndef _PBC_HPP_
#define _PBC_HPP_

#include "../utilities.hpp"
#include "boundary_conditions.hpp"

class PeriodicBoundaryConditions : public BoundaryConditions {
protected:
    arma::vec box_size;

public:
    void         setup(pt::ptree node) override;
    arma::mat    center_box(arma::mat const &x) override;
    arma::mat    wrap_coordinates(arma::mat const &x) override;
    arma::rowvec wrap_rowvector(arma::rowvec const &x) override;
    arma::vec    wrap_vector(arma::vec const &x) override;
};

REGISTER_TYPE_GENERAL(PeriodicBoundaryConditions, BoundaryConditions)

#endif
