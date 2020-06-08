#include "open_boundary_conditions.hpp"

void OpenBoundaryConditions::setup(pt::ptree node) {}

arma::rowvec OpenBoundaryConditions::wrap_rowvector(arma::rowvec const &x) { return x; }

arma::vec OpenBoundaryConditions::wrap_vector(arma::vec const &x) { return x; }

arma::mat OpenBoundaryConditions::wrap_coordinates(arma::mat const &x) { return x; }

arma::mat OpenBoundaryConditions::center_box(arma::mat const &x) { return x; }
