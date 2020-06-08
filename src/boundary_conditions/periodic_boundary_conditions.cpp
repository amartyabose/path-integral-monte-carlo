#include "periodic_boundary_conditions.hpp"

void PeriodicBoundaryConditions::setup(pt::ptree node) { box_size = arma::vec(node.get<std::string>("box_size")); }

arma::rowvec PeriodicBoundaryConditions::wrap_rowvector(arma::rowvec const &x) {
    arma::vec new_x = x;
    for (unsigned d = 0; d < x.n_cols; d++)
        if (new_x(d) > box_size(d) / 2)
            new_x(d) -= box_size(d);
        else if (new_x(d) < -box_size(d) / 2)
            new_x(d) += box_size(d);
    return new_x;
}

arma::vec PeriodicBoundaryConditions::wrap_vector(arma::vec const &x) {
    arma::vec new_x = x;
    for (unsigned d = 0; d < x.n_cols; d++)
        if (new_x(d) > box_size(d) / 2)
            new_x(d) -= box_size(d);
        else if (new_x(d) < -box_size(d) / 2)
            new_x(d) += box_size(d);
    return new_x;
}

arma::mat PeriodicBoundaryConditions::wrap_coordinates(arma::mat const &x) {
    auto new_x = x;
    auto com   = arma::mean(x, 0);
    for (unsigned a = 0; a < x.n_rows; a++)
        new_x.row(a) = wrap_vector(new_x.row(a) - com);
    return new_x;
}

arma::mat PeriodicBoundaryConditions::center_box(arma::mat const &x) {
    arma::mat new_x = x;
    auto      com   = arma::mean(x, 0);
    for (unsigned a = 0; a < x.n_rows; a++)
        new_x.row(a) = wrap_vector(new_x.row(a) - com) + box_size / 2;
    return new_x;
}