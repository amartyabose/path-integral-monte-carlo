#include <string>

#include "cv.hpp"

arma::mat CV::eval(std::shared_ptr<Configuration> const &x) {
    arma::mat ans = arma::zeros<arma::mat>(n_rows, n_cols);
    for (unsigned b = 0; b < x->num_beads(); b++) {
        double potential = (*pot)(x->time_slice(b));
        ans(0, 0) += potential * potential;
        ans(1, 0) += potential;
    }
    ans /= x->num_beads();
    ans(1, 0) *= beta;
    ans(0, 0) *= beta * beta;
    ans(0, 0) += x->num_dims() * x->num_atoms() / 2.0;
    return ans;
}
