#include "position_correlation.hpp"

arma::mat PositionCorrelation::eval(std::shared_ptr<Configuration> const &x) {
    arma::mat ans = arma::zeros<arma::mat>(n_rows, n_cols);
    double    tau = beta / (num_beads - 1);
    for (unsigned i = 0; i < num_beads; i++) {
        ans(i, 0) = i * tau;
        ans(i, 1) = arma::dot(x->time_slice(0), x->time_slice(i));
    }
    return ans;
}