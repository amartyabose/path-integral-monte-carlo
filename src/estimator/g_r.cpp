#include <cmath>
#include <string>

#include "../boundary_conditions/boundary_conditions.hpp"
#include "g_r.hpp"

arma::mat Gr::eval(std::shared_ptr<Configuration> const &x) {
    if (normalization.has_nan()) {
        if (x->num_dims() == 1)
            normalization = arma::ones<arma::vec>(arma::size(Rs));
        else if (x->num_dims() == 2)
            normalization = 2 * arma::datum::pi * bin_size * Rs;
        else if (x->num_dims() == 3)
            normalization = 4 * arma::datum::pi * bin_size * Rs % Rs;
    }

    arma::mat ans = arma::zeros<arma::mat>(n_rows, n_cols);
    ans.col(0)    = Rs;
    if (type == "pimc" || type == "pimd") {
        for (unsigned b = 0; b < x->num_beads(); b++) {
            auto time_slice = x->time_slice(b);
            for (unsigned a1 = 0; a1 < x->num_atoms(); a1++)
                for (unsigned a2 = 0; a2 < a1; a2++) {
                    auto   disp  = bc->wrap_vector(time_slice.row(a1) - time_slice.row(a2));
                    double dist  = arma::norm(disp);
                    int    index = std::round(dist / bin_size);
                    if (index < n_rows)
                        ans(index, 1) += 1;
                }
        }
        ans.col(1) /= x->num_beads() * normalization;
        ans(0, 1) = 0;
    } else
        spdlog::critical("Gr not yet defined for Wigner calculations.");
    return ans;
}
