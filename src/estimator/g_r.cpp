#include <cmath>
#include <string>

#include "../boundary_conditions/boundary_conditions.hpp"
#include "g_r.hpp"

void Gr::setup(std::string type_, arma::vec mass_, double beta_, unsigned num_beads_, pt::ptree node) {
    type   = type_;
    Rmax   = node.get<double>("<xmlattr>.rmax");
    n_cols = 2;
    n_rows = nbins = node.get<unsigned>("<xmlattr>.nbins", 20);

    Rs       = arma::linspace(0, Rmax, nbins);
    bin_size = Rs(1) - Rs(0);
    normalization =
        arma::datum::nan * arma::ones<arma::vec>(arma::size(Rs)); // To be constructed the first time eval() is called.

    spdlog::warn("\t\tYou are calculating g(r). The values printed would not be normalized with the number density. To "
                 "obtain the correct g(r), divide the output by the number density --- N/V.");
}

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
                        ans(index, 1) += 2; // one for a1-a2 interaction and one for a2-a1 interaction.
                }
        }
        ans.col(1) /= x->num_beads() * x->num_atoms() *
                      normalization; // The answer needs to be divided by the number density (N/V).
        ans(0, 1) = 0;
    } else {
        arma::mat pos = x->pos();
        for (unsigned a1 = 0; a1 < x->num_atoms(); a1++)
            for (unsigned a2 = 0; a2 < a1; a2++) {
                auto   disp  = bc->wrap_vector(pos.row(a1) - pos.row(a2));
                double dist  = arma::norm(disp);
                int    index = std::round(dist / bin_size);
                if (index < n_rows)
                    ans(index, 1) += 2; // one for a1-a2 interaction and one for a2-a1 interaction.
            }
        ans.col(1) /= x->num_atoms() * normalization; // The answer needs to be divided by the number density (N/V).
        ans(0, 1) = 0;
    }
    return ans;
}
