#include <string>

#include "potential_energy.hpp"

arma::mat PotentialEnergy::eval(std::shared_ptr<Configuration> const &x) {
    double ans = 0;
    if (type == "wigner")
        ans = (*pot)(x->pos());
    else {
        for (unsigned b = 0; b < x->num_beads(); b++)
            ans += (*pot)(x->time_slice(b));
        ans /= x->num_beads();
    }
    return ans * units.energy_to_non_base_units * arma::ones<arma::mat>(n_rows, n_cols);
}
