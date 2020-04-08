#include "virial_kinetic_energy.hpp"

arma::mat VirialKineticEnergy::eval(std::shared_ptr<Configuration> const &x) {
    double    ans = 0;
    double    val = 0;
    arma::mat com = arma::zeros<arma::mat>(x->num_atoms(), x->num_dims());
    for (unsigned i = 0; i < x->num_beads(); i++)
        com += x->time_slice(i);
    com /= x->num_beads();
    for (unsigned i = 0; i < x->num_beads(); i++) {
        arma::mat pos   = x->time_slice(i);
        arma::mat deriv = pot->derivative(pos);
        val += arma::accu((x->time_slice(i) - com) % deriv);
    }

    ans = x->num_dims() * x->num_atoms() / (2. * beta) + val / (2 * x->num_beads());
    return ans * units.energy_to_non_base_units * arma::ones<arma::mat>(n_rows, n_cols);
}
