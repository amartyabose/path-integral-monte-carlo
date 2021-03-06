#include "kinetic_energy.hpp"
#include "pi_wigner_config.hpp"

arma::mat KineticEnergy::eval(std::shared_ptr<Configuration> const &x) {
    double ans = 0;
    if (type == "wigner") {
        WignerConfiguration *conf = dynamic_cast<WignerConfiguration *>(x.get());
        ans                       = arma::accu(conf->get_momentum() % conf->get_momentum() / conf->get_mass()) / 2;
    } else {
        double val = 0;
        double tau = beta / x->num_beads();
        for (unsigned atom = 0; atom < x->num_atoms(); atom++)
            for (unsigned i = 0; i < x->num_beads(); i++) {
                arma::vec temp = x->bead_position(atom, i) - x->bead_position(atom, (i + 1) % x->num_beads());
                val += mass(atom) * arma::dot(temp, temp);
            }
        ans = (x->num_dims() * x->num_atoms() / (2. * tau) - val / (2 * x->num_beads() * tau * tau));
    }
    return ans * units.energy_to_non_base_units * arma::ones<arma::mat>(n_rows, n_cols);
}
