#include "potential.hpp"

std::map<std::string, PotentialFactory *> &Potential::get_factory() {
    static std::map<std::string, PotentialFactory *> factory;
    return factory;
}

void Potential::registerType(const std::string &name, PotentialFactory *factory) { get_factory()[name] = factory; }

std::shared_ptr<Potential> Potential::create(const std::string &name) {
    if (get_factory().find(name) == get_factory().end()) {
        std::cerr << "Not a valid potential" << std::endl;
        exit(1);
    }
    return get_factory()[name]->create();
}

double Potential::operator()(arma::mat const &x) {
    double pe = 0.0;
    for (unsigned atom = 0; atom < x.n_rows; atom++)
        pe += (*this)(x, atom);

    return pe;
}

arma::mat Potential::derivative(arma::mat const &x) {
    std::complex<double> epsilon(0, 1e-10);

    double    default_ans = (*this)(x);
    arma::mat ans         = arma::zeros<arma::mat>(arma::size(x));
    arma::mat zeros       = arma::zeros<arma::mat>(arma::size(x));
    for (unsigned r = 0; r < x.n_rows; r++)
        for (unsigned c = 0; c < x.n_cols; c++) {
            arma::cx_mat temp(x, zeros);
            temp(r, c) += epsilon;
            std::complex<double> temp_func = (*this)(temp);

            ans(r, c) = ((temp_func - default_ans) / epsilon).real();
        }
    return ans;
}

std::shared_ptr<Potential> pot;