#include "G2.hpp"
#include "G4.hpp"
#include "G6.hpp"
#include "G8.hpp"

double G2::operator()(const arma::cube &conf) {
    double total_amplitude = -tau / 2. *
                             ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                              (*V)(conf(arma::span::all, arma::span(1), arma::span::all)));
    return std::exp(total_amplitude);
}

arma::cube G2::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    ans(arma::span::all, arma::span(0), arma::span::all) =
        -tau / 2. * V->derivative(conf(arma::span::all, arma::span(0), arma::span::all));
    ans(arma::span::all, arma::span(1), arma::span::all) =
        -tau / 2. * V->derivative(conf(arma::span::all, arma::span(1), arma::span::all));
    return ans;
}

double G4::operator()(const arma::cube &conf) {
    double total_amplitude = 1;
    total_amplitude *= 4. / 3. *
                           std::exp(-tau / 2. *
                                    ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                     2. * (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) +
                                     (*V)(conf(arma::span::all, arma::span(2), arma::span::all)))) -
                       1. / 3. *
                           std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(2), arma::span::all))));
    return total_amplitude;
}

arma::cube G4::derivative(const arma::cube &conf) {
    arma::cube ans    = arma::zeros<arma::cube>(arma::size(conf));
    double     factor = 1. / (*this)(conf);
    std::cout << "factor = " << factor << std::endl;

    ans(arma::span::all, arma::span(0), arma::span::all) =
        4. / 3. *
            std::exp(-tau / 2. *
                     ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                      2. * (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) +
                      (*V)(conf(arma::span::all, arma::span(2), arma::span::all)))) *
            (-tau / 2) * V->derivative(conf(arma::span::all, arma::span(0), arma::span::all)) -
        1. / 3. *
            std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                             (*V)(conf(arma::span::all, arma::span(2), arma::span::all)))) *
            (-tau) * V->derivative(conf(arma::span::all, arma::span(0), arma::span::all));

    ans(arma::span::all, arma::span(0), arma::span::all) =
        4. / 3. *
        std::exp(-tau / 2. *
                 ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                  2. * (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) +
                  (*V)(conf(arma::span::all, arma::span(2), arma::span::all)))) *
        (-tau) * V->derivative(conf(arma::span::all, arma::span(0), arma::span::all));

    ans(arma::span::all, arma::span(2), arma::span::all) =
        4. / 3. *
            std::exp(-tau / 2. *
                     ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                      2. * (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) +
                      (*V)(conf(arma::span::all, arma::span(2), arma::span::all)))) *
            (-tau / 2) * V->derivative(conf(arma::span::all, arma::span(2), arma::span::all)) -
        1. / 3. *
            std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                             (*V)(conf(arma::span::all, arma::span(2), arma::span::all)))) *
            (-tau) * V->derivative(conf(arma::span::all, arma::span(2), arma::span::all));

    return ans * factor;
}

double G6::operator()(const arma::cube &conf) {
    double total_amplitude = 1;
    total_amplitude *= 64. / 45. *
                           std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) / 2. +
                                            (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(2), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(3), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(4), arma::span::all)) / 2.)) -
                       4. / 9. *
                           std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                            2. * (*V)(conf(arma::span::all, arma::span(2), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(4), arma::span::all)))) +
                       1. / 45. *
                           std::exp(-2 * tau *
                                    ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                     (*V)(conf(arma::span::all, arma::span(4), arma::span::all))));
    return total_amplitude;
}

arma::cube G6::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    return ans;
}

double G8::operator()(const arma::cube &conf) {
    double total_amplitude = 1;
    total_amplitude *= 54. / 35. *
                           std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) / 2. +
                                            (*V)(conf(arma::span::all, arma::span(1), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(2), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(3), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(4), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(5), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(6), arma::span::all)) / 2.)) -
                       27. / 40. *
                           std::exp(-tau * ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                            2. * (*V)(conf(arma::span::all, arma::span(2), arma::span::all)) +
                                            2. * (*V)(conf(arma::span::all, arma::span(4), arma::span::all)) +
                                            (*V)(conf(arma::span::all, arma::span(6), arma::span::all)))) +
                       2. / 15. *
                           std::exp(-tau * (1.5 * (*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                            3 * (*V)(conf(arma::span::all, arma::span(3), arma::span::all)) +
                                            1.5 * (*V)(conf(arma::span::all, arma::span(6), arma::span::all)))) -
                       1. / 840. *
                           std::exp(-3 * tau *
                                    ((*V)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                     (*V)(conf(arma::span::all, arma::span(6), arma::span::all))));
    return total_amplitude;
}

arma::cube G8::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    return ans;
}
