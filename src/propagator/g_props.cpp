#include "G2.hpp"
#include "G4.hpp"
#include "G6.hpp"
#include "G8.hpp"

double G2::operator()(const arma::cube &conf) {
    double total_amplitude = std::exp(-tau / 2. *
                                      ((*pot)(conf(arma::span::all, arma::span(0), arma::span::all)) +
                                       (*pot)(conf(arma::span::all, arma::span(1), arma::span::all))));
    return total_amplitude;
}

double G2::operator()(const arma::cube &conf, unsigned index) {
    double total_amplitude = std::exp(-tau / 2. *
                                      ((*pot)(conf(arma::span::all, arma::span(0), arma::span::all), index) +
                                       (*pot)(conf(arma::span::all, arma::span(1), arma::span::all), index)));
    return total_amplitude;
}

arma::cube G2::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    ans(arma::span::all, arma::span(0), arma::span::all) =
        -tau / 2. * pot->derivative(conf(arma::span::all, arma::span(0), arma::span::all));
    ans(arma::span::all, arma::span(1), arma::span::all) =
        -tau / 2. * pot->derivative(conf(arma::span::all, arma::span(1), arma::span::all));
    return ans;
}

double G4::operator()(const arma::cube &conf) {
    double V1 = (*pot)(conf(arma::span::all, arma::span(0), arma::span::all));
    double V2 = (*pot)(conf(arma::span::all, arma::span(1), arma::span::all));
    double V3 = (*pot)(conf(arma::span::all, arma::span(2), arma::span::all));

    double total_amplitude =
        4. / 3 * std::exp(-tau / 2 * (V1 + 2 * V2 + V3)) * (1 - 1. / 4 * std::exp(-tau / 2. * (V1 + 2 * V2 + V3)));
    return total_amplitude;
}

double G4::operator()(const arma::cube &conf, unsigned index) {
    double V1 = (*pot)(conf(arma::span::all, arma::span(0), arma::span::all), index);
    double V2 = (*pot)(conf(arma::span::all, arma::span(1), arma::span::all), index);
    double V3 = (*pot)(conf(arma::span::all, arma::span(2), arma::span::all), index);

    double total_amplitude =
        4. / 3 * std::exp(-tau / 2 * (V1 + 2 * V2 + V3)) * (1 - 1. / 4 * std::exp(-tau / 2. * (V1 + 2 * V2 + V3)));
    return total_amplitude;
}

arma::cube G4::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    return ans;
}

double G6::operator()(const arma::cube &conf) {
    double V1 = (*pot)(conf(arma::span::all, arma::span(0), arma::span::all));
    double V2 = (*pot)(conf(arma::span::all, arma::span(1), arma::span::all));
    double V3 = (*pot)(conf(arma::span::all, arma::span(2), arma::span::all));
    double V4 = (*pot)(conf(arma::span::all, arma::span(3), arma::span::all));
    double V5 = (*pot)(conf(arma::span::all, arma::span(4), arma::span::all));

    double total_amplitude = 64. / 45 * std::exp(-tau * (V1 / 2. + V2 + V3 + V4 + V5 / 2.)) -
                             4. / 9 * std::exp(-tau * (V1 + 2 * V3 + V5)) + 1. / 45 * std::exp(-2 * tau * (V1 + V5));
    return total_amplitude;
}

double G6::operator()(const arma::cube &conf, unsigned index) {
    double V1 = (*pot)(conf(arma::span::all, arma::span(0), arma::span::all), index);
    double V2 = (*pot)(conf(arma::span::all, arma::span(1), arma::span::all), index);
    double V3 = (*pot)(conf(arma::span::all, arma::span(2), arma::span::all), index);
    double V4 = (*pot)(conf(arma::span::all, arma::span(3), arma::span::all), index);
    double V5 = (*pot)(conf(arma::span::all, arma::span(4), arma::span::all), index);

    double total_amplitude = 64. / 45 * std::exp(-tau * (V1 / 2. + V2 + V3 + V4 + V5 / 2.)) -
                             4. / 9 * std::exp(-tau * (V1 + 2 * V3 + V5)) + 1. / 45 * std::exp(-2 * tau * (V1 + V5));
    return total_amplitude;
}

arma::cube G6::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    return ans;
}

double G8::operator()(const arma::cube &conf) {
    double V1 = (*pot)(conf(arma::span::all, arma::span(0), arma::span::all));
    double V2 = (*pot)(conf(arma::span::all, arma::span(1), arma::span::all));
    double V3 = (*pot)(conf(arma::span::all, arma::span(2), arma::span::all));
    double V4 = (*pot)(conf(arma::span::all, arma::span(3), arma::span::all));
    double V5 = (*pot)(conf(arma::span::all, arma::span(4), arma::span::all));
    double V6 = (*pot)(conf(arma::span::all, arma::span(5), arma::span::all));
    double V7 = (*pot)(conf(arma::span::all, arma::span(6), arma::span::all));

    double total_amplitude = 54. / 35 * std::exp(-tau * (V1 / 2 + V2 + V3 + V4 + V5 + V6 + V7 / 2)) -
                             27. / 40 * std::exp(-tau * (V1 + 2 * V3 + 2 * V5 + V7)) +
                             2. / 15 * std::exp(-3 * tau / 2 * (V1 + 2 * V4 + V7)) -
                             1. / 840 * std::exp(-3 * tau * (V1 + V7));
    return total_amplitude;
}

double G8::operator()(const arma::cube &conf, unsigned index) {
    double V1 = (*pot)(conf(arma::span::all, arma::span(0), arma::span::all), index);
    double V2 = (*pot)(conf(arma::span::all, arma::span(1), arma::span::all), index);
    double V3 = (*pot)(conf(arma::span::all, arma::span(2), arma::span::all), index);
    double V4 = (*pot)(conf(arma::span::all, arma::span(3), arma::span::all), index);
    double V5 = (*pot)(conf(arma::span::all, arma::span(4), arma::span::all), index);
    double V6 = (*pot)(conf(arma::span::all, arma::span(5), arma::span::all), index);
    double V7 = (*pot)(conf(arma::span::all, arma::span(6), arma::span::all), index);

    double total_amplitude = 54. / 35 * std::exp(-tau * (V1 / 2 + V2 + V3 + V4 + V5 + V6 + V7 / 2)) -
                             27. / 40 * std::exp(-tau * (V1 + 2 * V3 + 2 * V5 + V7)) +
                             2. / 15 * std::exp(-3 * tau / 2 * (V1 + 2 * V4 + V7)) -
                             1. / 840 * std::exp(-3 * tau * (V1 + V7));
    return total_amplitude;
}

arma::cube G8::derivative(const arma::cube &conf) {
    arma::cube ans = arma::zeros<arma::cube>(arma::size(conf));
    return ans;
}
