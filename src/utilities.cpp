#include "utilities.hpp"

double exp_series(double var, unsigned n) {
    double poly_expansion = 1;
    double term = var;
    for(unsigned long i=1; i<=n; i++) {
        poly_expansion += term;
        term *= var/(i+1);
    }
    return poly_expansion;
}
