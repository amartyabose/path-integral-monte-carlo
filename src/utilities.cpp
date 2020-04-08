#include "utilities.hpp"

int my_id, world_size;

namespace utilities {
bool boolparse(std::string const &st) {
    if (iequals(st, "yes") || iequals(st, "on") || iequals(st, "true"))
        return true;
    else
        return false;
}

double exp_series(double var, unsigned n) {
    double poly_expansion = 1;
    double term           = var;
    for (unsigned long i = 1; i <= n; i++) {
        poly_expansion += term;
        term *= var / (i + 1);
    }
    return poly_expansion;
}
} // namespace utilities
