#ifndef _G4_HPP_
#define _G4_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G4 : public Propagator {
public:
    double operator()(Configuration conf) {
        double total_amplitude = 1;
        for(unsigned t=0; t<conf.num_beads()-1; t+=2) {
            total_amplitude *=
                4./3. * std::exp(-tau/2. * ((*V)(conf.time_slice(t)) + 2. * (*V)(conf.time_slice(t+1)) + (*V)(conf.time_slice((t+2)%conf.num_beads()))))
                - 1./3. * std::exp(-tau * ((*V)(conf.time_slice(t)) + (*V)(conf.time_slice((t+2)%conf.num_beads()))));
        }
        return total_amplitude;
    }
};

REGISTER_TYPE_GENERAL(G4, Propagator)

#endif
