#ifndef _G8_HPP_
#define _G8_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G8 : public Propagator {
public:
    double operator()(Configuration conf) {
        double total_amplitude = 1;
        for(unsigned t=0; t<conf.num_beads()-5; t+=6) {
            total_amplitude *=
                54./35. * std::exp(-tau * ((*V)(conf.time_slice(t))/2. + (*V)(conf.time_slice(t+1)) + (*V)(conf.time_slice(t+2)) + (*V)(conf.time_slice(t+3)) + (*V)(conf.time_slice(t+4)) + (*V)(conf.time_slice(t+5)) + (*V)(conf.time_slice((t+6)%conf.num_beads()))/2.))
                - 27./40. * std::exp(-tau * ((*V)(conf.time_slice(t)) + 2.*(*V)(conf.time_slice(t+2)) + 2.*(*V)(conf.time_slice(t+4)) + (*V)(conf.time_slice((t+6)%conf.num_beads()))))
                + 2./15. * std::exp(-tau * (1.5*(*V)(conf.time_slice(t)) + 3*(*V)(conf.time_slice(t+3)) + 1.5*(*V)(conf.time_slice((t+6)%conf.num_beads()))))
                - 1./840. * std::exp(-3 * tau * ((*V)(conf.time_slice(t)) + (*V)(conf.time_slice((t+6)%conf.num_beads()))));
        }
        return total_amplitude;
    }
};

REGISTER_TYPE_GENERAL(G8, Propagator)

#endif
