#ifndef _G6_HPP_
#define _G6_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G6 : public Propagator {
public:
    double operator()(Configuration conf) {
        double total_amplitude = 1;
        for(unsigned t=0; t<conf.num_beads()-3; t+=4) {
            total_amplitude *=
                64./45. * std::exp(-tau * ((*V)(conf.time_slice(t))/2. + (*V)(conf.time_slice(t+1)) + (*V)(conf.time_slice(t+2)) + (*V)(conf.time_slice(t+3)) + (*V)(conf.time_slice((t+4)%conf.num_beads()))/2.))
                - 4./9. * std::exp(-tau * ((*V)(conf.time_slice(t)) + 2.*(*V)(conf.time_slice(t+2)) + (*V)(conf.time_slice((t+4)%conf.num_beads()))))
                + 1./45. * std::exp(-2*tau * ((*V)(conf.time_slice(t)) + (*V)(conf.time_slice((t+4)%conf.num_beads()))));
        }
        return total_amplitude;
    }
};

REGISTER_TYPE_GENERAL(G6, Propagator)

#endif
