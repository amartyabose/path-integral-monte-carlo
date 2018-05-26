#ifndef _G2_HPP_
#define _G2_HPP_

#include <boost/shared_ptr.hpp>

#include "../utilities.hpp"
#include "propagator.hpp"

class G2 : public Propagator {
public:
    double operator()(Configuration conf) {
        double total_pot = 0;
        for(unsigned t=0; t<conf.num_beads(); t++)
            total_pot += (*V)(conf.time_slice(t));
        return std::exp(-total_pot * tau);
    }
};

REGISTER_TYPE_GENERAL(G2, Propagator)

#endif
