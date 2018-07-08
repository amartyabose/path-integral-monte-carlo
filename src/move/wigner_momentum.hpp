#ifndef _WIGNER_MOMENTUM_HPP_
#define _WIGNER_MOMENTUM_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "../random.hpp"
#include "../utilities.hpp"
#include "move.hpp"

class WignerMomentum : public Move {
    double beta, mass, sigma;
public:
    void setup(pt::ptree::value_type node, pt::ptree params) {
        mass = params.get<double>("mass", 1);
        //beta = params.get<double>("beta");
        beta = node.second.get<double>("<xmlattr>.alpha");
        sigma = std::sqrt(mass/beta);
    }

    void operator()(boost::shared_ptr<Configuration> &conf, arma::uvec atom_nums) {
        unsigned dimensions = conf->num_dims();
        WignerConfiguration *wigner_newconf = dynamic_cast<WignerConfiguration*>(conf->duplicate());
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++) {
            for(unsigned dim=0; dim<dimensions; dim++)
                wigner_newconf->set_momentum(atom_nums(atom), dim, random_normal(0, sigma));
        }
        check_amplitude(conf, boost::shared_ptr<Configuration>(wigner_newconf));
    }

    virtual void check_amplitude(boost::shared_ptr<Configuration> &conf_old, boost::shared_ptr<Configuration> conf_new) {
        moves_tried++;

        WignerConfiguration *wigner_newconf = dynamic_cast<WignerConfiguration*>(conf_new.get());
        WignerConfiguration *wigner_oldconf = dynamic_cast<WignerConfiguration*>(conf_old.get());
        arma::mat temp_pold = wigner_oldconf->get_momentum();
        arma::mat temp_pnew = wigner_newconf->get_momentum();
        arma::cube pold = arma::cube(temp_pold.n_rows, temp_pold.n_cols, 1);
        pold.slice(0) = temp_pold;
        arma::cube pnew = arma::cube(temp_pnew.n_rows, temp_pnew.n_cols, 1);
        pnew.slice(0) = temp_pnew;

        double new_weight = 1, old_weight = 1;
        for(unsigned p=0; p<propagator.size(); p++) {
            new_weight *= (*propagator[p])(pnew);
            old_weight *= (*propagator[p])(pold);
        }
        //temp_pnew.print("New p");
        //temp_pold.print("Old p");
        //std::cout<<new_weight<<"\t"<<old_weight<<std::endl;
        if(new_weight/old_weight > random_float(0,1)) {
            conf_old.reset(conf_new->duplicate());
            moves_accepted++;
        }
    }
};

REGISTER_TYPE_GENERAL(WignerMomentum, Move)

#endif
