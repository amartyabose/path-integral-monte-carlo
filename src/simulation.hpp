#ifndef _SIM_H_
#define _SIM_H_

#include <memory>
#include <string>

#include <armadillo>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "boundary_conditions/boundary_conditions.hpp"
#include "configuration.hpp"
#include "estimator/estimator.hpp"
#include "move/move.hpp"
#include "output.hpp"
#include "pi_wigner_config.hpp"
#include "potential/potential.hpp"
#include "propagator/propagator.hpp"
#include "units.hpp"
#include "utilities.hpp"

class Simulation {
    int my_id, world_size;

    std::string type, initial_config;

    unsigned  ndimensions, natoms, nprops;
    double    beta;
    arma::vec mass, atom_specific_beta; // atom_specific_beta is different for
                                        // atoms only for Wigner calculations
    std::vector<double> dt_imaginary;

    // for PIMD simulations
    arma::mat   bead_specific_mass;
    double      dt;
    unsigned    burn, skip;
    std::string thermostat_order;

    std::vector<std::shared_ptr<Move>> moves;

    std::vector<std::string> move_names;

    unsigned nIC, nblocks;

    std::vector<unsigned> bead_nums;

    std::vector<std::shared_ptr<Propagator>> props;

    Output output;

    void initialize_output(pt::ptree output_node);
    void setup_output(pt::ptree output_node);
    void setup_logging(boost::optional<pt::ptree &> logging_node);
    void insert_position_momentum_moves(pt::ptree::value_type p);
    void get_propagator(pt::ptree prop);
    void get_dts(pt::ptree prop);
    void get_moves(pt::ptree move_params);
    void run_block(std::shared_ptr<Configuration> &conf);

    virtual void get_parameters(pt::ptree params);
    virtual void get_MC_params(pt::ptree params);
    virtual void get_potential(pt::ptree pot_tree);

public:
    Simulation(int myid, int ws) : my_id(myid), world_size(ws) {}
    void load(pt::ptree tree);
    void run();
};

#endif
