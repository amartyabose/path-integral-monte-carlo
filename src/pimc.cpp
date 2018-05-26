#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::vector;
using std::string;

#include <armadillo>

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>
namespace pt = boost::property_tree;
using boost::shared_ptr;

#include "random.hpp"
#include "configuration.hpp"
#include "potentials.hpp"
#include "propagators.hpp"
#include "moves.hpp"

struct Estimator {
    virtual double operator() (arma::mat x) = 0;
};

struct PE_estimator : public Estimator {
    shared_ptr<Potential> pot;
    PE_estimator() {}
    PE_estimator(shared_ptr<Potential> V) : pot(V) {}
    double operator()(arma::mat x) {
        return (*pot)(x);
    }
};

class Simulation {
    unsigned ndimensions, natoms, nbeads;
    double beta;

    vector<shared_ptr<Move> > moves;
    vector<string> move_names;

    unsigned nMC, nblocks;
    shared_ptr<Potential> pot;
    shared_ptr<Propagator> propagator;

    shared_ptr<Estimator> estimator;

    virtual void get_parameters(pt::ptree params) {
        beta = params.get<double>("beta");
        ndimensions = params.get<unsigned>("num_dimensions");
        natoms = params.get<unsigned>("num_atoms");
        nbeads = params.get<unsigned>("num_beads");
    }

    virtual void get_MC_params(pt::ptree params) {
        nMC = params.get<unsigned>("num_MC");
        nblocks = params.get<unsigned>("num_blocks");
    }

    virtual void get_system(pt::ptree sys) {
        get_potential(sys.get_child("potential"));
        if(sys.get<string>("estimator") == "potential")
            estimator = shared_ptr<Estimator>(new PE_estimator(pot));
    }

    virtual void get_potential(pt::ptree pot_tree) {
        pot = Potential::create(pot_tree.get<string>("<xmlattr>.name"));
        pot->setup(pot_tree);
    }

    void get_propagator(pt::ptree prop) {
        propagator = Propagator::create(prop.get<string>("<xmlattr>.name"));
        propagator->V = pot;
        propagator->tau = beta/nbeads;
    }

    void get_moves(pt::ptree move_params, pt::ptree params) {
        BOOST_FOREACH(const pt::ptree::value_type &m, move_params) {
            moves.push_back(Move::create(m.first));
            moves.back()->setup(m, params);
            moves.back()->propagator = propagator;
            move_names.push_back(m.first);
        }
    }

public:
    void load(const char *fname) {
        pt::ptree tree;
        read_xml(fname, tree);
        get_parameters(tree.get_child("pimc.parameters"));
        get_MC_params(tree.get_child("pimc.monte_carlo"));
        get_system(tree.get_child("pimc.system"));
        get_propagator(tree.get_child("pimc.propagator"));
        get_moves(tree.get_child("pimc.moves"), tree.get_child("pimc.parameters"));
    }

    void run_MC() {
        Configuration c(natoms, nbeads, ndimensions);
        std::ofstream ofs("test");
        arma::vec moves_accepted = arma::zeros<arma::vec>(moves.size());
        arma::vec moves_tried = arma::zeros<arma::vec>(moves.size());
        arma::vec est = arma::zeros<arma::vec>(nblocks);
        arma::vec xsq = arma::zeros<arma::vec>(nblocks);
        arma::uvec atoms = arma::regspace<arma::uvec>(0, natoms-1);
        for(unsigned b=0; b<nblocks; b++)
            for(unsigned i=0; i<nMC/nblocks; i++) {
                for(unsigned m=0; m<moves.size(); m++)
                    c = (*moves[m])(c, atoms);
                xsq(b) += c.time_slice(0)(0)*c.time_slice(0)(0);
                //ofs<<c.time_slice(0).st()<<std::endl;
                est(b) += (*estimator)(c.time_slice(0));
                ofs<<(*estimator)(c.time_slice(0))<<std::endl;
            }
        est /= nMC/nblocks;
        xsq /= nMC/nblocks;
        std::cout<<arma::mean(est)<<'\t'<<arma::stddev(est)<<std::endl;
        std::cout<<arma::mean(xsq)<<'\t'<<arma::stddev(xsq)<<std::endl;
        std::cout<<"\nFraction of moves accepted: \n";
        for(unsigned m=0; m<move_names.size(); m++)
            std::cout<<move_names[m]<<'\t'<<(double)moves[m]->moves_accepted/moves[m]->moves_tried<<std::endl;
    }
};


int main(int argc, char **argv) {
    if(argc<2) {
        std::cerr<<"Need an XML parameter file!"<<std::endl;
        exit(1);
    }
    Simulation sim;
    sim.load(argv[1]);
    sim.run_MC();
    return 0;
}
