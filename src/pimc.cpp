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

#include <mpi.h>

#include "random.hpp"
#include "potentials.hpp"
#include "moves.hpp"
#include "propagators.hpp"
#include "configuration.hpp"
#include "pi_wigner_config.hpp"

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
    int my_id, world_size;

    std::string type, output_name, initial_config;

    unsigned ndimensions, natoms, nprops;
    double beta, mass;

    vector<shared_ptr<Move> > moves;
    vector<string> move_names;

    unsigned nMC, nblocks;
    shared_ptr<Potential> pot;
    std::vector<boost::shared_ptr<Propagator> > props;
    std::vector<unsigned> bead_nums;

    shared_ptr<Estimator> estimator;

    virtual void get_parameters(pt::ptree params) {
        beta = params.get<double>("beta");
        ndimensions = params.get<unsigned>("num_dimensions");
        natoms = params.get<unsigned>("num_atoms");
        nprops = params.get<unsigned>("num_propagators");
        mass = params.get<double>("mass", 1);
        initial_config = params.get<std::string>("initial_configuration", "");
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

    std::vector<double> get_dts(pt::ptree prop) {
        //calculate the respective dt's
        std::vector<double> dt;
        double total_temp_accounted_for = 0;
        unsigned undetermined_dt = 0;
        BOOST_FOREACH(const pt::ptree::value_type &p, prop) {
            dt.push_back(p.second.get<double>("<xmlattr>.dt", 0));
            total_temp_accounted_for += dt.back();
            if (dt.back()<1e-10)
                undetermined_dt += p.second.get<double>("<xmlattr>.number");
        }
        if (std::abs(total_temp_accounted_for - beta)>1e-10 && undetermined_dt==0) {
            if (!my_id)
                std::cerr<<"Full beta not accounted for in propagators"<<std::endl;
            MPI_Finalize();
            exit(1);
        }
        double rem_dt = (beta - total_temp_accounted_for)/undetermined_dt;
        for(unsigned i=0; i<dt.size(); i++)
            if (dt[i]<1e-10)
                dt[i] = rem_dt;
        return dt;
    }

    void insert_position_momentum_moves(pt::ptree::value_type p, pt::ptree params) {
        moves.push_back(Move::create("WignerPosition"));
        moves.back()->setup(p, params);
        move_names.push_back("WignerPosition");

        boost::shared_ptr<Propagator> keprop(Propagator::create("WignerKE"));
        keprop->set_params(pot, 0., p, params);

        moves.push_back(Move::create("WignerMomentum"));
        moves.back()->setup(p, params);
        moves.back()->propagator.push_back(keprop);
        move_names.push_back("WignerMomentum");
    }

    void get_propagator(pt::ptree prop, pt::ptree param) {
        std::vector<double> dt = get_dts(prop);
        bead_nums.push_back(0);
        unsigned prop_num = 0;
        bool wigner_prop = false;
        BOOST_FOREACH(const pt::ptree::value_type &p, prop) {
            props.push_back(Propagator::create(p.second.get<string>("<xmlattr>.name")));
            props.back()->set_params(pot, dt[prop_num], p, param);
            if(p.second.get<string>("<xmlattr>.name") == "WignerG2") {
                if(type != "wigner") {
                    if(!my_id)
                        std::cerr<<"WignerG2 propagator allowed only in a Wigner calculation."<<std::endl;
                    MPI_Finalize();
                    exit(1);
                } else if (type == "wigner"){
                    if(!wigner_prop) {
                        if(p.second.get<unsigned>("<xmlattr>.number")>1) {
                            if(!my_id)
                                std::cerr<<"Only 1 WignerG2 propagator allowed in a Wigner calculation."<<std::endl;
                            MPI_Finalize();
                            exit(1);
                        }
                        wigner_prop = true;
                        insert_position_momentum_moves(p, param);
                        beta = beta - p.second.get<double>("<xmlattr>.dt") + props.back()->get_tau();
                    } else {
                        if(!my_id)
                            std::cerr<<"Only 1 WignerG2 propagator allowed in a Wigner calculation."<<std::endl;
                        MPI_Finalize();
                        exit(1);
                    }
                }
            }
            prop_num++;
            bead_nums.push_back(bead_nums.back() + props.back()->total_beads() - 1);
            unsigned n_props = p.second.get<unsigned>("<xmlattr>.number");
            for (unsigned i=1; i<n_props; i++) {
                props.push_back(props.back());
                bead_nums.push_back(bead_nums.back() + props.back()->total_beads() - 1);
            }
        }

        if(type == "wigner" && !wigner_prop) {
            if(!my_id)
                std::cerr<<"WignerG2 propagator necessary in a Wigner calculation."<<std::endl;
            MPI_Finalize();
            exit(1);
        }

        if(props.size() != nprops) {
            if (!my_id)
                std::cerr<<"Number of propagators in <parameters> and <propagators> do not match"<<std::endl;
            MPI_Finalize();
            exit(1);
        }

        if (!my_id)
            std::cout<<"Total number of beads: "<<bead_nums.back()+1<<std::endl;
    }

    void get_moves(pt::ptree move_params, pt::ptree params) {
        BOOST_FOREACH(const pt::ptree::value_type &m, move_params) {
            moves.push_back(Move::create(m.first));
            moves.back()->setup(m, params);
            moves.back()->propagator = props;
            move_names.push_back(m.first);

            if(type=="wigner" && m.first=="BrownianBridge") {
                moves.back()->set_beta(beta);
            }
        }
    }

public:
    Simulation(int myid, int ws) {
        my_id = myid;
        world_size = ws;
    }

    void load(const char *fname) {
        pt::ptree tree;
        read_xml(fname, tree);
        type = tree.get<string>("simulation.type");
        output_name = tree.get<string>("simulation.output_file");
        get_parameters(tree.get_child("simulation.parameters"));
        get_MC_params(tree.get_child("simulation.monte_carlo"));
        get_system(tree.get_child("simulation.system"));
        get_propagator(tree.get_child("simulation.propagators"), tree.get_child("simulation.parameters"));
        get_moves(tree.get_child("simulation.moves"), tree.get_child("simulation.parameters"));
    }

    void run_MC() {
        boost::shared_ptr<Configuration> c;
        if (type=="pimc")
            c.reset(new Configuration(natoms, ndimensions, bead_nums, mass));
        else if (type=="wigner")
            c.reset(new WignerConfiguration(natoms, ndimensions, bead_nums, mass));
        if (initial_config != "")
            c->load_config(initial_config);
        std::ofstream ofs(output_name + boost::lexical_cast<string>(my_id));
        arma::vec moves_accepted = arma::zeros<arma::vec>(moves.size());
        arma::vec moves_tried = arma::zeros<arma::vec>(moves.size());
        arma::cx_vec est = arma::zeros<arma::cx_vec>(nblocks);
        arma::uvec atoms = arma::regspace<arma::uvec>(0, natoms-1);
        arma::cx_vec total_weight = arma::zeros<arma::cx_vec>(nblocks);
        ofs<<c->header();
        for(unsigned b=0; b<nblocks; b++)
            for(unsigned i=0; i<nMC/(world_size * nblocks); i++) {
                for(unsigned m=0; m<moves.size(); m++)
                    (*moves[m])(c, atoms);

                ofs<<c->repr(pot);
                est(b) += c->weight() * (*estimator)(c->pos());
                total_weight(b) += c->weight();
            }
        est /= total_weight;

        arma::cx_vec est_global = arma::zeros<arma::cx_vec>(nblocks);
        MPI_Reduce(est.memptr(), est_global.memptr(), 2*nblocks, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        est_global /= world_size;

        std::vector<unsigned> global_moves_tried, global_moves_accepted;
        for(unsigned m=0; m<move_names.size(); m++) {
            unsigned tried = 0, accepted = 0;
            MPI_Reduce(&moves[m]->moves_tried, &tried, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&moves[m]->moves_accepted, &accepted, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
            global_moves_tried.push_back(tried);
            global_moves_accepted.push_back(accepted);
        }
        if(!my_id) {
            std::complex<double> mean_val = arma::mean(est_global);
            double real_std_val = arma::stddev(arma::real(est_global));
            double imag_std_val = arma::stddev(arma::imag(est_global));
            std::cout<<mean_val.real()<<'\t'<<mean_val.imag()<<'\t'<<real_std_val<<'\t'<<imag_std_val<<std::endl;
            for(unsigned m=0; m<move_names.size(); m++)
                std::cout<<move_names[m]<<'\t'<<(double)global_moves_accepted[m]/global_moves_tried[m]<<std::endl;
        }
    }
};

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int my_id, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    //generator.seed(static_cast<unsigned>(my_id) + 213);           // For reproducibility
    generator.seed(static_cast<unsigned>(my_id) + time(NULL));

    if(argc<2) {
        if (!my_id)
            std::cerr<<"Need an XML parameter file!"<<std::endl;
        MPI_Finalize();
        exit(1);
    }

    Simulation sim(my_id, world_size);
    sim.load(argv[1]);
    sim.run_MC();

    MPI_Finalize();
    return 0;
}
