#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using std::vector;
using std::string;

#include <armadillo>

#include <boost/random.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
namespace bRand = boost::random;
namespace pt = boost::property_tree;
using boost::shared_ptr;

bRand::mt19937 generator;

inline int random_integer(int low, int high) {
    bRand::uniform_int_distribution<> dist(low, high);
    return dist(generator);
}

inline double random_float(double low, double high) {
    bRand::uniform_real_distribution<double> dist(low, high);
    return dist(generator);
}

inline double random_normal(double mean, double sigma) {
    bRand::normal_distribution<double> dist(mean, sigma);
    return dist(generator);
}

struct Configuration {
    arma::cube positions;
    Configuration(unsigned natoms, unsigned nbeads, unsigned ndimensions) {
        positions = arma::zeros<arma::cube>(natoms, nbeads, ndimensions);
    }
    unsigned num_dims() const {
        return positions.n_slices;
    }
    unsigned num_beads() const {
        return positions.n_cols;
    }
    unsigned num_atoms() const {
        return positions.n_rows;
    }
    arma::vec bead_position(unsigned atom_num, unsigned time_ind) const {
        return positions.tube(atom_num, time_ind);
    }
    arma::mat necklace(unsigned atom_num) const {
        return positions(arma::span(atom_num), arma::span::all, arma::span::all);
    }
    arma::mat time_slice(unsigned time_ind) const {
        return positions(arma::span::all, arma::span(time_ind), arma::span::all);
    }
    arma::mat CoM() const {
        return arma::mean(positions, 1);
    }
};

struct Potential {
    virtual double operator()(arma::mat x) = 0;
};

struct HarmonicOscillator : public Potential {
    double omega;
    HarmonicOscillator() {}
    HarmonicOscillator(double w) : omega(w) {}
    double operator()(arma::mat x) {
        return 0.5 * omega * omega * arma::accu(x % x);
    }
};

struct Propagator {
    shared_ptr<Potential> V;
    double tau;
    Propagator(shared_ptr<Potential> pot, double Tau) {
        V = pot;
        tau = Tau;
    }
    virtual double operator()(Configuration conf) = 0;
};

struct G2 : public Propagator {
    G2(shared_ptr<Potential> pot, double tau) : Propagator(pot, tau) {}
    double operator()(Configuration conf) {
        double total_pot = 0;
        for(unsigned t=0; t<conf.num_beads(); t++)
            total_pot += (*V)(conf.time_slice(t));
        return std::exp(-total_pot * tau);
    }
};

struct G4 : public Propagator {
    G4(shared_ptr<Potential> pot, double tau) : Propagator(pot, tau) {}
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

struct G6 : public Propagator {
    G6(shared_ptr<Potential> pot, double tau) : Propagator(pot, tau) {}
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

struct G8 : public Propagator {
    G8(shared_ptr<Potential> pot, double tau) : Propagator(pot, tau) {}
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

struct Moves {
    unsigned moves_accepted, moves_tried;
    shared_ptr<Propagator> propagator;
    virtual Configuration operator()(Configuration conf, arma::uvec atom_nums) = 0;
    void check_amplitude(Configuration &conf, Configuration &conf_new) {
        moves_tried++;
        if((*propagator)(conf_new)/(*propagator)(conf)>random_float(0,1)) {
            conf = conf_new;
            moves_accepted++;
        }
    }
};

struct Translate : Moves {
    arma::mat max_step;
    Translate(double translation_length, unsigned natoms, unsigned ndims, shared_ptr<Propagator> prop) {
        max_step = translation_length * arma::ones<arma::mat>(natoms, ndims);
        moves_accepted = moves_tried = 0;
        propagator = prop;
    }
    virtual Configuration operator()(Configuration conf, arma::uvec atom_nums) {
        arma::vec step(max_step.n_cols);
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++) {
            for(unsigned i=0; i<max_step.n_cols; i++)
                step(i) = random_float(-1,1) * max_step(atom_nums(atom), i);

            Configuration conf_new = conf;
            for(unsigned time_ind=0; time_ind < conf.positions.n_cols; time_ind++)
                conf_new.positions.tube(atom_nums(atom), time_ind) += step;

            check_amplitude(conf, conf_new);
        }
        return conf;
    }
};

struct BrownianBridge : Moves {
    unsigned num_beads_moved;
    unsigned num_attempts;
    double tau;
    BrownianBridge(unsigned bb_length, unsigned bb_attempts, double Tau, shared_ptr<Propagator> prop) {
        num_beads_moved = bb_length;
        tau = Tau;
        num_attempts = bb_attempts;
        moves_accepted = moves_tried = 0;
        propagator = prop;
        /*****************************************************
        JUST ADDED THE VARIABLE, TAKE THIS ALSO INTO ACCOUNT WHEN REDESIGNING THE MOVES.
        YOU DONT JUST DO ONE B.B., YOU TRY AS MANY AS YOU NEED TO RECONFIGURE ON AVERAGE
        50% OF THE POLYMER
        *******************************************************/
    }
    virtual Configuration operator()(Configuration conf, arma::uvec atom_nums) {
        unsigned nbeads = conf.num_beads();
        for(unsigned atom=0; atom<atom_nums.n_rows; atom++)
            for(unsigned attempt=0; attempt<num_attempts; attempt++) {
                unsigned start = random_integer(0, nbeads-1);
                unsigned end = (start + num_beads_moved + 1) % nbeads;
                Configuration conf_new = conf;
                for(unsigned b = (start+1)%nbeads, nmoved=0; nmoved < num_beads_moved; b = (b+1)%nbeads, nmoved++) {
                    unsigned L = num_beads_moved - nmoved + 1;
                    for(unsigned d=0; d<conf.num_dims(); d++)
                        conf_new.positions(atom_nums(atom), b, d) = random_normal(((L-1.)*conf_new.positions(atom_nums(atom), (nbeads+b-1)%nbeads, d) + conf_new.positions(atom_nums(atom), end, d))/L, std::sqrt((L-1.)/L*tau));
                }
                check_amplitude(conf, conf_new);
            }
        return conf;
    }
};

class Simulation {
    unsigned ndimensions, natoms, nbeads;
    double beta;

    vector<shared_ptr<Moves> > moves;
    vector<string> move_names;

    unsigned nMC, nblocks;
    shared_ptr<Potential> pot;
    shared_ptr<Propagator> propagator;

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
    virtual void get_propagator(pt::ptree prop) {
        pot = shared_ptr<Potential>(new HarmonicOscillator(1.0));
        switch(prop.get<unsigned>("order")) {
        case 2:
            propagator = shared_ptr<Propagator>(new G2(pot, beta/nbeads));
            break;
        case 4:
            propagator = shared_ptr<Propagator>(new G4(pot, beta/nbeads));
            break;
        case 6:
            propagator = shared_ptr<Propagator>(new G6(pot, beta/nbeads));
            break;
        case 8:
            propagator = shared_ptr<Propagator>(new G8(pot, beta/nbeads));
            break;
        default:
            std::cerr<<"Not a valid propagator"<<std::endl;
            exit(1);
        }
    }
    void get_moves(pt::ptree move_params) {
        BOOST_FOREACH(const pt::ptree::value_type &m, move_params) {
            append_move(m);
        }
    }
    virtual void append_move(pt::ptree::value_type m) {
        if(m.first == "translation") {
            moves.push_back(shared_ptr<Moves>(new Translate(m.second.get<double>("length"), natoms, ndimensions, propagator)));
            move_names.push_back(m.first);
        } else if(m.first == "brownian_bridge") {
            moves.push_back(shared_ptr<Moves>(new BrownianBridge(m.second.get<unsigned>("length"), m.second.get<unsigned>("attempts"), beta/nbeads, propagator)));
            move_names.push_back(m.first);
        } else {
            std::cerr<<"Move not defined: "<<m.first<<std::endl;
            exit(1);
        }
    }
public:
    void load(const char *fname) {
        pt::ptree tree;
        read_xml(fname, tree);
        get_parameters(tree.get_child("pimc.parameters"));
        get_MC_params(tree.get_child("pimc.monte_carlo"));
        get_propagator(tree.get_child("pimc.propagator"));
        get_moves(tree.get_child("pimc.moves"));
    }
    void run_MC() {
        Configuration c(natoms, nbeads, ndimensions);
        std::ofstream ofs("test");
        arma::vec moves_accepted = arma::zeros<arma::vec>(moves.size());
        arma::vec moves_tried = arma::zeros<arma::vec>(moves.size());
        arma::vec xsq = arma::zeros<arma::vec>(nblocks);
        arma::uvec atoms = arma::regspace<arma::uvec>(0, natoms-1);
        for(unsigned b=0; b<nblocks; b++)
            for(unsigned i=0; i<nMC/nblocks; i++) {
                for(unsigned m=0; m<moves.size(); m++)
                    c = (*moves[m])(c, atoms);
                xsq(b) += c.time_slice(0)(0)*c.time_slice(0)(0);
                ofs<<c.time_slice(0).st()<<std::endl;
            }
        xsq /= nMC/nblocks;
        std::cout<<arma::mean(xsq)<<'\t'<<arma::stddev(xsq)<<std::endl;
        std::cout<<"Fraction of moves accepted: \n";
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
