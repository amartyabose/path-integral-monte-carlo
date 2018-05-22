#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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

struct Parameters {
    unsigned ndimensions, natoms, nbeads;
    double beta;
    unsigned nBrownianBridge, lBrownianBridge; 
    double nTranslateLength;
    
    unsigned nMC, nblocks;
    unsigned prop_level;
    void load(const char *fname) {
        pt::ptree tree;
        read_xml(fname, tree);
        prop_level = tree.get<unsigned>("pimc.propagator", 2);
        BOOST_FOREACH (pt::ptree::value_type const &v, tree.get_child("pimc")) {
            if(v.first == "parameters") {
                beta = v.second.get<double>("beta");
                ndimensions = v.second.get<unsigned>("num_dimensions");
                natoms = v.second.get<unsigned>("num_atoms");
                nbeads = v.second.get<unsigned>("num_beads");
                lBrownianBridge = v.second.get<unsigned>("brownian_bridge_length");
                nBrownianBridge = v.second.get<unsigned>("brownian_bridge_attempts");
                nTranslateLength = v.second.get<double>("num_translate_length");
                nMC = v.second.get<unsigned>("num_MC");
                nblocks = v.second.get<unsigned>("num_blocks");
            }
        }
    }
} params;

struct Configuration {
    arma::cube positions;
    Configuration(Parameters params) {
        positions.set_size(params.natoms, params.nbeads, params.ndimensions);
        for(unsigned r=0; r<positions.n_rows; r++)
            for(unsigned c=0; c<positions.n_cols; c++)
                for(unsigned s=0; s<positions.n_slices; s++)
                    positions(r, c, s) = random_float(-1, 1);
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

template <typename PotFunc>
struct Propagator {
    PotFunc V;
    double tau;
    Propagator(PotFunc pot, Parameters param) {
        V = pot;
        tau = param.beta/param.nbeads;
    }
    virtual double operator()(Configuration conf) = 0;
};

template <typename PotFunc>
struct G2 : public Propagator<PotFunc> {
    G2(PotFunc pot, Parameters param) : Propagator<PotFunc>(pot, param) {}
    double operator()(Configuration conf) {
        double total_pot = 0;
        for(unsigned t=0; t<conf.num_beads(); t++)
            total_pot += this->V(conf.time_slice(t));
        return std::exp(-total_pot * this->tau);
    }
};

template <typename PotFunc>
struct G4 : public Propagator<PotFunc> {
    G4(PotFunc pot, Parameters param) : Propagator<PotFunc>(pot, param) {}
    double operator()(Configuration conf) {
        double total_amplitude = 1;
        for(unsigned t=0; t<conf.num_beads()-2; t+=2) {
            total_amplitude *=
                4./3. * std::exp(-this->tau/2. * this->V(conf.time_slice(t) - this->tau * this->V(conf.time_slice(t+1)) - this->tau/2. * this->V(conf.time_slice(t+2))))
                - 1./3. * std::exp(-this->tau * (this->V(conf.time_slice(t)) + this->V(conf.time_slice(t+2))));
        }
        return total_amplitude;
    }
};

template <typename PotFunc>
struct G6 : public Propagator<PotFunc> {
    G6(PotFunc pot, Parameters param) : Propagator<PotFunc>(pot, param) {}
    double operator()(Configuration conf) {
        double total_amplitude = 1;
        for(unsigned t=0; t<conf.num_beads(); t+=4)
            total_amplitude *=
                64./45. * std::exp(-this->tau * (this->V(conf.time_slice(t))/2. + this->V(conf.time_slice(t+1)) + this->V(conf.time_slice(t+2)) + this->V(conf.time_slice(t+3)) + this->V(conf.time_slice(t+4))/2.))
                - 4./9. * std::exp(-this->tau * (this->V(conf.time_slice(t)) + 2.*this->V(conf.time_slice(t+2)) + this->V(conf.time_slice(t+4))))
                + 1./45. * std::exp(-2*this->tau * (this->V(conf.time_slice(t)) + this->V(conf.time_slice(t+4))));
        return total_amplitude;
    }
};

template <typename PotFunc>
struct G8 : public Propagator<PotFunc> {
    G8(PotFunc pot, Parameters param) : Propagator<PotFunc>(pot, param) {}
    double operator()(Configuration conf) {
        double total_amplitude = 1;
        for(unsigned t=0; t<conf.num_beads(); t+=6)
            total_amplitude *=
                54./35. * std::exp(-this->tau * (this->V(conf.time_slice(t))/2. + this->V(conf.time_slice(t+1)) + this->V(conf.time_slice(t+2)) + this->V(conf.time_slice(t+3)) + this->V(conf.time_slice(t+4)) + this->V(conf.time_slice(t+5)) + this->V(conf.time_slice(t+6))/2.))
                - 27./40. * std::exp(-this->tau * (this->V(conf.time_slice(t)) + 2.*this->V(conf.time_slice(t+2)) + 2.*this->V(conf.time_slice(t+4)) + this->V(conf.time_slice(t+6))))
                + 2./15. * std::exp(-this->tau * (1.5*this->V(conf.time_slice(t)) + 3*this->V(conf.time_slice(t+3)) + 1.5*this->V(conf.time_slice(t+6))))
                - 1./840. * std::exp(-3 * this->tau * (this->V(conf.time_slice(t)) + this->V(conf.time_slice(t+6))));
        return total_amplitude;
    }
};

struct HarmonicOscillator {
    double omega;
    HarmonicOscillator() {}
    HarmonicOscillator(double w) : omega(w) {}
    double operator()(arma::mat x) {
        return 0.5 * omega * omega * arma::accu(x % x);
    }
};

struct Moves {
    virtual Configuration operator()(Configuration conf, unsigned atom_num) = 0;
};

struct Translate : Moves {
    arma::mat max_step;
    Translate(Parameters param) {
        max_step = param.nTranslateLength*arma::ones<arma::mat>(param.natoms, param.ndimensions);
    }
    virtual Configuration operator()(Configuration conf, unsigned atom_num) {
        arma::vec step(max_step.n_cols);
        for(unsigned i=0; i<max_step.n_cols; i++)
            step(i) = random_float(-1,1) * max_step(atom_num, i);

        for(unsigned time_ind=0; time_ind < conf.positions.n_cols; time_ind++)
            conf.positions.tube(atom_num, time_ind) += step;
        return conf;
    }
};

struct BrownianBridge : Moves {
    unsigned num_beads_moved;
    unsigned num_attempts;
    double tau;
    BrownianBridge(Parameters param) {
        num_beads_moved = param.lBrownianBridge;
        tau = param.beta/param.nbeads;
	num_attempts = param.nBrownianBridge;
	/*****************************************************
		JUST ADDED THE VARIABLE, TAKE THIS ALSO INTO ACCOUNT WHEN REDESIGNING THE MOVES.
		YOU DONT JUST DO ONE B.B., YOU TRY AS MANY AS YOU NEED TO RECONFIGURE ON AVERAGE 
		50% OF THE POLYMER
	*******************************************************/
    }
    virtual Configuration operator()(Configuration conf, unsigned atom_num) {
        unsigned nbeads = conf.num_beads();
        unsigned start = random_integer(0, nbeads-1);
        unsigned end = (start + num_beads_moved + 1) % nbeads;
        for(unsigned b = (start+1)%nbeads, nmoved=0; nmoved < num_beads_moved; b = (b+1)%nbeads, nmoved++) {
            unsigned L = num_beads_moved - nmoved + 1;
            for(unsigned d=0; d<conf.num_dims(); d++) {
                //conf.positions(atom_num, b, d) = random_normal((conf.positions(atom_num, (nbeads+b-1)%nbeads, d) + (L-1.)*conf.positions(atom_num, end, d))/L, std::sqrt((L-1.)/L*tau));
		conf.positions(atom_num, b, d) = random_normal(((L-1.)*conf.positions(atom_num, (nbeads+b-1)%nbeads, d) + conf.positions(atom_num, end, d))/L, std::sqrt((L-1.)/L*tau));
            }
        }
        return conf;
    }
};

int main(int argc, char **argv) {
    if(argc<2) {
        std::cerr<<"Need an XML parameter file!"<<std::endl;
        exit(1);
    }
    params.load(argv[1]);
    Configuration c(params);
    HarmonicOscillator osc(1.0);
    shared_ptr<Propagator<HarmonicOscillator> > mep;
    switch(params.prop_level) {
    case 2:
        mep.reset(new G2<HarmonicOscillator>(osc, params));
        break;
    case 4:
        mep.reset(new G4<HarmonicOscillator>(osc, params));
        break;
    case 6:
        mep.reset(new G6<HarmonicOscillator>(osc, params));
        break;
    case 8:
        mep.reset(new G8<HarmonicOscillator>(osc, params));
        break;
    default:
        std::cout<<"Not a valid propagator"<<std::endl;
        exit(1);
    }
    std::vector<shared_ptr<Moves> > moves;
    moves.push_back(shared_ptr<Moves>(new Translate(params)));
    moves.push_back(shared_ptr<Moves>(new BrownianBridge(params)));
    std::ofstream ofs("test");
    arma::vec moves_accepted = arma::zeros<arma::vec>(moves.size());
    arma::vec moves_tried = arma::zeros<arma::vec>(moves.size());
    arma::vec xsq = arma::zeros<arma::vec>(params.nblocks);
    for(unsigned b=0; b<params.nblocks; b++)
        for(unsigned i=0; i<params.nMC/params.nblocks; i++) {
            for(unsigned m=0; m<moves.size(); m++)
                for(unsigned atom=0; atom<c.num_atoms(); atom++) {
                    
		    // OF COURSE THIS HORROR IS TEMPORARY
		    unsigned nmax=1;
		    if(m==1)
		        nmax=params.nBrownianBridge;
		    
		    for(unsigned nmv=0;nmv<nmax;nmv++)
		    {
		        Configuration cnew = (*moves[m])(c, atom);
			moves_tried(m)++;
                        if((*mep)(cnew)/(*mep)(c)>random_float(0,1)) {
                            c = cnew;
                            moves_accepted(m)++;
			}
                     }
		  }
                
            xsq(b) += c.time_slice(0)(0)*c.time_slice(0)(0);
            ofs<<c.time_slice(0).st()<<std::endl;
        }
    xsq /= params.nMC/params.nblocks;
    std::cout<<arma::mean(xsq)<<'\t'<<arma::stddev(xsq)<<std::endl;
    (moves_accepted/moves_tried).print("Fraction of moves accepted");  
    // BAD IDEA, BETTER TO KEEP TRACK OF THE TOTAL MOVES
    // SOMETIMES YOU DON'T KNOW EXACTLY HOW MANY MOVES YOU WILL TRY
    return 0;
}
