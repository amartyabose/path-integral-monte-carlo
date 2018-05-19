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
namespace bRand = boost::random;
namespace pt = boost::property_tree;

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
    void load(const char *fname) {
        pt::ptree tree;
        read_xml(fname, tree);
        BOOST_FOREACH (pt::ptree::value_type const &v, tree.get_child("pimc")) {
            if(v.first == "parameters") {
                beta = v.second.get<double>("beta");
                ndimensions = v.second.get<double>("num_dimensions");
                natoms = v.second.get<double>("num_atoms");
                nbeads = v.second.get<double>("num_beads");
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
    arma::mat time_view(unsigned time_ind) const {
        return positions(arma::span::all, arma::span(time_ind), arma::span::all);
    }
};

struct Moves {};

struct Translate : Moves {
    Configuration operator()(Configuration conf, unsigned atom_num, arma::vec dr) {
        for(unsigned time_ind=0; time_ind < conf.positions.n_cols; time_ind++)
            conf.positions.tube(atom_num, time_ind) += dr;
        return conf;
    }
};

struct BrownianBridge : Moves {
    unsigned num_beads_moved;
    double tau;
    BrownianBridge(unsigned nbeads_move, Parameters param) {
        num_beads_moved = nbeads_move;
        tau = param.beta/param.nbeads;
    }
    Configuration operator()(Configuration conf, unsigned atom_num) {
        unsigned nbeads = conf.num_beads();
        unsigned start = random_integer(0, nbeads);
        unsigned end = (start + num_beads_moved + 1) % nbeads;
        for(int b = start+1, nmoved=0; nmoved < num_beads_moved; b = (b+1)%nbeads, nmoved++) {
            unsigned L = num_beads_moved - nmoved;
            for(unsigned d=0; d<conf.num_dims(); d++)
                conf.positions(atom_num, b, d) = random_normal((conf.positions(atom_num, (b-1)%nbeads, d) + (L-1.)*conf.positions(atom_num, end, d))/L, std::sqrt((L-1.)/L)*tau);
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
    arma::mat necklace0 = c.necklace(0);
    BrownianBridge bb(5, params);
    arma::mat necklacebb = bb(c, 0).necklace(0);
    (necklace0 - necklacebb).print("change?");
    return 0;
}
