#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <armadillo>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include <mpi.h>

#include "spdlog/spdlog.h"

#include "simulation.hpp"

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    generator.seed(static_cast<unsigned>(my_id) + 213);

    pt::ptree tree;
    try {
        read_xml(argv[1], tree, pt::xml_parser::trim_whitespace);
    } catch (...) {
        spdlog::critical("xml input file is either misspelled or nonexistent");
        exit(1);
    }

    Simulation sim(my_id, world_size);
    try {
        sim.load(tree);
        sim.run();
    } catch (std::runtime_error &err) {
        spdlog::critical(err.what());
        MPI_Finalize();
        exit(1);
    }

    MPI_Finalize();
}