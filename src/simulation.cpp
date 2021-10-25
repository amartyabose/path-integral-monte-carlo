#include <mpi.h>

#include "spdlog/spdlog.h"

#include "random.hpp"
#include "simulation.hpp"
#include "utilities.hpp"

void Simulation::initialize_output(pt::ptree output_node) {
    output.each_node       = utilities::boolparse(output_node.get<std::string>("<xmlattr>.each_node", "NO"));
    output.each_block      = utilities::boolparse(output_node.get<std::string>("<xmlattr>.each_block", "NO"));
    output.out_progressive = utilities::boolparse(output_node.get<std::string>("<xmlattr>.progressive", "NO"));
    output.out_phasespace  = utilities::boolparse(output_node.get<std::string>("phasespace", "NO"));
    output.output_folder   = output_node.get<std::string>("output_folder", "output");
    output.type            = type;
}

void Simulation::setup_output(pt::ptree output_node) {
    spdlog::info("Creating estimators:");
    burn = output_node.get<unsigned>("<xmlattr>.equilibration", 1000);
    skip = output_node.get<unsigned>("<xmlattr>.skip", 20);
    spdlog::info("\tEstimators would be measured after an equilibration of " + std::to_string(burn) +
                 " steps, and would skip " + std::to_string(skip) + " MC steps in between.");
    for (auto const &o : output_node.get_child("estimators")) {
        std::string est_name    = o.first;
        std::string info_string = "\t" + est_name;

        output.ignor.push_back(utilities::boolparse(o.second.get<std::string>("<xmlattr>.IGNoR", "NO")));
        if (output.ignor.back()) {
            info_string += " with IGNOR.";
            output.estimator_names.push_back(est_name + ".IGNoR");
        } else {
            info_string += " without IGNOR.";
            output.estimator_names.push_back(est_name);
        }

        output.histogram.push_back(utilities::boolparse(o.second.get<std::string>("<xmlattr>.histogram", "NO")));
        info_string += output.histogram.back() ? " Histogramming enabled." : " Histogramming not enabled.";
        spdlog::info(info_string);

        output.estimators.push_back(Estimator::create(est_name));
        output.estimators.back()->setup(type, mass, beta, bead_nums.size(), o.second);
    }
    unsigned num_ests = output.estimator_names.size();
    output.est_vals.resize(num_ests);
    output.global_est_vals.resize(num_ests);
    output.est_re_vals_plus.resize(num_ests);
    output.est_im_vals_plus.resize(num_ests);
    output.est_re_vals_minus.resize(num_ests);
    output.est_im_vals_minus.resize(num_ests);
    output.weight_re_plus.resize(num_ests);
    output.weight_im_plus.resize(num_ests);
    output.weight_re_minus.resize(num_ests);
    output.weight_im_minus.resize(num_ests);
    for (unsigned i = 0; i < num_ests; i++) {
        output.est_vals[i].resize(nblocks);
        output.global_est_vals[i].resize(nblocks);
        output.est_re_vals_plus[i].resize(nblocks);
        output.est_im_vals_plus[i].resize(nblocks);
        output.est_re_vals_minus[i].resize(nblocks);
        output.est_im_vals_minus[i].resize(nblocks);
        output.weight_re_plus.resize(nblocks);
        output.weight_im_plus.resize(nblocks);
        output.weight_re_minus.resize(nblocks);
        output.weight_im_minus.resize(nblocks);

        for (unsigned b = 0; b < nblocks; b++) {
            output.est_vals[i][b] =
                arma::zeros<arma::cx_mat>(output.estimators[i]->n_rows, output.estimators[i]->n_cols);
            output.global_est_vals[i][b] =
                arma::zeros<arma::cx_mat>(output.estimators[i]->n_rows, output.estimators[i]->n_cols);
            output.est_re_vals_plus[i][b] =
                arma::zeros<arma::mat>(output.estimators[i]->n_rows, output.estimators[i]->n_cols);
            output.est_im_vals_plus[i][b] =
                arma::zeros<arma::mat>(output.estimators[i]->n_rows, output.estimators[i]->n_cols);
            output.est_re_vals_minus[i][b] =
                arma::zeros<arma::mat>(output.estimators[i]->n_rows, output.estimators[i]->n_cols);
            output.est_im_vals_minus[i][b] =
                arma::zeros<arma::mat>(output.estimators[i]->n_rows, output.estimators[i]->n_cols);
            output.weight_re_plus[b]  = 0;
            output.weight_im_plus[b]  = 0;
            output.weight_re_minus[b] = 0;
            output.weight_im_minus[b] = 0;
        }
    }
    output.setup(my_id, world_size, nblocks, bead_nums.size() - 1);
    spdlog::info("Estimators and output setup done.");
}

void Simulation::setup_logging(boost::optional<pt::ptree &> logging_node) {
    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::warn);
    if (my_id != 0)
        console_sink->set_level(spdlog::level::critical);

    boost::filesystem::create_directory(output.output_folder + "_logs");
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(output.output_folder + "_logs/nd." +
                                                                         std::to_string(my_id) + ".runlog");

    spdlog::set_default_logger(
        std::make_shared<spdlog::logger>("pimc", spdlog::sinks_init_list({console_sink, file_sink})));

    if (!logging_node) {
        spdlog::info("Using default logging settings");
        return;
    }

    boost::optional<pt::ptree &> fsettings = logging_node.get().get_child_optional("file");
    boost::optional<pt::ptree &> tsettings = logging_node.get().get_child_optional("terminal");

    if (!fsettings)
        spdlog::info("Logging to file set as default. You can change it in the "
                     "logging options");
    else {
        std::string ml = fsettings.get().get<std::string>("master.<xmlattr>.level", "trace");
        std::string ol = fsettings.get().get<std::string>("others.<xmlattr>.level", "critical");
        if (my_id == 0)
            utilities::setsink(file_sink, ml);
        else
            utilities::setsink(file_sink, ol);
    }

    if (!tsettings)
        spdlog::info("Logging to console set as default. You can change it in "
                     "the logging options");
    else {
        std::string ml = tsettings.get().get<std::string>("master.<xmlattr>.level", "warn");
        std::string ol = tsettings.get().get<std::string>("others.<xmlattr>.level", "critical");
        if (my_id == 0)
            utilities::setsink(console_sink, ml);
        else
            utilities::setsink(console_sink, ol);
    }
}

void Simulation::get_worm_data(boost::optional<pt::ptree &> worm_node) {
    if (worm_node) {
        pt::ptree wn    = worm_node.get();
        grand_canonical = wn.get<bool>("grand_canonical", false);
        auto mu         = wn.get_optional<double>("chemical_potential");
        if (grand_canonical && !mu)
            throw std::runtime_error("Chemical Potential required for grand canonical calculations.");
        chemical_potential = mu.get();
        max_worms          = wn.get<unsigned>("max_worms", 1);
    } else {
        grand_canonical = false;
        max_worms       = 0;
    }
}

void Simulation::get_parameters(pt::ptree params) {
    units.setup(params.get_child("units"));
    spdlog::info("System settings.");
    ndimensions = params.get<unsigned>("num_dimensions");
    natoms      = params.get<unsigned>("num_atoms");
    nprops      = params.get<unsigned>("num_propagators");
    get_worm_data(params.get_child_optional("worm"));

    auto beta_ = params.get_optional<double>("beta");
    auto T_    = params.get_optional<double>("temperature");
    if (!beta_ && !T_)
        throw std::runtime_error("Neither beta nor temperature provided.");
    else if (beta_ && T_)
        throw std::runtime_error("Both beta or temperature provided.");
    else if (beta_)
        beta = beta_.get();
    else if (T_)
        beta = 1.0 / (units.kB * T_.get());
    atom_specific_beta = beta * arma::ones<arma::vec>(natoms);

    std::string masses = params.get<std::string>("mass");
    mass               = arma::vec(masses);
    if (mass.n_rows == 1) {
        double m = mass(0);
        mass     = m * arma::ones<arma::vec>(natoms);
    } else if (mass.n_rows != natoms)
        throw std::runtime_error("Number of masses do not match the number of atoms.");

    get_potential(params.get_child("potential"));

    auto bound = params.get_child_optional("boundary_conditions");
    if (!bound) {
        bc = BoundaryConditions::create("OpenBoundaryConditions");
        spdlog::info("\tBoundary condition type not provided. Using open boundary conditions.");
    } else {
        auto bounds  = bound.get();
        auto bc_type = bounds.get<std::string>("<xmlattr>.name", "OpenBoundaryConditions");
        bc           = BoundaryConditions::create(bc_type);
        bc->setup(bounds);
        spdlog::info("\t" + bc_type + " created.");
        if (bc_type == "PeriodicBoundaryConditions") {
            arma::vec box = arma::vec(bounds.get<std::string>("box_size"));
            spdlog::info("\tThe particle number density is " + std::to_string(natoms / (box(0) * box(1) * box(2))) +
                         ".");
        }
    }

    auto open_chain_atoms = params.get_optional<std::string>("open_chains");
    if (open_chain_atoms) {
        std::string open_chain_vec = open_chain_atoms.get();
        spdlog::info("\tOpen chains on atoms: " + open_chain_vec);
        open_chains = arma::vec(open_chain_vec);
    }

    initial_config = params.get<std::string>("initial_configuration", "");
    spdlog::info("\tBasic system parameters setup done.");
}

void Simulation::get_MC_params(pt::ptree params) {
    nIC     = params.get<unsigned>("num_ICs");
    nblocks = params.get<unsigned>("num_blocks");
    spdlog::info("Monte Carlo parameters setup done.");
}

void Simulation::get_potential(pt::ptree pot_tree) {
    pot = Potential::create(pot_tree.get<std::string>("<xmlattr>.name"));
    pot->setup(pot_tree);
}

void Simulation::get_dts(pt::ptree prop) {
    double   total_temp_accounted_for = 0;
    unsigned undetermined_dt          = 0;
    for (auto const &p : prop) {
        dt_imaginary.push_back(p.second.get<double>("<xmlattr>.dt_frac", 0) * beta);
        total_temp_accounted_for += dt_imaginary.back();
        if (dt_imaginary.back() < 1e-10)
            undetermined_dt += p.second.get<double>("<xmlattr>.number");
    }
    if (arma::accu(total_temp_accounted_for - beta) > 1e-10 && undetermined_dt == 0)
        throw std::runtime_error("Full beta not accounted for in propagators");

    double rem_dt = (beta - total_temp_accounted_for) / undetermined_dt;
    for (unsigned i = 0; i < dt_imaginary.size(); i++)
        if (dt_imaginary[i] < 1e-10)
            dt_imaginary[i] = rem_dt;
}

void Simulation::insert_position_momentum_moves(pt::ptree::value_type p) {
    moves.push_back(Move::create("WignerPosition"));
    moves.back()->setup(p, beta, mass);
    move_names.push_back("WignerPosition");

    std::shared_ptr<Propagator> keprop(Propagator::create("WignerKE"));
    keprop->set_params(0., p, mass, beta);

    moves.push_back(Move::create("WignerMomentum"));
    moves.back()->setup(p, beta, mass);
    moves.back()->propagator.push_back(keprop);
    move_names.push_back("WignerMomentum");
}

void Simulation::get_propagator(pt::ptree prop) {
    get_dts(prop);
    bead_nums.push_back(0);
    unsigned prop_num    = 0;
    bool     wigner_prop = false;
    for (auto const &p : prop) {
        std::string prop_name = p.second.get<std::string>("<xmlattr>.name");
        if (type == "pimd" && prop_name != "G2")
            throw std::runtime_error("PIMD currently not possible with anything but G2 propagators.");
        unsigned nprop = p.second.get<unsigned>("<xmlattr>.number");
        props.push_back(Propagator::create(prop_name));
        props.back()->set_params(dt_imaginary[prop_num], p, mass, beta);
        if (prop_name == "WignerG2") {
            if (type != "wigner")
                throw std::runtime_error("WignerG2 propagator allowed only in a Wigner calculation.");
            else if (type == "wigner") {
                if (wigner_prop || prop_num > 1)
                    throw std::runtime_error("Only 1 WignerG2 propagator allowed in a Wigner calculation.");
                wigner_prop = true;
                insert_position_momentum_moves(p);
                for (unsigned atom = 0; atom < natoms; atom++)
                    atom_specific_beta(atom) =
                        atom_specific_beta(atom) * (1 - p.second.get<double>("<xmlattr>.dt_frac")) +
                        props.back()->get_tau(atom);
            }
        }
        prop_num++;
        bead_nums.push_back(bead_nums.back() + props.back()->total_beads() - 1);
        for (unsigned i = 1; i < nprop; i++) {
            props.push_back(props.back());
            bead_nums.push_back(bead_nums.back() + props.back()->total_beads() - 1);
        }
    }

    if (type == "wigner" && !wigner_prop)
        throw std::runtime_error("WignerG2 propagator allowed only in a Wigner calculation.");

    if (props.size() != nprops)
        throw std::runtime_error("number of propagators in <parameters> and <propagators> do not match");

    spdlog::info("Propagators setup done. Total number of beads: " + std::to_string(bead_nums.back() + 1) +
                 " counting x(0) and x(\\beta) separately.");
}

void Simulation::get_moves(pt::ptree move_params) {
    for (auto const &m : move_params) {
        moves.push_back(Move::create(m.first));
        moves.back()->setup(m, beta, mass);
        moves.back()->propagator = props;
        move_names.push_back(m.first);

        if (m.first == "BrownianBridge") {
            if (m.second.get<unsigned>("length") > bead_nums.back())
                throw std::runtime_error("Length of Brownian bridge chain is greater "
                                         "than the necklace length");
            moves.back()->set_beta(atom_specific_beta);
        }
    }
    spdlog::info("Moves setup done.");
}

void Simulation::load(pt::ptree tree) {
    type = tree.get<std::string>("simulation.type");
    initialize_output(tree.get_child("simulation.output"));
    auto logging_node = tree.get_child_optional("simulation.logging");
    setup_logging(logging_node);
    if (type != "pimc" && type != "wigner")
        throw std::runtime_error("Simulation type " + type + " is not supported.");
    get_parameters(tree.get_child("simulation.parameters"));
    get_propagator(tree.get_child("simulation.propagators"));
    get_moves(tree.get_child("simulation.moves"));
    get_MC_params(tree.get_child("simulation.monte_carlo"));
    setup_output(tree.get_child("simulation.output"));
    spdlog::info("Finished loading settings");
}

void Simulation::run() {
    std::shared_ptr<Configuration> conf;
    if (type == "pimc")
        conf.reset(new Configuration(natoms, ndimensions, bead_nums, mass, open_chains));
    else if (type == "wigner")
        conf.reset(new WignerConfiguration(natoms, ndimensions, bead_nums, mass));

    if (initial_config != "") {
        spdlog::info("Loading initial configuration: " + initial_config);
        conf->load_config(initial_config);
        spdlog::info("Initial configuration loaded.");
        spdlog::info("Potential at the initial config = " +
                     std::to_string(units.energy_to_non_base_units * (*pot)(conf->time_slice(0))));
    }

    spdlog::info("Starting equilibration.");
    arma::uvec atoms = arma::regspace<arma::uvec>(0, natoms - 1);
    for (unsigned i = 0; i < burn; i++)
        for (auto const &m : moves)
            (*m)(conf, atoms);

    spdlog::info("Starting actual runs.");
    for (unsigned i = 0; i < nblocks; i++) {
        run_block(conf);
        output.finalize_block();
        spdlog::info("Block " + std::to_string(i + 1) + " of " + std::to_string(nblocks) + " done.");
    }

    unsigned k = 0;
    for (auto const &m : moves) {
        unsigned tried = 0, accepted = 0;
        MPI_Reduce(&m->moves_tried, &tried, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&m->moves_accepted, &accepted, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
        spdlog::info(move_names[k] +
                     " acceptance ratio = " + std::to_string((double)m->moves_accepted / m->moves_tried));
        if (!my_id)
            spdlog::info("Global " + move_names[k++] +
                         " acceptance ratio = " + std::to_string((double)accepted / tried));
    }

    output.finalize_output();
}

void Simulation::run_block(std::shared_ptr<Configuration> &conf) {
    arma::uvec atoms = arma::regspace<arma::uvec>(0, natoms - 1);
    for (int i = 0; i < nIC; i++) {
        for (unsigned s = 0; s <= skip; s++)
            for (auto const &m : moves)
                (*m)(conf, atoms);
        output.add_config(conf);
    }
}
