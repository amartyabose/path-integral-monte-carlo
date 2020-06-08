#include "spdlog/spdlog.h"

#include "units.hpp"

Dimensions const Mass({{"mass", 1}, {"length", 0}, {"time", 0}, {"temperature", 0}});
Dimensions const Length({{"mass", 0}, {"length", 1}, {"time", 0}, {"temperature", 0}});
Dimensions const Time({{"mass", 0}, {"length", 0}, {"time", 1}, {"temperature", 0}});
Dimensions const Temperature({{"mass", 0}, {"length", 0}, {"time", 0}, {"temperature", 1}});
Dimensions const Force({{"mass", 1}, {"length", 1}, {"time", -2}, {"temperature", 0}});
Dimensions const Energy({{"mass", 1}, {"length", 2}, {"time", -2}, {"temperature", 0}});

bool operator<(Dimensions const &d1, Dimensions const &d2) {
    // This is an arbitrary ordering function to allow us to use Dimensions as
    // keys in a map. Change this with impunity if required --- doesn't affect
    // the code anywhere else.
    double val1 = d1.dim.at("mass") + 2 * d1.dim.at("length") + 4 * d1.dim.at("time") + 8 * d1.dim.at("temperature");
    double val2 = d2.dim.at("mass") + 2 * d2.dim.at("length") + 4 * d2.dim.at("time") + 8 * d2.dim.at("temperature");
    return val1 < val2;
}

void Units::setup(pt::ptree node) {
    spdlog::info("Units setup.");
    double mass_default = 1, length_default = 1, time_default = 1, temperature_default = 1, force_default = -1,
           energy_default = -1;
    std::string unit_type = node.get<std::string>("<xmlattr>.system", "");
    spdlog::info("\tusing " + unit_type + " unit system.");
    if (unit_type == "atomic") {
        mass_default        = 9.1093837e-31;
        length_default      = 5.291772109e-11;
        time_default        = 2.41888432658e-17;
        temperature_default = 3.1577464e5;
    } else if (unit_type == "real") {
        mass_default        = 1e-3 / 6.02214e23;
        length_default      = 1e-10;
        time_default        = 1e-15;
        temperature_default = 1;
        energy_default      = 4184 / 6.02214e23;
    }

    scaling[Mass]        = node.get<double>("mass", mass_default);
    scaling[Length]      = node.get<double>("length", length_default);
    scaling[Time]        = node.get<double>("time", time_default);
    scaling[Temperature] = node.get<double>("temperature", temperature_default);

    double energy = node.get<double>("energy", energy_default); //-1 means infer from base units

    if (energy > 0) {
        scaling[Energy] = energy;
        scaling[Force]  = energy / scaling[Length];
    }

    hbar = 1.0545718e-34 / (scale(Energy, true) * scale(Time, true));
    kB   = 1.38064852e-23 / scale(Energy, true) * scale(Temperature, true);
    spdlog::info("\thbar is " + std::to_string(hbar));
    spdlog::info("\tkB is " + std::to_string(kB));
    force_to_base_units      = non_specialify(Force, 1);
    energy_to_non_base_units = specialify(Energy, 1);
    spdlog::info("\tforce in base units is " + std::to_string(force_to_base_units));
    spdlog::info("\tenergy in default units is " + std::to_string(energy_to_non_base_units));
}

Units units;