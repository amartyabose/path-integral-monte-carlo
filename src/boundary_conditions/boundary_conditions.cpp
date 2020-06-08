#include "boundary_conditions.hpp"

std::map<std::string, BoundaryConditionsFactory *> &BoundaryConditions::get_factory() {
    static std::map<std::string, BoundaryConditionsFactory *> factory;
    return factory;
}

void BoundaryConditions::registerType(const std::string &name, BoundaryConditionsFactory *factory) {
    get_factory()[name] = factory;
}

std::shared_ptr<BoundaryConditions> BoundaryConditions::create(const std::string &name) {
    if (get_factory().find(name) == get_factory().end()) {
        std::cerr << "Not a valid boundary condition" << std::endl;
        exit(1);
    }
    return get_factory()[name]->create();
}

std::shared_ptr<BoundaryConditions> bc;