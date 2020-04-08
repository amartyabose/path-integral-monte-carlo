#include "propagator.hpp"

std::map<std::string, PropagatorFactory *> &Propagator::get_factory() {
    static std::map<std::string, PropagatorFactory *> factory;
    return factory;
}

void Propagator::registerType(const std::string &name, PropagatorFactory *factory) { get_factory()[name] = factory; }

std::shared_ptr<Propagator> Propagator::create(const std::string &name) {
    if (get_factory().find(name) == get_factory().end()) {
        std::cerr << "Not a valid propagator" << std::endl;
        exit(1);
    }
    return get_factory()[name]->create();
    // return std::shared_ptr<Propagator>(get_factory()[name]->create());
}
