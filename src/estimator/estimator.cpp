#include "estimator.hpp"

std::map<std::string, EstimatorFactory *> &Estimator::get_factory() {
    static std::map<std::string, EstimatorFactory *> factory;
    return factory;
}

void Estimator::registerType(const std::string &name, EstimatorFactory *factory) { get_factory()[name] = factory; }

std::shared_ptr<Estimator> Estimator::create(const std::string &name) {
    if (get_factory().find(name) == get_factory().end()) {
        std::cerr << "Not a valid estimator: " << name << std::endl;
        exit(1);
    }
    return get_factory()[name]->create();
}
