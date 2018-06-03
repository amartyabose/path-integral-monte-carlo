#include "potential.hpp"

std::map<std::string, PotentialFactory*>& Potential::get_factory() {
    static std::map<std::string, PotentialFactory*> factory;
    return factory;
}

void Potential::registerType(const std::string &name, PotentialFactory *factory) {
    get_factory()[name] = factory;
}

boost::shared_ptr<Potential> Potential::create(const std::string &name) {
    if(get_factory().find(name) == get_factory().end()) {
        std::cerr<<"Not a valid potential"<<std::endl;
        exit(1);
    }
    return boost::shared_ptr<Potential>(get_factory()[name]->create());
}
