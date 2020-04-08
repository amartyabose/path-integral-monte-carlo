#include "move.hpp"

std::map<std::string, MoveFactory *> &Move::get_factory() {
    static std::map<std::string, MoveFactory *> factory;
    return factory;
}

void Move::registerType(const std::string &name, MoveFactory *factory) { get_factory()[name] = factory; }

std::shared_ptr<Move> Move::create(const std::string &name) {
    if (get_factory().find(name) == get_factory().end()) {
        std::cerr << "Not a valid move: " << name << std::endl;
        exit(1);
    }
    return get_factory()[name]->create();
}
