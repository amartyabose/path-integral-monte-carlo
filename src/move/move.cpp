#include "move.hpp"

void Move::check_amplitude(std::shared_ptr<Configuration> &conf_old, std::shared_ptr<Configuration> conf_new) const {
    moves_tried++;
    double new_weight = 1, old_weight = 1;
    for (unsigned p = 0; p < propagator.size(); p++) {
        new_weight *= (*propagator[p])(conf_new->get_augmented_segment(p, p + 1));
        old_weight *= (*propagator[p])(conf_old->get_augmented_segment(p, p + 1));
    }
    if (new_weight / old_weight > random_float(0, 1)) {
        conf_old.reset(conf_new->duplicate());
        moves_accepted++;
    }
}

void Move::check_amplitude(std::shared_ptr<Configuration> &conf_old, std::shared_ptr<Configuration> conf_new,
                           unsigned index) const {
    moves_tried++;
    double new_weight = 1, old_weight = 1;
    for (unsigned p = 0; p < propagator.size(); p++) {
        new_weight *= (*propagator[p])(conf_new->get_augmented_segment(p, p + 1), index);
        old_weight *= (*propagator[p])(conf_old->get_augmented_segment(p, p + 1), index);
    }
    if (new_weight / old_weight > random_float(0, 1)) {
        conf_old.reset(conf_new->duplicate());
        moves_accepted++;
    }
}

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
