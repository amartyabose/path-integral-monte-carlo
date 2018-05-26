#ifndef _UTILITIES_HPP_
#define _UTILITIES_HPP_

#define REGISTER_TYPE_GENERAL(klass, type)                                       \
    class klass##Factory : public type##Factory {                        \
    public:                                                              \
        klass##Factory() {                                               \
            type::registerType(#klass, this);                            \
        }                                                                \
        virtual type* create() {                                         \
            return new klass();                                           \
        }                                                                \
    };                                                                   \
    static klass##Factory global_##klass##Factory;

#endif
