#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <memory>
#include <string>

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::iequals;

#include <boost/filesystem.hpp>

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

#define REGISTER_TYPE_GENERAL(klass, type)                                                                             \
    class klass##Factory : public type##Factory {                                                                      \
    public:                                                                                                            \
        klass##Factory() { type::registerType(#klass, this); }                                                         \
        virtual std::shared_ptr<type> create() { return std::make_shared<klass>(); }                                   \
    };                                                                                                                 \
    static klass##Factory global_##klass##Factory;

extern int my_id, world_size;

namespace utilities {
bool   boolparse(std::string const &st);
double exp_series(double var, unsigned n);

template <typename T> void setsink(std::shared_ptr<T> pfile, std::string clv) {
    if (iequals(clv, "all") || iequals(clv, "trace"))
        pfile->set_level(spdlog::level::trace);
    else if (iequals(clv, "debug"))
        pfile->set_level(spdlog::level::debug);
    else if (iequals(clv, "info"))
        pfile->set_level(spdlog::level::info);
    else if (iequals(clv, "warning"))
        pfile->set_level(spdlog::level::warn);
    else if (iequals(clv, "error"))
        pfile->set_level(spdlog::level::err);
    else if (iequals(clv, "critical"))
        pfile->set_level(spdlog::level::critical);
    else if (iequals(clv, "off"))
        pfile->set_level(spdlog::level::off);
    else
        spdlog::error("Log level " + clv +
                      " not found. Levels are: "
                      "all,debug,info,warning,error,critical,off");
}
} // namespace utilities

#endif
