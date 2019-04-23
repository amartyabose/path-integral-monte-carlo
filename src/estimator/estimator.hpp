#ifndef _ESTIMATOR_H_
#define _ESTIMATOR_H_

#include <boost/shared_ptr.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include <armadillo>

#include "../configuration.hpp"

class Estimator;

struct EstimatorFactory {
    virtual Estimator* create() = 0;
};

class Estimator {
    static std::map<std::string, EstimatorFactory*>& get_factory();
public:
    arma::vec final_values_real, final_values_imag;
    virtual void setup(std::string type, pt::ptree params, unsigned nblocks=10) = 0;
    virtual double eval(boost::shared_ptr<Configuration> x) = 0;
    virtual void operator() (boost::shared_ptr<Configuration> x, unsigned block_num) {
        std::complex<double> weight = x->weight();
        if (weight.real() >= 0) {
            i_real_plus(block_num) += weight.real();
            values_real_plus(block_num) += weight.real() * eval(x);
        } else {
            i_real_minus(block_num) += weight.real();
            values_real_minus(block_num) += weight.real() * eval(x);
        }
        if (weight.imag() >= 0) {
            i_imag_plus(block_num) += weight.imag();
            values_imag_plus(block_num) += weight.imag() * eval(x);
        } else {
            i_imag_minus(block_num) += weight.imag();
            values_imag_minus(block_num) += weight.imag() * eval(x);
        }
    }
    virtual void process() {
        arma::vec total_ampl_real = i_real_plus + i_real_minus;
        if (arma::any(i_real_minus)>0)
            final_values_real = 0.5 * (
                values_real_plus/total_ampl_real + values_real_minus/i_real_minus % (1. - i_real_plus/total_ampl_real)
                + values_real_minus/total_ampl_real + values_real_plus/i_real_plus % (1. - i_real_minus/total_ampl_real)
                );
        else
            final_values_real = values_real_plus/total_ampl_real;

        arma::vec total_ampl_imag = i_imag_plus + i_imag_minus;
        if (arma::any(i_imag_minus)>0)
            final_values_imag = 0.5/total_ampl_real % (
                values_imag_plus/total_ampl_imag + values_imag_minus/i_imag_minus % (1. - i_imag_plus/total_ampl_imag)
                + values_imag_minus/total_ampl_imag + values_imag_plus/i_imag_plus % (1. - i_imag_minus/total_ampl_imag)
                );
        else if (arma::any(total_ampl_imag) != 0)
            final_values_imag = values_imag_plus/total_ampl_imag;
        else
            final_values_imag = arma::zeros<arma::vec>(i_real_plus.n_rows);
    }

    static void registerType(const std::string &name, EstimatorFactory *factory);
    static boost::shared_ptr<Estimator> create(const std::string &name);

protected:
    arma::vec values_real_plus, values_real_minus;
    arma::vec values_imag_plus, values_imag_minus;
    arma::vec i_real_plus, i_real_minus;
    arma::vec i_imag_plus, i_imag_minus;
};

#endif
