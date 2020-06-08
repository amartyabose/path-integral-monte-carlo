#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <complex>
#include <string>
#include <vector>

#include "configuration.hpp"
#include "estimator/estimator.hpp"

struct Output {
    std::string output_folder;

    bool each_node, each_block;
    bool out_phasespace, phasespace_histogram;
    bool out_progressive;

    std::string type;

    int                  world_size, my_id;
    std::complex<double> weight_per_block;
    std::vector<bool>    ignor, histogram;

    std::vector<std::shared_ptr<Estimator>> estimators;
    std::vector<std::string>                estimator_names;
    std::vector<std::vector<arma::cx_mat>>  est_vals;
    std::vector<std::vector<arma::cx_mat>>  global_est_vals;

    std::vector<std::vector<arma::mat>> est_re_vals_plus, est_re_vals_minus;
    std::vector<std::vector<arma::mat>> est_im_vals_plus, est_im_vals_minus;

    std::vector<double> weight_re_plus, weight_re_minus, weight_im_plus, weight_im_minus;

    void setup(int rank, int world_size, int nblocks);
    void add_config(std::shared_ptr<Configuration> const &conf);
    void finalize_block();
    void finalize_output();

private:
    int           block_num;
    std::string   node_dir, progressive_dir;
    std::ofstream phase_space_stream;
    std::ofstream node_spec_est_stream;
    std::ofstream est_stream;
    std::ofstream progressive_stream;

    void write(std::ostream &os, arma::cx_mat const &val);
    void write(std::ostream &os, arma::cx_mat const &val, arma::mat const &err_re, arma::mat const &err_im);

    arma::cx_mat process_ignor(unsigned est_num, unsigned block_num, std::string obs_name);
};

#endif