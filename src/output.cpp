#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include <mpi.h>

#include "spdlog/spdlog.h"

#include "output.hpp"

void Output::setup(int rank, int size, int nblocks, int num_beads) {
    my_id            = rank;
    world_size       = size;
    block_num        = 0;
    weight_per_block = 0.0;
    spdlog::debug("Creating directories");
    if (!my_id) {
        if (boost::filesystem::exists(output_folder))
            spdlog::warn("Directory " + output_folder + " exists. Files may be rewritten.");
        else
            boost::filesystem::create_directory("./" + output_folder);

        if (out_progressive) {
            progressive_dir = output_folder + "/progressive";
            boost::filesystem::create_directory(progressive_dir);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (each_node) {
        node_dir = output_folder + "/node." + std::to_string(my_id);
        boost::filesystem::create_directory(node_dir);
    }

    bool any_hist = false;
    for (auto i = 0; i < histogram.size(); i++)
        if (histogram[i]) {
            if (!any_hist)
                boost::filesystem::create_directory(output_folder + "/histograms");

            std::ofstream ofs(output_folder + "/histograms/" + estimator_names[i] + ".node." + std::to_string(my_id));
            histogram_stream.push_back(std::move(ofs));
        }

    if (out_phasespace) {
        boost::filesystem::create_directory(output_folder + "/phasespace");
        boost::filesystem::create_directory(output_folder + "/phasespace/node-" + std::to_string(my_id));
        std::ofstream out;
        for (unsigned i = 0; i < num_beads; i++) {
            spdlog::info("Creating the phasespace file for bead number " + std::to_string(i));
            out.open(output_folder + "/phasespace/node-" + std::to_string(my_id) + "/bead" + std::to_string(i) +
                     ".configs.xyz");
            phase_space_stream.push_back(std::move(out));
        }
    }
}

void Output::add_config(std::shared_ptr<Configuration> const &conf) {
    static int frame_cnt = 0;
    if (out_phasespace)
        for (auto i = 0; i < phase_space_stream.size(); i++)
            phase_space_stream[i] << conf->repr(frame_cnt, i);

    std::complex<double> weight = conf->weight();
    if (weight.real() >= 0)
        weight_re_plus[block_num] += weight.real();
    else
        weight_re_minus[block_num] += weight.real();

    if (weight.imag() >= 0)
        weight_im_plus[block_num] += weight.imag();
    else
        weight_im_minus[block_num] += weight.imag();

    for (int i = 0; i < estimators.size(); i++) {
        arma::mat val = estimators[i]->eval(conf);
        est_vals[i][block_num] += conf->weight() * val;
        if (weight.real() >= 0)
            est_re_vals_plus[i][block_num] += weight.real() * val;
        else
            est_re_vals_minus[i][block_num] += weight.real() * val;

        if (weight.imag() >= 0)
            est_im_vals_plus[i][block_num] += weight.imag() * val;
        else
            est_im_vals_minus[i][block_num] += weight.imag() * val;

        if (histogram[i]) {
            histogram_stream[i] << frame_cnt << '\t' << conf->weight().real();
            for (unsigned r = 0; r < val.n_rows; r++)
                for (unsigned c = 0; c < val.n_cols; c++)
                    histogram_stream[i] << '\t' << val(r, c);
            histogram_stream[i] << std::endl;
        }
    }

    weight_per_block += conf->weight();
    frame_cnt++;
}

void Output::write(std::ostream &os, arma::cx_mat const &val) {
    if (type == "pimc")
        for (unsigned r = 0; r < val.n_rows; r++) {
            for (unsigned c = 0; c < val.n_cols; c++)
                os << boost::format("%.5e\t") % val(r, c).real();
            os << std::endl;
        }
    else
        for (unsigned r = 0; r < val.n_rows; r++) {
            for (unsigned c = 0; c < val.n_cols; c++)
                os << boost::format("%.5e\t%.5e\t") % val(r, c).real() % val(r, c).imag();
            os << std::endl;
        }
}

void Output::write(std::ostream &os, arma::cx_mat const &val, arma::mat const &err_re, arma::mat const &err_im) {
    if (type == "pimc")
        for (unsigned r = 0; r < val.n_rows; r++) {
            for (unsigned c = 0; c < val.n_cols; c++)
                os << boost::format("%.5e\t%.5e\t") % val(r, c).real() % err_re(r, c);
            os << std::endl;
        }
    else
        for (unsigned r = 0; r < val.n_rows; r++) {
            for (unsigned c = 0; c < val.n_cols; c++)
                os << boost::format("%.5e\t%.5e\t%.5e\t%.5e\t") % val(r, c).real() % err_re(r, c) % val(r, c).imag() %
                          err_im(r, c);
            os << std::endl;
        }
}

arma::cx_mat Output::process_ignor(unsigned est_num, unsigned block_num, std::string obs_name) {
    if (!ignor[est_num]) {
        return arma::cx_mat(est_re_vals_plus[est_num][block_num] + est_re_vals_minus[est_num][block_num],
                            est_im_vals_plus[est_num][block_num] + est_im_vals_minus[est_num][block_num]) /
               std::complex<double>(weight_re_plus[block_num] + weight_re_minus[block_num],
                                    weight_im_plus[block_num] + weight_im_minus[block_num]);
    } else {
        double               total_weight_re = weight_re_plus[block_num] + weight_re_minus[block_num];
        double               total_weight_im = weight_im_plus[block_num] + weight_im_minus[block_num];
        std::complex<double> total_weight(total_weight_re, total_weight_im);

        arma::mat re_part;
        if (std::abs(weight_re_minus[block_num]) == 0)
            re_part = est_re_vals_plus[est_num][block_num] / total_weight_re;
        else
            re_part = 0.5 * (est_re_vals_plus[est_num][block_num] / total_weight_re +
                             est_re_vals_minus[est_num][block_num] / weight_re_minus[block_num] *
                                 (1. - weight_re_plus[block_num] / total_weight_re) +
                             est_re_vals_minus[est_num][block_num] / total_weight_re +
                             est_re_vals_plus[est_num][block_num] / weight_re_plus[block_num] *
                                 (1. - weight_re_minus[block_num] / total_weight_re));

        arma::mat im_part;
        if (std::abs(weight_im_minus[block_num]) == 0) {
            if (std::abs(weight_im_plus[block_num]) == 0)
                im_part = arma::zeros<arma::mat>(arma::size(re_part));
            else
                im_part = est_im_vals_plus[est_num][block_num] / total_weight_re;
        } else
            im_part = 0.5 * (est_im_vals_plus[est_num][block_num] / total_weight_re +
                             est_im_vals_minus[est_num][block_num] / weight_im_minus[block_num] *
                                 (0. - weight_im_plus[block_num] / total_weight_re) +
                             est_im_vals_minus[est_num][block_num] / total_weight_re +
                             est_im_vals_plus[est_num][block_num] / weight_im_plus[block_num] *
                                 (0. - weight_im_minus[block_num] / total_weight_re));
        return arma::cx_mat(re_part, im_part);
    }
}

void Output::finalize_block() {
    for (int i = 0; i < estimators.size(); i++) {
        est_vals[i][block_num] = process_ignor(i, block_num, estimator_names[i]);
        arma::cx_mat val       = est_vals[i][block_num];
        if (estimator_names[i] == "CV") {
            val = arma::zeros<arma::cx_mat>(1, 1);
            val = est_vals[i][block_num](0, 0) - est_vals[i][block_num](1, 0) * est_vals[i][block_num](1, 0);
        }
        if (each_node) {
            node_spec_est_stream.open(node_dir + "/" + estimator_names[i] + ".block." + std::to_string(block_num));
            write(node_spec_est_stream, val);
            node_spec_est_stream.close();
        }

        // arma::cx_mat global_val = arma::zeros<arma::cx_mat>(estimators[i]->n_rows, estimators[i]->n_cols);
        arma::cx_mat global_val = arma::zeros<arma::cx_mat>(val.n_rows, val.n_cols);
        MPI_Reduce(val.memptr(), global_val.memptr(), 2 * val.n_rows * val.n_cols, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);

        if (!my_id) {
            global_val /= (double)world_size;
            global_est_vals[i][block_num] = global_val;
        }

        if (each_block && !my_id) {
            est_stream.open(output_folder + "/" + estimator_names[i] + ".block." + std::to_string(block_num));
            write(est_stream, global_val);
            est_stream.close();
        }

        if (out_progressive && !my_id) {
            val = arma::zeros<arma::cx_mat>(global_est_vals[i][0].n_rows, global_est_vals[i][0].n_cols);
            for (int j = 0; j <= block_num; j++)
                val += global_est_vals[i][j];
            val /= block_num + 1;
            arma::mat err_re = arma::zeros<arma::mat>(val.n_rows, val.n_cols);
            arma::mat err_im = arma::zeros<arma::mat>(val.n_rows, val.n_cols);
            for (int j = 0; j <= block_num; j++) {
                arma::cx_mat ind_obs = global_est_vals[i][j];
                for (unsigned r = 0; r < val.n_rows; r++)
                    for (unsigned c = 0; c < val.n_cols; c++) {
                        err_re(r, c) +=
                            (val(r, c).real() - ind_obs(r, c).real()) * (val(r, c).real() - ind_obs(r, c).real());
                        err_im(r, c) +=
                            (val(r, c).imag() - ind_obs(r, c).imag()) * (val(r, c).imag() - ind_obs(r, c).imag());
                    }
            }
            if (block_num) {
                err_re /= block_num * (block_num + 1);
                err_im /= block_num * (block_num + 1);
            }
            err_re = arma::sqrt(err_re);
            err_im = arma::sqrt(err_im);
            progressive_stream.open(progressive_dir + "/" + estimator_names[i] + ".block." + std::to_string(block_num));
            write(progressive_stream, val, err_re, err_im);
            progressive_stream.close();
        }
    }
    block_num++;
    weight_per_block = 0.0;
}

void Output::finalize_output() {
    if (!my_id)
        for (int i = 0; i < estimators.size(); i++) {
            arma::cx_mat val = arma::zeros<arma::cx_mat>(global_est_vals[i][0].n_rows, global_est_vals[i][0].n_cols);
            for (int j = 0; j < block_num; j++)
                val += global_est_vals[i][j];
            val /= block_num;
            arma::mat err_re = arma::zeros<arma::mat>(val.n_rows, val.n_cols);
            arma::mat err_im = arma::zeros<arma::mat>(val.n_rows, val.n_cols);
            for (int j = 0; j < block_num; j++) {
                arma::cx_mat ind_obs = global_est_vals[i][j];
                for (unsigned r = 0; r < val.n_rows; r++)
                    for (unsigned c = 0; c < val.n_cols; c++) {
                        err_re(r, c) +=
                            (val(r, c).real() - ind_obs(r, c).real()) * (val(r, c).real() - ind_obs(r, c).real());
                        err_im(r, c) +=
                            (val(r, c).imag() - ind_obs(r, c).imag()) * (val(r, c).imag() - ind_obs(r, c).imag());
                    }
            }
            if (block_num > 1) {
                err_re /= block_num * (block_num - 1);
                err_im /= block_num * (block_num - 1);
            }
            err_re = arma::sqrt(err_re);
            err_im = arma::sqrt(err_im);
            progressive_stream.open(output_folder + "/" + estimator_names[i]);
            write(progressive_stream, val, err_re, err_im);
            progressive_stream.close();
        }
}
