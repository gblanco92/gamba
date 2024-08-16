/*   GamBa: a Groebner basis engine
 *   Copyright (C) 2023 Guillem Blanco
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>. */

#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <CLI/CLI.hpp>

#include "debug.hpp"
#include "gamba.hpp"
#include "git_version.h"
#include "io.hpp"
#include "params.hpp"
#include "profile.hpp"

void add_command_line_options(CLI::App& app,
                              gamba::gamba_params& params,
                              std::string& input_file,
                              std::string& output_file)
{
    app.get_formatter()->column_width(40);

    app.add_option("-i,--input-file", input_file, "Input filename")
        ->required()
        ->check(CLI::ExistingFile.description(""));

    app.add_option("-o,--output-file", output_file, "Output filename");

    app.add_option("-e,--num-elim", params.num_elim_vars,
                   R"(Number of variables in the first elimination block
Must be a positive integer in the range [0, #variables))");

    app.add_option(
           "--max-spairs", params.max_spairs,
           R"(Maximum number of pairs with minimal degree selected in each round
Set max-spairs = 0 to select all pairs with minimal degree)")
        ->check(CLI::Range(0U, std::numeric_limits<uint32_t>::max())
                    .description(""));

    app.add_flag("--all-pairs", params.all_spairs,
                 "Select all pairs in the queue in each round");

    app.add_flag("--no-reduce", params.no_reduce,
                 "Do not compute a reduced Groebner basis")
        ->capture_default_str();

    app.add_option("-t,--threads", params.num_threads,
                   "Number of threads to be used")
        ->check(CLI::Range(0U, std::numeric_limits<uint32_t>::max())
                    .description(""))
        ->capture_default_str();

    app.add_option("-v,--verbose", params.verbose, "Verbosity level")
        ->check(CLI::Range(0U, std::numeric_limits<uint32_t>::max())
                    .description(""))
        ->capture_default_str();
}

void print_system_info()
{
    auto const* const processor = cpuinfo_get_processor(0);

    std::cout << "GamBa v" << GAMBA_VERSION;
    std::cout << " running on " << processor->package->name;  // NOLINT
}

int main(int argc, char** argv)
{
    std::string input_file;
    std::string output_file;
    gamba::gamba_params params{};

    cpuinfo_initialize();

    print_system_info();

#ifdef GAMBA_RELEASE
    std::random_device rnd;
    params.seed = rnd();
#endif
    std::cout << " [seed = " << params.seed << "]" << std::endl;

#ifndef GAMBA_RELEASE
    std::cout << "Git commit: " << git_hash << std::endl;
#endif

    GAMBA_DEBUG(gamba::print_compiler_version(gamba::debug);)
    GAMBA_DEBUG(gamba::print_libgmp_version(gamba::debug);)
    GAMBA_DEBUG(gamba::debug << std::endl;)

    GAMBA_PROFILE(gamba::perf.setup());

    CLI::App app{"GamBa: a Groebner basis engine"};

    add_command_line_options(app, params, input_file, output_file);

    CLI11_PARSE(app, argc, argv);

    try
    {
        /* check for AVX2 instrunction set at runtime */
        if (not cpuinfo_has_x86_avx2())
            throw std::runtime_error("CPU does not support AVX2 instructions.");

        /* read input file */
        std::ifstream infile{input_file};
        if (infile.fail())
            throw std::runtime_error("Error opening input file.");

        gamba::generators_data input_data;
        input_data.read(infile);

        GAMBA_DEBUG(input_data.write(gamba::debug);)
        GAMBA_DEBUG(gamba::debug << std::endl;)

        /* sanitize options that depend on input file */
        check_gamba_parameters(params, input_data);

        /* do the actual computation */
        gamba::generators_data const output_data =
            gamba::groebner_basis(input_data, params);

        /* write output file if available */
        if (not output_file.empty())
        {
            std::ofstream outfile{output_file};
            if (outfile.fail())
                throw std::runtime_error("Error opening output file.");

            output_data.write(outfile);
        }
    }
    catch (std::exception const& excep)
    {
        std::cerr << excep.what() << std::endl;

        return -1;
    }
    catch (...)
    {
        std::cerr << "Unexpected error." << std::endl;

        return -2;
    }

    return 0;
}
