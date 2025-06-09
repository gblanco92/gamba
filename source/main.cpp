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
#include <limits>
#include <stdexcept>
#include <string>

#include <CLI/CLI.hpp>

#include "config.hpp"
#include "cpuinfo.h"
#include "gamba.hpp"
#include "git_version.h"
#include "logger.hpp"
#include "utils.hpp"

#ifdef PROFILE_GAMBA
#    include "profile.hpp"
#endif

namespace
{

void add_command_line_options(CLI::App& app,
                              std::string& input_file,
                              std::string& output_file)
{
    app.get_formatter()->column_width(40);

    app.add_option("-i,--input-file", input_file, "Input filename")
        ->required()
        ->check(CLI::ExistingFile.description(""));

    app.add_option("-o,--output-file", output_file, "Output filename");

    app.add_option("-e,--num-elim", gamba::params::num_elim_vars,
                   R"(Number of variables in the first elimination block
Must be an integer in the range [0, #variables))");

    app.add_option(
           "--max-spairs", gamba::params::max_spairs,
           R"(Max. number of pairs with min. degree selected in each round
Set max-spairs = 0 to select all pairs with minimal degree)")
        ->check(CLI::Range(0U, std::numeric_limits<uint32_t>::max())
                    .description(""));

    app.add_flag("--all-spairs", gamba::params::all_spairs,
                 "Select all pairs in the queue in each round");

    app.add_flag("--no-reduce", gamba::params::no_reduce,
                 "Do not compute a reduced Groebner basis")
        ->capture_default_str();

    app.add_option("-s,--seed", gamba::params::seed,
                   "Seed to initialize internal RNG")
        ->check(CLI::Range(static_cast<uint64_t>(0),
                           std::numeric_limits<uint64_t>::max())
                    .description(""));

    app.add_option("-t,--threads", gamba::params::num_threads,
                   "Number of threads to be used")
        ->check(CLI::Range(0U, std::numeric_limits<uint32_t>::max())
                    .description(""))
        ->capture_default_str();

    app.add_option("-v,--verbose", gamba::params::verbose, "Verbosity level")
        ->check(CLI::Range(std::numeric_limits<int16_t>::min(),
                           std::numeric_limits<int16_t>::max())
                    .description(""))
        ->capture_default_str();
}

void print_system_info()
{
    using namespace gamba;  // NOLINT

    auto const* const processor = cpuinfo_get_processor(0);

    log::print(log::INFO0, "{}{}", colored<CYAN_BOLD>("GamBa v"),
               colored<CYAN_BOLD>(GAMBA_VERSION));
    log::print(log::INFO0, " running on {}", processor->package->name);
    log::print(log::INFO0, " [seed = {}]\n", gamba::params::seed);

    GAMBA_DEVELOP(log::print(log::INFO0, "Git commit: {}\n", git_hash));

    GAMBA_DEBUG(log::print(log::DEBG, "{}\n", compiler_version_string()));
    GAMBA_DEBUG(log::print(log::DEBG, "{}\n", libgmp_version_string()));
    GAMBA_DEBUG(log::print(log::DEBG, "{}\n", libflint_version_string()));
}

}  // namespace

int main(int argc, char** argv)
{
    int ret_val = 0;
    std::string input_file;
    std::string output_file;

    CLI::App app{"GamBa: a Groebner basis engine"};
    add_command_line_options(app, input_file, output_file);
    CLI11_PARSE(app, argc, argv);

    try
    {
        /* initialize subsystems */
        cpuinfo_initialize();
        GAMBA_PROFILE(gamba::profiler::init());

        /* print banner */
        print_system_info();

        /* check for AVX2 instrunction set at runtime */
        if (not cpuinfo_has_x86_avx2())
            throw std::runtime_error("CPU does not support AVX2 instructions.");

        /* set priority of the process */
        setpriority(PRIO_PROCESS, 0, GAMBA_PRIORITY);

        /* read input file */
        std::ifstream infile{input_file};

        if (infile.fail())
            throw std::runtime_error("Error opening input file.");

        gamba::generators_data input_data;
        input_data.read(infile);

        // GAMBA_DEBUG(input_data.write(std::cerr);)

        /* sanitize options that depend on input file */
        gamba::params::sanitize_input(input_data.num_vars);

        /* do the actual computation */
        gamba::generators_data const output_data =
            gamba::groebner_basis(input_data);

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
        gamba::log::print(gamba::log::ERROR, "\n{}\n", excep.what());
        ret_val = -1;
    }
    catch (...)
    {
        gamba::log::print(gamba::log::ERROR, "\nUnexpected error.\n");
        ret_val = -2;
    }

    /* clean-up subsystems */
    cpuinfo_deinitialize();
    GAMBA_PROFILE(gamba::profiler::destroy());

    return ret_val;
}
