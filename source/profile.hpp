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

#pragma once

#include <string>

#ifdef PROFILE
#    define GAMBA_PROFILE(code) code
#else
#    define GAMBA_PROFILE(code)
#endif

/*
Profile parts of the execution using the following 'perf' flags:
    --control=fifo:{$m_ctl_filename},{$m_ack_filename}
    -p {$PID}
*/

namespace gamba
{

class perf_linux_profiler
{
    static void create_fifo(std::string const& filename);

public:
    /* do not use constructor to setup due to the global singleton */
    void setup();

    ~perf_linux_profiler();

    /*  call to start recording events with 'perf'*/
    void enable_profiling();

    /* call to stop recording events with 'perf' */
    void disable_profiling();

private:
    /* 0 file discriptor should correspond to stdin */
    int m_ctl_fd{0};
    int m_ack_fd{0};

    /* only to read the 'ack' string */
    char m_buffer[10]{};

    constexpr static char const* m_ctl_filename = "/tmp/perf_ctl.fifo";
    constexpr static char const* m_ack_filename = "/tmp/perf_ack.fifo";
};

#ifdef PROFILE
extern perf_linux_profiler perf;
#endif

}  // namespace gamba
