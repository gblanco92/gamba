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

/*
Profile parts of the execution using the following 'perf' flags:
    --control=fifo:{$m_ctl_filename},{$m_ack_filename}
    -p {$PID}
*/

namespace gamba
{

class profiler
{
    static void create_fifo(std::string const& filename);

public:
    static void init();

    static void destroy() noexcept;

    /*  call to start recording events with 'perf'*/
    static void enable_profiling();

    /* call to stop recording events with 'perf' */
    static void disable_profiling();

private:
    /* 0 file discriptor should correspond to stdin */
    static int m_ctl_fd;
    static int m_ack_fd;

    /* only to read the 'ack' string */
    static char m_buffer[10];

    constexpr static char const* m_ctl_filename = "/tmp/perf_ctl.fifo";
    constexpr static char const* m_ack_filename = "/tmp/perf_ack.fifo";
};

}  // namespace gamba
