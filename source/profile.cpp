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

#include <cstring>
#include <iostream>

#include "profile.hpp"

#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace gamba
{

/* use a global singleton for controling the perf profiler in linux systems
 * not the best approach but it will work for now */
#ifdef PROFILE
perf_linux_profiler perf{};
#endif

void perf_linux_profiler::create_fifo(std::string const& filename)
{
    struct stat sb
    {};

    /* check if fifo file exist */
    if (stat(filename.c_str(), &sb) != -1)
    {
        /* if exists remove it; this forces executing perf on a standalone
         * process since executing via perf opens (if exists) the fifo file
         * before the execution begins; unlinking creates a deadlock when
         * opening the fifo later on */
        unlink(filename.c_str());
    }

    /* create fifo file with permisions 0750 */
    if (mkfifo(filename.c_str(), S_IRWXU | S_IRGRP | S_IXGRP) == -1)
    {
        throw std::runtime_error{std::string{"Failed to create FIFO file: "}
                                 + filename};
    }
}

void perf_linux_profiler::setup()
{
    /* create control fifo */
    create_fifo(m_ctl_filename);

    /* create acknowledge fifo */
    create_fifo(m_ack_filename);

    /* use stderr just like 'perf' does */
    std::cerr << "## PROFILE ## " << "PID: " << getpid()
              << "; FIFO control files: " << m_ctl_filename << ","
              << m_ack_filename << std::endl;

    /* open control fifo write only */
    m_ctl_fd = open(m_ctl_filename, O_WRONLY | O_CLOEXEC);

    if (m_ctl_fd == -1)
    {
        throw std::runtime_error{std::string{"Failed to open FIFO file: "}
                                 + m_ctl_filename};
    }

    /* open acknowledge fifo read only */
    m_ack_fd = open(m_ack_filename, O_RDONLY | O_CLOEXEC);

    if (m_ack_fd == -1)
    {
        throw std::runtime_error{std::string{"Failed to open FIFO file: "}
                                 + m_ack_filename};
    }

    /* start with events disabled */
    disable_profiling();
}

/* do the cleanup in the destructor so even in the event of an exception is
 * thrown we can make sure that the cleanup happens */
perf_linux_profiler::~perf_linux_profiler()
{
    if (m_ctl_fd)
    {
        close(m_ctl_fd);
        unlink(m_ctl_filename);
    }

    if (m_ack_fd)
    {
        close(m_ack_fd);
        unlink(m_ack_filename);
    }
}

void perf_linux_profiler::enable_profiling()
{
    static char const* enable = "enable";

    /* write 'enable' to control fifo file */
    if (write(m_ctl_fd, reinterpret_cast<void const*>(enable), strlen(enable))
        == -1)
    {
        throw std::runtime_error{"Failed to write to control FIFO file"};
    }

    static char const* ack = "ack";

    /* wait for acknoledgement before continue */
    while (read(m_ack_fd, reinterpret_cast<void*>(&m_buffer), 10))
    {
        if (strcmp(m_buffer, ack) != 0)  // NOLINT
            break;
    }
}

void perf_linux_profiler::disable_profiling()
{
    static char const* disable = "disable";

    /* write 'disable' to control fifo file */
    if (write(m_ctl_fd, reinterpret_cast<void const*>(disable), strlen(disable))
        == -1)
    {
        throw std::runtime_error{"Failed to write to control FIFO file"};
    }

    static char const* ack = "ack";

    /* wait for acknoledgement before continue */
    while (read(m_ack_fd, reinterpret_cast<void*>(&m_buffer), 10))
    {
        if (strcmp(m_buffer, ack) != 0)  // NOLINT
            break;
    }
}

}  // namespace gamba
