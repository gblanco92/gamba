/*   GamBa: a Groebner basis engine
 *   Copyright (C) 2025 Guillem Blanco
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

#include <unistd.h>

#include "params.hpp"

namespace gamba
{

#define CSI   "\x1B["  // Control Sequence Introducer (ANSI spec name)
#define CLEAR "\x1b[0m"

constexpr char const RED[]          = CSI "0;31m";
constexpr char const RED_BOLD[]     = CSI "1;31m";
constexpr char const GREEN[]        = CSI "0;32m";
constexpr char const GREEN_BOLD[]   = CSI "1;32m";
constexpr char const YELLOW[]       = CSI "0;33m";
constexpr char const YELLOW_BOLD[]  = CSI "1;33m";
constexpr char const BLUE[]         = CSI "0;34m";
constexpr char const BLUE_BOLD[]    = CSI "1;34m";
constexpr char const MAGENTA[]      = CSI "0;35m";
constexpr char const MAGENTA_BOLD[] = CSI "1;35m";
constexpr char const CYAN[]         = CSI "0;36m";
constexpr char const CYAN_BOLD[]    = CSI "1;36m";
constexpr char const WHITE[]        = CSI "0;37m";
constexpr char const WHITE_BOLD[]   = CSI "1;37m";

class log
{
public:
    enum level
    {
        ERROR = -3,
        WARN,
        DEBG,
        INFO0,
        INFO1,
        INFO2,
    };

    template <typename... T>
    static void print(level const lvl,
                      fmt::format_string<T...> fmt,
                      T&&... args);

    inline static bool colorize;
};

template <char const* code, typename T>
struct colored_type
{
    T const& value;
};

template <char const* code, typename T>
colored_type<code, T> colored(T const& value)
{
    colored_type<code, T> colored = {value};
    return colored;
}

}  // namespace gamba

namespace fmt
{

template <char const* code, typename T>
struct fmt::formatter<gamba::colored_type<code, T>> : formatter<T>
{
    auto format(gamba::colored_type<code, T> c, format_context& ctx) const
        -> format_context::iterator
    {
        if (gamba::log::colorize)
            fmt::format_to(ctx.out(), "{}", code);

        auto ctx_it = formatter<T>::format(c.value, ctx);

        if (gamba::log::colorize)
            fmt::format_to(ctx.out(), "{}", CLEAR);

        return ctx_it;
    }
};

}  // namespace fmt

namespace gamba
{

template <typename... T>
void log::print(level const lvl, fmt::format_string<T...> fmt, T&&... args)
{
    /* check verbosity level */
    if (lvl > params::verbose)
        return;

    std::FILE* const fp = lvl < 0 ? stderr : stdout;

    /* if output file is not a TTY do not colorize */
    log::colorize = ::isatty(::fileno(fp));

    switch (lvl)
    {
        case DEBG:
            fmt::print(fp, "{}", colored<BLUE_BOLD>("[DEBUG] "));
            break;
        case WARN:
            fmt::print(fp, "{}", colored<YELLOW_BOLD>("[WARNING] "));
            break;
        case ERROR:
            fmt::print(fp, "{}", colored<RED_BOLD>("[ERROR] "));
            break;
        default:
    }

    fmt::print(fp, fmt, std::forward<T>(args)...);
}

}  // namespace gamba
