/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#    include <psapi.h>
#    include <windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) \
    || (defined(__APPLE__) && defined(__MACH__))
#    include <sys/resource.h>
#    include <unistd.h>

#    if defined(__APPLE__) && defined(__MACH__)
#        include <mach/mach.h>

#    elif (defined(_AIX) || defined(__TOS__AIX__)) \
        || (defined(__sun__) || defined(__sun)     \
            || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#        include <fcntl.h>
#        include <procfs.h>

#    elif defined(__linux__) || defined(__linux) || defined(linux) \
        || defined(__gnu_linux__)
#        include <cstdio>

#    endif

#else
#    error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

size_t getPeakRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) \
    || (defined(__sun__) || defined(__sun)     \
        || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return 0UL; /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return 0UL; /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1'024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) \
    || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    rusage rusage{};

    getrusage(RUSAGE_SELF, &rusage);
#    if defined(__APPLE__) && defined(__MACH__)
    return static_cast<size_t>(rusage.ru_maxrss);
#    else
    return static_cast<size_t>(rusage.ru_maxrss * 1'024L);  // NOLINT
#    endif

#else
    /* Unknown OS ----------------------------------------------- */
    return 0UL; /* Unsupported. */
#endif
}

size_t getCurrentRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    mach_task_basic_info info{};
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;

    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info),
                  &infoCount)
        != KERN_SUCCESS)
        return 0UL; /* Can't access? */

    return static_cast<size_t>(info.resident_size);

#elif defined(__linux__) || defined(__linux) || defined(linux) \
    || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    int64_t rss = 0L;

    FILE* fp = nullptr;
    // Read VmRSS
    if ((fp = fopen("/proc/self/statm", "r")) == nullptr)  // NOLINT
        return 0UL;                                        /* Can't open? */

    if (fscanf(fp, "%*s%ld", &rss) != 1)  // NOLINT
    {
        fclose(fp);  // NOLINT
        return 0UL;  /* Can't read? */
    }

    fclose(fp);  // NOLINT
    return static_cast<size_t>(rss * sysconf(_SC_PAGESIZE));

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return 0UL; /* Unsupported. */
#endif
}
