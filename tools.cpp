// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// various utilities.
#include "tools.hpp"
#include <cstdlib>
#include <signal.h>
#include <fstream>
#include <stdexcept>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>

uint64_t gclock ()
{
    static uint64_t base = time (0);
    struct timeval tv;
    gettimeofday(&tv, 0);
    return (tv.tv_sec-base) * 1000000ull + tv.tv_usec;
}

static bool sigterm_installed = false, sigterm_caught_flag;

static void sigterm_handler (int)
{
    sigterm_caught_flag = true;
}

bool sigterm_caught ()
{
    if (!sigterm_installed)
    {
        sigterm_caught_flag = false;
        sigterm_installed = true;
        signal (SIGTERM, sigterm_handler);
    }

    return sigterm_caught_flag;
}

void shell (string_ref command)
{
    if (std::system (command.c_str ()) != 0)
        rt_error ("command failed: " + command);
}

bool file_is_readable (string_ref path)
{
    struct stat stat_result;
    if (stat (path.c_str (), &stat_result) != 0)
        return false;
    return !! (stat_result.st_mode & S_IRUSR);
}

void rt_error (string_ref msg)
{
    throw std::runtime_error (msg);
}

std::ostream &ABORT (std::ostream &os)
{
    os << std::endl;
    std::abort ();
}

static std::ofstream cout_logfile;

void redirect_cout (string_ref filename, bool append)
{
    std::ios::openmode mode = std::ios::out;
    if (append)
        mode |= std::ios::app;
    cout_logfile.open (filename, mode);
    if (!cout_logfile)
        rt_error ("error opening logfile" + filename);
    std::cerr.rdbuf (cout_logfile.rdbuf ());
}


// argument helpers

std::pair <string, string> int_pop_arg (const char **argv)
{
    assert (*argv);
    string kw = *argv++;
    if (!*argv)
        rt_error ("missing argument for kw " + kw);
    return std::make_pair (kw, string (*argv));
}

template <>
string read_arg (const char **argv)
{
    string kw, value;
    std::tie (kw, value) = int_pop_arg (argv);
    return value;
}

template <>
unsigned long read_arg (const char **argv)
{
    string kw, value;
    std::tie (kw, value) = int_pop_arg (argv);
    try
    {
        return std::stoul (value);
    }
    catch (...)
    {
        rt_error ("cannot convert value " + value + " for keyword " + kw);
    }
}

template <>
double read_arg (const char **argv)
{
    string kw, value;
    std::tie (kw, value) = int_pop_arg (argv);
    try
    {
        return std::stod (value);
    }
    catch (...)
    {
        rt_error ("cannot convert value " + value + " for keyword " + kw);
    }
}


// RandomContext

void RandomContext::seed (unsigned seed)
{
    e_.seed (seed);
}

// FIXME RandomContext::seed_from

double RandomContext::real ()
{
    std::uniform_real_distribution<> dis (0, 1);
    return dis (e_);
}

double RandomContext::real (double high)
{
    return high * real ();
}

double RandomContext::real (double low, double high)
{
    assert (low < high);
    std::uniform_real_distribution<> dis (low, high);
    return dis (e_);
}

unsigned RandomContext::uint (unsigned max_excluded)
{
    std::uniform_int_distribution <> dis (0, max_excluded-1u);
    return dis (e_);
}

// standard exponential (lambda = 1)
double RandomContext::exponential ()
{
    return -std::log (1-real ());
}

double RandomContext::exponential (double scale)
{
    assert (scale > 0.);
    return exponential () / scale;
}

double RandomContext::pareto (double paralpha, double minimum)
{
    assert (paralpha > 0.);
    assert (minimum > 0.);
    double u = 1 - real ();
    return pow (u, -1/paralpha) * minimum;
}

template <>
vector <2> RandomContext::unitbox ()
{
    return {{ real (-.5, .5), real (-.5, .5) }};
}

template <>
vector <3> RandomContext::unitbox ()
{
    return {{ real (-.5, .5), real (-.5, .5), real (-.5, .5) }};
}

template <>
vector <2> RandomContext::fromsphere <2> (double r)
{
    double phi = real (2*M_PI);
    return {{ r*std::cos (phi), r*std::sin (phi) }};
}

template <>
vector <3> RandomContext::fromsphere <3> (double r)
{
    double x = real (-1, 1);
    double y = real (-1, 1);
    double z = real (-1, 1);
    double s = pow (x*x + y*y + z*z, -0.5) * r;
    return {{ s*x, s*y, s*z }};
}

template <>
vector <2> RandomContext::fromball <2> (double r)
{
    double x = real () + real ();
    if (x > 1.) x = 2-x;
    return fromsphere <2> (x*r);
}

template <>
vector <3> RandomContext::fromball <3> (double r)
{
    double s = std::cbrt (real ());
    return fromsphere <3> (s*r);
}

template <size_t DIM>
vector <DIM> RandomContext::fromshell (double rmin, double rmax)
{
    assert (rmin >= 0.);
    assert (rmax > rmin);
    double dDIM = DIM;
    double u = real (pow (rmin, dDIM), pow (rmax, dDIM));
    return fromsphere <DIM> (pow (u, 1/dDIM));
}

// instantiate
template vector <2> RandomContext::fromshell (double, double);
template vector <3> RandomContext::fromshell (double, double);

Histogram::Histogram (double bin_wid, double low, double high)
    : data_ ((high-low)/bin_wid, 0), bin_wid_ (bin_wid), low_ (low), num_sampl_ (0)
{
}

Histogram::~Histogram ()
{
}

void Histogram::savetxt (string_ref filename) const
{
    std::ofstream ofs (filename);
    ofs.exceptions (std::ios::failbit | std::ios::badbit);
    savetxt (ofs);
}

void Histogram::savetxt (std::ostream &os) const
{
    for (unsigned i = 0; i != size (); ++i)
        os << bin_center (i) << ' ' << bin_density (i) << '\n';
}
