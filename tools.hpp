// (c) 2015-2018 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// various utilities.
#pragma once

#include <array>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <climits>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include "vector.hpp"
#include "ts.mk.hpp"

#define RANDOM_EXPONENTIAL_BUFFER_SIZE 64

using std::string;
typedef const string &string_ref;

template <size_t DIM>
using vector = std::array <double, DIM>;
// import utility operators
using namespace vector_math_for_std_array;
using namespace vector_printing_for_std_array;

template <size_t DIM>
std::array <double, DIM> zero_vector ()
{
    std::array <double, DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = 0.;
    return ret;
}

template <size_t DIM>
std::array <double, DIM> unit_vector (size_t nz)
{
    std::array <double, DIM> ret = zero_vector <DIM> ();
    ret[nz] = 1.;
    return ret;
}

// return a microsecond-precision timestamp
uint64_t gclock ();

bool sigterm_caught ();

// run shell command.  throws if exitcode != 0
void shell (string_ref command);

// return true if the file exists and could probably be read
bool file_is_readable (string_ref path);
void rename_file (string_ref old_path, string_ref new_path);
bool remove_file (string_ref path, bool ignore_failure);

string hostname ();

// split string into space-separated words
std::vector <string> string_split (string_ref str);

// throw a std::runtime_error
void rt_error (string_ref msg) __attribute__((noreturn));

struct AbortObject { AbortObject () {} } const ABORT;
std::ostream &operator<< (std::ostream &, const AbortObject &) __attribute__((noreturn));

// machinery to redirect std::cerr to a logfile.
// no error checking against misuse is done.
// call this in the beginning
void buffer_cerr ();
// and one of these once you know whether you want to redirect
void redirect_cerr (string_ref filename, bool append);
void dont_redirect_cerr ();

template <typename TYPE>
TYPE read_arg (const char **argv);


// float utilities
using std::fmax;
using std::fmin;
using std::pow;
using std::sqrt;

// shorthand: cast to float and divide
inline
double fdivide (double x, double y)
{
    return x / y;
}

template <typename TYPE, size_t DIM, typename VEC2>
auto elementwise_fdivide (const std::array <TYPE, DIM> &lhs, const VEC2 &rhs)
    -> std::array <double, DIM>
{
    std::array <double, DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = double (lhs[n]) / double (rhs[n]);
    return ret;
}

inline
double sq (double x)
{
    return x*x;
}

// new name, same thing
inline
double fsq (double x)
{
    return x*x;
}

inline
double fmin (double a, double b, double c)
{
    return fmin (a, fmin (b, c));
}

template <typename ITERATOR>
double fproduct (ITERATOR begin, const ITERATOR &end)
{
    double prod = 1.;
    while (begin != end)
        prod *= *begin++;
    return prod;
}

template <size_t DIM>
// return a vector q1...qDIM such that
// (1) alpha = q1:p1 = q2:p2 = ... = qDIM:pDIM
// and
// (2) q1 * q2 * ... * qDIM = q.
vector <DIM> factor_quantity (double q, const vector <DIM> &p)
{
    double volume = fproduct (p.begin (), p.begin () + DIM);
    double alpha = pow (q / volume, 1./DIM);
    return alpha * p;
}

inline
void update_max (double *m, double val)
{
    *m = fmax (*m, val);
}

inline
void update_min (double *m, double val)
{
    *m = fmin (*m, val);
}


// hypersphere volume, ball volume
inline constexpr
double sphere_surface (unsigned DIM)
{
    return
        DIM==3 ? 4*M_PI :
        DIM==2 ? 2*M_PI :
        (std::abort (), NAN);
}

inline constexpr
double sphere_surface (unsigned DIM, double radius)
{
    return sphere_surface (DIM) * pow (radius, DIM-1);
}

inline constexpr
double ball_volume (unsigned DIM)
{
    return sphere_surface (DIM) / DIM;
}

inline constexpr
double ball_volume (unsigned DIM, double radius)
{
    return ball_volume (DIM) * pow (radius, DIM);
}

// factories

template <typename ABSTRACT>
class Factory
{
public:
    Factory (string_ref typecode)
    {
        registry () [typecode] = this;
    }

    virtual ABSTRACT *make () = 0;

    static ABSTRACT *make (string typecode)
    {
        Factory *f = registry () [typecode];
        if (!f)
            rt_error ("no factory to make product " + typecode);
        ABSTRACT *prod = f->make ();
        prod->typecode_ = typecode;
        return prod;
    }

private:
    static std::map <string, Factory *> &registry ()
    {
        static std::map <string, Factory *> themap;
        return themap;
    }
};

template <typename ABSTRACT>
class FactoryProduct
{
public:
    typedef ABSTRACT abstract_t;
    string_ref typecode () const { return typecode_; }
private:
    string typecode_;
    friend class Factory <ABSTRACT>;
};

template <typename CONCRETE>
class Register : Factory <typename CONCRETE::abstract_t>
{
    typedef typename CONCRETE::abstract_t abstract_t;
public:
    Register (string_ref typecode)
        : Factory <abstract_t> (typecode)
    {
        // base constructor registers the factory
    }

    virtual abstract_t *make ()
    {
        return new CONCRETE;
    }
};

// random numbers

// hack to make profiles more readable.
struct MersenneTwister : std::mt19937_64
{
};

struct RandomContext
{
    RandomContext ();
    void seed (unsigned);
    void seed_from (RandomContext *);
    double real (/* 0., 1. */);
    double real (/* 0., */ double high);
    double real (double low, double high);
    unsigned uint (unsigned max_excluded);
    double exponential ();
    double exponential (double scale);
    double pareto (double alpha, double min);

    template <size_t DIM>
    vector <DIM> unitbox ();

    template <size_t DIM>
    vector <DIM> fromsphere (double radius = 1.);
    template <size_t DIM>
    vector <DIM> fromball (double radius = 1.);
    template <size_t DIM>
    vector <DIM> fromshell (double rmin, double rmax);

private:
    MersenneTwister e_;
#if RANDOM_EXPONENTIAL_BUFFER_SIZE != 0
    void refill_exp_buffer_ ();
    double exp_buffer_[RANDOM_EXPONENTIAL_BUFFER_SIZE];
    unsigned num_exp_used_;
#endif
};

inline
double RandomContext::exponential (double scale)
{
    assert (scale > 0.);
#if RANDOM_EXPONENTIAL_BUFFER_SIZE != 0
    if (num_exp_used_ == RANDOM_EXPONENTIAL_BUFFER_SIZE)
        refill_exp_buffer_ ();
    return exp_buffer_[num_exp_used_++] / scale;
#else
    return exponential () / scale;
#endif
}

class Histogram
{
public:
    typedef size_t count_t;

    Histogram (double bin_wid, double low, double high);
    ~Histogram ();

    size_t size () const
    {
        return data_.size ();
    }

    void add (double value)
    {
        size_t i = (value-low_) / bin_wid_;
        if (i < size ())
            ++data_[i];
        ++num_sampl_;
    }

    double bin_center (size_t i) const
    {
        return bin_wid_ * (i+0.5) + low_;
    }

    double bin_density (size_t i) const
    {
        double norm = 1. / bin_wid_ / num_sampl_;
        return norm * data_.at (i);
    }

    void savetxt (string_ref ) const;
    void savetxt (std::ostream &) const;

private:
    std::vector <count_t> data_;
    const double bin_wid_, low_;
    count_t num_sampl_;
};

template <typename TYPE>
class Extrema
{
    TYPE min_, max_;
    size_t num_sampl_;
public:
    Extrema ()
    {
        reset ();
    }

    void reset ()
    {
        min_ = std::numeric_limits <TYPE>::max ();
        max_ = std::numeric_limits <TYPE>::lowest ();
        num_sampl_ = 0;
    }

    const TYPE &min () const { return min_; }
    const TYPE &max () const { return max_; }

    void add (const TYPE &value)
    {
        min_ = fmin (min_, value);
        max_ = fmax (max_, value);
        ++num_sampl_;
    }

    size_t num_samples () const
    {
        return num_sampl_;
    }
};
