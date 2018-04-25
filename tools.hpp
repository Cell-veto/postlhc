// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// various utilities.
#ifndef TOOLS_HPP_INCLUDED 
#define TOOLS_HPP_INCLUDED 

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
#include "ts.mk.hpp"

using std::string;
typedef const string &string_ref;

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

// redirect std::cout to a file
void redirect_cout (string_ref filename, bool append);

template <typename TYPE>
TYPE read_arg (const char **argv);


// float utilities

using std::pow;
using std::sqrt;

inline
double sq (double x)
{
    return x*x;
}

inline
double pow6 (double x)
{
    x *= x;
    return x*x*x;
}

inline
double fdivide (double x, double y)
{
    return x / y;
}


// vector math

template <size_t DIM>
using vector = std::array <double, DIM>;

template <size_t DIM>
vector <DIM> zero_vector ()
{
    vector <DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = 0.;
    return ret;
}

template <size_t DIM>
vector <DIM> unit_vector (size_t nz)
{
    vector <DIM> ret = zero_vector <DIM> ();
    ret[nz] = 1.;
    return ret;
}

template <size_t DIM>
vector <DIM> operator+ (const vector <DIM> &lhs, const vector <DIM> &rhs)
{
    vector <DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] + rhs[n];
    return ret;
}

template <size_t DIM>
vector <DIM> operator- (const vector <DIM> &lhs, const vector <DIM> &rhs)
{
    vector <DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] - rhs[n];
    return ret;
}

template <size_t DIM, typename VECTOR2>
vector <DIM> &operator*= (vector <DIM> &lhs, const VECTOR2 &rhs)
{
    for (size_t n = 0; n != DIM; ++n)
        lhs[n] *= rhs[n];
    return lhs;
}

template <size_t DIM>
vector <DIM> &operator*= (vector <DIM> &lhs, const double &rhs)
{
    for (size_t n = 0; n != DIM; ++n)
        lhs[n] *= rhs;
    return lhs;
}

template <size_t DIM, typename VECTOR2>
vector <DIM> operator* (const vector <DIM> &lhs, const VECTOR2 &rhs)
{
    vector <DIM> ret = lhs;
    return ret *= rhs;
}

template <size_t DIM>
vector <DIM> operator* (const double &lhs, const vector <DIM> &rhs)
{
    vector <DIM> ret = rhs;
    return ret *= lhs;
}

template <size_t DIM>
vector <DIM> operator/ (const vector <DIM> &lhs, double rhs)
{
    vector <DIM> ret;
    double scale = 1./rhs;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n] * scale;
    return ret;
}

template <size_t DIM>
double inner (const vector <DIM> &lhs, const vector <DIM> &rhs)
{
    double ret = lhs[0]*rhs[0];
    for (size_t n = 1; n != DIM; ++n)
        ret += lhs[n]*rhs[n];
    return ret;
}

template <size_t DIM, typename TYPE2>
vector <DIM> elementwise_prod (const vector <DIM> &lhs, const TYPE2 &rhs)
{
    vector <DIM> ret;
    for (size_t n = 0; n != DIM; ++n)
        ret[n] = lhs[n]*rhs[n];
    return ret;
}

template <size_t DIM>
double norm_sq (const vector <DIM> &v)
{
    return inner (v, v);
}

template <size_t DIM>
double norm (const vector <DIM> &v)
{
    return sqrt (norm_sq (v));
}

template <size_t DIM>
std::ostream &operator<< (std::ostream &os, const vector <DIM> &rhs)
{
    for (unsigned n = 0; n != DIM-1; ++n)
        os << rhs[n] << ' ';
    return os << rhs[DIM-1];
}

// misc small utilities
// FIXME these have namesakes in std:: which do THE OPPOSITE
// (swallow NaN silently)

inline
double fmin (double a, double b)
{
    if (a < b)
        return a;
    // NAN propagation
    return a==a ? b : a;
}

inline
double fmin (double a, double b, double c)
{
    return fmin (a, fmin (b, c));
}

inline
double fmax (double a, double b)
{
    if (a > b)
        return a;
    // NAN propagation
    return a==a ? b : a;
}

template <typename ITERATOR>
double fproduct (ITERATOR begin, const ITERATOR &end)
{
    double prod = 1.;
    while (begin != end)
        prod *= *begin++;
    return prod;
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

class RandomContext
{
public:
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
    std::mt19937_64 e_;
};

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
        min_ = std::numeric_limits <TYPE>::max ();
        max_ = std::numeric_limits <TYPE>::lowest ();
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

#endif /* TOOLS_HPP_INCLUDED */
