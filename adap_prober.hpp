// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef ADAPTIVE_PROBER_HPP_INCLUDED 
#define ADAPTIVE_PROBER_HPP_INCLUDED 

#include "tools.hpp"

#include <boost/format.hpp>
using boost::format;

class AdaptiveProber
{
public:
    void init (unsigned DIM, double /* rtol */,
        double rinner, double inner_paralpha,
        double rcut, double tail_paralpha)
    {
        binrate_.fill (0.);

        // set up radial bins
        binrad_[N_-1] = INFINITY;
        binexp_[N_-1] = tail_paralpha;
        unsigned i = N_-2;
        for (; i > 0 && rcut > rinner; --i)
        {
            binrad_[i] = rcut;
            binexp_[i] = inner_paralpha;
            rcut *= .9;
        }
        for (; i > 0; --i)
        {
            binrad_[i] = rcut;
            binexp_[i] = -1.; // const. rate
            rcut *= .9;
        }
        binrad_[0] = rcut;
        binexp_[0] = NAN;

        // fix bad exponents
        for (unsigned i = 1; i != N_; ++i)
        {
            double iexpon = DIM - (1.+binexp_[i]);
            if (iexpon == 0.)
            {
                std::cerr << "adap_prober: fixing singular exponent " << binexp_[i];
                binexp_[i] += .1;
                std::cerr << " -> " << binexp_[i] << "\n";
            }
        }

        assert (bin_from_radius (-1) == 0);
        assert (bin_from_radius (0) == 0);
        assert (bin_from_radius (1e-19) == 0);
        assert (bin_from_radius (1e19) == N_-1);
    }

    void calib_add (double r, double rate, double rtol = 0.)
    {
        if (rate == 0.) return;

        if (! (r >= 0.))
            std::cerr << "calib_add: negative distance" << ABORT;

        if (!std::isfinite (rate))
            std::cerr << "calib_add: nonfinite lifting rate" << ABORT;

        // range of bins affected (including imax)
        unsigned i =    bin_from_radius (r-rtol);
        unsigned imax = bin_from_radius (r+rtol);
        assert (imax < N_);

        // zero bin (const probe rate)
        if (i == 0)
        {
            binrate_[i] = fmax (binrate_[i], rate);
            ++i;
        }

        // outer bins (correct for appropriate exponent)
        for (; i <= imax; ++i)
        {
            double ir = fmin (binrad_[i], r+rtol); // FIXME??
            double irate = rate / pow (ir, -(1.+binexp_[i]));
            binrate_[i] = fmax (binrate_[i], irate);
        }
    }

    void calib_finish (unsigned DIM)
    {
        // fill in cumsum_
        cumsum_[0] = binrate_[0] * ball_volume (DIM, binrad_[0]);
        auto stat = format ("adap_prober_stat %4i %10.5f %10.5f || % 5g || %10.5f\n");
        std::cerr << boost::str (format (stat) % 0 % 0. % binrad_[0] % binexp_[0] % cumsum_[0]);

        for (unsigned i = 1; i != N_; ++i)
        {
            double iexpon = DIM - (1.+binexp_[i]);
            double w = pow (binrad_[i]/binrad_[i-1], iexpon);
            double tot = (w - 1) * pow (binrad_[i-1], iexpon);
            tot *= sphere_surface (DIM);
            tot /= iexpon;
            tot *= binrate_[i];

            std::cerr << boost::str (stat % i % binrad_[i-1] % binrad_[i] % binexp_[i] % tot);
            cumsum_[i] = cumsum_[i-1] + tot;
        }

        std::cerr << "adap_prober_stat total " << cumsum_[N_-1] << '\n';
    }

    double total_probe_rate () const
    {
        return cumsum_[N_-1];
    }

    template <size_t DIM>
    double probe_rate (const vector <DIM> &r_)
    {
        double r = norm (r_);
        unsigned i = bin_from_radius (r);
        if (i == 0)
            return binrate_[0];
        else if (i < N_)
            return binrate_[i] * pow (r, -(1.+binexp_[i]));

        std::cerr << "probe_rate: impossible" << ABORT;
        return 0.;
    }

    template <size_t DIM>
    double random_probe (vector <DIM> *ret, RandomContext *random)
    {
        double x = random->real (total_probe_rate ());
        unsigned i = bin_from_quantile (x);

        if (i == 0)
        {
            *ret = random->fromball <DIM> (binrad_[0]);

            // return probe rate
            return binrate_[0];
        }
        else if (i < N_)
        {
            // sample radial distance
            double iexpon = DIM - (1.+binexp_[i]);
            double w_min = pow (binrad_[i]/binrad_[i-1], iexpon);
            double w = 1 - random->real () * (1-w_min);
            double r = binrad_[i-1] * pow (w, 1. / iexpon);

            // add direction
            *ret = random->fromsphere <DIM> (r);

            // return probe rate
            return binrate_[i] * pow (r, -(1.+binexp_[i]));
        }

        std::cerr << "random_probe: impossible" << ABORT;
        return 0.;
    }

private:
    static const int N_ = 64;
    vector <N_> cumsum_;
    vector <N_> binrad_; // upper bounds of bins
    vector <N_> binexp_; // == pareto alpha for that bin
    vector <N_> binrate_;

    unsigned bin_from_radius (double r) const
    {
        auto low = binrad_.begin ();
        return std::lower_bound (low, binrad_.end (), r) - low;
    }

    unsigned bin_from_quantile (double x) const
    {
        assert (x >= 0.);
        auto low = cumsum_.begin ();
        return std::lower_bound (low, cumsum_.end (), x) - low;
    }
};

#endif /* ADAPTIVE_PROBER_HPP_INCLUDED */
