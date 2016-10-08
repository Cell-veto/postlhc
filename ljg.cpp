// (c) 2016 Sebastian Kapfer <sebastian.kapfer@fau.de>
//          Miriam Martinsons <miriam.martinsons@fau.de>, FAU Erlangen
#include "ecmc.hpp"

// LennardJonesGauss
//
// Lennard-Jones-Gauss potential
// beta U(r) = inv_temperature * f(min (r, cutoff))
// f(r) = r^{-12} - 2 r^{-6} - gauss_epsilon exp (-(r - gauss_r0)^2 / 2 gauss_sigma_sq)
//
// no long-range events are included, the potential is truncated at "cutoff"

struct LennardJonesGauss : Interaction
{
    // minimum of LJ part of the potential
    static constexpr double LJ_MINIMUM = 1.;

    LennardJonesGauss ()
    {
        inv_temperature = 5;
        gauss_epsilon = 1.8;
        gauss_r0 = 1.52;
        gauss_sigma_sq = .02;
        set_parameter ("cutoff", 3.6);
    }

    double inv_temperature;
    double gauss_epsilon, gauss_r0, gauss_sigma_sq;
    double cutoff;
    double pref_g;
    double g_at_cutoff, g_at_origin, lj_at_cutoff;

    void set_parameter (string name, double value)
    {
        if (! (value > 0.))
            rt_error ("Invalid " + name + " value");

        if (name == "temperature")
            inv_temperature = 1./value;
        else if (name == "inv_temperature")
            inv_temperature = value;
        else if (name == "cutoff")
            cutoff = value;
        else if (name == "gauss_epsilon")
            gauss_epsilon = value;
        else if (name == "gauss_r0")
            gauss_r0 = value;
        else if (name == "gauss_sigma_sq")
            gauss_sigma_sq = value;
        else
            rt_error ("Invalid parameter: " + name);

        if (gauss_r0 < LJ_MINIMUM || cutoff < gauss_r0)
            std::cerr << "Code is not correct in this regime" << ABORT;

        // precompute some interesting quantities
        lj_at_cutoff = evaluate_lj (sq (cutoff));
        pref_g = -.5 / gauss_sigma_sq;
        g_at_cutoff = evaluate_g (sq (cutoff));
        g_at_origin = evaluate_g (sq (0.));
    }

    double sr_repulsion_range () const
    {
        return gauss_r0;
    }

    double sr_attraction_range () const
    {
        return cutoff;
    }

    static
    double evaluate_lj (double rsq)
    {
        double rsix = pow (rsq, -3.);
        double e = rsix * (rsix - 2.);
        return e;
    }

    double evaluate_g (double rsq)
    {
        double r = sqrt (rsq);
        return -gauss_epsilon * exp (pref_g * sq (r-gauss_r0));
    }

    double random_repulsive_lift_lj (double rsq, RandomContext *random)
    {
        rsq = fmin (rsq, sq (LJ_MINIMUM));

        // add thermal energy increase on current interaction energy
        double e_now = evaluate_lj (rsq);
        double e_star = random->exponential (inv_temperature);
        double e_evt = e_now + e_star;

        // invert the LJ potential, picking the solution in repulsive region
        // FIXME this is naive and not sufficient for very small rsq.
        double rsix = 1. + sqrt (1. + e_evt);
        double rsq_evt = pow (rsix, -1./3);

        assert (rsq_evt <= sq (LJ_MINIMUM));
        return rsq_evt;
    }

    double random_repulsive_lift_g (double rsq, RandomContext *random)
    {
        rsq = fmin (rsq, sq (gauss_r0));

        // add thermal energy increase on current interaction energy
        double e_now = evaluate_g (rsq);
        double e_star = random->exponential (inv_temperature);
        double e_evt = e_now + e_star;

        if (e_evt >= g_at_origin)
        {
            // no event (thermal energy large enough)
            return -rsq;
        }
        else
        {
            // invert the Gauss, picking the solution in repulsive region
            double r_evt = gauss_r0 - sqrt (log (e_evt / (-gauss_epsilon)) / pref_g);
            assert (r_evt > 0.);
            return sq (r_evt);
        }
    }

    double random_sr_repulsion (double rsq, RandomContext *random)
    {
        return fmax (random_repulsive_lift_lj (rsq, random),
                     random_repulsive_lift_g  (rsq, random)) - rsq;
    }

    double random_attractive_lift_lj (double rsq, RandomContext *random)
    {
        rsq = fmax (rsq, sq (LJ_MINIMUM));

        // add thermal energy increase on current interaction energy
        double e_now = evaluate_lj (rsq);
        double e_star = random->exponential (inv_temperature);
        double e_evt = e_now + e_star;

        if (e_evt >= lj_at_cutoff)
        {
            // thermal energy too large, particle escapes to infinity
            return ESCAPES_TO_INFINITY ();
        }
        else
        {
            // invert the LJ potential, picking the solution in attractive region
            double rsix = 1. - sqrt (1. + e_evt);
            double rsq_evt = pow (rsix, -1./3);
            assert (rsq_evt <= sq (cutoff) && rsq_evt >= sq (LJ_MINIMUM));
            return rsq_evt;
        }
    }

    double random_attractive_lift_g (double rsq, RandomContext *random)
    {
        rsq = fmax (rsq, sq (gauss_r0));

        // add thermal energy increase on current interaction energy
        double e_now = evaluate_g (rsq);
        double e_star = random->exponential (inv_temperature);
        double e_evt = e_now + e_star;

        if (e_evt >= g_at_cutoff)
        {
            // thermal energy too large, particle escapes to infinity
            return ESCAPES_TO_INFINITY ();
        }
        else
        {
            // invert the Gauss, picking the solution in attractive region
            double r_evt = gauss_r0 + sqrt (log (e_evt / (-gauss_epsilon)) / pref_g);
            return sq (r_evt);
        }
    }

    double random_sr_attraction (double rsq, RandomContext *random)
    {
        return fmin (random_attractive_lift_lj (rsq, random),
                     random_attractive_lift_g  (rsq, random)) - rsq;
    }
};

static Register <ChainRunner <LennardJonesGauss, Monodisperse2D>>
    one ("ljg/mono2d");
