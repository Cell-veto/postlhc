// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// LJ interactions.
#ifndef LENNARD_JONES_INCLUDED 
#define LENNARD_JONES_INCLUDED 
    
#include "ecmc.hpp"
#include "pareto_prober.hpp"

// LennardJones
//
// Simulates a Lennard-Jones-type potential
// beta U(r) = strength * (rsix^2 - rsix)
// where rsix = (scale/rmod)^6
// and   rmod = min (r, cutoff * scale)
// and   beta = 1 / (k_B T)
//
// To obtain the conventional LJ form
//      U(r) = 4 epsilon ((sigma/r)^12 - (sigma/r)^6)
// at temperature T, set the parameters as
//      cutoff = infinity
//      strength = 4 epsilon / (k_B T)
//      scale = sigma
//
// The remaining parameter sr_lr_split sets the distance at which
// the long-range code takes over, in units of "scale".  It merely
// affects the computational efficiency of the simulation.


// minimum of LJ potential
static constexpr double LJ_MINIMUM = 1.122462048309373;

struct LennardJones : public Interaction
{
    LennardJones ()
    {
        strength = 1.;
        sr_lr_split = 3.8;
        scale = 2.;
        set_parameter ("cutoff", INFINITY);
    }

    double strength;
    double sr_lr_split, scale, cutoff;
    double rep_cutoff, attr_cutoff, rsix_attr_min;

    void set_parameter (string name, double value)
    {
        if (! (value > 0.))
            rt_error ("Invalid " + name + " value");

        if (name == "strength")
            strength = value;
        else if (name == "scale")
            scale = value;
        else if (name == "cutoff")
            cutoff = value;
        else if (name == "sr_lr_split")
            sr_lr_split = value;
        else
            rt_error ("Invalid parameter: " + name);

        // precompute some interesting quantities
        rep_cutoff  = fmin (LJ_MINIMUM, cutoff, sr_lr_split);
        attr_cutoff = fmin (cutoff, sr_lr_split);
        rsix_attr_min = pow (attr_cutoff, -6.);
    }

    double sr_repulsion_range () const
    {
        return rep_cutoff * scale;
    }

    double sr_attraction_range () const
    {
        return attr_cutoff < LJ_MINIMUM
            ? 0.
            : attr_cutoff * scale;
    }

    static
    double evaluate_lj (double rsq)
    {
        double rsix = pow (rsq, -3.);
        double e = rsix * (rsix - 1.);
        return e;
    }

    double random_sr_repulsion (double rsq_, RandomContext *random)
    {
        double rsq = fmin (rsq_ / sq (scale), sq (rep_cutoff));

        // add thermal energy increase on current interaction energy
        double e_now = evaluate_lj (rsq);
        double e_star = random->exponential (strength);
        double e_evt = e_now + e_star;

        // invert the LJ potential, picking the solution in repulsive region
        // FIXME look at corner cases
        double rsix = .5 + sqrt (.25 + e_evt);
        double rsq_evt = pow (rsix, -1./3);
        assert (rsq_evt <= sq (rep_cutoff));

        return rsq_evt * sq (scale) - rsq_;
    }

    double random_sr_attraction (double rsq_, RandomContext *random)
    {
        double rsq = fmax (rsq_ / sq (scale), sq (LJ_MINIMUM));

        // add thermal energy increase on current interaction energy
        double e_now = evaluate_lj (rsq);
        double e_star = random->exponential (strength);
        double e_evt = e_now + e_star;

        // invert the LJ potential, picking the solution in attractive region
        // FIXME look at corner cases
        double rsix = .5 - sqrt (.25 + e_evt);
        if (rsix <= rsix_attr_min)
            return -rsq_;  // escapes to infinity (no event)
        double rsq_evt = pow (rsix, -1./3);
        assert (rsq_evt <= sq (attr_cutoff) && rsq_evt >= sq (LJ_MINIMUM));

        return rsq_evt * sq (scale) - rsq_;
    }

    ParetoProber prober;

    void notify_error_bound (const AbstractStorage *,
        double err_bound, size_t DIM)
    {
        if (prober.error_bound != err_bound)
        {
            prober.error_bound = err_bound;
            prober.setup (sr_lr_split * scale, 6., DIM);
            std::cerr << "total_probe_rate " << total_probe_rate (0) << "\n";
        }
    }

    double probe_pref () const
    {
        return strength * 6. * pow6 (scale) * 5.;
    }

    double total_probe_rate (unsigned /* direction */) const
    {
        if (sr_lr_split > cutoff)
            return 0.;
        else
            return probe_pref () * prober.total_probe_rate ();
    }

    template <typename VECTOR>
    double probe_rate (const VECTOR &r, size_t direction)
    {
        return probe_pref () * prober.probe_rate (r, direction);
    }

    template <typename VECTOR>
    double lr_event_rate (const VECTOR &r, size_t direction)
    {
        double rsq = norm_sq (r) / sq (scale);

        if (rsq < sq (sr_lr_split) || rsq > sq (cutoff))
            return 0.;

        double rsix = pow (rsq, -3.);
        double rate = (12.*rsix - 6.) * r[direction];

        if (rate <= 0.)
            return 0.;
        else
            return rate * strength * rsix / (scale * scale * rsq);
    }

    template <typename VECTOR>
    double random_probe (VECTOR *ret, unsigned, RandomContext *random)
    {
        return probe_pref () * prober.random_probe (ret, random);
    }
};

#endif /* LENNARD_JONES_INCLUDED */
