// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef IPL_HPP_INCLUDED
#define IPL_HPP_INCLUDED

// FIXME
// U = stren / exponent / r**exponent
// U = stren * log(-r)      (for exponent = 0)

// note that Gamma... FIXME

struct PowerLawInteraction : Interaction
{
    PowerLawInteraction ()
    {
        // we set beta = 1.  scale density accordingly!
        gamma = 1.;
        set_parameter ("exponent", 32.);
        // sr_lr_split = cutoff for a truncated IPL.  otherwise, it's the
        // handshake distance between probing events and direct enumeration.
        sr_lr_split = 3.8;
    }

    double gamma;
    double stren, exponent;
    double sr_lr_split;

    void set_parameter (string_ref name, double value)
    {
        if (name == "strength")
        {
            if (! (value >= 0.))
                rt_error ("Invalid strength for power law");
            gamma = value;
        }
        else if (name == "exponent")
        {
            exponent = value;
        }
        else if (name == "sr_lr_split")
        {
            if (! (value > 0.))
                rt_error ("Invalid sr_lr_split for power law");
            sr_lr_split = value;
        }
        else
        {
            rt_error ("Invalid parameter: " + name);
        }

        if (exponent == 0)
            stren = gamma;
        else
            stren = gamma * fabs (exponent);
    }

    // raw potential, without the "stren" prefactor
    // input value rsq is r^2.
    double unit_potential (double rsq)
    {
        if (exponent == 0)
            return -.5 * std::log (rsq);
        else
            return pow (rsq, -.5 * exponent) / exponent;
    }

    template <typename VECTOR>
    double unit_potential (const VECTOR &vec)
    {
        return unit_potential (norm_sq (vec));
    }

    double invert_unit_potential (double e)
    {
        if (exponent == 0)
            return std::exp (-2.*e);
        else
        {
            e *= exponent;
            if (e < 0.)
                return -1.;
            return pow (e, -2./exponent);
        }
    }

    /*
    double potential (double rsq)
    {
        return stren * unit_potential (rsq);
    }
    */

    double unit_directional_derivative (double x, double rsq)
    {
        return x * pow (rsq, -.5 * (exponent+2));
    }

    double sr_repulsion_range () const
    {
        return sr_lr_split;
    }

    double random_repulsive_lift (double rsq, RandomContext *random)
    {
        rsq = fmin (rsq, sq (sr_lr_split));

        double u_now = unit_potential (rsq);
        double u_star = random->exponential (stren);
        double u_evt = u_now + u_star;
        double rsq_evt = invert_unit_potential (u_evt);

        if (! (rsq_evt <= rsq))
            std::cerr << rsq << " " << rsq_evt << ABORT;

        return rsq_evt;
    }
};

struct JelliumInteraction : PowerLawInteraction
{
    static constexpr
    bool poison_xdisp_pressure ()
    {
        // screening charge invalidates ECMC pressure formula
        return true;
    }
};

#endif /* IPL_HPP_INCLUDED */
