// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef IPL_HPP_INCLUDED
#define IPL_HPP_INCLUDED

// U = stren / exponent / r**exponent
// note that the parameter 'strength' for historical reasons sets the
// conventional 'dimensionless interaction strength' Gamma == stren/exponent
//
// for exponent == 0, it takes the special form
// U = stren * log(1/r)      (for exponent = 0)

struct PowerLawInteraction : Interaction
{
    PowerLawInteraction ()
    {
        // we set beta = 1.  scale density accordingly!
        gamma = 1.;
        // sr_lr_split = cutoff for a truncated IPL.  otherwise, it's the
        // handshake distance between probing events and direct enumeration.
        sr_lr_split = 3.8;

        // built-in self-test, will leave exponent at 32 at the end
        double rsq_list[] = { 1e-6, 1e-2, 1., 1e2, 1e6 };
        double ex_list[] = { 1024., 6., 0., .5, /*-.5,*/ 1., 32. };
        for (double rsq : rsq_list)
        for (double ex : ex_list)
        {
            set_parameter ("exponent", ex);
            double delta = repulsive_lift (rsq, 1e-2);
            if (! (delta <= 0.))
            {
                std::cerr << "ipl: self-test error\n" << ABORT;
            }
        }
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

        // for compatibility with conventional normalization
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

    double random_sr_repulsion (double rsq, RandomContext *random)
    {
        double u_star = random->exponential (stren);
        return repulsive_lift (rsq, u_star);
    }
    
    double repulsive_lift (double rsq, double u_star)
    {
        // handle truncation of short-range potential
        if (rsq > sq (sr_lr_split))
        {
            double delta = repulsive_lift (sq (sr_lr_split), u_star);
            return (sq (sr_lr_split) - rsq) + delta;
        }

        if (exponent > 0.)
        {
            double x = u_star * exponent * pow (rsq, .5*exponent);
            assert (x >= 0.);

            if (x < 1e-8)
                // for very close particles, avoid loss of precision
                return -2./exponent * x * rsq;
            else
                return (pow (1+x, -2./exponent) - 1.) * rsq;
        }
        else if (exponent == 0.)
        {
            double x = 2 * u_star;
            assert (x >= 0.);

            if (x < 1e-8)
                // for very close particles, avoid loss of precision
                return -x * rsq;
            else
                return (exp (-x) - 1.) * rsq;
        }
        else
        {
            double x = u_star * exponent * pow (rsq, .5*exponent);
            assert (x <= 0.);

            if (x > -1e-8)
                // for very close particles, avoid loss of precision
                return -2./exponent * x * rsq;
            else if (x <= -1.)
                // active particle penetrates the target potential
                return -rsq;
            else
                return (pow (1+x, -2./exponent) - 1.) * rsq;
        }
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
