// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef PARETO_PROBER_HPP_INCLUDED
#define PARETO_PROBER_HPP_INCLUDED
    
struct ParetoProber
{
    ParetoProber ()
    {
        error_bound = 0.;
    }

    double error_bound;
private:
    double marginal_rate, zeroshell_radius;
    double total_pr, paralpha0;
    double towerp[3];

public:
    // configure the ParetoProber:  produce a flat probe rate in
    // a sphere of radius zeroshell_radius; outwards it decays smoothly
    // as a powerlaw.  everything is shifted a distance error_bound
    // outwards with respect to the actual cutoff to compensate for
    // discretization of the Storage.
    // returns the total probe rate.
    void setup (double int_dist, double exponent, unsigned DIM)
    {
        if (exponent <= DIM-1.)
            rt_error ("IPL exponent too small");

        // inner radius of zero shell would be rep_int_distance - error_bound
        // FIXME not implemented
        paralpha0 = exponent - (DIM-1);
        marginal_rate = pow (int_dist, -exponent-1);
        zeroshell_radius = int_dist + error_bound;

        double longrange_pr[3];
        if (DIM == 2)
        {
            longrange_pr[0] = 1.;
            longrange_pr[1] = error_bound;
            longrange_pr[2] = NAN;
        }
        else if (DIM == 3)
        {
            longrange_pr[0] = 1.;
            longrange_pr[1] = 2.*error_bound;
            longrange_pr[2] = sq (error_bound);
        }
        else
        {
            std::cerr << "DIM > 3 not implemented in pareto prober" << ABORT;
        }

        total_pr = 0.;
        static bool once = true;
        if (once) std::cerr << "total_pr ";
        for (unsigned n = 0; n != DIM; ++n)
        {
            // integrate out pareto distribution
            longrange_pr[n] *= pow (int_dist, -paralpha0-n) / (paralpha0+n);
            longrange_pr[n] *= sphere_surface (DIM);
            assert (longrange_pr[n] >= 0.);
            total_pr += longrange_pr[n];
            towerp[n] = total_pr;
            if (once) std::cerr << total_pr << ' ';
        }

        // zeroshell
        total_pr += marginal_rate * ball_volume (DIM, zeroshell_radius);
        if (once) std::cerr << total_pr << '\n';
        once = false;
    }

    double total_probe_rate () const
    {
        return total_pr;
    }

    // sample a LR probe. returns the probe rate with which this *ret would fire.
    template <size_t DIM>
    double random_probe (vector <DIM> *ret, RandomContext *random)
    {
        double tow = random->real (total_pr);

        // the long-range part is split up into DIM pareto's since we need
        // to expand the (r+error_bound)^(DIM-1) in the numerator.
        // tower sample: which exponent to use
        for (unsigned n = 0; n != DIM; ++n)
            if (tow < towerp[n])
            {
                double r = random->pareto (paralpha0+n,
                    zeroshell_radius-error_bound) + error_bound;
                *ret = random->fromsphere <DIM> (r);
                return pow (r-error_bound, -paralpha0-(DIM-1)-1);
            }

        // zeroshell event
        *ret = random->fromball <DIM> (zeroshell_radius);
        return marginal_rate;
    }

    template <size_t DIM>
    double probe_rate (const vector <DIM> &dist, unsigned direction)
    {
        (void)direction;  // not used, we sample azimuth uniformly
        double r = norm (dist);
        double s = fmax (zeroshell_radius, r);
        double ret = pow (s - error_bound, -paralpha0-(DIM-1)-1);
        assert (std::isfinite (ret) || ret == 0);
        return ret;
    }
};

#endif /* PARETO_PROBER_HPP_INCLUDED */
