// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef ECMC_HPP_INCLUDED
#define ECMC_HPP_INCLUDED

#include "storage.hpp"

// each interaction is encoded by a class derived from Interaction.
// you need to override at least the methods sr_repulsion_range
// and random_sr_repulsion.
struct Interaction
{
    // REPULSIVE SHORT-RANGE INTERACTIONS
    // see IPL for an example.
    static constexpr
    double sr_repulsion_range ()
    {
        return 0;
    }

    // compute next event distance
    // specifically compute the change in rsq (distance-squared)
    // which is negative for repulsive potentials.
    double random_sr_repulsion (double /* rsq_now */, RandomContext *)
    {
        std::cerr << "random_sr_repulsion dummy called" << ABORT;
    }

    // ATTRACTIVE SHORT-RANGE INTERACTIONS
    // default: no attractive interactions.  see LennardJones for an example.
    static constexpr
    unsigned sr_attraction_range ()
    {
        return 0;
    }

    // compute next event distance
    // specifically compute the change in rsq (distance-squared)
    // which is positive for repulsive potentials.
    double random_sr_attraction (double /* rsq_now */, RandomContext *)
    {
        std::cerr << "random_sr_attraction dummy called" << ABORT;
    }

    // large distance, signals escape to infinity
    static constexpr
    double ESCAPES_TO_INFINITY ()
    {
        return 1e300;
    }

    // LONG-RANGE INTERACTIONS
    // rate of probe events, see IPL for an example
    // default: no probes, just short-range events
    static constexpr
    unsigned total_probe_rate (unsigned /* direction */)
    {
        return 0;
    }

    template <typename VECTOR>
    double random_probe (VECTOR *r, unsigned direction, RandomContext *)
    {
        (void)r;
        (void)direction;
        std::cerr << "random_probe dummy called" << ABORT;
    }

    template <typename VECTOR>
    double lr_event_rate (const VECTOR &r, unsigned direction)
    {
        (void)r;
        (void)direction;
        return 0;
    }

    template <typename VECTOR>
    unsigned probe_rate (const VECTOR &r, unsigned direction)
    {
        (void)r;
        (void)direction;
        return 0;
    }

    // return true if interaction is not simple pair potential (ECMC pressure
    // formula does not hold)
    static constexpr
    bool poison_xdisp_pressure ()
    {
        return false;
    }

    // default: do nothing
    void notify_random_context (RandomContext *)
    {
    }

    // default: do nothing
    void notify_error_bound (const AbstractStorage *, double /* error_bound */)
    {
    }
};

template <typename INTERACT, typename ENCODING>
struct ChainRunner : AbstractChainRunner
{
    typedef INTERACT inter_t;
    typedef ENCODING encoding_t;
    typedef Storage <encoding_t> stor_t;
    static const unsigned DIM = encoding_t::DIM;
    typedef typename stor_t::vector_t vector_t;
    typedef typename stor_t::key_t key_t;

    double calib_constant;
    inter_t inter;
    RandomContext random;

    // motion data
    key_t active;
    unsigned direction;
    double disp_left;

    // currently scheduled event
    key_t planned_next;
    double planned_disp, planned_xdisp;
#ifdef XDISP_HISTO
    double planned_rsq_event;
#endif
    std::array <double, MAX_DIM> reported_probe_rate;

    ChainRunner ()
    {
        calib_constant = 1.;
        reported_probe_rate.fill (-1.);
    }

    static
    stor_t *downcast (AbstractStorage *stor_)
    {
        // convert to reference to force bad_cast exceptions
        return &dynamic_cast <stor_t &> (*stor_);
    }

    virtual
    void seed_random (unsigned seed)
    {
        random.seed (seed);
    }

    virtual
    void set_parameter (string_ref name, double value)
    {
        inter.set_parameter (name, value);
        AbstractChainRunner::set_parameter (name, value);
    }

    virtual
    void calibrate (AbstractStorage *stor_)
    {
        stor_t *stor = downcast (stor_);
        double chexp = measure_chain_expansion (stor, 0);
        double Lmax = stor->periods ().max ();
        calib_constant = Lmax / chexp / 2;
    }

    virtual
    void do_collide (AbstractStorage *stor_, double disp_per_particle)
    {
        stor_t *stor = downcast (stor_);
        unsigned num_chains = 1 + stor->num_particles () * disp_per_particle / calib_constant;
        this->run (stor, num_chains, calib_constant);
    }

    // norm-squared, skipping the component 'skip'
    template <size_t DIM>
    static
    double norm_ortho_sq (const vector <DIM> &v, size_t skip)
    {
        assert (skip < DIM);
        double ret = 0.;
        for (size_t i = 0; i != DIM; ++i)
        {
            if (i != skip)
                ret += v[i]*v[i];
        }
        return ret;
    }

    // after the planned_* fields have been initialized to a nil event,
    // advancing the active particle by planned_disp without any collisions,
    // factor in the short-range (including hard-core) events. updates planned_*.
    void find_sr_events (stor_t *stor)
    {
        // FIXME these should, in principle, depend on the active particle.
        const double attract_range = inter.sr_attraction_range ();
        const double repuls_range = inter.sr_repulsion_range ();
        const double strip_width = fmax (repuls_range, attract_range);

        // due to the design of the Storage, we cannot allow the enumerate to wrap
        double max_disp = stor->strip_max_extent (direction) - strip_width - attract_range;
        if (max_disp < planned_disp)
            planned_disp = max_disp;

        // loop over all particles in strip
        auto g = stor->enumerate_strip (active, direction, strip_width,
            -attract_range, planned_disp + strip_width);
        for (; g.not_done (); g.next ())
        {
            vector_t r_now = stor->distance_vector (g.key (), active);
            double x_now = r_now[direction];
            double xsq = sq (x_now);
            double ortho_rsq = norm_ortho_sq (r_now, direction);
            double x_min = x_now - planned_disp;

            double disp;

            // repulsive interactions
            double rsq_min = ortho_rsq + (x_min > 0. ? sq (x_min) : 0.);
            if (x_now > 0. && rsq_min < sq (repuls_range))
            {
                double rsq = ortho_rsq + xsq;
                double delta = inter.random_sr_repulsion (rsq, &random);
                assert (delta <= 0.);
                // delta <= -rsq signals 'no event' (will not enter the next if clause)

                if (xsq + delta > 0.)
                {
                    disp = -delta / (sqrt (xsq + delta) + x_now);
                    goto have_event;
                }
            }

            // attractive short-range interactions
            // these always occur at displacements later than repulsive ones.
            // FIXME more aggressive pruning
            if (x_min < 0. && ortho_rsq < sq (attract_range))
            {
                // note that even if the other particle is ahead _now_, we might
                // schedule an event after passing it.
                double x_virt = (x_now > 0.) ? 0. : x_now;
                double xsq_virt = sq (x_virt);
                double rsq_virt = ortho_rsq + xsq_virt;
                double delta = inter.random_sr_attraction (rsq_virt, &random);
                assert (delta >= 0.);
                disp = (x_now-x_virt) + delta / (sqrt (xsq_virt + delta) - x_virt);
                goto have_event;
            }

            // no event, try next particle
            continue;

        have_event:
            // when jumping here, disp will have been initialized.
            // enter new event into bookkeeping
            ++shortrange_predicts;

            if (! (disp >= 0.))
            {
                std::cerr << "computed event is not in the future" << ABORT;
            }

            // see if this event preempts the previous one
            if (disp < planned_disp)
            {
                planned_disp = disp;
                planned_xdisp = x_now - disp;
                planned_next = g.key ();
#ifdef XDISP_HISTO
                planned_rsq_event = sq (planned_xdisp) + ortho_rsq;
#endif
                g.clip (active, disp + strip_width);
            }
        }
    }

    // after all short-range events have been factored into the planned_* fields,
    // find any long-range events.
    // updates planned_*, and returns true if a long-range event takes precendence,
    // false otherwise.
    bool find_lr_events (stor_t *stor)
    {
        // find the probe rate
        auto probe_rate = inter.total_probe_rate (direction);
        if (probe_rate == 0)
            return false;

        probe_rate *= stor->cell_max_density ();
        if (probe_rate != reported_probe_rate[direction])
            std::cerr << "total_probe_rate "
                << (reported_probe_rate[direction] = probe_rate) << " "
                << direction << "\n";

        double lr_disp = 0.;

        for (;;)
        {
            lr_disp += random.exponential (probe_rate);
            if (lr_disp >= planned_disp)
                return false;

            vector_t r;
            double qselect = inter.random_probe (&r, direction, &random);
            r[direction] += lr_disp;
#ifdef DEBUG
            vector_t r_before = r;
#endif
            key_t k;
            if (stor->probe (&k, active, &r, &random) && k != active)
            {
                ++longrange_predicts;
#ifdef DEBUG
                double error = norm (r - r_before);
                assert (error < stor->cell_diagonal ());
#endif
                // somebody's home
                r[direction] -= lr_disp;
                double qevent = inter.lr_event_rate (r, direction);
                if (! (qevent <= qselect))
                {
                    std::cerr << "fatal: event rate is higher than probe rate\n"
                        << "qevent " << qevent << " qselect " << qselect << "\n"
                        << "direction " << direction << "\n"
                        << "distance " << r << " = " << norm (r) << "\n" << ABORT;
                }

                double paccept = qevent/qselect;
                probe_paccept.add (paccept);
                if (random.real () < paccept)
                {
                    // LR event
                    planned_disp = lr_disp;
                    planned_xdisp = r[direction];
                    planned_next = k;
#ifdef XDISP_HISTO
                    planned_rsq_event = norm_sq (r);
#endif
                    return true;
                }
            }
        }
    }

    void run (AbstractStorage *stor_,
        unsigned num_chains, double chain_disp)
    {
        assert (chain_disp > 0.);
        stor_t *stor = downcast (stor_);
        inter.notify_random_context (&random);

        for (unsigned c = 0; c != num_chains; ++c)
        {
            active = stor->random_particle (&random);
            direction = random.uint (DIM);
            disp_left = chain_disp;
#ifdef XDISP_HISTO
            double saved_disp = 0.;
#endif

            while (disp_left > 0.)
            {
                planned_next = active;
                planned_disp = disp_left;
                planned_xdisp = 0.;

                inter.notify_error_bound (stor, stor->cell_diagonal ());

                find_sr_events (stor);

                if (find_lr_events (stor))
                {
                    // there can be no events preempting a LR event
                    ++longrange_lifts;
                    goto handle_event;
                }

                if (active != planned_next)
                {
                handle_event:
                    ++total_lifts;
                    total_xdisp += planned_xdisp;
#ifdef XDISP_HISTO
                    xdisp_histo.add (planned_xdisp);
                    disp_histo.add (saved_disp + planned_disp);
                    saved_disp = 0.;
                    log_revent_histo.add (.5 * std::log (planned_rsq_event));
                    revent_histo.add (sqrt (planned_rsq_event));
#endif
                }
                else
                {
#ifdef XDISP_HISTO
                    saved_disp += planned_disp;
#endif
                }

                // displace actice particle
                disp_left -= planned_disp;
                active = stor->displace (active, direction, planned_disp, planned_next);
            }

            // end-of-chain event
            ++total_lifts;
            ++total_chains;
#ifdef XDISP_HISTO
            disp_histo.add (saved_disp);
#endif
        }

        report_xdisp_pressure &= !inter.poison_xdisp_pressure ();
        total_disp += num_chains * chain_disp;
    }

    double measure_chain_expansion (stor_t *stor, unsigned direction_, unsigned num_samples = 100)
    {
        inter.notify_random_context (&random);
        inter.notify_error_bound (stor, stor->cell_diagonal ());

        direction = direction_;
        disp_left = stor->periods () [direction_];

        // estimate median of disp and xdisp
        std::vector <double> ds, xs;

        for (unsigned n = 0; n != num_samples; ++n)
        {
            planned_next = active = stor->random_particle (&random);
            planned_disp = disp_left;
            planned_xdisp = 0.;

            find_sr_events (stor);
            find_lr_events (stor);

            ds.push_back (planned_disp);
            xs.push_back (planned_xdisp);
        }

        std::sort (ds.begin (), ds.end ());
        std::cerr << "calib_stat disp " << ds[0] << " " << ds[num_samples/2]
            << " " << ds.back () << "\n";
        std::sort (xs.begin (), xs.end ());
        std::cerr << "calib_stat xdisp " << xs[0] << " " << xs[num_samples/2]
            << " " << xs.back () << "\n";

        double ret = 1. + xs[num_samples/2]/ds[num_samples/2];
        std::cerr << "calib_stat chexp " << ret << "\n";
        return (ret < 0.) ? 1. : ret;
    }

    /*
set format "%L"; set log xy; hypot (x,y) = sqrt (x*x+y*y); p "pt.out" index 4 u (hypot($1,$2)):3 t "actual probes" , "" index 0:1 u (hypot($1,$2)):3 t "probe rate" , "" index 2:3 u (hypot($1,$2)):3 t "LR event rate"
    */

    void pt_dump (const vector_t &r, double value)
    {
        std::cout << r[0] << ' ' << r[1] << ' '
            << value * sphere_surface (DIM, norm (r)) << '\n';
    }

    virtual
    void probe_test_pattern (AbstractStorage *stor_, unsigned direction)
    {
        stor_t *stor = downcast (stor_);
        RandomContext random;
        inter.notify_random_context (&random);
        inter.notify_error_bound (stor, stor->cell_diagonal ());

        vector <DIM> r = zero_vector <DIM> ();

        std::cout << std::scientific;

        // initialize
        const double box_size = 2*stor->periods ().max ();
        const double tpr = inter.total_probe_rate (direction);
        const double big_box_size = 1e20 * box_size;
        const double step = box_size / 100;

        // reference curve
        std::cerr << "probe rate\n";
        std::cout << "# probe rate\n";
        for (r[0] = -box_size; r[0] < box_size; r[0] += step)
        {
            for (r[1] = -box_size; r[1] < box_size; r[1] += step)
                pt_dump (r, inter.probe_rate (r, direction));
            std::cout << '\n';
        }
        std::cout << "\n";
        for (r[0] = r[1] = .0001; r[0] < big_box_size; r[0] *= -1.01)
            pt_dump (r, inter.probe_rate (r, direction));
        std::cout << "\n\n";

        // what should be
        std::cerr << "LR event rate\n";
        std::cout << "# LR event rate\n";
        for (r[0] = -box_size; r[0] < box_size; r[0] += step)
        {
            for (r[1] = -box_size; r[1] < box_size; r[1] += step)
                pt_dump (r, inter.lr_event_rate (r, direction));
            std::cout << '\n';
        }
        std::cout << "\n";
        for (r[0] = r[1] = .0001; r[0] < big_box_size; r[0] *= -1.01)
            pt_dump (r, inter.lr_event_rate (r, direction));
        std::cout << "\n\n";

        // make a histogram
        std::cerr << "probe samples\n";
        std::cout << "# histogram of sampled probes\n";
        Histogram rprobe_histo (step/100, 0., 1e3 * box_size);
        static int num_warnings_left = 10;
        for (unsigned i = 0; i != 1000000; ++i)
        {
            double reported = inter.random_probe (&r, direction, &random);
            double reference = inter.probe_rate (r, direction);

            (void)reported;
            assert (reported > 0.);
            assert (fabs (reported-reference) / reported < 1e-6);

            if (num_warnings_left && inter.lr_event_rate (r, direction) > reference)
            {
                --num_warnings_left;
                std::cerr << "warn: probe rate too small\n"
                    << reference << " vs "
                    << inter.lr_event_rate (r, direction)
                    << " for r = (" << r << "), rnorm = " << norm (r) << "\n";
            }
            rprobe_histo.add (norm (r));
        }

        // dump it
        r = zero_vector <DIM> ();
        for (unsigned i = 0; i != rprobe_histo.size (); ++i)
        {
            double rate = rprobe_histo.bin_density (i) * tpr;
            if (rate == 0.) continue;
            r[0] = rprobe_histo.bin_center (i);
            rate /= sphere_surface (DIM, r[0]);
            pt_dump (r, rate);
        }
    }
};

#endif /* ECMC_HPP_INCLUDED */
