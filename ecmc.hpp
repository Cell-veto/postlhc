// Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#pragma once

#include "storage.hpp"
#include <boost/format.hpp>

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
    using vector_t = vector <DIM>;
    using uivector_t = std::array <uint64_t, DIM>;
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
    vector_t reported_probe_rate;
    vector_t target_epc, target_tdisppc;

    // cumulative statistics
    uint64_t walltime;
    vector_t cumulative_disp, cumulative_xdisp;
    uivector_t all_events, longrange_events, chainend_events;
    uint64_t longrange_predicts, shortrange_predicts;
#ifdef XDISP_HISTO
    Histogram xdisp_histo;
    Histogram disp_histo;
    Histogram log_revent_histo;
    Histogram revent_histo;
#endif

    ChainRunner ()
#ifdef XDISP_HISTO
        : xdisp_histo (.02, -10., 40.)
        , disp_histo (.001, 0., 3.)
        , log_revent_histo (.02, -5., 40.)
        , revent_histo (.02, 0., 8.)
#endif
    {
        calib_constant = 1.;
        reported_probe_rate.fill (-1.);
        target_epc.fill (1.);
        target_tdisppc.fill (1.);
        reset_statistics ();
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
    }

    virtual
    void reset_statistics ()
    {
        walltime = 0;
        all_events.fill (0);
        longrange_events.fill (0);
        chainend_events.fill (0);
        longrange_predicts = shortrange_predicts = 0;
        cumulative_disp.fill (0.);
        cumulative_xdisp.fill (0.);
    }

    static
    uint64_t sum_uivector (const uivector_t &v)
    {
        uint64_t ret = v[0];
        for (size_t i = 1; i != v.size (); ++i)
            ret += v[i];
        return ret;
    }

    virtual
    void dump_statistics (std::ostream &os)
    {
        using boost::format;

        // some of these are physics data, so increase precision
        os.precision (15);
        // dump performance metrics
        uint64_t num_events = sum_uivector (all_events);
        os
            << "# new-style logfile -- episode dump\n"
            << "perfstat walltime "
            << walltime*1e-6 << " seconds\n"
            << "perfstat num_events "
            << num_events << "\n"
            << "perfstat events_per_hour "
            << format ("%.3e\n") % (3600. * num_events / walltime * 1e6)
            << "perfstat ratio_lr_events "
            << fdivide (sum_uivector (longrange_events), num_events) << "\n"
            << "perfstat predicts_per_event "
            << fdivide (shortrange_predicts, num_events) << " sr "
            << fdivide (longrange_predicts, num_events) << " lr\n";

        // EC dynamics
        os
            << "ecstat cumulative_disp "
            << cumulative_disp << "\n"
            << "ecstat cumulative_xdisp "
            << cumulative_xdisp << "\n"
            << "ecstat mean_free_path "
            << elementwise_fdivide (cumulative_disp, all_events) << "\n"
            << "ecstat mean_xdisp "
            << elementwise_fdivide (cumulative_xdisp, all_events) << "\n"
            << "ecstat mean_events_per_chain "
            << elementwise_fdivide (all_events, chainend_events) << "\n"
            << "ecstat scaled_mean_events_per_chain "
            << elementwise_fdivide (all_events,
                elementwise_product (chainend_events, target_epc)) << "\n"
            << "ecstat mean_tdisp_per_chain "
            << elementwise_fdivide (cumulative_disp + cumulative_xdisp, chainend_events) << "\n"
            << "ecstat scaled_mean_tdisp_per_chain "
            << elementwise_fdivide (cumulative_disp + cumulative_xdisp,
                elementwise_product (chainend_events, target_tdisppc)) << "\n";

        // compressibility factor: beta pressure / rho
        if (!inter.poison_xdisp_pressure ())
        {
            os
                << "compressib_factor instantaneous "
                << elementwise_fdivide (cumulative_disp + cumulative_xdisp, cumulative_disp) << "\n";
        }

        // zero
        os << "# new-style logfile -- clearing counters\n";
        reset_statistics ();
    }

    virtual
    void save_histograms (string_ref prefix)
    {
#ifdef XDISP_HISTO
        xdisp_histo.savetxt (prefix + "xdisp.histo");
        disp_histo.savetxt (prefix + "disp.histo");
        revent_histo.savetxt (prefix + "revent.histo");
#else
        (void)prefix;
#endif
    }

    uint64_t time_for_collisions (stor_t *stor, double sr_lr_split, int num_samples)
    {
        set_parameter ("sr_lr_split", sr_lr_split);
        inter.notify_random_context (&random);
        inter.notify_error_bound (stor, stor->cell_diagonal ());

        uint64_t start_time = gclock ();

        direction = 0;
        disp_left = stor->periods ()[0];

        for (int n = 0; n != num_samples; ++n)
        {
            planned_next = active = stor->random_particle (&random);
            planned_disp = disp_left;
            planned_xdisp = 0.;

            find_sr_events (stor);
            find_lr_events (stor);
        }

        uint64_t end_time = gclock ();
        std::cerr << "info optimizing sr_lr_split " << sr_lr_split << " " << (end_time-start_time) << " usec\n";
        return end_time - start_time;
    }

    virtual
    void optimize_sr_lr_split (AbstractStorage *stor_)
    {
        stor_t *stor = downcast (stor_);

        const int num_samples = 1000;
        double try_splits[] = {  2.1,  1.6,  1.7,  2.3,  2.9,  1.8,  1.5,  3.6,
            0.8,  3.9,  3.5, 0.9,  3.4,  2. ,  1.2,  1.4,  1.9,  3.1,  3.7,  3.,
            3.3,  2.4, 2.2,  2.5,  2.8,  2.7,  1. ,  0.6,  0.7,  3.2,  1.3,  3.8,
            2.6, 1.1 };

        uint64_t smallest_time_so_far = time_for_collisions (stor, try_splits[0], num_samples);
        int best_split_so_far = 0;

        const int num_splits = sizeof(try_splits) / sizeof (try_splits[0]);

        for (int i = 1; i != num_splits; ++i)
        {
            uint64_t time_now = time_for_collisions (stor, try_splits[i], num_samples);
            if (time_now < smallest_time_so_far)
            {
                smallest_time_so_far = time_now;
                best_split_so_far = i;
            }
        }

        std::cerr << "info best sr_lr_split " << try_splits[best_split_so_far] <<"\n";
        set_parameter ("sr_lr_split", try_splits[best_split_so_far]);
    }

    static
    vector_t trunc_periods (stor_t *stor)
    {
        vector_t ret;
        for (size_t n = 0; n != DIM; ++n)
            ret[n] = stor->period (n);
        for (size_t n = DIM; n != MAX_DIM; ++n)
            assert (stor->period (n) == 0.);
        return ret;
    }

    virtual
    void calibrate (AbstractStorage *stor_)
    {
        stor_t *stor = downcast (stor_);
        double chexp = measure_chain_expansion (stor, 0);
        double Lmax = stor->periods ().max ();
        calib_constant = Lmax / chexp / 2;

        target_epc = factor_quantity (stor->num_particles (), trunc_periods (stor));
        target_tdisppc = trunc_periods (stor);
    }

    virtual
    void collide (AbstractStorage *stor_, double disp_per_particle)
    {
        walltime -= gclock ();
        stor_t *stor = downcast (stor_);
        unsigned num_chains = 1 + stor->num_particles () * disp_per_particle / calib_constant;
        this->run (stor, num_chains, calib_constant);
        walltime += gclock ();
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
            ++shortrange_predicts;

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

                if (! (delta <= 0.))
                {
                    std::cerr << "computed short-range repulsive event is not in the future" << ABORT;
                }

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
                if (! (delta >= 0.))
                {
                    std::cerr << "computed short-range attractive event is not in the future" << ABORT;
                }

                disp = (x_now-x_virt) + delta / (sqrt (xsq_virt + delta) - x_virt);
                goto have_event;
            }

            // no event, try next particle
            continue;

        have_event:
            // when jumping here, disp will have been initialized.
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
                    ++longrange_events[direction];
                    goto handle_event;
                }

                if (active != planned_next)
                {
                handle_event:
                    ++all_events[direction];
                    cumulative_xdisp[direction] += planned_xdisp;
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
            ++all_events[direction];
            cumulative_disp[direction] += chain_disp;
            ++chainend_events[direction];
#ifdef XDISP_HISTO
            disp_histo.add (saved_disp);
#endif
        }
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
