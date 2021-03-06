// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// a simple test program which exercises some of the long-range potentials available.
#include "storage.hpp"
#include <fstream>
#include <boost/format.hpp>

using boost::format;

static
string format_out (string_ref prefix, unsigned snap, string_ref type = "coords")
{
    return prefix + type + std::to_string (snap) + ".dat";
}

int main (int, const char **argv)
{
    buffer_cerr ();
    AbstractStorage *stor = 0;
    AbstractChainRunner *cr = 0;
    AbstractCorrelator *gofr = 0;

    string inter, prefix = "out/";
    string in_prefix;
    bool probe_test = false;
    bool skip_calib = false;
    bool recover = false;
    bool redirect_log = false;
    bool write_gofr = false;
    double snap_disp = 100;
    unsigned long random_seed = 42;
    unsigned long snap_target = 100000;
    long recover_from = -1;
    bool dont_optimize_sr_lr_split = false;

    // parse command line, initializing our tools along the way
    for (++argv; *argv; argv += 2)
    {
        string kw = *argv;

        if (kw == "load")
        {
            in_prefix = read_arg <string> (argv);
            in_prefix += '/';
        }
        else if (kw == "addprefix")
        {
            prefix += boost::str (format ("%s/") % read_arg <string> (argv));
        }
        else if (kw == "setprefix")
        {
            prefix = read_arg <string> (argv);
        }
        else if (kw == "log")
        {
            redirect_log = true;
            --argv; // no argument
        }
        else if (kw == "seed")
        {
            random_seed = read_arg <unsigned long> (argv);
        }
        else if (kw == "snap")
        {
            snap_target = read_arg <unsigned long> (argv);
        }
        else if (kw == "snap_disp")
        {
            snap_disp = read_arg <double> (argv);
        }
        else if (kw == "recover")
        {
            recover = true;
            --argv; // no argument
        }
        else if (kw == "recover_from")
        {
            recover_from = static_cast <long> (read_arg <unsigned long> (argv));
            recover = true;
        }
        else if (kw == "skip_calib")
        {
            skip_calib = true;
            --argv; // no argument
        }
        else if (kw == "probe_test")
        {
            probe_test = true;
            --argv; // no argument
        }
        else if (!stor)
        {
            // create storage
            if (kw == "stor")
                stor = make_storage (read_arg <string> (argv));
            else
                rt_error ("First set storage mechanism with keyword stor");

            prefix += boost::str (format ("%s/") % read_arg <string> (argv));
        }
        else if (kw == "gofr")
        {
            // we will write a density/density pair correlator
            write_gofr = true;
            --argv; // no argument
        }
        else if (!cr)
        {
            // as soon as we encounter 'inter', instantiate the ECMC algorithm
            if (kw == "inter")
            {
                inter = read_arg <string> (argv);
                cr = make_chainrunner (inter, stor);
                prefix += inter + '/';
            }
            else
                rt_error ("First set interaction type with keyword inter");
        }
        else if (kw == "sr_lr_split")
        {
            // we don't want a subdirectory for this parameter.
            double value = read_arg <double> (argv);
            cr->set_parameter (kw, value);
            dont_optimize_sr_lr_split = true;
        }
        else
        {
            // unrecognized argument -- assume it's an interaction parameter
            double value = read_arg <double> (argv);
            prefix += boost::str (format ("%s%f/") % kw % value);
            std::cerr << "inter_param " << kw << " " << value << "\n";
            cr->set_parameter (kw, value);
        }
    }

    if (!stor)
        rt_error ("Missing stor");
    if (!cr)
        rt_error ("Missing inter");
    if (in_prefix == "")
        rt_error ("Missing load");

    // make the output directory
    shell ("mkdir -p " + prefix);

    prefix += boost::str (format ("seed%lu_") % random_seed);

    // redirect std::cerr to logfile
    // (append if recovering)
    if (redirect_log)
        redirect_cerr (prefix + "log", recover);
    else
        dont_redirect_cerr ();

    std::cerr << "hostname " << hostname () << '\n';

    // load data
    stor->load_periods (in_prefix + "periods");
    unsigned snap = 0;
    string in_filename = in_prefix + "coords.dat";
    if (recover_from != -1)
    {
        // check if final snapshot exists already
        if (file_is_readable (format_out (prefix, snap_target)))
        {
            std::cerr << "recover: the final target snapshot already exists.  exiting.\n";
            return 0;
        }

        if (file_is_readable (format_out (prefix, recover_from)))
        {
            snap = recover_from + 1u;
            std::cerr << "recover: attempting to load snapshot " << (snap-1u) << "\n";
            in_filename = format_out (prefix, snap-1u);
        }
        else
        {
            std::cerr << "recover snapshot not found; aborting\n";
            return -1;
        }
    }
    else if (recover)
    {
        // check if final snapshot exists already
        if (file_is_readable (format_out (prefix, snap_target)))
        {
            std::cerr << "recover: the final target snapshot already exists.  exiting.\n";
            return 0;
        }

        // find the first snapshot which is _not_ on disk
        if (file_is_readable (format_out (prefix, 0)))
        {
            // bisection-style algorithm
            // loop invariant is: low is present, high is missing
            unsigned low = 0, high = snap_target;
            while (high - low > 1)
            {
                unsigned mid = (low+high) / 2;
                std::cerr << "recover: checking " << low << " "
                    << mid << " " << high << "\n";
                if (file_is_readable (format_out (prefix, mid)))
                    low = mid;
                else
                    high = mid;
            }

            snap = high;
            std::cerr << "recover: attempting to load snapshot " << (snap-1u) << "\n";
            in_filename = format_out (prefix, snap-1u);
        }
        else
        {
            std::cerr << "recover: no snapshots found, starting clean slate\n";
        }
    }
    stor->add_data (in_filename);
    stor->dump_statistics (std::cerr);

    if (write_gofr)
        gofr = make_correlator ("density", stor, 0.01, 0.);

    if (!recover)
    {
        stor->periods ().savetxt (prefix + "periods");
        stor->save_data (prefix + "init-coords.dat");
    }

    // init chainrunner
    std::cerr << "inter " << inter << '\n';
    cr->seed_random (random_seed);

    if (probe_test)
    {
        cr->probe_test_pattern (stor);
        return 0;
    }

    if (inter == "lj" && !dont_optimize_sr_lr_split)
    {
        cr->optimize_sr_lr_split (stor);
    }

    if (!skip_calib)
    {
        cr->calibrate (stor);
    }

    // warm-up
    cr->collide (stor, 1.);
    cr->dump_statistics (std::cerr);
    unsigned jmany = 1, gofr_samples = 1;
    for (; snap != snap_target; ++snap)
    {
    redo_snap:
        for (unsigned j = 0; !sigterm_caught () && j != jmany; ++j)
        {
            cr->collide (stor, 1.);
            if (gofr) gofr->sample (stor, gofr_samples);
        }

        stor->dump_statistics (std::cerr);
        cr->dump_statistics (std::cerr);
        cr->save_histograms (prefix);
        if (gofr) gofr->savetxt (format_out (prefix, snap, "gofr"));
        stor->save_data (format_out (prefix, snap));

        if (sigterm_caught ())
        {
            std::cerr << "exit caught SIGTERM, exiting\n";
            return 1;
        }

        if (jmany == 1)
        {
            // now for real
            jmany = snap_disp;
            gofr_samples = 1000;
            goto redo_snap;
        }
    }

    std::cerr << "exit number of snapshots reached\n";
    return 0;
}
