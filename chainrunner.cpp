// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "ecmc.hpp"

AbstractChainRunner::AbstractChainRunner ()
#ifdef XDISP_HISTO
    : xdisp_histo (.02, -10., 40.)
    , disp_histo (.001, 0., 3.)
    , log_revent_histo (.02, -5., 40.)
    , revent_histo (.02, 0., 8.)
#endif
{
    reset_statistics ();
}

void AbstractChainRunner::set_parameter (string_ref name, double value)
{
    inter_params_[name] = value;
}

void AbstractChainRunner::reset_statistics ()
{
    walltime = 0;
    total_lifts = longrange_lifts = total_chains = 0;
    longrange_predicts = shortrange_predicts = 0;
    total_disp = total_xdisp = 0.;
    report_xdisp_pressure = true;
}

void AbstractChainRunner::dump_report (std::ostream &os) const
{
    os << "walltime " << walltime*1e-6 << " seconds\n";
    os << "num_events total " << total_lifts
        << " lr " << longrange_lifts << "\n";
    os << "perf_stat events_per_hour "
        << fdivide (total_lifts, walltime*1e-6) * 3600. << "\n";
    os << "predict lr " << longrange_predicts
        << " sr " << shortrange_predicts << "\n";
    os << "disp " << total_disp
        << " xdisp " << total_xdisp << "\n";
    os << "mean_chain_tdisp " << (total_disp+total_xdisp) / total_chains << "\n";
    if (probe_paccept.num_samples () != 0u)
        os << "probe_paccept min " << probe_paccept.min ()
            << " max " << probe_paccept.max () << "\n";
    // compressibility factor: beta pressure / rho
    if (report_xdisp_pressure)
        os << "compressib_factor " << ((total_disp+total_xdisp) / total_disp) << "\n";
    // interaction parameters
    for (auto kv : inter_params_)
        os << "inter_param " << kv.first << " " << kv.second << "\n";
    inter_params_.clear ();  // output only once
}

void AbstractChainRunner::collide (AbstractStorage *stor,
    double disp_per_particle)
{
    walltime -= gclock ();
    do_collide (stor, disp_per_particle);
    walltime += gclock ();
}

// construct chain runners

AbstractChainRunner *make_chainrunner (string_ref interaction_type, AbstractStorage *stor)
{
    return Factory <AbstractChainRunner>::make (
        interaction_type + "/" + stor->typecode ());
}
