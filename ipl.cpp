// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "ecmc.hpp"
#include "ipl.hpp"
#include "pareto_prober.hpp"

// an interaction which does not implement the long-range part
// of the power law, i.e. truncates the potential.
struct TruncatedIPL : PowerLawInteraction
{
    void set_parameter (string_ref name, double value)
    {
        if (name == "cutoff")
        {
            assert (value > 0.);
            sr_lr_split = value;
        }
        else if (name == "sr_lr_split")
        {
            rt_error ("Invalid parameter sr_lr_split for truncated IPL");
        }
        else
        {
            PowerLawInteraction::set_parameter (name, value);
        }
    }
};

struct IPL : PowerLawInteraction
{
    ParetoProber prober;

    void notify_error_bound (const AbstractStorage *,
        double err_bound, size_t DIM)
    {
        if (prober.error_bound != err_bound)
        {
            prober.error_bound = err_bound;
            prober.setup (sr_lr_split, exponent, DIM);
        }
    }

    double total_probe_rate (unsigned /* direction */) const
    {
        return stren * prober.total_probe_rate ();
    }

    template <typename VECTOR>
    double probe_rate (const VECTOR &dist, unsigned direction)
    {
        return stren * prober.probe_rate (dist, direction);
    }

    template <typename VECTOR>
    double lr_event_rate (const VECTOR &dist, unsigned direction)
    {
        double rsq = norm_sq (dist);
        if (rsq < sq (sr_lr_split))
            return 0.;
        double rate = unit_directional_derivative (dist[direction], rsq);
        if (rate <= 0.)
            return 0.;
        else
            return stren * rate;
    }

    template <size_t DIM>
    double random_probe (vector <DIM> *ret, unsigned, RandomContext *random)
    {
        return stren * prober.random_probe (ret, random);
    }
};

static Register <ChainRunner <IPL, Monodisperse2D>>
    one ("ipl/mono2d");
static Register <ChainRunner <IPL, Tagged <Monodisperse2D>>>
    one_a ("ipl/tagged_mono2d");
static Register <ChainRunner <TruncatedIPL, Monodisperse2D>>
    two ("truncipl/mono2d");
static Register <ChainRunner <IPL, Monodisperse3D>>
    three ("ipl/mono3d");
