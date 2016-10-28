// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// IPL potentials in the screened flavor (no lattice sum, but polar pairs)
// FIXME documentation
#include "ecmc.hpp"
#include "ipl.hpp"
#include "adap_prober.hpp"

static
const double SAFETY_CUSHION = 2.;

template <size_t DIM, typename CALLABLE>
static
void for_all_angles (unsigned discretization, const CALLABLE &yield)
{
    // FIXME explicit 3D missing (strictly speaking, we might underestimate
    // the LR if the box is very different from cubic.  ecmc.hpp would
    // violently complain, however, and this would not go undetected)
    double incr = 2*M_PI / discretization;
    vector <DIM> d = zero_vector <DIM> ();
    for (unsigned i = 0; i != discretization; ++i)
    {
        d[0] = cos (incr*i);
        d[1] = sin (incr*i);
        yield (d);
    }
}

struct Jellium3 : JelliumInteraction
{
    Periods period;
    double error_bound;
    AdaptiveProber prober[MAX_DIM];

    Jellium3 ()
        : period (MAX_DIM, NAN)
    {
        set_parameter ("exponent", 1.);
    }

    void set_parameter (string_ref name, double value)
    {
        error_bound = 0.;
        JelliumInteraction::set_parameter (name, value);
    }

    void add_sample (unsigned direction, double r, double rate)
    {
        rate *= SAFETY_CUSHION;
        prober[direction].calib_add (r, rate, error_bound);
    }

    template <unsigned DIM>
    void calibrate (unsigned direction)
    {
        typedef vector <DIM> vector_t;
        const double L = period[direction];
        const double r_decay = sr_lr_split; // first decaying inner shell
        const double r_tail  = 2*L;  // last (decaying) inner shell, begin tail
        std::cerr << "jellium3_init shell_radii " << r_decay << " " << r_tail << "\n";
        std::cerr << "jellium3_init error_bound " << error_bound << "\n";

        prober[direction].calib_begin (DIM,
            r_decay, exponent+0.5,               // inner paralpha
            r_tail,  exponent+3);                // tail paralpha (= exponent+4 decay)

        // scan along the sr_lr_split
        for_all_angles <DIM> (100, [&] (const vector_t &unit_vec)
        {
            const double r = sr_lr_split;
            add_sample (direction, r, lr_event_rate (r*unit_vec, direction));
        });

        // in extremely small systems, the first few copies can cause trouble
        for (int i = -10; i <= 10; ++i)
        {
            vector_t ridge = zero_vector <DIM> ();
            ridge[direction] = i*L;
            for (; ridge[direction] < 3*sr_lr_split; ridge[direction] += 1e-3*sr_lr_split)
                add_sample (direction, norm (ridge), lr_event_rate (ridge, direction));
        }

        for_all_angles <DIM> (30, [&] (const vector_t &unit_vec)
        {
            for (double r = 0; r < 2*L; r += 1e-3*sr_lr_split)
                add_sample (direction, r, lr_event_rate (r*unit_vec, direction));

            for (double r = 1.; r < 1e15*L; r *= 1.0001)
                add_sample (direction, r, lr_event_rate (r*unit_vec, direction));
        });

        prober[direction].calib_finish (DIM);
    }

    void notify_error_bound (const AbstractStorage *stor,
        double err_bound)
    {
        const unsigned DIM = stor->dimension ();

        if (error_bound != err_bound || period != stor->periods ())
        {
            error_bound = err_bound;
            period = stor->periods ();

            if (DIM >= exponent + 3)
                std::cerr << "jellium3: exponent " << exponent
                    << " too small for dimension " << DIM << ABORT;

            switch (DIM)
            {
            case 2:
                for (unsigned n = 0; n != DIM; ++n)
                    calibrate <2> (n);
                return;
            case 3:
                for (unsigned n = 0; n != DIM; ++n)
                    calibrate <3> (n);
                return;
            default:
                std::cerr << "jellium3: DIM = " << DIM
                    << " not implemented" << ABORT;
            }
        }
    }

    double total_probe_rate (unsigned direction) const
    {
        return prober[direction].total_probe_rate ();
    }

    template <typename VECTOR>
    double probe_rate (const VECTOR &r, unsigned direction)
    {
        return prober[direction].probe_rate (r);
    }

    template <size_t DIM>
    double lr_event_rate (const vector <DIM> &r_, unsigned direction)
    {
        const double L = period[direction];
        const vector <DIM> shift = unit_vector <DIM> (direction) * (L/2);

        // assign polar pair partner
        const double x_prim = r_[direction];
        const double x_pair = x_prim - 2*L * floor (x_prim/L) - L;

        bool first_copy = fabs (x_prim - x_pair) < 2*L;
        vector <DIM> r[2] = { r_, r_ };
        if (x_prim >= 0.)
            r[0][direction] = x_pair;
        else
            r[1][direction] = x_pair;
        double rsq[2] = { norm_sq (r[0]), norm_sq (r[1]) };

        assert (r[0][direction] < 0.);
        assert (r[1][direction] >= 0.);

        // compute the rate for the polar pair
        double rate = 0.;

        // for large distances, the naive formula for the event rate can be
        // numerically problematic.
        if (fmax (rsq[0], rsq[1]) < 1e6*L*L)
        {
            for (int i = 0; i != 2; ++i)
            {
                double x = r[i][direction];

                double u_rp = (i==0 && first_copy)
                    ? 0.
                    : unit_potential (r[i] + shift);
                double u_rm = (i==1 && first_copy)
                    ? 0.
                    : unit_potential (r[i] - shift);

                rate += (u_rp-u_rm) / L;

                // main charge -- parts smaller than sr_lr_split are handled by the SR code
                if (rsq[i] > sq (sr_lr_split))
                    rate += unit_directional_derivative (x, rsq[i]);
            }
        }
        else if (fmax (rsq[0], rsq[1]) < 1e20*L*L)
        {
            for (int i = 0; i != 2; ++i)
            {
                double x = r[i][direction];
                // avoid catastrophic cancellation in the far field
                // (valid also for exponent = 0)
                rate += x * pow (rsq[i], -.5*6. -.5*exponent) * L*L
                    * (2+exponent) * (1./24) * (3.*rsq[i] - (4+exponent)*x*x);
            }
        }
        else
        {
            // rate < numerical precision
        }

        // asymptotic decay of particle event rate is ~ 1 / r^(n+4).
        double w = (r[0][direction] == r_[direction])
            ? norm_sq (r[1]) / norm_sq (r[0])
            : norm_sq (r[0]) / norm_sq (r[1]);

        rate /= 1 + pow (w, -.5 * (exponent+4.));

        if (rate <= 0.)
            return 0.;
        else
            return stren * rate;
    }

    template <typename VECTOR>
    double random_probe (VECTOR *ret, unsigned direction, RandomContext *random)
    {
        return prober[direction].random_probe (ret, random);
    }
};

static Register <ChainRunner <Jellium3, Monodisperse2D>>
    one ("jellium3/mono2d");
static Register <ChainRunner <Jellium3, Tagged <Monodisperse2D>>>
    one_tagged ("jellium3/tagged_mono2d");
static Register <ChainRunner <Jellium3, Monodisperse3D>>
    two ("jellium3/mono3d");
