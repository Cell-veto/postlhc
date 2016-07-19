// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// IPL potentials in the screened flavor (no lattice sum)
#include "ecmc.hpp"
#include "ipl.hpp"
#include "adap_prober.hpp"

static
const double SAFETY_CUSHION = 2.;

static
double zmin = 1e99;

static
double signum (double x)
{
    if (x < 0)
        return -1;
    else
        return 1;
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
        error_bound = 0.;
    }

    template <size_t DIM>
    void calib_add (const vector <DIM> &r)
    {
        for (unsigned direction = 0; direction != DIM; ++direction)
        {
            double rate = lr_event_rate (r, direction);
            prober[direction].calib_add (norm (r), SAFETY_CUSHION * rate, error_bound);
        }
    }

    template <unsigned DIM>
    void setup_calib ()
    {
        double r_max = 20 * period.max ();

        for (unsigned i = 0; i != 10; ++i)
        {
            vector <DIM> u_vec = zero_vector <DIM> ();
            u_vec[0] = cos (2*M_PI/10 * i);
            u_vec[1] = sin (2*M_PI/10 * i);

            for (double r = 0; r < r_max; r += 1e-3)
                calib_add (r * u_vec);
            for (double r = 1e-10 * r_max; r < 1e10 * r_max; r *= 1.1)
                calib_add (r * u_vec);
        }
    }

    void notify_error_bound (const AbstractStorage *stor,
        double err_bound, size_t DIM)
    {
        if (error_bound != err_bound || period != stor->periods ())
        {
            error_bound = err_bound;
            period = stor->periods ();

            if (DIM >= exponent + 3)
                std::cerr << "jellium3: exponent " << exponent
                    << " too small for dimension " << DIM << ABORT;

            for (unsigned direction = 0; direction != DIM; ++direction)
                prober[direction].init (DIM, error_bound,
                    .25*period[direction] + error_bound /* tail radius */,
                    exponent+1 /* tail paralpha (= exponent+2 decay) */,
                    3*period[direction] + error_bound /* tail radius */,
                    exponent+2 /* tail paralpha (= exponent+3 decay) */);

            switch (DIM)
            {
            case 2:
                setup_calib <2> ();
                break;
            case 3:
                setup_calib <3> ();
                break;
            default:
                std::cerr << "jellium3: DIM = " << DIM
                    << " not implemented" << ABORT;
            }

            for (unsigned direction = 0; direction != DIM; ++direction)
                prober[direction].calib_finish (DIM);
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
    double compute_shift (const vector <DIM> &r, unsigned direction, double rhs)
    {
        rhs *= .5;
        double rsq = invert_unit_potential (rhs);

        if (rsq <= 0.)
            return 0.;

        double x = r[direction];
        double c = norm_sq (r) - rsq;
        if (x*x - c <= 0)
            return 0.;
        double D = sqrt (x*x - c);
        double q = (x<0) ? (D-x) : (-D-x);  // q has opposite sign of x
        double x1 = q, x2 = c/q;
        if (x1 < 0 && x2 > 0)
            return x1;
        else if (x2 < 0 && x1 < 0)
            return x2;
        std::cerr << "hihihi" << norm_sq (r) << " "<<  rhs << " " << (x*x-c) << " " << x1 << " " << x2 << ABORT;
        return 0.;
    }

    template <size_t DIM>
    double lr_event_rate (const vector <DIM> &r, unsigned direction)
    {
        typedef vector <DIM> vector_t;

        const double rsq = norm_sq (r);
        const double x = r[direction];
        const double L = period[direction];

        double rate;

        /*
        if (rsq == 0.)
            return 0;

        if (x*x < L*L && (rsq-x*x) < L*L)
        {
            // screening charges (discarded in the first periodic cell)
            vector_t rpp  = r + ( 1.5*L) * e;
            double u_rpp  = unit_potential (rpp);
            vector_t imp  = r + (     L) * e;
            double lr_imp = unit_directional_derivative (imp[direction], norm_sq (imp));
            vector_t rp   = r + (  .5*L) * e;
            double u_rp   = unit_potential (rp);

            double lr_im0 = unit_directional_derivative (x, rsq);

            vector_t rm   = r + ( -.5*L) * e;
            double u_rm   = unit_potential (rm);
            vector_t imm  = r + (    -L) * e;
            double lr_imm = unit_directional_derivative (imm[direction], norm_sq (imm));
            vector_t rmm  = r + (-1.5*L) * e;
            double u_rmm  = unit_potential (rmm);

            std::cerr << x << " " << rsq << "\n";
            std::cerr << u_rpp << " " << u_rp << " " << u_rm << " " << u_rmm << "\n";
            std::cerr << lr_imp << " " << lr_im0 << " " << lr_imm << "\n";

            rp[direction] += compute_shift (rp, direction, (u_rpp-u_rm) / L + lr_imp - lr_im0);
            rm[direction] += compute_shift (rm, direction, (u_rp-u_rmm) / L + lr_im0 - lr_imm);

            u_rp   = unit_potential (rp);
            u_rm   = unit_potential (rm);

            rate = (u_rp-u_rm) / L;

            if (rsq > sq (sr_lr_split))
                rate += lr_im0;
        }
        else */
        if (rsq < 1e6 * L*L)
        {
            vector_t rp = r;
            double &xp  = rp[direction];
            xp         += .5*L;
            xp         += signum (xp) * sqrt (fmax (0., .25*L*L - norm_sq (rp)));
            double u_rp = unit_potential (rp);
            vector_t rm = r;
            double &xm  = rm[direction];
            xm         += -.5*L;
            xm         += signum (xm) * sqrt (fmax (0., .25*L*L - norm_sq (rm)));
            double u_rm = unit_potential (rm);

            rate  = (u_rp-u_rm) / L;

            if (sqrt (fmin (norm_sq (rm), norm_sq (rp))) < zmin)
            {
                zmin = sqrt (fmin (norm_sq (rm), norm_sq (rp)));
                std::cerr << "ZMIN " << zmin << "\n";
            }

            // main charge -- parts smaller than sr_lr_split are handled by the SR code
            if (rsq > sq (sr_lr_split))
                rate += unit_directional_derivative (x, rsq);
        }
        else
        {
            // avoid catastrophic cancellation in the far field
            // (valid also for exponent = 0)
            rate = x * pow (rsq, -.5*6. -.5*exponent) * L*L
                * (2+exponent) * (1./24) * (3.*rsq - (4+exponent)*x*x);
        }

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
static Register <ChainRunner <Jellium3, Monodisperse3D>>
    two ("jellium3/mono3d");
