// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// IPL potentials in the screened-lattice flavor
#include "ecmc.hpp"
#include "ipl.hpp"
#include "marsaglia.hpp"

static
const int LATTICE_ORDER = 4;

static
const double SAFETY_CUSHION = 1.5;

template <int DIM>
struct GridIndices
{
    typedef std::array <int, DIM> index_tuple_t;

    // constructors
    GridIndices (int low, int high)
    {
        low_.fill (low);
        high_.fill (high);
    }

    GridIndices (const index_tuple_t &low, const index_tuple_t &high)
        : low_ (low), high_ (high)
    {
    }

    // iteration
    struct iterator;

    iterator begin () const
    {
        iterator ret (*this);
        for (size_t n = 0; n != DIM; ++n)
            ret.pos_[n] = low_[n];
        return ret;
    }

    iterator end () const
    {
        iterator ret (*this);
        ret.pos_[0] = high_[0];
        for (size_t n = 1; n != DIM; ++n)
            ret.pos_[n] = low_[n];
        return ret;
    }

    struct iterator
    {
        const index_tuple_t &operator* () const
        {
            return pos_;
        }

        iterator &operator++ ()
        {
            for (size_t n = DIM-1; n != 0; --n)
            {
                if (++pos_[n] == high_[n])
                    pos_[n] = low_[n];
                else
                    return *this;
            }
            ++pos_[0];
            return *this;
        }

        friend
        std::ostream &operator<< (std::ostream &os, const iterator &me)
        {
            os << me.pos_[0];
            for (size_t n = 1; n != DIM; ++n)
                os << " " << me.pos_[n];
            return os;
        }

        friend
        bool operator!= (const iterator &lhs, const iterator &rhs)
        {
            return lhs.pos_ != rhs.pos_;
        }

    private:
        iterator (const GridIndices &grid)
            : low_ (grid.low_), high_ (grid.high_)
        {
        }

        index_tuple_t pos_;
        const index_tuple_t low_, high_;

        friend struct GridIndices;
    };

private:
    index_tuple_t low_, high_;
};


template <unsigned DIM, typename ENCODING>
struct Jellium4 : JelliumInteraction
{
    typedef vector <DIM> vector_t;
    typedef MtwDiscreteSampler <vector_t> prober_t;
    typedef ENCODING encoding_t;
    Periods periods;
    double cell_volume;
    double error_bound;
    prober_t prober[DIM];

    Jellium4 ()
        : periods (MAX_DIM, NAN)
    {
        set_parameter ("exponent", 1.);
    }

    void set_parameter (string_ref name, double value)
    {
        error_bound = 0.;
        JelliumInteraction::set_parameter (name, value);
    }

    double lr_event_rate (const vector_t &r, unsigned direction, bool clip = true)
    {
        // FIXME 3D missing
        if (DIM != 2)
            std::cerr << "not impl" << ABORT;

        double pr = 0.;

        for (int i = -LATTICE_ORDER; i <= LATTICE_ORDER; ++i)
        {
            double x = r[direction] + i*periods[direction];
            for (int j = -LATTICE_ORDER; j <= LATTICE_ORDER; ++j)
            {
                double y = r[!direction] + j*periods[!direction];
                double rsq = x*x + y*y;
                if (rsq >= sr_lr_split*sr_lr_split)
                    pr += unit_directional_derivative (x, rsq);
            }
        }

        // screening charges
        double spr = 0;
        double xf = r[direction] + (LATTICE_ORDER+.5) * periods[direction];
        double xb = r[direction] - (LATTICE_ORDER+.5) * periods[direction];
        for (int j = -LATTICE_ORDER; j <= LATTICE_ORDER; ++j)
        {
            double y = r[!direction] + j*periods[!direction];
            spr += unit_potential (xf*xf + y*y);
            spr -= unit_potential (xb*xb + y*y);
        }

        pr += spr / periods[direction];

        if (clip && pr < 0.)
            return 0.;
        else
            return stren * pr;
    }

    template <typename STORAGE>
    static
    GridIndices <DIM> get_cell_grid (const STORAGE &stor)
    {
        std::array <int, DIM> low, high;
        for (size_t n = 0; n != DIM; ++n)
        {
            int c = stor.cell_count (n);
            low[n] = 0 - c/2;
            high[n] = c - c/2;
        }
        return GridIndices <DIM> (low, high);
    }

    void setup_probe_distr (const AbstractStorage *stor_)
    {
        // FIXME introduce CellStorage
        auto stor = dynamic_cast <const Storage <ENCODING> *> (stor_);
        if (!stor)
            std::cerr << "jellium4: incompatible Storage passed in" << ABORT;

        vector_t cwid = stor->cell_widths ();
        cell_volume = fproduct (cwid.begin (), cwid.end ());

        GridIndices <DIM> cell_grid = this->get_cell_grid (*stor);

        std::cerr << "jellium4_init cell_volume " << cwid << " = " << cell_volume << "\n";
        std::cerr << "jellium4_init cell_count(DIM) " << stor->cell_count (DIM) << "\n";

        for (unsigned direction = 0; direction != DIM; ++direction)
        {
            prober_t *pr = &prober[direction];
            pr->clear ();

            auto it = cell_grid.begin (), it_end = cell_grid.end ();
            for (; it != it_end; ++it)
            {
                vector_t cell_vector = elementwise_prod (cwid, *it);
                int discret = (norm (cell_vector) < 2 * sr_lr_split) ? 50 : 6;
                GridIndices <DIM> discret_grid = GridIndices <DIM> (-discret, discret+1);

                double max_event_rate_found = 0.;

                auto iit = discret_grid.begin (), iit_end = discret_grid.end ();
                for (; iit != iit_end; ++iit)
                {
                    vector_t dr = elementwise_prod (cwid, *iit) / discret;
                    double er = lr_event_rate (cell_vector + dr, direction);
                    update_max (&max_event_rate_found, er);
                }

                // report maximum to prober
                pr->add (max_event_rate_found * SAFETY_CUSHION, cell_vector);
            }

            pr->finish ();
        }
    }

    double total_probe_rate (unsigned direction)
    {
        // hack: ecmc will multiply this by cell_max_density, effectively
        // producing cell_count(DIM) * total_weights() which is what we want.
        return cell_volume * prober[direction].total_weights ();
    }

    void notify_error_bound (const AbstractStorage *stor,
        double err_bound)
    {
        if (stor->dimension () != DIM)
            std::cerr << "Wrong dimension in Jellium4" << ABORT;
        if (error_bound == err_bound && periods == stor->periods ())
            return;

        error_bound = err_bound;
        periods = stor->periods ();

        setup_probe_distr (stor);
    }

    unsigned probe_rate (const vector_t &r, unsigned direction)
    {
        (void)r;
        (void)direction;
        // not implemented
        std::cerr << "jellium4 probe_rate called" << ABORT;
    }

    template <typename VECTOR>
    double random_probe (VECTOR *ret, unsigned direction, RandomContext *random)
    {
        return prober[direction].random_sample (ret, random);
    }
};

static Register <ChainRunner <Jellium4 <2, Monodisperse2D>, Monodisperse2D>>
    one ("jellium4/mono2d");
static Register <ChainRunner <Jellium4 <2, Tagged <Monodisperse2D>>, Tagged <Monodisperse2D>>>
    one_tagged ("jellium4/tagged_mono2d");
