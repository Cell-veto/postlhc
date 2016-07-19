// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef STORAGE_HPP_INCLUDED
#define STORAGE_HPP_INCLUDED

#include "tools.hpp"

const unsigned MAX_DIM = 3;

struct MostGeneralParticle
{
    vector <MAX_DIM> coords;
    double radius;
    double tag;
};

struct Periods : vector <MAX_DIM>
{
    Periods ()
    {
        fill (0.);
    }

    Periods (unsigned MD, double value)
    {
        (void)MD;
        assert (MD == MAX_DIM);
        fill (value);
    }

    static Periods from_file (string_ref filename);
    void savetxt (string_ref filename);

    double min (unsigned MAX) const
    {
        assert (MAX <= MAX_DIM);
        return *std::min_element (begin (), begin () + MAX);
    }

    double max () const
    {
        return *std::max_element (begin (), end ());
    }

    double volume (unsigned MAX) const
    {
        assert (MAX <= MAX_DIM);
        return fproduct (begin (), begin () + MAX);
    }
};

struct AbstractParticleGenerator
{
    virtual ~AbstractParticleGenerator () {}
    virtual bool get (MostGeneralParticle *) = 0;
};

struct AbstractParticleSink
{
    virtual ~AbstractParticleSink () {}
    virtual void put (const MostGeneralParticle &) = 0;
};

struct AbstractStorage : AbstractParticleSink, FactoryProduct <AbstractStorage>
{
    virtual ~AbstractStorage () {}
    virtual void dump_report (std::ostream &) const = 0;

    virtual unsigned dimension () const = 0;
    virtual Periods periods () const = 0;

    virtual size_t num_particles () const = 0;
    virtual AbstractParticleGenerator *all_particles () const = 0;
            void save_data (string_ref filename);

    virtual void set_periods (const Periods &) = 0;
            void load_periods (string_ref filename);
    virtual void put (const MostGeneralParticle &) = 0;
            void add_data (string_ref filename);
};

// in-memory encodings
struct Monodisperse2D
{
    static const unsigned DIM = 2;
    struct Payload {};

    static
    void decode_props (MostGeneralParticle *most, const Payload &)
    {
        most->radius = 1.;
        most->tag = 0.;
    }

    static
    void encode_props (Payload *, const MostGeneralParticle &)
    {
    }
};

struct Monodisperse3D : Monodisperse2D
{
    static const unsigned DIM = 3;
};

// exception: no more free cells, need to subdivide
struct StorageFull {};

template <typename ENCODING>
struct Storage : ENCODING, AbstractStorage
{
public:
    typedef uint64_t key_t;
    typedef ENCODING encoding_t;
    using typename ENCODING::Payload;
    static const unsigned DIM = ENCODING::DIM;
    typedef vector <DIM> vector_t;

    struct CellData : Payload
    {
        uint64_t coord[DIM];
        bool in_use () const { return coord[1] & 1ull; }
        void clear () { coord[1] = 0; }
    };

    void reduce_cell_width (double target)
    {
        while (cell_width (0) > target || cell_width (1) > target)
            subdivide ();
    }

private:
    // MEMORY MANAGEMENT -- messy internals
    // rely on Linux' late allocation strategy
    // we don't touch all this memory until we need it
    static const size_t MAX_KEY = 100000000;

    size_t num_cells_, num_;

    unsigned nextsub_;
    unsigned bc_[DIM+1];
    Periods peri_;
    double real2frac[DIM], frac2real[DIM];

public:
    size_t cell_count (unsigned n) const
    {
        assert (n <= DIM);
        return 1ull << bc_[n];
    }

private:
    double cell_width (unsigned n) const
    {
        assert (n < DIM);
        return period (n) / (1ull<<bc_[n]);
    }

    uint64_t cell_stride (unsigned n) const
    {
        assert (n < DIM);
        return 1ull << (64 - bc_[n]);
    }

    uint64_t cell_floor (uint64_t frac, unsigned n) const
    {
        assert (n < DIM);
        frac >>= 64 - bc_[n];
        frac <<= 64 - bc_[n];
        return frac;
    }

    uint64_t cell_ceil (uint64_t frac, unsigned n) const
    {
        assert (n < DIM);
        frac >>= 64 - bc_[n];
        ++frac;
        frac <<= 64 - bc_[n];
        return frac;
    }

    // collect values of an accessor function and return as a vector datatype
    template <typename VECTOR, typename PTR_TO_MEMBERFUNC>
    VECTOR to_vector (PTR_TO_MEMBERFUNC func) const
    {
        VECTOR ret;
        for (unsigned n = 0; n != ret.size (); ++n)
            ret[n] = (this->*func) (n);
        return ret;
    }

    // compute base key where a particle with those fractional
    // coordinates should go (feed to find_spot)
    key_t make_base_key (const uint64_t frac[DIM]) const
    {
        uint64_t ret = 0;
        for (unsigned n = 0; n != DIM; ++n)
        {
            assert (bc_[n] > 0);
            ret <<= bc_[n];
            ret |= frac[n] >> (64-bc_[n]);
        }
        return ret << bc_[DIM];
    }

    // find an unused cell in the array
    key_t find_spot (key_t base_key)
    {
        key_t k = base_key;
        key_t k_end = base_key + cell_count (DIM);
        for (; k != k_end; ++k)
            if (!data_[k].in_use ())
                return k;
        throw StorageFull ();
    }

    // increase capacity of cell array.
    // INVALIDATES existing keys.
    void subdivide (bool silent = false)
    {
        // back up data
        std::vector <CellData> backup;
        backup.reserve (num_);
        for (size_t k = 0; k != num_cells_; ++k)
            if (data_[k].in_use ())
                backup.push_back (data_[k]);
        // at most two particles might be in hiatus in displace()
        assert (num_ <= backup.size ()+2u);
        assert (num_ >= backup.size ());

        // subdivide one spatial direction
        ++bc_[nextsub_++];
        nextsub_ %= DIM;

        // compute new number of cells
        num_cells_ = 1ull << bc_[DIM];
        for (unsigned n = 0; n != DIM; ++n)
            num_cells_ *= 1ull << bc_[n];
        if (num_cells_ >= MAX_KEY)
            rt_error ("excessive subdivision - particles seem to be clustering");

        // init new cells
        std::memset (data_, 0, num_cells_ * sizeof (CellData));
        // restore particles
        for (size_t n = 0; n != backup.size (); ++n)
            put_encoded (backup[n]); // cannot throw

        // emit updated parameters to logfile
        if (silent) return;
        std::cerr << "cell_widths " << cell_widths () << "\n";
        std::cerr << "cell_max_density " << cell_max_density () << "\n";
    }

    void encode_coords (uint64_t frac[DIM], const vector <MAX_DIM> &real) const
    {
        unsigned n = 0;
        for (; n != DIM; ++n)
        {
            assert (real[n] >= 0.);
            assert (real[n] < peri_[n]);
            frac[n] = real2frac[n] * real[n];
        }
        for (; n != MAX_DIM; ++n)
            assert (real[n] == 0.);
        // in_use flag
        frac[1] |= 1;
    }

    void decode_coords (vector <MAX_DIM> *real, const uint64_t frac[DIM]) const
    {
        unsigned n = 0;
        for (; n != DIM; ++n)
            (*real)[n] = frac2real[n] * frac[n];
        for (; n != MAX_DIM; ++n)
            (*real)[n] = 0.;
    }

    // put CellData in a place where it fits; returns the key.
    // if array is full, this throws StorageFull.
    key_t put_encoded (const CellData &cd)
    {
        key_t k = make_base_key (cd.coord);
        k = find_spot (k);
        data_[k] = cd;
        return k;
    }

    // put CellData in a place where it fits; returns the key.
    // if array is full, reallocates.
    key_t put_encoded_possibly_subdivide (const CellData &cd)
    {
        try
        {
            return put_encoded (cd);
        }
        catch (StorageFull)
        {
            subdivide ();
            return put_encoded_possibly_subdivide (cd);
        }
    }

    // remove a CellData from the array, and clear the cell.
    CellData pull_encoded (key_t k)
    {
        assert (k < num_cells_);
        assert (data_[k].in_use ());
        CellData cd = data_[k];
        data_[k].clear ();
        return cd;
    }

public:
    typedef Storage <ENCODING> this_t;

    void *operator new (size_t /* ignored */)
    {
        size_t len = MAX_KEY * sizeof (CellData) + sizeof (Storage);
        return ::new char[len];
    }

    void operator delete (void *mem)
    {
        return ::delete[] static_cast <char *> (mem);
    }

    Storage ()
    {
        num_cells_ = num_ = 0;
        nextsub_ = 0;
        for (unsigned n = 0; n != DIM; ++n)
            bc_[n] = 4;
        bc_[0] = 3;
        bc_[DIM] = 3;
        subdivide (true);
    }

    virtual
    unsigned dimension () const
    {
        return DIM;
    }

    size_t num_cells () const
    {
        return num_cells_;
    }

    virtual
    size_t num_particles () const
    {
        return num_;
    }

    double particle_density () const
    {
        return num_particles () / peri_.volume (DIM);
    }

    double period (unsigned compo) const
    {
        assert (compo < MAX_DIM);
        return peri_[compo];
    }

    virtual
    Periods periods () const
    {
        return peri_;
    }

    virtual
    void set_periods (const Periods &periods)
    {
        if (num_particles () != 0)
            rt_error ("don't set_periods with particles");

        double MAX = 1. + UINT64_MAX;
        unsigned n = 0;
        for (; n != DIM; ++n)
        {
            peri_[n] = periods[n];
            real2frac[n] = periods[n]==0. ? 0. : MAX/periods[n];
            frac2real[n] = periods[n]==0. ? 0. : 1./real2frac[n];
        }
        for (; n != MAX_DIM; ++n)
        {
            if (periods[n] != 0.)
                rt_error ("periods invalid for this Storage");
        }

        // emit updated parameters to logfile
        std::cerr << "periods " << peri_ << '\n';
    }

    double cell_diagonal () const
    {
        double nsq = 0.;
        for (unsigned n = 0; n != DIM; ++n)
            nsq += sq (cell_width (n));
        return sqrt (nsq);
    }

    vector_t cell_widths () const
    {
        return to_vector <vector_t> (&this_t::cell_width);
    }

    double cell_max_density () const
    {
        double md = 1.;
        for (unsigned n = 0; n != DIM; ++n)
            md /= cell_width (n);
        return md * cell_count (DIM);
    }

    virtual
    void dump_report (std::ostream &os) const
    {
        os << "fillfrac " << fdivide (num_, num_cells_) << "\n";

        // cell occupation histogram
        size_t max_occup = cell_count (DIM);
        std::vector <size_t> occ_hist (max_occup + 1, 0);
        for (size_t k = 0; k != num_cells_; k += max_occup)
        {
            size_t c = 0;
            for (size_t j = 0; j != max_occup; ++j)
                c += data_[k+j].in_use ();
            ++occ_hist[c];
        }
        os << "cell_occup";
        for (size_t occ : occ_hist)
            os << ' ' << fdivide (occ, num_cells_/max_occup);
        os << '\n';
    }

    virtual
    void put (const MostGeneralParticle &part)
    {
        (void)insert (part);
    }

    key_t insert (const MostGeneralParticle &part)
    {
        CellData cd;
        this->encode_coords (cd.coord, part.coords);
        this->encode_props (static_cast <Payload *> (&cd), part);
        key_t ret = put_encoded_possibly_subdivide (cd);
        ++num_;
        return ret;
    }

    void get (MostGeneralParticle *ret, key_t k) const
    {
        assert (k < num_cells_);
        assert (data_[k].in_use ());
        this->decode_coords (&ret->coords, data_[k].coord);
        this->decode_props (ret,
            static_cast <const Payload &> (data_[k]));
    }

    void extract (MostGeneralParticle *ret, key_t k)
    {
        get (ret, k);
        remove (k);
    }

    void remove (key_t k)
    {
        data_[k].clear ();
        --num_;
    }

    key_t random_particle (RandomContext *random)
    {
        assert (num_ > 0);
        key_t ret;
        do {
            ret = random->uint (num_cells_);
        } while (!data_[ret].in_use ());
        return ret;
    }

    // choose an arbitrary cell, return true if there's a particle
    bool random_cell (key_t *ret, RandomContext *random) const
    {
        *ret = random->uint (num_cells_);
        return data_[*ret].in_use ();
    }

    key_t displace (key_t whom, unsigned direction, double distance,
        key_t preserve)
    {
        assert (distance >= 0.);
        CellData moving = pull_encoded (whom);

        moving.coord[direction] += uint64_t (real2frac[direction] * distance);
        
        // maintain in_use flag
        moving.coord[1] |= 1;

        if (preserve == whom)
            return put_encoded_possibly_subdivide (moving);

        try
        {
            put_encoded (moving);
            return preserve;
        }
        catch (StorageFull)
        {
            // damage control - we need to resize the cell array,
            // but our caller also needs a valid reference to
            // inext.  pull it out.
            CellData next = pull_encoded (preserve);
            put_encoded_possibly_subdivide (moving);
            return put_encoded_possibly_subdivide (next);
        }
    }

    bool probe (key_t *k, key_t ref, vector_t *r, RandomContext *rng)
    {
        const double MAX = 1. + UINT64_MAX;
        uint64_t kv[DIM];
        assert (data_[ref].in_use ());
        // find fractional coordinates of probe location
        for (unsigned n = 0; n != DIM; ++n)
        {
            (*r)[n] *= real2frac[n];
            double abs = (*r)[n] + data_[ref].coord[n];
            double rem = std::remainder (abs, MAX);
            if (rem < 0.) // FIXME use std::remquo here
                rem += MAX;
            assert (rem >= 0. && rem < MAX);
            kv[n] = rem;
            (*r)[n] -= kv[n] << bc_[n] >> bc_[n];
        }
        *k = make_base_key (kv);
        *k += rng->uint (1ull << bc_[DIM]);
        // correct r for actual position of particle
        // (if there is no particle at *k, the result is bollocks,
        // but it won't be used.)
        for (unsigned n = 0; n != DIM; ++n)
        {
            (*r)[n] += data_[*k].coord[n] << bc_[n] >> bc_[n];
            (*r)[n] *= frac2real[n];
        }
        return data_[*k].in_use ();
    }

    vector_t distance_vector (key_t j, key_t i)
    {
        vector_t r;
        for (unsigned n = 0; n != DIM; ++n)
            r[n] = frac2real[n] *
                int64_t (data_[j].coord[n] - data_[i].coord[n]);
        return r;
    }

    // generators

    struct AllNewGenerator;
    struct StripGenerator;

    struct GeneratorBase : AbstractParticleGenerator
    {
        GeneratorBase (const this_t &stor)
        {
            stor_ = &stor;
        }

        virtual bool not_done () = 0;

        void next ()
        {
            do
                next_cell ();
            while (not_done () && !cell ()->in_use ());
        }

        key_t key () const
        {
            return k_;
        }

        const CellData *cell ()
        {
            assert (k_ < stor_->num_cells_);
            return stor_->data_ + k_;
        }

        // AbstractParticleGenerator interface
        virtual bool get (MostGeneralParticle *dst)
        {
            if (not_done ())
            {
                stor_->get (dst, k_);
                next ();
                return true;
            }

            return false;
        }

    protected:
        virtual void next_cell () = 0;

        void init (key_t k)
        {
            k_ = k;
            while (not_done () && !cell ()->in_use ())
                next_cell ();
        }

        key_t k_;
        const this_t *stor_;
    };

    virtual
    AllNewGenerator *all_particles () const
    {
        return new AllNewGenerator (*this);
    }

    AllNewGenerator enumerate_all ()
    {
        return AllNewGenerator (*this);
    }

    StripGenerator enumerate_strip (key_t whence, unsigned direction,
        double strip_width, double begin, double end)
    {
        if (2*cell_diagonal () + 2*strip_width >= period (direction))
            rt_error ("sr_lr_split too large / periods too small");
        if (fabs (end) > strip_max_extent (direction))
            rt_error ("strip too long (this should not happen)");

        return StripGenerator (*this, data_[whence].coord, direction,
            strip_width, begin, end);
    }

    StripGenerator enumerate_box (key_t whence, double a)
    {
        if (! (a > 0))
            rt_error ("box width negative");
        for (unsigned n = 0; n != DIM; ++n)
            if (a > strip_max_extent (n))
                rt_error ("enumeration box too large / periods too small");
        return StripGenerator (*this, data_[whence].coord, 0, a, -a, a);
    }

    double strip_max_extent (unsigned direction) const
    {
        double ml = .5 * period (direction) - 3 * cell_width (direction);
        assert (ml > 0.);
        return ml;
    }

    struct StripGenerator : GeneratorBase
    {
        using GeneratorBase::k_;
        using GeneratorBase::stor_;

        // this assumes a sample of sufficient size.
        // self-overlaps of the strip are not handled gracefully.
        // for EC dynamics, this generator traverses particles in
        // direction-of-motion order (roughly), and can be aborted early.
        StripGenerator (const this_t &stor, uint64_t base[DIM],
            unsigned direction,
            double strip_width, double begin, double end)
            : GeneratorBase (stor)
        {
            direction_ = direction;
            for (unsigned n = 0; n != DIM; ++n)
            {
                uint64_t b, e;
                if (n == direction)
                {
                    b = int64_t (stor_->real2frac[n] * begin);
                    e = stor_->real2frac[n] * end;
                }
                else
                {
                    uint64_t rep_frac = stor_->real2frac[n] * strip_width;
                    b = 0ull - rep_frac;
                    e = rep_frac;
                }

                lo_bound_[n] = stor_->cell_floor (base[n] + b, n);
                hi_bound_[n] = stor_->cell_ceil (base[n] + e, n);
                assert (lo_bound_[n] != hi_bound_[n]);
                position_[n] = lo_bound_[n];
            }

            position_[DIM] = 0;

            GeneratorBase::init (stor_->make_base_key (position_));
        }

        // clip the enumeration longitudinally
        void clip (key_t active, double end)
        {
            uint64_t whence = stor_->data_[active].coord[direction_];
            uint64_t end_frac = stor_->real2frac[direction_] * end;
            uint64_t new_hi_bound = stor_->cell_ceil (whence + end_frac, direction_);
#ifdef DEBUG
            // check we're not jumping over the end marker.
            for (uint64_t t = new_hi_bound; t != hi_bound_[direction_];
                    t += stor_->cell_stride (direction_))
                assert (t != position_[direction_]);
            assert (lo_bound_[direction_] != new_hi_bound);
#endif
            hi_bound_[direction_] = new_hi_bound;
        }

        virtual
        bool not_done ()
        {
            return position_[direction_] != hi_bound_[direction_];
        }

    private:
        virtual
        void next_cell ()
        {
            ++k_;
            if (++position_[DIM] < stor_->cell_count (DIM))
                return;
            position_[DIM] = 0;

            unsigned n;
            for (n = DIM-1; n+1u != 0u; --n)
            {
                if (n == direction_)
                    continue;
                position_[n] += stor_->cell_stride (n);
                if (position_[n] != hi_bound_[n])
                    goto recompute;
                position_[n] = lo_bound_[n];
            }

            n = direction_;
            position_[n] += stor_->cell_stride (n);

        recompute:
            k_ = stor_->make_base_key (position_);
        }

        uint64_t position_[DIM+1];
        uint64_t lo_bound_[DIM], hi_bound_[DIM];
        unsigned direction_;
    };

    struct AllNewGenerator : GeneratorBase
    {
        using GeneratorBase::k_;
        using GeneratorBase::stor_;

        AllNewGenerator (const this_t &stor) : GeneratorBase (stor)
        {
            GeneratorBase::init (0);
        }

        virtual bool not_done ()
        {
            return k_ != stor_->num_cells ();
        }

    private:
        virtual void next_cell ()
        {
            ++k_;
        }
    };

private:
    CellData data_[1]; // must be last data member
};

struct AbstractChainRunner : FactoryProduct <AbstractChainRunner>
{
    AbstractChainRunner ();

    // cumulative statistics
    uint64_t walltime;
    uint64_t total_lifts, longrange_lifts, total_chains;
    uint64_t longrange_predicts, shortrange_predicts;
    double total_disp, total_xdisp;
#ifdef XDISP_HISTO
    Histogram xdisp_histo;
    Histogram disp_histo;
    Histogram log_revent_histo;
    Histogram revent_histo;
#endif
    Extrema <double> probe_paccept;
    bool report_xdisp_pressure;

    // virtual interface
    virtual ~AbstractChainRunner () {}
    virtual void seed_random (unsigned) = 0;
    virtual void set_parameter (string_ref name, double value);
    virtual void reset_statistics ();
    virtual void dump_report (std::ostream &) const;
    virtual void calibrate (AbstractStorage *) = 0;
            void collide (AbstractStorage *, double disp_per_particle);
    virtual void probe_test_pattern (AbstractStorage *, unsigned direction = 0) = 0;
protected:
    mutable std::map <string, double> inter_params_;
    virtual void do_collide (AbstractStorage *, double disp_per_particle) = 0;
};

struct AbstractCorrelator : FactoryProduct <AbstractCorrelator>
{
    virtual ~AbstractCorrelator ();
    virtual void init (AbstractStorage *stor, double bin_width, double rmax) = 0;
    virtual void sample (AbstractStorage *stor, size_t multiplier) = 0;
    virtual void savetxt (std::ostream &os) = 0;
            void savetxt (string_ref filename);
};

AbstractStorage *make_storage (string_ref typecode);
AbstractChainRunner *make_chainrunner (string_ref typecode, AbstractStorage *stor);
AbstractCorrelator *make_correlator (string_ref attribute, AbstractStorage *stor,
    double bin_width, double rmax);
void add_data (AbstractParticleSink *, string_ref filename);
void save_data (string_ref filename, AbstractParticleGenerator *);

#endif /* STORAGE_HPP_INCLUDED */
