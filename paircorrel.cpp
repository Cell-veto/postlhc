// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
// compute pair correlators.
#include "storage.hpp"
#include <fstream>
#include <iomanip>

// FIXME implement attributes

template <typename ENCODING>
struct Correlator : AbstractCorrelator
{
    typedef Storage <ENCODING> stor_t;
    typedef typename stor_t::key_t key_t;
    static const unsigned DIM = ENCODING::DIM;

    static stor_t *downcast (AbstractStorage *stor_)
    {
        // convert to reference to force bad_cast exceptions
        return &dynamic_cast <stor_t &> (*stor_);
    }

    virtual void init (AbstractStorage *stor_, double bin_width, double rmax)
    {
        stor_t *stor = downcast (stor_);

        binwid_ = bin_width;

        if (rmax == 0.)
        {
            method_ = "full";
            rmax_ = stor->periods ().min (DIM) / 2;
        }
        else
        {
            method_ = "window";
            assert (rmax > 0.);
            rmax_ = rmax;
        }

        hist_.resize (bin (rmax_) - 1., 0);
        ndata_ = 0;
        walltime_ = 0;
    }

    virtual void sample (AbstractStorage *stor_, size_t multiplier)
    {
        stor_t *stor = downcast (stor_);

        // initialized late since stor might be empty at ctor time
        rho_ = stor->particle_density ();
        ndata_ += multiplier;
        walltime_ -= gclock ();

        for (size_t n1 = 0; n1 != multiplier; ++n1)
        {
            key_t k1 = stor->random_particle (&random);

            if (method_ == "full")
                do_sample (stor, k1, stor->enumerate_all ());
            else // window
                do_sample (stor, k1, stor->enumerate_box (k1, rmax_));
        }

        walltime_ += gclock ();
    }

    virtual void savetxt (std::ostream &os)
    {
        const double exp_sample_dens = ndata_ * rho_;
        hist_[0] = 0;

        double rv1 = 0., rv2;
        for (unsigned i = 0; i != hist_.size (); ++i)
        {
            double r = (i+.5) * binwid_;
            rv2 = ball_volume (DIM, (i+1) * binwid_);
            double expected = exp_sample_dens * (rv2-rv1);
            double gofr = hist_[i] / expected;
            os << std::setprecision (3) << std::fixed
                << r << ' '
                << std::setprecision (12) << std::scientific
                << gofr << '\n';
            rv1 = rv2;
        }

        os << "# gofr_stat " << ndata_ << " fin\n";
        std::cerr << "gofr_stat " << ndata_ << " wrote histo\n";
        std::cerr << "gofr_walltime " << walltime_*1e-6 << " seconds\n";
    }

private:
    template <typename GENERATOR>
    void do_sample (stor_t *stor, key_t k1, GENERATOR g)
    {
        for (; g.not_done (); g.next ())
        {
            key_t k2 = g.key ();
            double r = norm (stor->distance_vector (k1, k2));
            unsigned i = bin (r);
            if (i < hist_.size ())
                ++hist_[i];
        }
    }

    size_t bin (double r) const
    {
        return std::floor (r/binwid_);
    }

    RandomContext random; 
    double binwid_, rho_;
    string method_;
    double rmax_;
    std::vector <size_t> hist_;
    size_t ndata_;
    uint64_t walltime_;
};

AbstractCorrelator::~AbstractCorrelator ()
{
}

void AbstractCorrelator::savetxt (string_ref filename)
{
    std::ofstream ofs (filename);
    if (!ofs)
        rt_error ("cannot open file: " + filename);
    savetxt (ofs);
    if (!ofs)
        rt_error ("error writing data: " + filename);
    ofs.close ();
    if (!ofs)
        rt_error ("error writing data (2): " + filename);
}

AbstractCorrelator *make_correlator (string_ref attribute, AbstractStorage *stor,
    double bin_width, double rmax)
{
    AbstractCorrelator *ret = Factory <AbstractCorrelator>::make (
        attribute + "/" + stor->typecode ());
    ret->init (stor, bin_width, rmax);
    return ret;
}

static Register <Correlator <Monodisperse2D>> one ("density/mono2d");
static Register <Correlator <Monodisperse3D>> two ("density/mono3d");
