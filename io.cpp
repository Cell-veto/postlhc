// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "storage.hpp"
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::ends_with;

Periods Periods::from_file (string_ref filename)
{
    Periods periods (MAX_DIM, 0.);
    std::ifstream ifs (filename);
    if (!ifs)
        rt_error ("error opening file: " + filename);
    ifs >> periods[0] >> periods[1];
    if (!ifs)
        goto error;
    ifs >> periods[2] >> std::ws;
    if (ifs.eof ())
        return periods;
error:
    rt_error ("error reading file: " + filename);
}

void Periods::savetxt (string_ref filename)
{
    std::ofstream ofs (filename);
    if (!ofs)
        rt_error ("error opening file: " + filename);
    ofs.precision (15);

    for (double l : *this)
    {
        ofs << l << '\n';
    }

    if (!ofs)
        rt_error ("error writing file: " + filename);
}

static
bool more_data_on_line (std::istream &is)
{
    while (is)
    {
        switch (is.peek ())
        {
        case EOF:
            // ignore EOF caused by our peeking
            is.clear (std::ios::eofbit);
        case '\n':
            return false;
        case ' ':
        case '\t':
            // skip over whitespace
            is.get ();
            assert (!!is);
            continue;
        default:
            return true;
        }
    }

    // stream was in error state even before our checking
    // leave the error bits set so caller can detect this
    return false;
}

template <size_t DIM> static
std::istream &parse_line (std::istream &is, vector <DIM> &d)
{
    unsigned n = 0;
    while (is)
    {
        if (n == DIM)
            rt_error ("line too long in parse_line");
        if (!more_data_on_line (is >> d[n++]))
            return is;
    }

    return is;
}

void add_data (AbstractParticleSink *sink, string_ref filename)
{
    Extrema <double> particle_coords[MAX_DIM];

    if (ends_with (filename, ".dat"))
    {
        MostGeneralParticle most;
        std::memset (&most, 0, sizeof most);
        most.radius = 1.;

        std::ifstream ifs (filename);
        unsigned N = 0;
        while (parse_line (ifs, most.coords))
        {
            for (unsigned n = 0; n != MAX_DIM; ++n)
                particle_coords[n].add (most.coords[n]);
            sink->put (most);
            ++N;
        }

        if (ifs.eof ())
        {
            std::cerr << "npart " << N << " particles read successfully\n";
            std::cerr << "coord_min";
            for (unsigned n = 0; n != MAX_DIM; ++n)
                std::cerr << ' ' << particle_coords[n].min ();
            std::cerr << "\ncoord_max";
            for (unsigned n = 0; n != MAX_DIM; ++n)
                std::cerr << ' ' << particle_coords[n].max ();
            std::cerr << "\n";
            return;
        }

        if (!ifs)
        {
            rt_error ("error reading input file " + filename);
        }
    }
    else
    {
        rt_error ("Not implemented");
    }
}

void save_data (string_ref filename, AbstractParticleGenerator *source)
{
    if (ends_with (filename, ".dat"))
    {
        MostGeneralParticle most;

        std::ofstream ofs (filename);
        ofs.precision (15);

        while (source->get (&most))
        {
            if (!ofs)
                rt_error ("error writing to file " + filename);
            ofs << most.coords[0] << ' ' << most.coords[1] << ' ' << most.coords[2] << '\n';
        }

        ofs.close ();

        if (!ofs)
            rt_error ("error writing to file " + filename);
    }
    else
    {
        rt_error ("Not implemented");
    }
}
