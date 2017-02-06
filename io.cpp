// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "storage.hpp"
#include <fstream>
#include <memory>
#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::starts_with;
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

struct IoError {};

// basic I/O for monodisperse particles.
static
MostGeneralParticle parse_line_from_file (unsigned no, string_ref line)
{
    if (MAX_DIM > 3u)
        std::cerr << "DIM>3 is not implemented" << ABORT;

    // read columns
    std::vector <string> cs = string_split (line);
    if (cs.size () > MAX_DIM)
        rt_error ("too many columns");
    double cd[6] = { 0 };
    for (size_t n = 0; n != cs.size (); ++n)
    {
        size_t consumed;
        cd[n] = std::stod (cs[n], &consumed);
        if (! (cd[n] >= 0.))
            rt_error ("negative coordinate: \"" + cs[n] + "\"");
        if (consumed != cs[n].size ())
            rt_error ("field not converted completely: \"" + cs[n] + "\"");
    }

    // put into MostGeneralParticle
    MostGeneralParticle part;
    std::memset (&part, 0, sizeof part);
    part.coords[0] = cd[0];
    part.coords[1] = cd[1];
    part.coords[2] = cd[2];
    part.tag = no;
    part.disp[0] = cd[3];
    part.disp[1] = cd[4];
    part.disp[2] = cd[5];
    return part;
}

void add_data (AbstractParticleSink *sink, string_ref filename)
{
    Extrema <double> particle_coords[MAX_DIM];

    if (ends_with (filename, ".dat"))
    {
        std::ifstream ifs (filename);
        unsigned N = 0;
        string line;
        while (std::getline (ifs, line))
        {
            MostGeneralParticle part = parse_line_from_file (N, line);
            for (unsigned n = 0; n != MAX_DIM; ++n)
                particle_coords[n].add (part.coords[n]);
            sink->put (part);
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

static
void write_mono_particle_to_stream (std::ostream &os, const MostGeneralParticle &part, bool tagged)
{
    if (MAX_DIM > 3u)
        std::cerr << "DIM>3 is not implemented" << ABORT;

    os << part.coords[0] << ' ' << part.coords[1] << ' ' << part.coords[2];
    
    if (tagged)
        os << ' ' << part.disp[0] << ' ' << part.disp[1] << ' ' << part.disp[2];
        
    os << '\n';

    if (!os)
        throw IoError ();
}

static
void write_mono_particles_to_stream (std::ostream &os, AbstractParticleGenerator &source, bool tagged)
{
    if (tagged)
    {
        // sort particles by tag to keep original order
        // also, we'll output displacements
        std::vector <MostGeneralParticle> list;
        MostGeneralParticle part;
        while (source.get (&part))
            list.push_back (part);
        std::sort (list.begin (), list.end (), MostGeneralParticle::by_tag);
        size_t tag = 0;
        for (const auto &p : list)
        {
            if (p.tag != tag++)
                std::cerr << "missing particle with tag " << (tag-1) << ABORT;
            write_mono_particle_to_stream (os, p, true);
        }
    }
    else
    {
        // just write them out in generator order
        MostGeneralParticle part;
        while (source.get (&part))
            write_mono_particle_to_stream (os, part, false);
    }
}

static
void save_data (string_ref filename, AbstractParticleGenerator &source, bool tagged)
{
    if (ends_with (filename, ".dat"))
    {
        std::ofstream ofs (filename + ".tmp");
        ofs.precision (15);

        try
        {
            if (!ofs)
                throw IoError ();

            write_mono_particles_to_stream (ofs, source, tagged);
            ofs.close ();

            if (!ofs)
                throw IoError ();
        }
        catch (IoError)
        {
            // try to unlink temporary file if writing failed
            remove_file (filename + ".tmp", true);
            rt_error ("error writing to file " + filename + ".tmp");
        }

        rename_file (filename + ".tmp", filename);
    }
    else
    {
        rt_error ("Not implemented");
    }
}

void save_data (string_ref filename, AbstractStorage *stor)
{
    std::unique_ptr <AbstractParticleGenerator> gen (stor->all_particles ());
    bool tagged = starts_with (stor->typecode (), "tagged_");
    save_data (filename, *gen, tagged);
}
