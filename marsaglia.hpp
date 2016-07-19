// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef MTW_DISCR_DISTR_HPP_INCLUDED 
#define MTW_DISCR_DISTR_HPP_INCLUDED 

#include "tools.hpp"

template <typename VALUE>
struct MtwDiscreteSampler
{
private:
    typedef VALUE value_t;
    typedef unsigned long index_t;
    std::vector <double> weights;
    double weights_total;
    std::vector <value_t> values;
    static const unsigned NUM_LEVELS = 5;
    static const index_t P_SCALE = 1ull << 30;
    std::vector <index_t> table[NUM_LEVELS];
    index_t tbase[NUM_LEVELS+1];

    index_t sample (index_t u) const
    {
        for (unsigned n = 0; n != NUM_LEVELS; ++n)
        {
            if (u < tbase[n+1])
            {
                u -= tbase[n];
                u >>= 6 * (NUM_LEVELS-1-n);
                return table[n].at (u);
            }
        }

        std::cerr << "Broken lookup table in MtwDiscreteSampler" << ABORT;
        return 0; // squelch warning
    }

    void init_table ()
    {
        // convert to scaled representation
        std::vector <index_t> scaled_ws;
        for (double w : weights)
            // + .5 would be fair rounding
            // + 1. to over-allocate a bit (compensates for rounding)
            scaled_ws.push_back (P_SCALE / weights_total * w + 1.);

        // set up lookup tables
        index_t entry_weight = P_SCALE;
        index_t total_alloc = 0;
        for (unsigned n = 0; n != NUM_LEVELS; ++n)
        {
            table[n].clear ();
            tbase[n] = total_alloc;
            entry_weight /= 64;
            for (index_t i = 0; i != weights.size (); ++i)
            {
                while (scaled_ws[i] >= entry_weight)
                {
                    scaled_ws[i] -= entry_weight;
                    table[n].push_back (i);
                    total_alloc += entry_weight;
                }
            }
        }
        tbase[NUM_LEVELS] = total_alloc;

        std::cerr << "mtw_init alloci " << total_alloc << ' ' << P_SCALE
            << "\nmtw_init rawrate " << weights_total << '\n';

        if (total_alloc < P_SCALE)
            std::cerr << "Initializing MtwDiscreteSampler: rounding error" << ABORT;
    }

public:
    void clear ()
    {
        weights.clear ();
        values.clear ();
    }

    void add (double weight, const value_t &value)
    {
        if (weight == 0)
            return;
        weights.push_back (weight);
        values.push_back (value);
    }

    void finish ()
    {
        // normalize weights
        weights_total = 0.;
        for (double w : weights)
            weights_total += w;
        if (weights.size () == 0)
            return;

        init_table ();
    }

    double total_weights () const
    {
        return weights_total;
    }

    const double &random_sample (value_t *out, RandomContext *random) const
    {
        index_t i = sample (random->uint (P_SCALE));
        *out = values.at (i);
        return weights.at (i);
    }
};

#endif /* MTW_DISCR_DISTR_HPP_INCLUDED */
