// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#ifndef HARDDISK_HPP_INCLUDED 
#define HARDDISK_HPP_INCLUDED 

#include "ecmc.hpp"

struct HardInteraction : Interaction
{
    double diameter;

    HardInteraction ()
    {
        diameter = 2;
    }

    double sr_repulsion_range () const
    {
        return diameter;
    }

    void set_parameter (string_ref name, double value)
    {
        if (name == "diameter")
        {
            if (! (value > 0.))
                rt_error ("Invalid particle diameter");
            diameter = value;
        }
        else
        {
            rt_error ("Invalid parameter: " + name);
        }
    }

    double random_repulsive_lift (double, RandomContext *)
    {
        // hard disks/spheres always collide at distance 'diameter'
        return sq (diameter);
    }
};

#endif /* HARDDISK_HPP_INCLUDED */
