// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "lj.hpp"

static Register <ChainRunner <LennardJones, Monodisperse2D>>
    one ("lj/mono2d");
static Register <ChainRunner <LennardJones, Monodisperse3D>>
    two ("lj/mono3d");
