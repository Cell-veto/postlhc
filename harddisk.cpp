// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "harddisk.hpp"

static Register <ChainRunner <HardInteraction, Monodisperse2D>>
    one ("harddisk/mono2d");
static Register <ChainRunner <HardInteraction, Monodisperse3D>>
    two ("hardsphere/mono3d");
