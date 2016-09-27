// (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen
#include "storage.hpp"
#include <memory>

void AbstractStorage::load_periods (string_ref filename)
{
    set_periods (Periods::from_file (filename));
}

void AbstractStorage::add_data (string_ref filename)
{
    ::add_data (this, filename);
}

void AbstractStorage::save_data (string_ref filename)
{
    ::save_data (filename, this);
}

AbstractStorage *make_storage (string_ref encoding_typecode)
{
    return Factory <AbstractStorage>::make (encoding_typecode);
}

static Register <Storage <Monodisperse2D>> one ("mono2d");
static Register <Storage <Tagged <Monodisperse2D>>> three ("tagged_mono2d");
static Register <Storage <Monodisperse3D>> two ("mono3d");
