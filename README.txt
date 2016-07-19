

`postlhc': C++ implementation of Cell-veto Monte Carlo


This is a complete implementation of Cell-veto Monte Carlo (CVMC) for a number
of typical model potentials, in 2D and 3D.  The code has been tested in
particular in the 2D / IPL (inverse power law) case.  Some more obscure
features may still have bugs.  Also, performance tuning is still missing.  If
you run into problems, feel free to contact <sebastian.kapfer@fau.de>.

If you're using this, please cite the original publication,
Sebastian C. Kapfer, Werner Krauth:
"Cell-veto Monte Carlo algorithm for long-range systems"
arXiv:1606.06780 [cond-mat.stat-mech]


Building

In order to build the program, type the following:

touch features.mk
make


Demo program

The demo program is intended to be modified according to your needs.
As is, it loads some starting configuration (you can find examples
in test/data), performs CVMC simulations, and writes periodic
snapshots.  Usage:

./postlhc stor mono2d load test/data/hex-18144/ inter jellium4 exponent 1 strength 140.000000 seed 752 gofr

stor:  select storage mechanism; mono2d for 2D, mono3d for 3D.
load:  load a dataset, which is a directory containing two files:
       periods  gives the dimensions of the periodic box
       coords.dat  gives the particle's coordinates
inter: set interaction type.
       jellium4: lattice-screened IPL
       jellium3: screened IPL
       ipl: bare IPL (exponent > D-1)
       lj: Lennard-Jones
       harddisk: hard disks
       truncipl:  truncated IPL
exponent, strength: set parameters of the interaction
seed:  for the RNG.
gofr:  enable writing radially averaged pair correlators along with the snapshots.


Implementing new interactions

See lj.cpp for a simple example with long-range interactions, and both
attractive and repulsive forces.  If you have questions, feel free to mail
<sebastian.kapfer@fau.de>.


Order parameters

order.py computes some of the observables required to assess the
phase transition behavior of two-dimensional systems.  If you
use this, please cite http://dx.doi.org/10.1103/PhysRevLett.114.035702
as well.
