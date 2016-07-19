#!/usr/bin/env python
import sys
import cPickle as pickle
import bz2

if sys.argv[1].endswith ('.bz2'):
    infile = bz2.BZ2File (sys.argv[1], 'r')
else:
    infile = open (sys.argv[1], 'r')
outfile = open ('coords.dat', 'w')
data = pickle.load (infile)
nzero = 0

peri = data.size_

periods = open ('periods', 'w')
print >>periods, peri[0]
print >>periods, peri[1]
print >>periods, peri[2]
periods.close ()

if data.__dict__.has_key ('extra_particles'):
    from array import array
    N = data.extra_particles
    for _ in range (N):
        d = array ('d')
        d.fromstring (infile.read (5*8))
        # extra attributes is not supported
        outfile.write ('%.15e %.15e\n' % (d[0], d[1]))
else:
    rt_error ('old format not supported yet')

print >>sys.stderr, '%i particles converted' % N

outfile.close ()
