#!/usr/bin/env python
# coding=utf8
# compute various derived quantities from a particle configuration in D=2.
# (c) 2015-2016 Sebastian Kapfer <sebastian.kapfer@fau.de>, FAU Erlangen

# usage:
#
# load data from directory dir/, compute pair correlator with bins of .2
#   order.py rho load dir/ corr rho2.dat discret .2
#
# load data explicitly, compute pair correlator for Minkowski tensor psi6
# (please cite http://dx.doi.org/10.1063/1.4774084 if you're using this)
#   order.py psi6 periods dir/periods coords dir/coords.dat corr g6.dat discret .2
# (this corresponds to Fig4a in http://dx.doi.org/10.1103/PhysRevLett.114.035702)
#
# make a color-coded image showing the psi6 order parameter
#   order.py psi6 load dir/ png psi6.png
# (Fig 2b-g in http://dx.doi.org/10.1103/PhysRevLett.114.035702)
#
# compute normalized scan through g(x,y), suitable for coherent averaging
# (see http://dx.doi.org/10.1103/PhysRevLett.107.155704)
#   order.py rho load dir/ cohgofr coherent.dat
# (Fig 3a in http://dx.doi.org/10.1103/PhysRevLett.114.035702)

import os
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import qhull
import matplotlib.mlab
import matplotlib.colors
import math

FFT_MEM_LIMIT = 1e9   # 1 GB
FFT_DATATYPE = np.complex128
FFT_MAX_CELLS = FFT_MEM_LIMIT / 8

def die (why):
    raise RuntimeError (why)

def read_flag (which):
    try:
        idx = sys.argv.index (which)
    except ValueError:
        return False
    del sys.argv[idx]
    return True

def read_arg (which, default):
    try:
        idx = sys.argv.index (which)
    except ValueError:
        return default
    del sys.argv[idx]
    try:
        ret = sys.argv[idx]
        del sys.argv[idx]
        return ret
    except IndexError:
        die ('missing argument for keyword %s' % which)

def load_periods (filename):
    """ load a periods file """
    peri = np.loadtxt (filename)
    if peri.shape == (3,):
        # drop third dimension (we don't care)
        return peri[:2]
    return peri

def load_coords (filename):
    """ load a coords.dat file """
    try:
        import pandas
        data = pandas.read_csv (filename, sep=' ')
        return data[:,:2]
    except:
        # maybe we don't have pandas. try again with plain numpy.
        pass
    try:
        data = np.loadtxt (filename)
    except:
        print >>sys.stderr, 'error loading %s (pwd = %s)' % (filename, os.getcwd ())
        raise
    return data[:,:2]

def normals_and_weights ():
    """compute normals and the associated weights from a Voronoï diagram"""
    diagr = setup_voronoi ()
    # represent vertices by complex numbers
    vert = diagr.vertices[:,0] + 1.j*diagr.vertices[:,1]
    rot90 = -1.j
    for i in xrange (len (diagr.points)):
        reg = diagr.regions[diagr.point_region[i]]
        pnt = vert[reg]
        avg = np.mean (pnt)
        pnt = pnt[np.argsort (np.angle (pnt-avg))]
        edge = pnt - np.roll (pnt, 1)
        w = np.abs (edge) + 1e-99
        n = edge / w * rot90
        w /= np.sum (w)
        yield (n, w)

def polygon_area (x,y):
    # https://stackoverflow.com/a/30408825
    return 0.5 * np.abs (np.dot (x, np.roll (y,1)) - np.dot (y, np.roll (x,1)))

def cell_areas ():
    """compute cell areas in the Voronoï diagram"""
    diagr = setup_voronoi ()
    # represent vertices by complex numbers
    vert = diagr.vertices[:,0] + 1.j*diagr.vertices[:,1]
    for i in xrange (len (diagr.points)):
        reg = diagr.regions[diagr.point_region[i]]
        pnt = vert[reg]
        yield polygon_area (pnt.real, pnt.imag)

def pad_points (coords, periods, overlap):
    """ add periodic padding to point set """
    imax = int (overlap/periods[0]) + 1
    jmax = int (overlap/periods[1]) + 1

    with_copies = list (coords)

    # pad with periodic copies
    for i in range (-imax, imax+1):
        for j in range (-jmax, jmax+1):
            if i == 0 and j == 0:
                continue
            shifted = coords + (i*periods[0], j*periods[1])
            select = np.all ((
                    shifted[:,0] > -overlap,
                    shifted[:,0] < periods[0]+overlap,
                    shifted[:,1] > -overlap,
                    shifted[:,1] < periods[1]+overlap
                ), axis=0)
            with_copies.extend (shifted[select])

    #print >>sys.stderr, 'padded', N, '->', len (with_copies)

    return with_copies

def cossin (angle):
    return np.array ([math.cos (angle), math.sin (angle)])

def round_1024 (x):
    """ round up to next multiple of 1024 """
    x = np.array (x)
    ret = np.ones_like (x, dtype=np.int32)
    for i, x in enumerate (x):
        ret[i] = int (math.ceil (x / 1024)) * 1024
    return ret

# idea from http://stackoverflow.com/questions/21242011/most-efficient-way-to-calculate-radial-profile
def radial_profile (x, y, data, lencell, discret):
    xx, yy = np.meshgrid (x, y)
    r = np.hypot (xx, yy) / discret
    r = r.astype (np.int)
    tbin = np.bincount (r.ravel (), data.ravel ())
    nr = np.bincount (r.ravel()) + 1e-99
    rapro = tbin / nr
    rmax = (data.shape * lencell).max () * .5
    rapro = rapro[:int (rmax / discret)]
    r = (np.arange (len (rapro)) + .5) * discret
    return r, rapro 

def pair_correlator (coords, periods, data, discret = .2):
    N = len (coords)
    assert N == len (data)
    # pad the density function
    numcell = round_1024 (periods[:2] / discret + .5)
    lencell = periods[:2] / numcell
    if 0:
        print 'N', N
        print 'numcell', numcell
        print 'lencell', lencell
    if np.product (numcell) > FFT_MAX_CELLS:
        die ("excessive memory use")
    # density function
    da = np.zeros (numcell, dtype=FFT_DATATYPE)
    idx0 = np.floor (coords[:,0] / lencell[0]).astype (np.int)
    idx1 = np.floor (coords[:,1] / lencell[1]).astype (np.int)
    assert np.all (idx0 >= 0)
    assert np.all (idx1 >= 0)
    da[idx0,idx1] += data
    # FFT convolution
    fda = np.fft.fft2 (da)
    corr = np.fft.ifft2 (fda.real**2 + fda.imag**2)
    corr = np.real (corr)
    corr /= N*N
    corr *= np.product (numcell)
    corr[0,0] = 0
    return corr, lencell

def complex_number_to_rgb (c):
    """
    represent a complex number by a RGB color.
    magnitude 1 is full saturation, magnitude 0 gives gray.
    complex phase is mapped to hue.
    """
    angle = np.angle (c) + math.pi
    angle /= 2 * math.pi
    satur = np.clip (np.abs (c), 0., 1.)
    value = np.ones_like (satur) * .8
    hsv = np.transpose ([angle, satur, value])
    return matplotlib.colors.hsv_to_rgb (hsv)

def RGB_tuples_to_PIL_image (rgb):
    assert len (rgb.shape) == 3
    from PIL import Image
    # some assembly required: flip y
    rgb = rgb[::-1,:,:]

    # work around a PIL bug: http://stackoverflow.com/questions/18325154/numpy-error-when-creating-gif-using-images2gif-py
    class WrapThing:
        def __init__ (self, other):
            self.other = other
            self.__array_interface__ = other.__array_interface__

        def tobytes (self):
            try:
                return self.other.tostring ()
            except:
                return self.other.tobytes ()

    return Image.fromarray (WrapThing (rgb))

def save_real_field_png (filename, data):
    assert len (data.shape) == 2
    # put data into green channel
    stuff = data.clip (0., 1.) * 255
    stuff = stuff.astype (np.uint8)
    zero  = np.zeros_like (stuff, dtype=np.uint8)
    stuff = np.transpose ([zero, stuff, zero], (1, 2, 0))
    # write to .png file
    img = RGB_tuples_to_PIL_image (stuff)
    img.save (filename, 'PNG')

def save_complex_field_png (filename, data):
    assert len (data.shape) == 2
    # convert to HSV color
    stuff = complex_number_to_rgb (data)
    stuff = (255*stuff).astype (np.uint8)
    # write to .png file
    img = RGB_tuples_to_PIL_image (stuff)
    img.save (filename, 'PNG')

def complex_number_to_gnuplot (data):
    """ represent complex numbers by HSV color code and put into GNUPLOT color format """
    stuff = complex_number_to_rgb (data)
    stuff = (255*stuff).astype (np.uint32)
    return 65536 * stuff[:,0] + 256 * stuff[:,1] + stuff[:,2]

# START
padded_coords = None
current_overlap = 0.
periods = None
coords = None
voro = None

def init ():
    """ peel input data from command line, load it, compute some simple statistics """
    global periods_filename, coords_filename
    coords_filename = read_arg ('coords', '')
    periods_filename = read_arg ('periods', '')
    load_filename = read_arg ('load', '')
    if load_filename != '':
        coords_filename = '%s/coords.dat' % load_filename
        periods_filename = '%s/periods' % load_filename
    global periods, coords
    periods = load_periods (periods_filename)
    coords = load_coords (coords_filename)
    for compo in (0, 1):
        print 'coord_minmax', compo, np.min (coords[:,compo]), np.max (coords[:,compo])
    print 'periods', periods
    global N, A, rho, dee
    N = len (coords)
    A = periods[0] * periods[1]
    rho = N/A
    print 'rho', rho
    dee = (math.pi*rho) ** (-.5)
    global padded_coords, current_overlap, voro
    padded_coords = None
    current_overlap = 0.
    voro = None

def setup_padding (min_overlap, rot_angle = 0.):
    global padded_coords, current_overlap, voro
    if min_overlap > current_overlap or rot_angle != 0.:
        padded_coords = pad_points (coords, periods, min_overlap)
        if rot_angle != 0.:
            c, s = cossin (rot_angle)
            def rotate (p):
                return np.array ([c*p[0] - s*p[1], s*p[0] + c*p[1]])
            center = periods/2
            padded_coords = [ center + rotate (p-center) for p in padded_coords ]
        padded_coords = np.array (padded_coords)
        current_overlap = min_overlap
        voro = None

def setup_voronoi ():
    global voro
    if voro is not None:
        return voro
    setup_padding (10*dee)
    while True:
        voro = qhull.Voronoi (padded_coords)
        # FIXME check no boundary vertices are touched
        return voro

def grid_and_save_field_png (filename, values, resolu):
    pixX = np.arange (.5, resolu+.5) / resolu * periods[0]
    pixY = np.arange (.5, resolu+.5) / resolu * periods[1]
    dataX = padded_coords[:,0]
    dataY = padded_coords[:,1]

    togrid = lambda (dataF): matplotlib.mlab.griddata (dataX, dataY, dataF, pixX, pixY, interp='linear')

    values = np.array (values)

    assert len (dataX) == len (values)
    assert len (dataY) == len (values)

    if values.dtype == np.complex128:
        # griddata doesn't do complex
        gridded = togrid (np.real (values)) + 1j * togrid (np.imag (values))
        save_complex_field_png (filename, gridded)
    elif values.dtype == np.int:
        die ('integer field not implemented for PNG output')
    else:
        gridded = togrid (values)
        save_real_field_png (filename, gridded)

def compute (what):
    if what[0:3] == 'psi':
        # Minkowski structure metrics, see http://dx.doi.org/10.1063/1.4774084
        angular = int (what[3:])
        MT = [ np.dot (np.power (n, angular), w) for (n, w) in normals_and_weights () ]
        return np.asarray (MT)
    elif what == 'rho':
        return np.ones (N)
    elif what == 'Z':
        # (Voronoi) coordination number.
        # Z is an integer, and most outputs don't now what to do with that.
        # "dat" output works, though.
        Z = [ len (n) for (n, w) in normals_and_weights () ]
        return np.asarray (Z)
    elif what == 'voronoi_volume':
        # obserable pinned to each particle is the Voronoï cell volume.
        Vc = np.array ([ a for a in cell_areas () ])
        return np.asarray (Vc)
    elif what == 'voronoi_density':
        # obserable pinned to each particle is the local density, defined as
        # 1/Vc where Vc is the Voronoï cell volume.
        Vc = np.array ([ a for a in cell_areas () ])
        return np.power (Vc, -1.)
    else:
        die ('unknown order parameter: %s' % what)

# command-line driver (if not used as a module)
if __name__ == "__main__":
    init ()
    dat_out = read_arg ('dat', '')
    png_out = read_arg ('png', '')
    gnuplot_out = read_arg ('gnu', '')
    corr_field_out = read_arg ('corrfield', '')
    corr_field_png_out = read_arg ('corrpng', '')
    corr_radial_out = read_arg ('corr', '')
    cohgofr_out = read_arg ('cohgofr', '')
    cohgofr_angular = int (read_arg ('cohgofr_angular', '6'))
    debug_cohgofr = read_flag ('debug_cohgofr')

    if dat_out + png_out + gnuplot_out + corr_field_out + corr_radial_out + corr_field_png_out + cohgofr_out == '':
        die ('no output, master?')

    if png_out + gnuplot_out != '':
        # we'll need padding to produce sensible images
        setup_padding (10*dee)

    order_name = sys.argv[1]
    del sys.argv[1]
    order = compute (order_name)

    # save as a two-column text file, col1 = real part, col2 = imag part
    # we don't output data for the padding points.
    # we don't output the second column for real data.
    if dat_out != '':
        if order.dtype == np.complex128:
            re = np.real (order[:N])
            im = np.imag (order[:N])
            np.savetxt (dat_out, zip (re, im))
        elif order.dtype == np.int:
            np.savetxt (dat_out, order[:N], fmt='%d')
        else:
            np.savetxt (dat_out, order[:N])

    # save as either a color-coded (HSV) .png or as a greenscale .png
    if png_out != '':
        resolu = int (read_arg ('resolu', 450))
        grid_and_save_field_png (png_out, order, resolu)

    # output for gnuplot or rdisk
    # FIXME this doesn't know what to do with real (as opposed to complex) data
    if gnuplot_out != '':
        colors = complex_number_to_gnuplot (order)
        outfile = open (gnuplot_out, 'w')
        outfile.write ('# %s: gnuplot me like\n' % order_name)
        outfile.write ('# p "%s"  u 1:2:(1/sqrt(pi)):3 w circl lc rgbcolor vari fs solid\n' % gnuplot_out)
        dataX = padded_coords[:,0]
        dataY = padded_coords[:,1]
        for xyc in zip (dataX, dataY, colors):
            outfile.write ('%.8e %.8e %i\n' % xyc)
        outfile.close ()

    # various correlation functions
    if corr_field_out + corr_radial_out + corr_field_png_out + cohgofr_out != '':
        # FIXME all uses of pair_correlator seem to call for the shifted version
        discret = float (read_arg ('discret', .01))
        corr, lencell = pair_correlator (coords, periods, order[:N])
        c_corr = np.fft.fftshift (corr)
        del corr
        c_center = np.array (c_corr.shape) // 2
        c_rmax = np.min (periods-8*lencell) / 2
        x = lencell[0] * (np.arange (c_corr.shape[0]) - c_center[0])
        y = lencell[1] * (np.arange (c_corr.shape[1]) - c_center[1])
        if corr_field_out != '':
            l,h = c_center-100, c_center+100
            export = c_corr[l[0]:h[0],l[1]:h[1]]
            np.savetxt (corr_field_out, export)
        print 'shapes', x.shape, y.shape, c_corr.shape
        r, rad_corr  = radial_profile (x, y, c_corr, lencell, discret)
        if corr_radial_out != '':
            np.savetxt (corr_radial_out, zip (r, rad_corr))
        # coherent g(r), see http://dx.doi.org/10.1103/PhysRevLett.107.155704
        if cohgofr_out != '':
            mark_scanline = 5 * np.max (np.abs (c_corr))
            num_scan_angles = cohgofr_angular // 2
            scan_angle_incr = 2 * math.pi / cohgofr_angular
            scan_angle_incr_deg = 360. / cohgofr_angular
            interp = RectBivariateSpline (x, y, c_corr)
            order_par = np.mean (compute ('psi%i' % cohgofr_angular)[:N])
            detected_ori = np.angle (order_par) / cohgofr_angular
            detected_ori_deg = 180. / math.pi * detected_ori
            print 'detected_ori_deg', detected_ori_deg
            r = np.arange (0., c_rmax, discret)
            r_dotted = r[::100]
            cohgofr = [ r ]
            for scan_angle in detected_ori + scan_angle_incr * np.arange (num_scan_angles):
                 c, s = cossin (scan_angle)
                 cg = interp.ev (c*r, s*r)
                 cohgofr.append (cg)
                 if debug_cohgofr:
                     for ii,jj in zip (c*r_dotted, s*r_dotted):
                         ii = int (ii/lencell[0] + c_center[0] + .5)
                         jj = int (jj/lencell[1] + c_center[1] + .5)
                         if 0 <= ii < c_corr.shape[0]:
                             if 0 <= jj < c_corr.shape[1]:
                                 c_corr[ii,jj] = mark_scanline
            fp = open (cohgofr_out, 'w')
            fp.write ('# coherent g(x,y), aligned with psi%(cohgofr_angular)i, along %(num_scan_angles)i scan directions in %(scan_angle_incr_deg).1f degree increments\n# sample orient %(detected_ori).4f %(detected_ori_deg).4f\n' % locals ())
            np.savetxt (fp, np.transpose (cohgofr))
            fp.close ()
        if corr_field_png_out != '':
            field = c_corr / np.max (np.abs (c_corr))
            field = np.clip (field, 0., 1.)
            field = np.power (field, .4)
            save_real_field_png (corr_field_png_out, field)
