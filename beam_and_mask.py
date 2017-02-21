import healpy
from enlib import powspec,curvedsky,enmap,sharp
import numpy
import mcm_code
import utils
import pylab
import time
import dict_utils
import sys
import os


def makeBeam(beamSize,lmax):
    rad2min = 180./numpy.pi * 60.
    ell = numpy.arange(lmax)
    bl = numpy.exp(-ell*(ell+1)*(beamSize/rad2min)**2 / (16*numpy.log(2)) )
    return(ell,bl)

def makemask(nside,fracsky, apod):
    npix=12*nside**2
    mask=numpy.ones(npix)
    min=int(npix/2-(1-fracsky)*npix/2)
    max=int(npix/2+(1-fracsky)*npix/2)
    mask[min:max]=0
    mask=healpy.smoothing(mask,0.1)
    return(mask)


print "Reading dict file"
p = dict_utils.flipperDict()
p.read_from_file(sys.argv[1])

nside=p['nside']
beamsize=p['beamsize']
mask_fsky=p['mask_fsky']
lmax=p['lmax']

l,bl=makeBeam(beamsize,2*lmax)
pylab.plot(l,bl)
pylab.show()

g = open('data/beamfile_%0.2f.dat'%beamsize,mode="w")

for k in xrange(len(l)):
    g.write("%f %e \n"%(l[k],bl[k]))
g.close()

mask=makemask(nside,mask_fsky, 0.1)

healpy.mollview(mask)
pylab.show()

healpy.write_map('data/mask_%d_%0.2f.fits'%(nside,mask_fsky),mask)

