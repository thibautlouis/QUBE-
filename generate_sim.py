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

print "Reading dict file"
p = dict_utils.flipperDict()
p.read_from_file(sys.argv[1])

nside=p['nside']
ncomp=3
plot=p['plot']
nd=p['nd']
lmax=p['lmax']
noise=p['noise']
nSims=p['nSims']



minfo = sharp.map_info_healpix(nside)
ainfo = sharp.alm_info(lmax)
sht = sharp.sht(minfo, ainfo)
map = numpy.zeros([ncomp,12*nside**2])
ps=powspec.read_spectrum(p['theoryPS'])[:ncomp,:ncomp]
l,bl=numpy.loadtxt(p['beamfile'],unpack=True)

ps=utils.apply_beam(ps,bl)

try:
    os.makedirs('results')
except:
    pass


for iii in range(nSims):
    
    patchDir='results/patches_%03d'%iii
    try:
        os.makedirs(patchDir)
    except:
        pass

    alm =curvedsky.rand_alm(ps,ainfo)
    sht.alm2map(alm[0], map[0])
    sht.alm2map(alm[1:],map[1:],spin=2)

    for j in range(nd):

        if noise['apply']==True:
            rmsArcmin= float(noise['rmsArcmin'])
            print "Add White Noise to the realization %s uK.arcmin"%(rmsArcmin)
            T=utils.addWhiteNoise(map[0].copy(),rmsArcmin*numpy.sqrt(nd),nside)
            Q=utils.addWhiteNoise(map[1].copy(),rmsArcmin*numpy.sqrt(nd)*numpy.sqrt(2),nside)
            U=utils.addWhiteNoise(map[2].copy(),rmsArcmin*numpy.sqrt(nd)*numpy.sqrt(2),nside)

        healpy.write_map(patchDir+'/map_%d.fits'%j, (T,Q,U))

