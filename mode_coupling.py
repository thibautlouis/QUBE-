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

lmax=p['lmax']
binsize=p['binsize']
nbin=lmax/binsize-1

patchDir='results/patches_000'
mcmDir='results/mcm'
try:
    os.makedirs(mcmDir)
except:
    pass

mask= healpy.read_map(p['maskfile'])


healpy.mollview(mask)
pylab.savefig(patchDir+'/mask.png')
pylab.clf()
pylab.close()

ell,bl=numpy.loadtxt(p['beamfile'],unpack=True)

# get the mask and beam ready to be sent to the fortran code
Pl=healpy.anafast(mask)
l=numpy.arange(len(Pl))

bl=bl[:len(l)]
Pl*=(2*l+1)

# There is 4 mode coupling matrix to compute
mcm=numpy.zeros((lmax,lmax))
mcm_p=numpy.zeros((lmax,lmax))
mcm_pp=numpy.zeros((lmax,lmax))
mcm_mm=numpy.zeros((lmax,lmax))
# and the binned one
mbb=numpy.zeros((nbin,nbin))
mbb_p=numpy.zeros((nbin,nbin))
mbb_pp=numpy.zeros((nbin,nbin))
mbb_mm=numpy.zeros((nbin,nbin))

########### FORTRAN COMPUTATION ###########
mcm_code.calc_mcm(Pl[:],bl[:],mcm.T,mcm_p.T,mcm_pp.T,mcm_mm.T)
############################################
mcm_code.bin_mcm(mcm.T, binsize, mbb.T)
mcm_code.bin_mcm(mcm_p.T, binsize, mbb_p.T)
mcm_code.bin_mcm(mcm_pp.T, binsize, mbb_pp.T)
mcm_code.bin_mcm(mcm_mm.T, binsize, mbb_mm.T)
############################################

numpy.savetxt(mcmDir+'/mbb.dat',mbb)
numpy.savetxt(mcmDir+'/mbb_p.dat',mbb_p)
numpy.savetxt(mcmDir+'/mbb_pp.dat',mbb_pp)
numpy.savetxt(mcmDir+'/mbb_mm.dat',mbb_mm)