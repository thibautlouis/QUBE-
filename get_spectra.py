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
import spec_utils

print "Reading dict file"
p = dict_utils.flipperDict()
p.read_from_file(sys.argv[1])

nside=p['nside']
ncomp=3
nd=p['nd']
lmax=p['lmax']
binsize=p['binsize']
nbin=lmax/binsize-1


patchDir='results/patches'
plotDir='results/plot'
specDir='results/spectra'
if len(sys.argv)> 2:
    patchDir='results/patches_%s'%sys.argv[2]
    plotDir='results/plot_%s'%sys.argv[2]
    specDir='results/spectra_%s'%sys.argv[2]


try:
    os.makedirs(specDir)
except:
    pass

try:
    os.makedirs(plotDir)
except:
    pass

minfo = sharp.map_info_healpix(nside)
ainfo = sharp.alm_info(lmax)
sht = sharp.sht(minfo, ainfo)
# Initialize alm but should not be needed ####################
ps=powspec.read_spectrum(p['theoryPS'])[:ncomp,:ncomp]
alm =curvedsky.rand_alm(ps,ainfo)
##########################################################

mask= healpy.read_map(p['maskfile'])



fields=['TT','EE','BB','TE','EB','TB']
alm_all=[]
for j in range(nd):
    maps= healpy.read_map(patchDir+'/map_%d.fits'%j, field=[0,1,2])
    T=maps[0]*mask
    Q=maps[1]*mask
    U=maps[2]*mask

    sht.map2alm( T, alm[0]); sht.map2alm((Q,U),alm[1:],spin=2)
    del maps
    alm_all+=[alm.copy()]

cl_auto_mean={}
for l1 in fields:
    cl_auto_mean[l1]=[]

for j in range(nd):
    Cl_array=healpy.alm2cl(alm_all[j])
    l,Cl= spec_utils.getClDict(Cl_array,fields)
    lBin,clbin= spec_utils.binCl(l,Cl,nbin,binsize,fields)
    filename='/autospec_windowed_%d%d'%(j,j)
    #spec_utils.plot(fields,plotDir+filename,lBin,clbin,theoryfile=p['theoryPS'])
    spec_utils.writeCl(fields,specDir+filename,lBin,clbin)
    for l1 in fields:
        cl_auto_mean[l1]+=[clbin[l1]]

for l1 in fields:
    cl_auto_mean[l1]=numpy.mean(cl_auto_mean[l1],axis=0)

filename='/mean_autospec_windowed'
#spec_utils.plot(fields,plotDir+filename,lBin,cl_auto_mean,theoryfile=p['theoryPS'])
spec_utils.writeCl(fields,specDir+filename,lBin,cl_auto_mean)


cl_cross_mean={}
for l1 in fields:
    cl_cross_mean[l1]=[]

for i in range(nd):
    for j in range(nd):
        if i >=j: continue
        print i,j
        Cl_array=healpy.alm2cl(alm_all[i],alm_all[j])
        l,Cl= spec_utils.getClDict(Cl_array,fields)
        lBin,clbin= spec_utils.binCl(l,Cl,nbin,binsize,fields)
        filename='/crossspec_windowed_%d%d'%(i,j)
        #spec_utils.plot(fields,plotDir+filename,lBin,clbin,theoryfile=p['theoryPS'])
        spec_utils.writeCl(fields,specDir+filename,lBin,clbin)
        for l1 in fields:
            cl_cross_mean[l1]+=[clbin[l1]]


for l1 in fields:
    print len(cl_cross_mean[l1])
    cl_cross_mean[l1]=numpy.mean(cl_cross_mean[l1],axis=0)

filename='/mean_crossspec_windowed'
#spec_utils.plot(fields,plotDir+filename,lBin,cl_cross_mean,theoryfile=p['theoryPS'])
spec_utils.writeCl(fields,specDir+filename,lBin,cl_cross_mean)



