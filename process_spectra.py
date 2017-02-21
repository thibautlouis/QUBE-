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
nS=p['nS']
lmax=p['lmax']
binsize=p['binsize']
nbin=lmax/binsize-1


patchDir='results/patches'
mcmDir='results/mcm'
specDir='results/spectra'
plotDir='results/plot'

if len(sys.argv)> 2:
    patchDir='results/patches_%s'%sys.argv[2]
    plotDir='results/plot_%s'%sys.argv[2]
    specDir='results/spectra_%s'%sys.argv[2]



mbb=numpy.loadtxt(mcmDir+'/mbb.dat')
mbb_p=numpy.loadtxt(mcmDir+'/mbb_p.dat')
mbb_pp=numpy.loadtxt(mcmDir+'/mbb_pp.dat')
mbb_mm=numpy.loadtxt(mcmDir+'/mbb_mm.dat')

mode_coupling=numpy.zeros((6*nbin,6*nbin))
mode_coupling[:nbin,:nbin]=mbb
mode_coupling[nbin:2*nbin,nbin:2*nbin]=mbb_p
mode_coupling[2*nbin:3*nbin,2*nbin:3*nbin]=mbb_p
mode_coupling[3*nbin:4*nbin,3*nbin:4*nbin]=mbb_pp
mode_coupling[4*nbin:5*nbin,4*nbin:5*nbin]=mbb_pp-mbb_mm
mode_coupling[5*nbin:6*nbin,5*nbin:6*nbin]=mbb_pp
mode_coupling[3*nbin:4*nbin,5*nbin:6*nbin]=mbb_mm
mode_coupling[5*nbin:6*nbin,3*nbin:4*nbin]=mbb_mm

fields=['TT','TE','TB','EE','EB','BB']


Cb_cross_vec=[]
Cb_auto_vec=[]

for l1 in fields:
    lbin,clbin_cross_mean=numpy.loadtxt(specDir+'/mean_crossspec_windowed'+'_%s.dat'%l1,unpack=True)
    Cb_cross_vec=numpy.append(Cb_cross_vec,clbin_cross_mean)
    lbin,clbin_auto_mean=numpy.loadtxt(specDir+'/mean_autospec_windowed'+'_%s.dat'%l1,unpack=True)
    Cb_auto_vec=numpy.append(Cb_auto_vec,clbin_auto_mean)


t=time.time()
mbb_inv=numpy.linalg.inv(mode_coupling)
print time.time()-t
Cb_vec_deconvolve_cross=numpy.dot(mbb_inv,Cb_cross_vec)
Cb_vec_deconvolve_auto=numpy.dot(mbb_inv,Cb_auto_vec)

Cb_auto={}
Cb_cross={}
count=0
for l1 in fields:
    Cb_auto[l1]=Cb_vec_deconvolve_auto[count*nbin:(count+1)*nbin]
    Cb_cross[l1]=Cb_vec_deconvolve_cross[count*nbin:(count+1)*nbin]
    count+=1

Cb_auto['ET']=Cb_auto['TE']
Cb_auto['BE']=Cb_auto['EB']
Cb_auto['BT']=Cb_auto['TB']


mask= healpy.read_map(p['maskfile'])

#Cheap way to get fsky, will do something better in the futur
s1=numpy.sum(mask)
s2=numpy.sum(mask*0+1)
fsky=s1/s2
print fsky

cov={}
counts1=0
for spec1 in fields:
    counts2=0
    for spec2 in fields:
        a,b,c,d=spec1[0],spec1[1],spec2[0],spec2[1]
        cov[spec1+spec2]=Cb_auto['%s%s'%(a,c)]*Cb_auto['%s%s'%(b,d)]+Cb_auto['%s%s'%(a,d)]*Cb_auto['%s%s'%(b,c)]
        fac=1./(fsky*(2*lbin[:]+1)*binsize)
        cov[spec1+spec2]*=fac
        counts2+=1
    counts1+=1


Cl_th={}
lth,Cl_th['TT'],Cl_th['EE'],Cl_th['BB'],Cl_th['TE']=numpy.loadtxt('data/bode_almost_wmap5_lmax_1e4_lensedCls.dat',unpack=True)
Cl_th['TB']=Cl_th['TT']*0
Cl_th['EB']=Cl_th['TT']*0

filename='/spectrum'
spec_utils.plot(fields,plotDir+filename,lbin,Cb_cross,cov=cov, theoryfile=p['theoryPS'])
spec_utils.writeCl(fields,specDir+filename,lbin,Cb_cross,cov=cov)



