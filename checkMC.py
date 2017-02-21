import dict_utils
import sys
import os
from enlib import powspec
import numpy
import pylab
import spec_utils


print "Reading dict file"
p = dict_utils.flipperDict()
p.read_from_file(sys.argv[1])

nSims=p['nSims']
lmax=p['lmax']
binsize=p['binsize']
nbin=lmax/binsize-1


fields=['TT','EE','BB','TE','EB','TB']


ps=powspec.read_spectrum(p['theoryPS'])[:3,:3]
clth={}
clth['TT']=ps[0,0,:][:lmax+100]
clth['EE']=ps[1,1,:][:lmax+100]
clth['BB']=ps[2,2,:][:lmax+100]
clth['TE']=ps[0,1,:][:lmax+100]
clth['TB']=ps[0,1,:][:lmax+100]*0
clth['EB']=ps[0,1,:][:lmax+100]*0
lth=numpy.arange(len(clth['TT']))


lb,Cb= spec_utils.binCl(lth,clth,nbin,binsize,fields)

print lb.shape



for l1 in fields:
    cl_all=[]
    for iii in range(nSims):
        l,cl,error=numpy.loadtxt('results/spectra_%03d/spectrum_%s.dat'%(iii,l1),unpack=True)
        cl_all+=[cl]
    mean=numpy.mean(cl_all,axis=0)
    std=numpy.std(cl_all,axis=0)

    pylab.plot(lb,Cb[l1])
    pylab.errorbar(l,mean,std,fmt='.',color='red')
    pylab.show()