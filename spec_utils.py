import numpy
import healpy
from enlib import powspec
import pylab

def plot(fields,filename,lbin,Cb,cov=None,theoryfile=None):
    lmax=int(numpy.max(lbin))
    
    if theoryfile != None:
        
        ps=powspec.read_spectrum(theoryfile)[:3,:3]
        clth={}
        clth['TT']=ps[0,0,:][:lmax+100]
        clth['EE']=ps[1,1,:][:lmax+100]
        clth['BB']=ps[2,2,:][:lmax+100]
        clth['TE']=ps[0,1,:][:lmax+100]
        clth['TB']=ps[0,1,:][:lmax+100]*0
        clth['EB']=ps[0,1,:][:lmax+100]*0
        lth=numpy.arange(len(clth['TT']))
        fth=lth*(lth+1)/(2*numpy.pi)
    for l1 in fields:
        if theoryfile != None:
            pylab.plot(lth,clth[l1]*fth,color='gray')
        if cov != None:
            pylab.errorbar(lbin,Cb[l1],numpy.sqrt(cov[l1+l1]),fmt='.')
        else:
            pylab.errorbar(lbin,Cb[l1],fmt='.')
        pylab.xlabel(r'$\ell$',fontsize=22)
        pylab.ylabel(r'$\ell(\ell+1) C^{%s}_\ell/(2 \pi)$'%l1,fontsize=22)
        pylab.savefig(filename+'_%s.png'%l1)
        pylab.clf()
        pylab.close()


def writeCl(fields,filename,lbin,Cb,cov=None):
    for l1 in fields:

        if cov != None:
            g = open(filename+'_%s.dat'%l1,mode="w")
            for k in xrange(len(lbin)):
                g.write("%f %e %e\n"%(lbin[k],Cb[l1][k],numpy.sqrt(cov[l1+l1])[k]))
            g.close()
        else:
            g = open(filename+'_%s.dat'%l1,mode="w")
            for k in xrange(len(lbin)):
                g.write("%f %e \n"%(lbin[k],Cb[l1][k]))
            g.close()


def getClDict(Cl_array,fields):
    l=numpy.arange(len(Cl_array[0]))
    Cl={}
    count=0
    for l1 in fields:
        Cl[l1]=Cl_array[count]
        count+=1
    return(l,Cl)


def getbinfile(nbin,binsize):
    binLower=numpy.zeros(nbin)
    binCenter=  numpy.zeros(nbin)
    binUpper=numpy.zeros(nbin)
    for i in range(nbin):
        binLower[i]=2+i*binsize
        binUpper[i]=2+(i+1)*binsize
        binCenter[i]=(binLower[i]+binUpper[i])/2
    
    return(binLower,binUpper,binCenter)

def binCl(l,cl,nbin,binsize,fields):
    
    binLower,binUpper,binCenter = getbinfile(nbin,binsize)
    nBins = len(binCenter)
    clBin={}
    f=l*(l+1)/(2*numpy.pi)
    
    for l1 in fields:
    
        lBin = numpy.zeros(nBins)
        clBin[l1] =  numpy.zeros(nBins)
        for i in xrange(len(binCenter)):
            idx = numpy.where((l<binUpper[i]) & (l>=binLower[i]))
            lBin[i] = (l[idx]).mean()
            clBin[l1][i] = (cl[l1][idx]*f[idx]).mean()
    
    return lBin,clBin




