import numpy
import healpy


def apply_beam(ps,bl):
    nc=ps.shape[0]
    lhigh=numpy.minimum(bl.shape,ps.shape[2])[0]
    ps_conv=numpy.zeros((nc,nc,lhigh))
    for i in range(nc):
        for j in range(nc):
            ps_conv[i,j,:]=ps[i,j,:lhigh]*bl[:lhigh]**2
    print 'convolved PS up to: %d'%lhigh
    return(ps_conv)


def addWhiteNoise(map,rmsArcmin,nside):
    
    radToMin = 180/numpy.pi*60
    pixArea = radToMin**2 * healpy.nside2pixarea(nside)
    rms = rmsArcmin/numpy.sqrt(pixArea)
    noise = numpy.random.normal( scale = rms, size = map.shape[0] )
    map += noise
    return map
