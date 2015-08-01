import numpy as np
from numpy import radians, degrees, sin, cos, arctan2, hypot, tan

def spheredist(ra1, dec1, ra2, dec2):
    """Returns great circle distance (and position angle).  Inputs in degrees.

    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.  From E. Tollerud, with position angle calculation added."""

    from numpy import radians, degrees, sin, cos, arctan2, hypot, tan

    # terminology from the Vicenty formula - lambda and phi and
    # "standpoint" and "forepoint"
    lambs = radians(ra1)
    phis = radians(dec1)
    lambf = radians(ra2)
    phif = radians(dec2)

    dlamb = lambf - lambs

    numera = cos(phif) * sin(dlamb)
    numerb = cos(phis) * sin(phif) - sin(phis) * cos(phif) * cos(dlamb)
    numer = hypot(numera, numerb)
    denom = sin(phis) * sin(phif) + cos(phis) * cos(phif) * cos(dlamb)

    theta  = arctan2(sin(dlamb), cos(phis) * tan(phif) - sin(phis) * cos(dlamb))
    
    return degrees(arctan2(numer, denom)), degrees(theta)



def deproject_coords(center = None, ra = None, dec = None, pa  = 0, inclination = 0., **kwargs):
    """Adapted from J. Moustakas' im_hii_region_deproject.pro.
    Uses simple geometry, should be changed to use spherical
    distances as in my IDL code"""   

    from numpy import radians, degrees, sin, cos, arctan2


    dra = (ra - center[0]) * np.cos(radians(dec)) * 3600
    ddec  = (dec -center[1]) * 3600
    
    x1 = 0-dra # east is usually positive; switch the coordinate axis
    y1 =  ddec

    theta = radians(pa - 90.0)  #because PA is measured east of north (y-axis) instead of east of west.
    xp =  x1*cos(theta) + y1*sin(theta)
    yp = -x1*sin(theta) + y1*cos(theta)

    ypp = yp/cos(radians(inclination))   # de-project
    radii = hypot(xp,ypp) # [arcsec]

    if (x1 == 0.0):
        phi = 0.0
    else :
        phi = degrees(arctan2(ypp, xp))

    return radii, phi

def read_kappa(kappaname, wavelength):
    wave, kappa = np.loadtxt(kappaname, usecols=(0, 4), unpack = True)
    oo  = np.argsort(wave)
    wave = wave[oo]*1e4
    kappa = kappa[oo]
    new_kappa = np.interp(wave,kappa,wavelength,new_kappa)
    return new_kappa/new_kappa.max()

def read_isrf(isrfname, wavelength):
    wave, isrf_nu = np.loadtxt(isrfname,unpack = True, usecols = (0,1))
    wave = wave*1e4
    isrf_nu=isrf_nu/isrf_nu.max()
    isrf = np.interp(wave,isrf_nu/(wave**2),wavelength)
    return isrf/isrf.max()  #I_lambda

def read_mp(self, galaxy = '', component = 'tot', mpdir = None):
    fname = mpdir+galaxy+component+'.sed'
    logwave, logflux, logflux_intrinsic = np.loadtxt(fname, usecols = (0,1,2), unpack = True)
    spec = {}
    spec['wavelength'] = 10**logwave
    spec['f_lambda'] = 10**logflux
    spec['f_lambda_int'] = 10**logflux_intrinsic
    
    pass




