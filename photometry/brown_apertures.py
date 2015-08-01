import glob
import numpy as np
from sedpy import photometer, ds9region
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

def read_brown_apertures(cat={}, aperturefile='data/brown_apertures.txt'):
    with open(aperturefile) as f:
        lines = f.readlines()

    # dtype = [('galaxy_name', 'S100'), ('height', np.float64),
    #          ('width', np.float64), ('PA', np.float64)]
    # data = np.zeros(len(lines), dtype=np.dtype(dtype))

    hdr_line = len(lines)
    for i, line in enumerate(lines):
        if line.startswith('(arcsec)'):
            hdr_line = i
        if i > hdr_line:
            n, h, w, p, s = process_aperture_line(line)
            aperture = {'height':h, 'width':w, 'PA':p, 'specref':s}
            # data[i] = tuple(dat)
            if n in cat:
                cat[n].update(aperture)
            else:
                cat[n] = aperture
    return cat  # data[(hdr_line+1):]

def process_aperture_line(line):
    cols = line.split('\t')
    name = cols[0].lower().replace(' ','')
    height, width = [float(c.strip()) for c in cols[1].split('x')]
    pa = float(cols[2])
    specsource = cols[3]
    return name, height, width, pa, specsource


def read_brown_coordinates(cat={}, coordfile='data/brown_coordinates.txt'):
    with open(coordfile) as f:
        lines = f.readlines()
    hdr_line = len(lines)
    for i, line in enumerate(lines):
        if line.startswith('Class'):
            hdr_line = i
        if i > hdr_line:
            if line.startswith('\n'):
                break
            n, ra, dec = process_coordinate_line(line)
            # data[i] = tuple(dat)
            if n in cat:
                cat[n].update({'ra':ra, 'dec':dec})
            else:
                cat[n] = {'ra':ra, 'dec':dec}
    return cat


def process_coordinate_line(line):
    cols = line.split('\t')
    name = cols[0].lower().replace(' ','')
    coords = [float(c) for c in cols[1].split()]
    ra = 15 * (coords[0] + coords[1] / 60. + coords[2] / 3600.)
    dec = (coords[3] + coords[4] / 60. + coords[5] / 3600.)
    if '-' in cols[1] and dec > 0:
        dec *= -1
    
    return name, ra, dec


def make_ds9_region(info, **extras):
    ra, dec = dumb_corners(**info)
    defstring = [str(v) for pair in zip(ra, dec) for v in pair]
    defstring = ','.join(defstring)
    reg = ds9region.Polygon(defstring)
    return reg


def dumb_corners(ra=None, dec=None, PA=None,
                 height=None, width=None, **extras):
    """Approximate using 2-D plane"""
    size = np.array([width / 3600., height / 3600.])
    theta = np.deg2rad(PA)
    T = np.array([[np.cos(theta), np.sin(theta)],
                  [-np.sin(theta), np.cos(theta)]])
    corners = np.array([[-0.5, -0.5],
                        [-0.5, 0.5],
                        [0.5, 0.5],
                        [0.5, -0.5]])
    newcorners = np.dot(corners * size[None, :], T) 
    outra = ra + newcorners[:,0] / np.cos(np.deg2rad(dec))
    outdec = dec + newcorners[:,1]
    return outra, outdec


def smart_corners(ra=None, dec=None, PA=None,
                 height=None, width=None):
    #convert rotation axis and corners to cartesian coordinates

    #build rotation matrix

    pass


def photometer(image, header, regions, pad=[0,0], mef=False):
    """Given an image name (including path) and a set of regions,
    produce total fluxes within those regions
    """
    wcs = pywcs.WCS(header)
    try:
        cd = wcs.wcs.cd
    except:
        cd = np.zeros([2, 2])
        cd[np.diag_indices_from(cd)] = wcs.wcs.cdelt[0:2]
    ps = np.hypot(*cd*3600.)
    yy, xx = np.indices(image.shape)
    if mef:
        ra, dec, _ = wcs.wcs_pix2world(xx.flatten(), yy.flatten(), np.array([0]), 0)
    else:
        ra, dec = wcs.wcs_pix2world(xx.flatten(), yy.flatten(), 0)
    ra = ra.flatten()
    dec = dec.flatten()
    points = np.vstack((ra,dec)).T
    flatim = image.flatten()
    flux = []
    for region in regions:
        sel = region.contains(points=points, fast=False, pad=pad)
        flux.append(np.nansum(flatim[sel]))
    return np.array(flux), ps, header.get('BUNIT', None)

   
def measure_pacs_flux(imagenames, gal_info):
    pacsbands = ['pacs70', 'pacs100', 'pacs160']
    reg = make_ds9_region(gal_info)
    fluxes, uncertainties, bands = [], [], []
    for imname in imagenames:
        image = pyfits.getdata(imname)
        hdr = pyfits.getheader(imname)
        im, unc = image[0:2,:,:]
        flux, ps, units =  photometer(im, hdr, [reg], mef=True)
        var, _, _ = photometer(unc**2, hdr, [reg], mef=True)
        # get flux in Jy, assuming image data in MJy/sr
        # flux = (flux * 1e6) * (ps.prod() * 2.35e-11)
        unc = np.sqrt(var) # * 1e6 * (ps.prod() * 2.35e-11)
        fluxes.append(flux)
        uncertainties.append(unc)
        bands.append([b for b in pacsbands if b in imname][0])
    return fluxes, uncertainties, bands


def measure_spire_flux(imagenames, gal_info):
    spirebands = ['spire250', 'spire350', 'spire500']
    reg = make_ds9_region(gal_info)
    fluxes, uncertainties, bands = [], [], []
    for imname in imagenames:
        uncname = imname.replace('scan.fits', 'scan.unc.fits')
        im = pyfits.getdata(imname)
        hdr = pyfits.getheader(imname)
        unc = pyfits.getdata(uncname)
        flux, ps, units =  photometer(im, hdr, [reg])
        var, _, _ = photometer(unc**2, hdr, [reg])
        # get flux in Jy, assuming image data in MJy/sr
        flux = (flux * 1e6) * (ps.prod() * 2.35e-11)
        unc = (np.sqrt(var) * 1e6) * (ps.prod() * 2.35e-11)
        fluxes.append(flux)
        uncertainties.append(unc)
        bands.append([b for b in spirebands if b in imname][0])
    return fluxes, uncertainties, bands


def find_images(imnames, galaxyname):
    possible = []
    for p in imnames:
        if galaxyname.upper() in p.upper():
            possible.append(p)
    return possible
    
if __name__ == "__main__":
    cat = read_brown_coordinates()
    cat = read_brown_apertures(cat=cat)

    pacs = glob.glob('../imaging/kingfish_pacs_scanam_v17/*fits')
    spire = glob.glob('../imaging/KINGFISH_SPIRE_v3.0_updated/*scan.fits')

    fields = ['ra','dec','height', 'width', 'PA', 'flag',
              'pacs70', 'pacs70_unc', 'pacs100', 'pacs100_unc',
              'pacs160', 'pacs160_unc', 'spire250', 'spire250_unc',
              'spire350', 'spire350_unc', 'spire500', 'spire500_unc']
    dt = [('name', 'S20')] + [(f, np.float) for f in fields]
        
    galaxy_names = cat.keys()
    # galaxy_names = ['ngc3190']
    fcat = np.zeros(len(galaxy_names), dtype=np.dtype(dt))
    for i, name in enumerate(galaxy_names):
        fcat[i]['name'] = name
        for k in cat[name]:
            try:
                fcat[i][k] = cat[name][k]
            except:
                pass
        thispacs = find_images(pacs, name)
        thisspire = find_images(spire, name)
        pflux, punc, pband = measure_pacs_flux(thispacs, cat[name])
        sflux, sunc, sband = measure_spire_flux(thisspire, cat[name])
        if (len(pflux) == 0) and (len(sflux) ==0):
            fcat[i]['flag'] = 1
        for f, u, b in zip(pflux+sflux, punc+sunc, pband+sband):
            fcat[i][b] = f
            fcat[i][b + '_unc'] = u


    h = pyfits.hdu.image.BinTableHDU(fcat)
    h.header['FlUXUNIT'] = 'Jy'
    h.header['PAUNIT'] = 'Degrees E of N'
    h.header['APUNIT'] = 'Arcseconds'
    pyfits.writeto('kingfish.brownapertures.flux.fits', fcat, h.header, clobber=True)


        
