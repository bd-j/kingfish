import numpy as np
import kingfish
import isrf
import utils
import matplotlib.pyplot as pl

rp = {'galaxy':'NGC0628',
      'maxrad': 380,
      'dataset':'S250_100_SSS_100',
      'utype':'UBAR',
      'ulim': 0.2,
      'snrcut':3,
      'bands':['']}

### Initialization
umodel = isrf.Umodel()

## this reads the dust parameter images, the convolved flux images,
## and some catalog information, and derives some basic quantities.
## these are all stored as attributes of the galaxy object
galaxy = kingfish.KingfishGalaxy(**rp)


### Pixel quantities
## define the good pixels and store them in the object
galaxy.goodpixels = np.where( (galaxy.dust['Mdust'] > 0.0) &
                              (galaxy.dust['Mdust']/galaxy.dust['Mdust_unc'] > rp['snrcut']) &
                              (np.isfinite(galaxy.dust['Mdust']/galaxy.dust['Mdust_unc'])) & 
                              (galaxy.dust[rp['utype']] > rp['ulim'] ) &
                             )
## get the dust parameters for the good pixels
pixel_udust = galaxy.dust[rp['utype']][galaxy.goodpixels]
pixel_mdust = galaxy.dust['Mdust'][galaxy.goodpixels]
## get the deprojected radii
pixel_radii = galaxy.pixel_deprojected_radii(galaxy.goodpixels)
## get the scaled model profiles for this galaxy at the pixl locations
pixel_ubulge, pixel_udisk = umodel.profiles(galaxy.info, radii = pixel_radii, angular = True)
pixel_utot = (pixel_ubulge*galaxy.spectrum['bulge']['eta']+
              pixel_udisk*galaxy.spectrum['disk']['eta'])

## get a smooth model profile, convolved with a gaussian PSF
model_radii = np.arange(rp['maxrad']-1) +1 #1 arcsec intervals up to maxrad, starting at 1 arcsecond
model_ubulge, model_udisk = umodel.profiles(galaxy.info, radii = model_radii, angular = True, convolved = True)
model_utot = (model_ubulge*galaxy.spectrum['bulge']['eta']+
              model_udisk*galaxy.spectrum['disk']['eta'])

### Do it all for another resolution
#galaxy.dataset = 'S350_110_SSS_110'
#galaxy.load_data()
#galaxy.goodpixels =
#pixel_radii = galaxy.pixel_deprojected_radii(galaxy.goodpixels)
#pixel_utot, pixel_ubulge, pixel_udisk = umodel.physical_profile(galaxy.info, radii = pixel_radii)

#galaxy.write_table()

pl.figure(1)
pl.plot(pixel_radii, pixel_udust, 'ro', alpha = 0.5)
pl.plot(model_radii, model_utot, '-k')
pl.plot(model_radii, model_ubulge*galaxy.spectrum['bulge']['eta'], ':k')
pl.plot(model_radii, model_udisk*galaxy.spectrum['disk']['eta'], ':k')
pl.savefig(figname1)

pl.figure(2)
pl.plot(pixel_radii,pixel_udust/pixel_utot, 'ro', alpha = 0.5)
pl.plot([0, rp['maxrad']], [1,1], 'k')
pl.plot([0, rp['maxrad']], [avg_ratio, avg_ratio], ':k')
pl.savefig(figname2)

#plotter.plot_uofr()


