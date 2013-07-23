import os
import numpy as np
import pyfits, pywcs, pyregion
import observate
import utils

class AnianoResults(object):

    kdir = os.path.expanduser('~/SFR_FIELDS/Nearby/KINGFISH/')
    datadir = self.kdir+'data/'
    

    partype = 'All'
    parnames = ['Mdust','Ldust','LPDR','U_min','f_PDR','Gamaa','U_bar','qpah','Alpha']
    imprefix = ['PerPixel_', 'PerPixel_', 'PerPixel_', '','','','','','']
    image_bands = ['GALEX_FUV','GALEX_NUV',
                   'Optical_U','Optical_B','Optical_V','Optical_R','Optical_I',
                   'Optical_Ha','IRAC_3.6','IRAC_4.5','IRAC_5.8','IRAC_8.0',
                   'MIPS_24','MIPS_160',
                   'PACSS_70','PACSS_100','PACSS_160',
                   'SPIRE_250','SPIRE_350', 'SPIRE_500']
    
    def __init__(self,galaxy = None, dataset = None, bands = None, **kwargs):
        self.galaxy = galaxy
        self.aname, self.kname = self.sanitize_name(galaxy)
        self.dataset = dataset

        self.read_dustmaps()
        self.read_imaging(bands = bands)
        

    def read_dustmaps(self):
        self.dust = {}
        imname, uncname, maskname = self.dust_map_names()
        self.dust_header = pyfits.getheader(imname[0])
        for i, par in enumerate(self.parnames):
            im = pyfits.open(imname[i])
            unc = pyfits.open(uncname[i])
            self.dust[par] = im[0].data
            self.dust[par+'_unc'] = unc[0].data
            im.close()
            unc.close()
        mask = pyfits.open(maskname)
        self.dustpars['mask'] = mask[0].data

    def read_imaging(self, bands = None):
        if bands is None : bands = self.image_bands
        self.image = {}
        self.image_headers = {}
        for i, b in enumerate(bands):
            imname = 
            im = pyfits.open(imname)
            self.image[b] = im[0].data
            self.image_header[b] = im[0].header
            ##!!!!!!!!!FIX THIS!!!!!!!!
            recal = False  
            if recal is True:
                pixel_size = 1.0 #in arcsec, depending on telescope
                self.image[b] =self.image[b]*(206265.0/pixel_size)**2
            unc = pyfits.open(imname.replace('.fits.gz','_unc.fits.gz'))
            self.image[b+'_unc'] = unc[0].data
            unc.close()
            im.close()
            
    def read_aniano_results(self):
        pass

    def dust_map_names(self):
        dirname = self.datadir+self.aname.lower()+'/herschel/dustmaps/'+self.dataset+'/'
        subdir = +'/Results/'
        imname = [dirname+'Results/' + self.aname+'_'+self.dataset+'_Model_'+self.partype+
                  '{0}{1}.fits.gz'.format(self.imprefix[i], self.parnames[i])
                  for i in xrange(len(self.parnames))]
        uncname = imname.replace('.fits.gz','_unc.fits.gz')
        maskname = dirname+'Convolved_Images/'+self.aname+'_'+self.dataset+'_Gal_mask.fits.gz'
        return imname, uncname, maskname
        
    def self.sanitize_name(self, galaxy):
        kname = galaxy.upper()
        #kname = kname.replace('N(0-9)', 'NGC') #blerg
        kname = kname.replace('HOL1','HOLI')
        kname = kname.replace('HOL2','HOLII')
        kname = kname.replace('M81DB','M81dB')
        aname = name.replace('HOLII','Hol2')
        aname = aname.replace('HOLI','Hol1')
        aname = aname.replace('M81DB','M81dB')
        
        return aname, kname

class KingfishGalaxy(AnianoResults):
    
    infodir = self.kdir+'info/'
    info = {}

    def __init__(self,galaxy = None, dataset = None, bands = None, **kwargs):
        self.galaxy = galaxy
        self.aname, self.kname = self.sanitize_name(galaxy)
        self.dataset = dataset

        self.load_data(bands = bands)

    def load_data(self,bands = None):
        self.read_dustmaps()
        self.read_imaging(bands = bands)
        self.set_basic_info()
        self.load_stellar_population()
        self.dust['Radius'] = np.zeros(self.dust['Mdust'].shape)
        
    def set_basic_info(self):
        """compile and rename parameters from several catalogs"""
        
        q0 = 0.2
        ks = pyfits.open(infodir+'kingfish_sample.fits')[1].data
        dk = pyfits.open(infodir+'dalephot_kingfish.fits')[1].data
        fkh = pyfits.open(infodir+'fisher_Hband_kingfish.fits')[1].data
        sk = pyfits.open(infodir+'skibba_kingfish.fits')[1].data
        ak = pyfits.open(infodir+'aniano_kingfish_S250_100_SSS_100_vOct2012.fits')[1].data
        reg = pyregion.open(self.kdir+'photometry/dale_region_files'+self.kname+'.deg.reg')
        
        ind = (np.where(ks['NAME'] == self.kname))[0]

        self.info['ra'] = ks['RA'][ind]
        self.info['dec'] = ks['DEC'][ind]
        self.info['dist'] = ak['DIST'][ind]
        self.arcsec_to_kpc = self.info['dist']*1e3/206265.0
        self.info['a'] = ks['SEMIMAJOR'][ind]
        self.info['b'] = ks['SEMIMINOR'][ind]
        self.info['pa'] = reg[0].coord_list[2]
        q = self.info['b']/self.info['a']
        self.info['inclination'] = np.arccos(np.sqrt( (q^2- q0^2)/(1 - q0^2)))*(180./np.pi) + 3.0
        self.info['R_e'] = {'obs':fkh['R_EFF_BULGE'][ind]} #half-light radius in arcsec
        self.info['R_e'] = {'phys':fkh['R_EFF_BULGE'][ind]*self.arcsec_to_kpc} #half-light radius in kpc
        self.info['mu_e'] = {'obs':fkh['MU_EFF'][ind]} #in AB mag/arcsec**2 at H band
        self.info['n_bulge'] = fkh['N_BULGE'][ind] #sersic
        self.info['H_disk'] = {'obs':fkh['R_DISK'][ind]} #exponential scale length in arcsec
        self.info['H_disk'] = {'phys':fkh['R_DISK'][ind]*self.arcsec_to_kpc} #exponential scale length in arcsec
        self.info['mu0'] = {'obs':fkh['MU_0_DISK'][ind]} #in AB mag/arcsec**2 at H band
        self.info['morph_band'] = 'twomass_H'
        self.info['type'] = sk['MORPH'][ind]
        #self.info['units'] = 'observed (", mag, etc)'


    def load_stellar_population(self, complist = None):
        """read all the stellar pops info from Magphys outputs"""

        mpdir = os.expanduser('~/WRITING/KINGFISH/Uexample/plots/eta/magphys_out/')
        if complist is None:
            complist = ['tot','bulge','disk']
        self.spectrum = {}
        for i,c in enumerate(complist):
            spec = utils.read_mp(galaxy = self.aname, component = c, mpdir = mpdir)
            lbol = observate.Lbol(spec['wavelength'], spec['f_lambda_int'], wave_max = 1e5)
            self.spectrum[c] = spec
            self.spectrum[c]['eta'] =  self.get_eta(spec['wavelength'], spec['f_lambda_int'])
            self.spectrum[c]['log_Lbol']  = np.log10(observate.Lbol(spec['wavelength'],
                                                                    spec['f_lambda_int'],wave_max = 1e5))

        self.info['mu_e']['phys']=
        
    def bolometric_correction(self, component = 'tot', filternamelist = None):
        lightspeed = 2.998e18 #AA/s
        filterlist = observate.load_filters(filternamelist)
        mags = observate.getSED(self.spectrum[component]['wavelength'],
                                self.spectrum[component]['f_lambda_int'],
                                filterlist)
        nueff = np.array([lightspeed/f.wave_effective for f in filterlist])
        
        log_bc = self.spectrum[component]['log_Lbol'] - log_nuLnu

        return bc

    def get_eta(self, wavelength, f_lambda):
        kappaname='$SPECFIT_DIR/data/kext_albedo_WD_MW_3.1_60_D03.all.dat'
        isrfname = '$SPECFIT_DIR/data/ISRF_MATHIS.DAT'
        kappa = utils.read_kappa(kappaname, wavelength)
        isrf = utils.read_isrf(isrfname, wavelength)
        u_mmp83 = observate.Lbol(wavelength,isrf, wave_min = 1e3, wave_max = 1e5)
        u_abs_mmp83 = observate.Lbol(wavelength,isrf*kappa, wave_min = 1e3, wave_max = 1e5)
        u_star = observate.Lbol(wavelength,f_lambda, wave_min = 1e3, wave_max = 1e5)
        u_abs_star = observate.Lbol(wavelength,f_lambda*kappa, wave_min = 1e3, wave_max = 1e5)

        return u_abs_star/u_star*u_mmp83/u_abs_mmp83

    def pixel_deprojected_radii(self, pixels = None):
        wcs = pywcs.WCS(header = self.dust_header)
        pcoords = wcs.wcs_pix2sky(pixels, origin = 0, ra_dec_order = True)
        radii, phi = utils.deproject_coords(center = (self.info['ra'], self.info['dec']),
                                            ra = pcoords[:,0], dec = pcoords[:,1],
                                            **self.info)


        
        return radii


        
