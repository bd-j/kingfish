import numpy as np

def get_b_sersic(n):
    return 2.0*n-1/3.+(4./405.)/n+(46.0/25515.0)/(n**2)+(131.0d/1148175.0)/(n**3)

class Umodel(object):

    def __init__(self):
        self.load_bulge_models()
        self.load_disk_models()

    def profiles(self, info, radii = None, angular = True, convolved = False):

        if angular is True: ptype = 'obs' else : ptype = 'phys'
            
        scale_disk = info['H_disk'][ptype]
        scale_bulge = info['R_e'][ptype]
        
        u_bulge, u_disk = self._profiles(mu_e_phys = info['mu_e']['phys'],
                                         mu0_phys = info['mu0']['phys'],
                                         n_bulge = info['n_bulge'], k = 8)

        all_radii = np.concatenate([self.radii_bulge,self.radii_disk])
        if radii is None:
            radii = np.sort(all_radii)

        rr = self.radii_bulge.append(all_radii.max())
        uu = u_bulge.append(5e-4*u_bulge.min())
        log_u_bulge = np.interp(np.log10(rr*scale_bulge), np.log10(uu), np.log10(radii))

        rr = self.radii_disk
        uu = u_disk
        log_u_disk = np.interp(np.log10(rr*scale_disk), np.log10(uu), np.log10(radii))

        return 10**log_u_bulge, 10**log_u_disk

            
    def _profiles(self, mu_e_phys = 1e9, mu0_phys = 1e9, n_bulge = 4, k = 8):

        moduub, n_model =  self.moduub, self.n_model
        moduud, k_model =  self.moduub, self.k_model
        
        ind = np.searchsorted(n_model, n_bulge)#nearest(n_model,n,/below)
        weight = ( np.array([(n_bulge-n_model[ind]),(n_model[ind+1]-n_bulge)])/
                  (n_model[ind+1]-n_model[ind]) )
        u_bulge = (moduub[ind,:]*weight[0]+moduub[ind+1,:]*weight[1])*(mu_e_phys/1e9)
        
        #add a disk with k = k_sub
        ind = np.where(k_model == k)
        u_disk = moduud[ind,:]*(mu0_phys/1e9) 
        return u_bulge, u_disk
 
    def load_disk_models(self):
        
        kdir = '~/WRITING/KINGFISH/Uexample/plots/'
        self.radii_disk, self.moduud, self.k_model = self.load_model(kdir+'model_uofr/udist_disk.fits')

    def load_bulge_models(self):

        kdir = '~/WRITING/KINGFISH/Uexample/plots/'
        self.radii_bulge, self.moduub, self.n_model = self.load_model(kdir+'model_uofr/udist_sersic.fits')

    def load_model(self, fname):
        f = pyfits.open(fname)
        darray = f[0].data
        unit = f[0].header['BUNIT']
        crp = f[0].header['CRPIX1']
        crv = f[0].header['CRVAL1']
        cd = f[0].header['CDELT1']
        radius = (np.arange(darray.shape[0])-(crp-2))*cd+crv
        
        uu = darray[:,crp-1:]
        mpar = darray[:,crp-2]
        f.close()
        return radius, uu, mpar
    

def convolve_profile():
    pass
