linename =['hb', 'ha', 'oiii_88', 'oiii_5007', 'oiii_4959','sii_','sii_']
atomname = ['H', 'H', 'OIII', 'OIII', 'SII', 'SII']
lineflux = []
line_unc = []

atom = [getattr(fivel, a) for a in atomname]

attenuator = dustcurves.MilkyWay()

def prior(abundance = abundance, n = 100, Te = 1e4, Av = 0, Rv = 3.1, apcorr = 1):
    pass

def model(abundance = abundance, n = 100, Te = 1e4, Av = 0, Rv = 3.1, apcorr = 1):
    #get fluxes from fivel and hummer and story
    for i,line in enumerate(linename):
        modelflux[i] = abundance[i]*atom[i].get_emissivity(n,Te,line)
    #attenuate by dust
    tau = attenuator.acurve(linewave, A_v = A_v, R_v = R_v)/1.086
    modelflux = modelflux*np.exp(-tau)
    #add an aperture correction to ppak derived lines
    modelflux[ppak] = modelflux[ppak]*apcorr

    return modelflux

def lnprob(theta, obsflux = obsflux, prior = prior):
    
    Ne = theta[0]
    NO = theta[1]
    NS = theta[2]
    abundance = [Ne, Ne, NO, NO, NO, NS, NS]
    n = theta[3]
    Te = theta[4]
    Av = theta[5]
    Rv = theta[6]
    apcorr = theta[7]

    lnp_prior = prior(abundance = abundance, n = n, Te = Te, Av = Av, Rv = Rv, apcorr = apcorr)
    mflux = model(abundance = abundance, n = n, Te = Te, Av = Av, Rv = Rv, apcorr = apcorr)
    lnp = -((mflux-obsflux)/obsunc).sum()/2
    #lnp = lnp + lnpprior
    
    return lnp
