u_mmp83=2.17E-2                                ;erg/cm^2/s
conv1=10^(alog10(!lsun)-2.*alog10(!pc2cm*1E3)) ;L_sun/kpc^2 to erg/s/cm^2

n=findgen(28)*0.2+0.6
nn=n_elements(n)

I_e=1D9
R_e=1

nr=60
uu=fltarr(nn,nr)
for i=0,nn-1 do begin
   b=get_b_sersic(n[i])
   p=(1.0d)-0.6097/n[i]+0.05563/n[i]^2
   norm=I_e/n[i]*b^(n[i]*(3-p))/gamma(n[i]*(3-p))/4.0 ;rho_0 * R_e/2
   norm=norm*conv1/u_mmp83
   radii=(findgen(nr)*0.1+0.10)*R_e
   uu[i,*]=norm*sersic_to_tdist(n[i],R_e,radii=radii,int1=int1,int2=int2)
endfor

array=[[n],[uu]]
mwrfits,array,'udist_sersic.fits',/create
array=mrdfits('udist_sersic.fits',0,hdr)
sxaddpar,hdr,'BUNIT','(U/u_mmp83) (10^9L_sun kpc^-2/\mu_e)'
sxaddpar,hdr,'CRPIX1',2
sxaddpar,hdr,'CRVAL1',0.10,'r/R_e'
sxaddpar,hdr,'CDELT1',0.1
mwrfits,array,'udist_sersic.fits',hdr,/create

;to get the correct output (in units of u_mmp83)
;  1)determine mu_e (from I_0 or from luminosity within some radius)
;  2)interpolate array in n
;  3)multiply array by (I_e)

end



