
u_mmp83=2.17E-2                                ;erg/cm^2/s
conv1=10^(alog10(!lsun)-2.*alog10(!pc2cm*1E3)) ;L_sun/kpc^2 to erg/s/cm^2
k=[5,8,10,12,16]

nk=n_elements(k)

H=1 ;kpc
I_0=1D9 ;L_sun/kpc^2

nr=60
uu=fltarr(nk,nr)
for i=0,nk-1 do begin
   radii=(findgen(nr)*0.1+0.1)*H
   uu[i,*]=disk_to_tdist_fixedy(k[i],H,radii=radii,eps=1E-5)
   norm=k[i]*I_0*conv1/u_mmp83
   uu[i,*]=uu[i,*]*norm;*(1E10/rho_0)
   aa=radii/H
endfor

array=[[k],[uu]]
mwrfits,array,'udist_disk.fits',/create
array=mrdfits('udist_disk.fits',0,hdr)
sxaddpar,hdr,'BUNIT','(U/u_mmp83) (10^9L_sun kpc^-2/\mu_0)'
sxaddpar,hdr,'CRPIX1',2
sxaddpar,hdr,'CRVAL1',0.10,'r/H'
sxaddpar,hdr,'CDELT1',0.1
mwrfits,array,'udist_disk.fits',hdr,/create

;to get the correct output (in units of u_mmp83)
;  1)determine mu_e (from I_0 or from luminosity within some radius)
;  2)interpolate array in n
;  3)multiply array by (mu_e)

end



