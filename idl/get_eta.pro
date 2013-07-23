;calculate the ratio of the fractional absorbed flux in the MMP83 ISRF
;to the fractional absorbed flux in an arbitrary spectrum

FUNCTION get_eta,wave,flux,lambda_max=lambda_max,wave_range=wave_range

sz=size(flux)
nmod=sz[2]

;;----------------
;;GET DUST KAPPA
kappafile='$SPECFIT_DIR/data/kext_albedo_WD_MW_3.1_60_D03.all.dat'
readcol,kappafile,lambda,albedo,cos,ext,kappa,cos2,/silent
lambda=reverse(lambda*1E4)
oo=sort(lambda)
lambda=lambda[oo]
kappa=(reverse(kappa))[oo]
linterp,lambda,kappa,wave,new_kappa,missing=0
kappa=new_kappa/max(new_kappa)

;;-------------
;; GET ISRF
;;------------
isrffile='$SPECFIT_DIR/data/ISRF_MATHIS.DAT'
readcol,isrffile,lambda,isrf_nu,/silent
lambda=lambda*1E4
isrf_nu=isrf_nu/max(isrf_nu)
linterp,lambda,isrf_nu/lambda^2,wave,isrf
isrf_lambda=isrf/max(isrf)

;;-----------------
;;CALCULATE RATIO
;;-----------------
if keyword_set(wave_range) EQ 0 then wave_range=[1E3,8E4]
region=where(wave GT wave_range[0] and wave LT wave_range[1],nr)

u_mmp83=tsum(wave[region],isrf_lambda[region])
u_abs_mmp83=tsum(wave[region],isrf_lambda[region]*kappa[region])

u_star=tsum_2d(wave[region],flux[region,*])
u_abs_star=tsum_2d(wave[region],flux[region,*]*(kappa[region]#(fltarr(nmod)+1)))
j1=u_abs_star/u_star
j2=u_abs_mmp83/u_mmp83
;we=k_lambda_to_edges(wave[region])
;dlambda=we[1:nr]-we[0:nr-1]
jj=max(flux[region,*]*((kappa[region]*(wave[region]))#(fltarr(nmod)+1)),sub_max,dimension=1)
lambda_max=(wave[region]#(fltarr(nmod)+1))[sub_max]
;stop

return,j1/j2

end
