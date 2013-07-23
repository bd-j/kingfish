FUNCTION ufunc,x

common ufunc_pars,a_radius, n_sersic ;sets n and a
;a=0.1
;n=4.0
a=a_radius
n=n_sersic
b=get_b_sersic(n)
p=(1.0d)-0.6097/n+0.05563/n^2
power=2.0d*n-n*p-1

expression=(x/b)^(power)*exp(0-x)*alog((a+(x/b)^n)/abs(a-(x/b)^n))

return,expression

end



FUNCTION sersic_to_tdist,n,R_e,radii=radii,$
                         b=b,p=p,$
                         int1=int1,int2=int2,norm=norm,$
                         LofR=LofR,rofL=rofL,I_0=I_0,$
                         cgs=cgs,mmp83=mmp83,rho_0=rho_0

;;input units are solar luminosities and kpc.  output is L_sun/kpc^3 unless
;;\cgs or \mmp83 is set

  u_mmp83=2.17E-2 ;erg/cm^2/s
  conv1=10^(alog10(!lsun)-3.*alog10(!pc2cm*1E3) ) ;L_sun*s/kpc^3 to erg/cm^3
  conv2=conv1/(u_mmp83/(!lightspeed*1E-8)) ;from L_sun*s/kpc^3 to units of U_mmp/c
  conv=1
  if keyword_set(cgs) then conv=conv1
  if keyword_set(mmp83) then conv=conv2

;stop

  common ufunc_pars,a_radius, n_sersic
  n_sersic=n                    ;add to the common block

;sersic model
;I(R)=I_0*exp(-b*(R/R_e)^(1/n))

;------------
;;get b for this n
;----------
;\Gamma(2n)=2\igamma(2n,b)
;approximate
  ;b=2.0d*n-1./3.+0.009876/n
  b=get_b_sersic(n)

;----------
;luminosity density
;----------
  p=(1.0d)-0.6097/n+0.05563/n^2

;Lrho(r)=Lrho_0*(r/R_e)^(-p)*exp(-b*(r/R_e)^(1/n)
;Lrho_0=(L/L_lambda)*I_{0,\lambda}*b^(n*(1-p))*

;;---------------
;;normalize
;----------------
;L(r)=\int_0^r I(R)2\pi R dR=   ;;use x=b(R/R_e)^(1/n)
;L(r)=2.0d*!PI*I_0*R_e^2*n*b^(0-2*n)*\gamma(2n,(r/R_e)^(1/n))

  c=!lightspeed*1E-8/!pc2cm/1E3 ;kpc/s
  if n_elements(I_0) EQ 0 then $
     I_0=LofR/(2*!PI*R_e^2*n*b^(0.-2*n)*igamma(2.*n,b*(RofL/R_e)^(1./n))*Gamma(2.*n)) ;L_sun/kpc^2
  rho_0=I_0*b^(n*(1.0-p))*Gamma(2.0d*n)/(2.*R_e*gamma(n*(3.-p))) ;L_sun/kpc^3
  norm=rho_0*R_e*n/(2*c*b)*conv
;norm=1

;----------
;U profile
;--------------
  if keyword_set(radii) EQ 0 then radii=R_e*(findgen(20)*0.1+0.1)
  a=radii/R_e

;U(r')=\int_0^\infty rho(r) r/(2r'c) ln[(r'+r)/(|r'-r|)] dr
;a=r'/R_e
;u=b*(r/R_e)^(1/n)
;
;U(a*R_e)=rho_0*R_e*n/(2*a*c)*$
;          integral_0^\infty{
;(u/b)^(2*n-n*p-1)*exp(0-u)*alog((a+(u/b)^n)/abs(a-(u/b)^n))*du}

  eps=[1E-5,5E-4]
  neps=n_elements(eps)
  nr=n_elements(a)
  int1=fltarr(nr,neps)
  int2=fltarr(nr,neps)
  
  for iradius=0,nr-1 do begin
     a_radius=a[iradius] ;add to common block
     singularity=(b*a_radius^(1.0d/n))
     int1[iradius,*]=qromo('ufunc',0,singularity-5*eps*(a_radius/0.1));,jmax=50)
     int2[iradius,*]=qromo('ufunc',singularity+5*eps*(a_radius/0.1),/midexp) ;,jmax=50)
  endfor


Uprofile=norm/a*(int1[*,0]+int2[*,0])

;expression=(u/b)^(2*n-n*p-1)*exp(0-u)*alog((a[i]+(u/b)^n)/abs(a[i]-(u/b)^n))
;plot,u,expression,/xlog,/ylog,xrange=[0,20],yrange=[0,0.4]
;for i=0,19 do oplot,u,(u/b)^(2*n-n*p-1)*exp(0-u)*alog((a[i]+(u/b)^n)/abs(a[i]-(u/b)^n))

;stop

return,Uprofile

end
