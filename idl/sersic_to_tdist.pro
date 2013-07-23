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
                         int1=int1,int2=int2
 
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
  p=(1.0d)-0.6097/n+0.05563/n^2

;----------
;U profile
;--------------
  if keyword_set(radii) EQ 0 then radii=(findgen(20)*0.1+0.1)
  a=radii

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

  Uprofile=1.0d/a*(int1[*,0]+int2[*,0])
  return,Uprofile

end
