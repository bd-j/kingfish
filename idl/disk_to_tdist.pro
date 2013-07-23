FUNCTION ufunc_disk,x

  ;eps=[1E-6,5E-6,1E-5,5E-5,1E-4,5E-4,1E-3,5E-3]
  ;;sets ratio of scale length to scale height, a, and x, and chooses the
  ;;accuracy to return
  common ufunc_pars, a_radius,k_disk,x_integral,epsilon_integral 
  a=a_radius
  k=k_disk
  x_integral=x
  eps=epsilon_integral

  ;;the returned expression should be the integral over y at the given x
  singularity=sqrt((a-x)^2.0)
  int1=qromo('ufunc_disk_fixedx',0,singularity-eps,k=8)
  int2=qromo('ufunc_disk_fixedx',singularity+eps,/midexp,k=8)

  return,int1+int2
end

FUNCTION ufunc_disk_fixedx,y

  ;;sets ratio of scale length to scale height, a, and x
  common ufunc_pars, a_radius,k_disk,x_integral,epsilon_integral 
  a=a_radius
  k=k_disk
  x=x_integral

  denom=sqrt((a-x)^2+y^2)*sqrt((a+x)^2+y^2)
  expression=x*exp(0-x)*exp(0-k*y)/denom

  return,expression
end


FUNCTION disk_to_tdist,k,H,radii=radii,eps=eps
  ;eps=[1E-6,5E-6,1E-5,5E-5,1E-4,5E-4]
  common ufunc_pars,a_radius,k_disk,x_integral,epsilon_integral
  k_disk=k*1.0d                 ;add to the common block
  if n_elements(eps) GT 0 then $
     epsilon_integral=eps else $
        epsilon_integral=1E-5

  a=radii/H*1.0d
  nr=n_elements(a)

  int=fltarr(nr)
  for iradius=0,nr-1 do begin
     a_radius=a[iradius]        ;add to common block
     
     int[iradius]=qromo('ufunc_disk',0,/midexp,k=8)
  endfor


return,int

end
