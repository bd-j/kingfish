FUNCTION convolve_profile,radius, value, fwhm

maxr = max(radius)+5*fwhm

ss = findgen(floor(maxr))+1

linterp, alog(radius), alog(value), alog(ss), lognewvalue

nv = exp(lognewvalue)
rr = [0-reverse(ss), ss]
oo = [reverse(nv), nv]

out = gaus_convol(rr, oo, fwhm/2.35)
good = where(rr GT 0)
return, [[rr[good]],[out[good]]]

end
