FUNCTION get_b_sersic,n,fisher=fisher

n=n*1.0d
b=2.0*n-1/3.+(4./405.)/n+(46.0d/25515.0d)/(n^2)+(131.0d/1148175.0d)/(n^3)

if keyword_set(fisher) then $
   b=2.17*n-0.355

return,b

end
