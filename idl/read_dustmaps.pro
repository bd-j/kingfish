;;procedure to read in data for a given galaxy


PRO READ_DUSTMAPS,gal,dust_map,dust_header,dataset=dataset,dust_unc=dust_unc,$
                  dustpar=dustpar,unc_dustpar=unc_dustpar,partype=partype,$
                  chisq=chisq,dopars=dopars,results=results,mask=mask


  g2=strupcase(gal)
  g2=repstr(g2,'HOL','Hol')

  dir='~/SFR_FIELDS/Nearby/KINGFISH/data/'
  subdir=gal+'/herschel/dustmaps/'+dataset+'/Results/'
  imname=g2+'_'+dataset+'_Model_All_PerPixel_Mdust.fits.gz'
  uncname=g2+'_'+dataset+'_Model_All_PerPixel_Mdust_unc.fits.gz'
  dust_map=mrdfits(dir+subdir+imname,0,dust_header,/silent)
  dust_unc=mrdfits(dir+subdir+uncname,0,/silent)
  mask=mrdfits(dir+subdir+'../Convolved_Images/'+g2+'_'+dataset+'_Gal_mask.fits.gz',0,/silent)

  if keyword_set(dopars) then begin
     sz=size(dust_map)
     if keyword_set(partype) EQ 0 then partype='All_' ;'All_' or 'Gal_', latter has a mask applied
     par=partype+['PerPixel_Mdust','PerPixel_Ldust','U_min',$
                  'f_PDR','PerPixel_LPDR','Gamaa','U_bar','qpah','Alpha']
     par_unc=par+'_unc'
     npar=n_elements(par)
     
     dustpar=fltarr(sz[1],sz[2],npar)
     unc_dustpar=dustpar
     imroot=g2+'_'+dataset+'_Model_'
     for ipar=0,npar-1 do begin
        dustpar[*,*,ipar]=mrdfits(dir+subdir+imroot+par[ipar]+$
                                  '.fits.gz',/silent)
        unc_dustpar[*,*,ipar]=mrdfits(dir+subdir+imroot+par_unc[ipar]+$
                                      '.fits.gz',/silent)
     endfor
     chisq=mrdfits(dir+subdir+imroot+partype+'chi_sq.fits.gz',/silent)
  endif

  results=read_aniano(gal,dataset)

end
