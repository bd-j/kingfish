PRO read_imaging,gal,images,bands=bands,dataset=dataset,$
                 headers=headers,images_unc=images_unc,mask=mask

dir='~/SFR_FIELDS/Nearby/KINGFISH/data/'

if n_elements(bands) EQ 0 then $
   bands=['GALEX_FUV','GALEX_NUV','Optical_U','Optical_B','Optical_V','Optical_R','Optical_I','Optical_Ha','IRAC_3.6','IRAC_4.5','IRAC_5.8','IRAC_8.0','MIPS_24','PACSS_70','PACSS_100','PACSS_160','SPIRE_250','SPIRE_350']
;recal=[0,0,1,1,1,1,1,2,fltarr(10)]
nband=n_elements(bands)

subdir=gal+'/herschel/dustmaps/'+dataset+'/Convolved_Images/'
maskname=gal+'_'+dataset+'_Gal_mask.fits.gz'

mask=mrdfits(dir+subdir+maskname,/silent)
sz=size(mask)
images=fltarr(sz[1],sz[2],nband)
images_unc=images
headers=strarr(nband,4000)


for ib=0,nband-1 do begin
   imname=gal+'_'+dataset+'_'+bands[ib]+'.fits.gz'
   uncname=repstr(imname,'.fits.gz','_unc.fits.gz')

   if file_test(dir+subdir+imname) then begin
      tmp=mrdfits(dir+subdir+imname,0,hdr,/silent)
      nhl=n_elements(hdr)
      headers[ib,0:nhl-1]=hdr
      
      cal=where(strpos(hdr,'DN/pixel/s') GE 0,nc)
      if nc GT 0 then begin 
         cal=hdr[cal[0]-2]
         cal=strpos(cal,'1.0000E+00') GE 0
      endif
      if cal EQ 1 then begin
;         photflam=sxpar(hdr,'PHOTFLAM')
         obs=sxpar(hdr,'OBSERVAT')
         
         ps=1.0
         if strpos(obs,'KPNO') GE 0 then ps=0.305 else begin
            if strpos(obs,'CTIO') GE 0 then begin
               if strpos(sxpar(hdr,'TELESCOP'),'4.0') GE 0 then ps=0.305 else $;ps=sxpar(hdr,'PIXSCAL1') else $
                  ps=0.433
            endif
           
         endelse
;        print,'No pixscale info for '+imname
;        tmp=tmp*photflam/1E6/ps^2 ;MJy/arcsec^2 in input image
         tmp=tmp*206265.0d^2/ps^2 ;convert from DN/pixel/s to DN/sr/s
print,gal+' ',bands[ib],' '+obs,ps,' '+sxpar(hdr,'TELESCOP')
      endif
      images[*,*,ib]=tmp
   endif else print,"couldn't find "+dir+subdir+imname

   if file_test(dir+subdir+uncname) then begin
      tmp=(mrdfits(dir+subdir+uncname,/silent))[*,*,3]
  ;    if recal[ii] EQ 1 then tmp=tmp*photflam/1E6/ps^2*206265.0d^2
      images_unc[*,*,ib]=tmp
   endif
endfor


end
