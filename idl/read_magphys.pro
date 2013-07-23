;+
; NAME:
;   READ_MAGPHYS
;
; VERSION:
;   1.0 (Nov, 2011)
;
; PURPOSE:
;   to read the outputs of the sed fitting code MAGPHYS (daCunha et
;   al. 2008) and put them in a structure
;
; REFERENCE:
;   Da Cunha et al. 2008 ApJ
;
; CATEGORY:
;   File reading
;
; CALLING SEQUENCE:
;   structure=read_magphys(name[,DIR=DIR])
;
; INPUTS:
;    name - the root name (i.e. without file extensions) of the
;           magphys output file (string)
;
; KEYWORD PARAMETERS:
;    dir - string containing the path to the directory where the
;          magphys output is located

;
; OUTPUT:
;    A 1-element structure containing the results of the sed fit,
;    including the input filters and fluxes, pdfs for each returned
;    parameter, and the spectrum of the best fit
;
;
; COMMENTS:
;    See the MAGPHYS documentations for description of the returned
;    parameters, their PDFs, percentiles of the PDF, etc..  Right now
;    the location of data in the magphys out files is hard coded; if
;    the magphys output format changes this code will break.
;                                                                  
;
; REVISION HISTORY:
;    Nov 2011 - written, B. Johnson (IAP)
;
;--------------------------------------------------------------------



FUNCTION read_magphys,name,dir=dir

gal=name
if keyword_set(dir) EQ 0 then dir=''
;dir=''
j=file_search(dir+gal+'.fit',count=count)

nl=(file_lines(j))[0]
if count EQ 0 then begin
   print, 'Error finding file: '+gal
   return,-99
endif
c=strarr(nl)
d=''
get_lun,lun
openr,lun,j
for i=0, nl-1 do begin 
   readf,lun,d,format='(A)'
   c[i]=d
endfor
close,lun
free_lun,lun

readfast,repstr(j[0],'.fit','.sed'),sed,skipline=10
nw=(size(sed))[2]

parnames=['f_mu_sfh','f_mu_ir','mu','tau_v','ssfr','m_star','l_dust','t_ism','t_bc','xi_cold','xi_pah','xi_mir','xi_warm','mutau_v','m_dust']
npar=n_elements(parnames)
;could do this also by string searching for parnames, in case output format changes
parline=[17,40,63,86,137,210,273,336,349,382,405,428,451,474,557]-1 
npdf_par=[(parline[1:*]-parline[0:npar-2])-3,60.] ;number of elements in the pdf for each parameter
npdf_par_string=strcompress(string(fix(npdf_par),format='(I3.0)'),/remove_all) ;convert to string

;pull out some data from the string array
filters=(strsplit(c[1],' ',/extract))[1:*]
nfilt=n_elements(filters)
obs=float(strsplit(c[2],' ',/extract))
unc=float(strsplit(c[3],' ',/extract))
model=float(strsplit(c[12],' ',/extract))
bestfit=float(strsplit(c[10],' ',/extract))
info=float(strsplit(c[8],' ',/extract))

;;CREATE STRUCTURE

mp={name:'',date:'',redshift:0.,filters:strarr(nfilt),$
    obs_flux:fltarr(nfilt),obs_flux_unc:fltarr(nfilt),model_flux_bestfit:fltarr(nfilt),$
    spec_wave:fltarr(nw),spec_flux:fltarr(nw),spec_flux_intrinsic:fltarr(nw),$
    chibest:0.,optical_index:0.,ir_index:0.}


;add percentile tags
start=n_tags(mp)
mp=struct_addtags(mp,[parnames+'_median',parnames+'_bestfit',$
                      parnames+'_2p5',parnames+'_16',parnames+'_84',$
                      parnames+'_97p5'],strarr(npar*6)+'0.')
;add pdf tags
pdfstart=n_tags(mp)
mp=struct_addtags(mp,[parnames+'_pdf'],$
                  [strarr(npar)+'fltarr(2,'+npdf_par_string+')'])


;;FILL THE STRUCTURE

mp.name=gal
mp.filters=filters
mp.obs_flux=obs & mp.obs_flux_unc=unc
mp.model_flux_bestfit=model
mp.redshift=info[3] & mp.chibest=info[2]
mp.optical_index=info[0] & mp.ir_index=info[1]

for ipar=0,npar-1 do begin

   percentiles=float(strsplit(c[parline[ipar]+npdf_par[ipar]+1],' ',/extract))
   mp.(start+ipar)=percentiles[2] ;median
   mp.(start+npar+ipar)=bestfit[ipar] ;besfit
   mp.(start+npar*2+ipar)=percentiles[0] ;2.5
   mp.(start+npar*3+ipar)=percentiles[1] ;16
   mp.(start+npar*4+ipar)=percentiles[3] ;84
   mp.(start+npar*5+ipar)=percentiles[4] ;97.5

   pdf=fltarr(2,npdf_par[ipar])
   for j=0,npdf_par[ipar]-1 do $
      pdf[*,j]=float(strsplit(c[parline[ipar]+j],' ',/extract))
   mp.(pdfstart+ipar)=pdf
endfor

mp.spec_wave=reform(sed[0,*])
mp.spec_flux=reform(sed[1,*])
mp.spec_flux_intrinsic=reform(sed[2,*])

return,mp


end
