FUNCTION read_aniano,gal,dataset,root=root

npar=12

g2=strupcase(gal)
g2=repstr(g2,'HOL','Hol')


if keyword_set(root) EQ 0 then root='~/SFR_FIELDS/Nearby/KINGFISH/data/'
j=file_search(root+gal+'/herschel/dustmaps/'+dataset+'/'+g2+'_'+dataset+'_results.dat',count=count)

nmethod=2*2
nphot=35
as={name:'',dataset:'',date:'',dist:0.,ra:0.,dec:0.,global_phot:fltarr(nphot),global_phot_err:fltarr(nphot),$
    Temp:fltarr(nmethod),beta:fltarr(nmethod),M_dust_bb:fltarr(nmethod),$
    m_dust:fltarr(nmethod),L_dust:fltarr(nmethod),L_star:fltarr(nmethod),$
    L_pdr:fltarr(nmethod),q_pah:fltarr(nmethod),f_cold:fltarr(nmethod),$
    f_pdr:fltarr(nmethod),u_min:fltarr(nmethod),u_bar:fltarr(nmethod),$
    u_cold:fltarr(nmethod),alpha:fltarr(nmethod),gamma:fltarr(nmethod),$
    chi2:fltarr(nmethod)}


;stop
as.name=gal
as.dataset=dataset

if count EQ 0 then begin
   print, 'Error finding file: '+gal+'_'+dataset
   return,as
endif


nl=(file_lines(j))[0]
get_lun,lun
;print,lun
openr,lun,j

c=strarr(nl)
d=''
for i=0,nl-1 do begin
   readf,lun,d,format='(A)'
   c[i]=d 
endfor

c=strcompress(c,/remove_all)
c=repstr(c,'chi^2','chi2')

;;BASIC DATA
sub=where(strpos(c,'AssumedDistance') GE 0,count)
if count GT 0 then begin
   dd=(strsplit(c[sub],':',/extract))
   ;p=strpos(dd[1],'Mpc')
   as.dist=float(dd[1])
endif

sub=where(strpos(c,'generated') GE 0,count)
if count GT 0 then begin
   dd=c[sub]
   p=strpos(dd,'201')
   date=strmid(dd,p,10)
   as.date=date
endif

;;GLOBAL PHOTOMETRY
sub=where(strpos(c,'Globalphotometry') GE 0, count)
if count GT 0 then begin
   for k=0,nphot-1 do begin
      valstring=strsplit(c[sub+k+1],'[()]',/extract,/regex)
      if n_elements(valstring) GT 1 then begin
         dd=strsplit(valstring[1],'(\+/-)',/extract,/regex , count=count)
         as.global_phot[k]=float(dd[0])
         if count GT 1 then as.global_phot_err[k]=float(dd[1])

      endif
   endfor
endif

;;DERIVED PARAMETERS, from single pixel and sum over pixels.  single
;;pixel result goes in the first element, sum in the second element
si=(where(strpos(c,'onebigpixel') GE 0,count))[0]
if count GT 1 then print,'problem'
sum=(where(strpos(c,'Sumovertheindividual') GE 0,count))[0]
if count GT 1 then print,'problem'
sbb=(where(strpos(c,'ModifiedBlackBody') GE 0,hasbb))
if count LT 2 then start=11 else start=8
tags=tag_names(as)
ntag=n_elements(tags)

for i=8,ntag-1 do begin
   hack=1 & hack2=0
   if tags[i] EQ 'M_DUST' and hasbb GT 0 then hack=2 ;deal with the fact that M_dust parameter is repeated for bb fits
   if tags[i] EQ 'M_DUST_BB' then begin
      sub=where(strpos(strupcase(c),'M_DUST') GE 0) 
      hack2=1 ;deal with the fact that M_dust parameter is repeated for bb fits
   endif else $
      sub=where(strpos(strupcase(c),tags[i]) GE 0)
   ss=where(sub GT si and sub LT sum,fc)

   if fc EQ 0 then continue

;stop


   for k=0,1 do begin
      valstring=strsplit(c[sub[ss[0]+k*hack+hack2]],'[()]',/extract,/regex)
      if n_elements(valstring) GT 1 then begin
         dd=strsplit(valstring[1],'(\+/-)',/extract,/regex , count=count)
         ddval=float(dd[0])
         as.(i)[k]=ddval
         if count GT 1 then begin
            ddunc=float(dd[1])
            as.(i)[2+k]=ddunc
         endif
      endif else begin
         valstring=(strsplit(c[sub[ss+k*hack+hack2]],'=',/extract,/regex))[1]
         ddval=float((stregex(valstring,'[0-9.]*',/extract))[0])
         as.(i)[k]=ddval
      endelse
   endfor
endfor

free_lun,lun

return,as
;mwrfits,as,root+gal+'_'+dataset+'.strct.fits',/create


end
