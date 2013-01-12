pro read_run_all,file=file,lat0=lat0,lon0=lon0,$
	conc_all,emis_all,drydep_all,wetdep_all,chem_all,$
	longitude,latitude,convert,prof_all,gridinfo,transp_all,hg2frac=hg2frac

; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;  This program reads in concentrations, deposition, and emissions
;  for the htap Hg runs
; 
;  add ocean parameters as well.
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

;  use the correct diaginfo/tracerinfo files
; ctm_diaginfo,diagn,diagstru, filename='~/4-geos/1-GEOS-Chem.v8-02-03.yanxu/run.05x0666/dataout/diaginfo.dat'
; ctm_tracerinfo,tracer,tracerstru, filename='~/4-geos/1-GEOS-Chem.v8-02-03.yanxu/run.05x0666/dataout/tracerinfo.dat'

 ctm_diaginfo,diagn,diagstru, filename='./diaginfo.dat'
 ctm_tracerinfo,tracer,tracerstru,filename='./tracerinfo.dat'

ctm_get_data,datainfo, 'ij-avg-$',tracer=1,file=file
idate=tau2yymmdd(datainfo.tau0)
idate=idate.year*10000L+idate.month*100L+idate.day

;if not keyword_set(year) then year=2004
;idate=long(year)*10000L+1+(1L+indgen(12))*100L
imonth=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
nmonth=n_elements(idate)

;  Avogadro's number
Avo = 6.023d23
MolWt = 201d ;  molecular weight of Hg (g/mole)

;  read mixed layer depth
; if datainfo(0).dim(0) eq 144 then $
; 	file_mld='/data/ctm/GEOS_2x2.5/mercury_200501/mld.geos.2x25' else $
; 	file_mld='/data/ctm/GEOS_4x5/mercury_200501/mld.geos.4x5' 
    file_mld='/home/disk/oreo/nobackup/data/ctm/GEOS_05x0666/mercury_201007/MLD_DReqDT.geos5.05x0666' 

 for k=0,12-1 do begin
  success=ctm_get_datablock(z,'BXHGHT-$',tracer=5,file=file_mld,$
  	tau0=nymd2tau(ymd2date(1985,k+1,1)),lat=lat0,lon=lon0,non=non)
  if k eq 0 then mld=replicate({z:z*0.},12)
 ;  ocean mixed layer depth in meters
  mld(k).z=z
 endfor


;  Loop over months
; file2='~yanxuz/4-geos/Run.05x0666/dataout/ctm.bpch.08-16-10.1.2nd.part'
for k=0,nmonth-1 do begin
 date=idate(k)

 tau0=nymd2tau(date)
 month=imonth(where(date eq idate))

;  if tau0 ge 183336.0 then file=file2

 print,'Reading file...',file

 ;  extract Hg emissions
 success=ctm_get_datablock(anthr,'HG-SRCE',tracer=1,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,y=latitude,x=longitude,$
	thisdatainfo=datainfo,gridinfo=gridinfo,non=non)
 success=ctm_get_datablock(ocean,'HG-SRCE',tracer=3,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(land,'HG-SRCE',tracer=4,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(land_natural,'HG-SRCE',tracer=5,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(anthr_rgm,'HG-SRCE',tracer=6,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(anthr_hgp,'HG-SRCE',tracer=9,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(bb,'HG-SRCE',tracer=13,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(land_veg,'HG-SRCE',tracer=14,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(land_soil,'HG-SRCE',tracer=15,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)

;  add ocean aqueous mass of Hg0, Hg2 and HgC (particulate Hg)
 success=ctm_get_datablock(hg0aq,'HG-SRCE',tracer=2,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(hg2aq,'HG-SRCE',tracer=7,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)
 success=ctm_get_datablock(hgcaq,'HG-SRCE',tracer=11,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)

;  add Hgnr sink in ocean (particulate sinking)
 success=ctm_get_datablock(hgsink,'HG-SRCE',tracer=8,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,lev=level,non=non)

 ;  gross ocean flux out
  success=ctm_get_datablock(up,'HG-SRCE',tracer=16,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,lev=level,non=non)
  if success eq 0 then begin
     print,'UP not found...'
     up=0.
  endif

 ;  dry dep of Hg0 to ocean
  success=ctm_get_datablock(down,'HG-SRCE',tracer=17,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,lev=level,non=non)
  if success eq 0 then begin
     print,'down not found...'
     down=0.
  endif

 ;  extract Hg surface concentrations
 success=ctm_get_datablock(hg0, 'IJ-AVG-$',tracer=1,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,non=non,z=altitude)
 success=ctm_get_datablock(rgm, 'IJ-AVG-$',tracer=2,file=file,tau0=tau0,$
	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(hgp, 'IJ-AVG-$',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)

 success=ctm_get_datablock(newrgm, 'PL-HG2-$',tracer=10,file=file,tau0=tau0,$  
	lat=lat0,lon=lon0,non=non)  ; fg
  success=ctm_get_datablock(newrpm, 'PL-HG2-$',tracer=9,file=file,tau0=tau0,$  
 	lat=lat0,lon=lon0,non=non)  ; fp

 ; calculate rgm/rpm, note: reference model stores fraction, while new dataoutput stores concentrations
 if keyword_set(hg2frac) then begin
    newrgm*=(rgm+hgp)/(newrgm+newrpm)
    newrpm*=(rgm+hgp)/(newrgm+newrpm)

    temp=newrgm
    newrgm=newrpm
    newrpm=temp
 endif

 ;  extract dry dep loss for HG (molec/cm2/s)
  success=ctm_get_datablock(dry_hg0,'DRYD-FLX',tracer=1,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(dry_rgm,'DRYD-FLX',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(dry_hgp,'DRYD-FLX',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)

 ;  extract LS wet dep loss for HG (kg/s)
  success=ctm_get_datablock(wetls_rgm,'WETDLS-$',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,total=4,non=non,thisdatainfo=datainfo,y=latitude,x=longitude,gridinfo=gridinfo)
  success=ctm_get_datablock(wetls_hgp,'WETDLS-$',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,total=4,non=non)
 ;  extract CV wet dep loss for HG (kg/s)
  success=ctm_get_datablock(wetcv_rgm,'WETDCV-$',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,total=4,non=non)		
  success=ctm_get_datablock(wetcv_hgp,'WETDCV-$',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,total=4,non=non)

 ;  extract sea-salt loss of Hg2 (kg) - for the new Holmes simulation, this loss
 ;  is treated separately and is not included in the dry dep above
  success=ctm_get_datablock(ssloss_rgm,'PL-HG2-$',tracer=4,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)

 ; extract Br and BrO filed (molec/cm3)
  success=ctm_get_datablock(net_ox, 'PL-HG2-$',tracer=1,file=file,tau0=tau0,$  
	lat=lat0,lon=lon0,non=non)  
  success=ctm_get_datablock(br, 'PL-HG2-$',tracer=7,file=file,tau0=tau0,$  
	lat=lat0,lon=lon0,non=non)  
  success=ctm_get_datablock(bro, 'PL-HG2-$',tracer=8,file=file,tau0=tau0,$  
	lat=lat0,lon=lon0,non=non)  

 ;  extract EW flux of Hg 
  success=ctm_get_datablock(ew_hg0,'EW-FLX-$',tracer=1,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(ew_rgm,'EW-FLX-$',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(ew_hgp,'EW-FLX-$',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  if  success eq 0 then begin
       ew_hg0=0. & ew_rgm=0. & ew_hgp=0.
  endif

 ;  extract NS flux of Hg 
  success=ctm_get_datablock(ns_hg0,'NS-FLX-$',tracer=1,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(ns_rgm,'NS-FLX-$',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(ns_hgp,'NS-FLX-$',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  if  success eq 0 then begin
       ns_hg0=0. & ns_rgm=0. & ns_hgp=0.
  endif

 ;  extract UP  flux of Hg 
  success=ctm_get_datablock(up_hg0,'UP-FLX-$',tracer=1,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(up_rgm,'UP-FLX-$',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  success=ctm_get_datablock(up_hgp,'UP-FLX-$',tracer=3,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  if  success eq 0 then begin
       up_hg0=0. & up_rgm=0. & up_hgp=0.
  endif

 ;  read mass of air in grid box (in kg)
  success=ctm_get_datablock(airmass,'BXHGHT-$',tracer=2,file=file,tau0=tau0,$
 	lat=lat0,lon=lon0,non=non)
  if success eq 0 then begin
     print,'Airmass not found...'
     airmass=0.
  endif

 ;   calculate the number of seconds per month
  Seconds = Total( ( datainfo.Tau1 - datainfo.Tau0 ) * 3600d0 )

 ;  surface area
  A_Cm2 = CTM_BoxSize( GridInfo, /Cm2,/GEOS)
  A_m2 = CTM_BoxSize( GridInfo, /m2,/GEOS)
  a_cm2=a_cm2(where(gridinfo.xmid ge min(longitude) and $
  	           gridinfo.xmid le max(longitude) ),*)
  a_cm2=a_cm2(*,where(gridinfo.ymid ge min(latitude) and $
  	             gridinfo.ymid le max(latitude) ))
  a_m2=a_m2(where(gridinfo.xmid ge min(longitude) and $
  	           gridinfo.xmid le max(longitude) ),*)
  a_m2=a_m2(*,where(gridinfo.ymid ge min(latitude) and $
  	             gridinfo.ymid le max(latitude) ))

 if k eq 0 then begin
      tmp=hg0(*,*,0)*0.
      tmp1=hg0*0.
      conc_all=replicate({hg0:tmp,rgm:tmp,hgp:tmp,$
             newrgm:tmp, newrpm:tmp, $
			 hg:tmp,burden:0.,burden_rgm:0.,burden_hgp:0.,$
			 hg0aq:tmp,hg2aq:tmp,hgcaq:tmp,$
			 burden_hg0aq:0.,burden_hg2aq:0.,burden_hgcaq:0.},12)
      prof_all=replicate({hg0:tmp1,rgm:tmp1,hgp:tmp1,$
 			 hg:tmp1,alt:altitude},12)
      emis_all=replicate({anthr:tmp,anthr_rgm:tmp,anthr_hgp:tmp,$
      			 ocean:tmp,up:tmp,down:tmp,$
 			 land:tmp,$
      			 land_natural:tmp,bb:tmp,$
 			 veg:tmp,soil:tmp,hgsink:tmp,total:tmp},12)

      drydep_all=replicate({hg0:tmp,rgm:tmp,hgp:tmp,rgm_ss:tmp},12)
      wetdep_ls_all=replicate({rgm:tmp,hgp:tmp},12)
      wetdep_cv_all=replicate({rgm:tmp,hgp:tmp},12)
      wetdep_all=replicate({rgm:tmp,hgp:tmp},12)
      chem_all=replicate({net_ox:tmp1,br:tmp1,bro:tmp1},12)
      transp_all=replicate({ew_hg0:tmp1,ew_rgm:tmp1,ew_hgp:tmp1,$
                            ns_hg0:tmp1,ns_rgm:tmp1,ns_hgp:tmp1,$
                            up_hg0:tmp1,up_rgm:tmp1,up_hgp:tmp1 },12)
      tmp=0.
 endif 
 
 pptv_ngm3 = 8.93d0
 conc_all(k).hg0=hg0(*,*,0)*pptv_ngm3
 conc_all(k).rgm=rgm(*,*,0)*pptv_ngm3
 conc_all(k).hgp=hgp(*,*,0)*pptv_ngm3
 conc_all(k).newrgm=newrgm(*,*,0)*pptv_ngm3
 conc_all(k).newrpm=newrpm(*,*,0)*pptv_ngm3
 
  conc_all(k).hg=conc_all(k).hg0+conc_all(k).rgm+conc_all(k).hgp

 ;  aqueous concentration (convert from kg to ng/L)
  conc_convert= 1d11   / ( a_m2 ) / (mld(k).z*1e2)
  zero=where(mld(k).z eq 0)
  minzero=min(zero)
  if minzero ne -1 then conc_convert(zero)=0.

  conc_all(k).hg0aq=hg0aq*conc_convert
  conc_all(k).hg2aq=hg2aq*conc_convert
  conc_all(k).hgcaq=hgcaq*conc_convert
 ;  atmospheric profiles
  prof_all(k).hg0=hg0
  prof_all(k).rgm=rgm
  prof_all(k).hgp=hgp
  prof_all(k).hg=hg0+rgm+hgp
 ;  calculate burden in Mg
  conc_all(k).burden=total( (hg0+rgm+hgp) * 1d-12 * MolWt/28.97 * airmass * 1d-3)
  conc_all(k).burden_rgm=total( (rgm) * 1d-12 * MolWt/28.97 * airmass * 1d-3)
  conc_all(k).burden_hgp=total( (hgp) * 1d-12 * MolWt/28.97 * airmass * 1d-3)
 ;  same for the ocean in Mg
  conc_all(k).burden_hg0aq=total( hg0aq )*1e-3
  conc_all(k).burden_hg2aq=total( hg2aq )*1e-3
  conc_all(k).burden_hgcaq=total( hgcaq )*1e-3

 ;  emissions in Mg/month
  emis_all(k).anthr=anthr / 1e3
  emis_all(k).anthr_rgm=anthr_rgm / 1e3
  emis_all(k).anthr_hgp=anthr_hgp / 1e3
  emis_all(k).bb=bb / 1e3
  emis_all(k).land=land / 1e3
  emis_all(k).ocean=ocean / 1e3
  emis_all(k).land_natural=land_natural / 1e3
  emis_all(k).veg=land_veg / 1e3
  emis_all(k).soil=land_soil / 1e3
  emis_all(k).total=(anthr+bb+land+land_natural+land_veg+land_soil+$
  	anthr_rgm+anthr_hgp+ocean)/1e3
  emis_all(k).up=up / 1e3
  emis_all(k).down=down / 1e3
  emis_all(k).hgsink=hgsink / 1e3

 ;  convert from kg/s to Mg/month
  wetdep_ls_all(k).rgm=wetls_rgm * Seconds / 1e3
  wetdep_ls_all(k).hgp=wetls_hgp * Seconds / 1e3
  wetdep_cv_all(k).rgm=wetcv_rgm * Seconds / 1e3
  wetdep_cv_all(k).hgp=wetcv_hgp * Seconds / 1e3
  wetdep_all(k).rgm= wetdep_ls_all(k).rgm + wetdep_cv_all(k).rgm
  wetdep_all(k).hgp= wetdep_ls_all(k).hgp + wetdep_cv_all(k).hgp

 ;  convert from molec/cm2/s to Mg/month
  drydep_all(k).hg0=dry_hg0 * ( A_Cm2 / Avo ) * (MolWt * 1e-3  * Seconds ) / 1e3
  drydep_all(k).hgp=dry_hgp * ( A_Cm2 / Avo ) * (MolWt * 1e-3  * Seconds ) / 1e3
  drydep_all(k).rgm=dry_rgm * ( A_Cm2 / Avo ) * (MolWt * 1e-3  * Seconds ) / 1e3
  drydep_all(k).rgm_ss=ssloss_rgm / 1e3

 ; br
  chem_all(k).net_ox=net_ox / 1e3  ; Mg/month
  chem_all(k).br=br                ; molec/cm3
  chem_all(k).bro=bro

 ;  transport flux, convert from kg/s to Mg/month
  transp_all(k).ew_hg0=ew_hg0 * Seconds / 1e3
  transp_all(k).ew_rgm=ew_rgm * Seconds / 1e3
  transp_all(k).ew_hgp=ew_hgp * Seconds / 1e3
  transp_all(k).ns_hg0=ns_hg0 * Seconds / 1e3
  transp_all(k).ns_rgm=ns_rgm * Seconds / 1e3
  transp_all(k).ns_hgp=ns_hgp * Seconds / 1e3
  transp_all(k).up_hg0=up_hg0 * Seconds / 1e3
  transp_all(k).up_rgm=up_rgm * Seconds / 1e3
  transp_all(k).up_hgp=up_hgp * Seconds / 1e3

endfor
 ; conversion factor from Mg/year to ug/m2/yr
 ;  1 Mg = 1e3 kg = 1e6 g = 1e6*1e6 ug
  convert=1e6*1e6/a_m2

end
