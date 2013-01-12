pro plot_mdn,file,resolution,lon0=lon0,lat0=lat0,select_year=select_year,$
               year=year,debug=debug,rgm=rgm,hgp=hgp,region=region,bound=bound,regrid=regrid,sites=sites,$
	    	   monthly=monthly,seasonal=seasonal,interannual=interannual,snow=snow, sn0w=sn0w, anomaly=anomaly,$
	           version2=version2,resolution2=resolution2,$
	           version3=version3,resolution3=resolution3,$		   
	    	   version4=version4,resolution4=resolution4, $
               mettype=mettype, metresolution=metresolution, $
               vertical=vertical
; try to use ug/m2 as the unit
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;  Create plot of annual deposition from model overlayed with annual deposition
;  from the MDN sites. Also plot seasonal cycle of monthly mean deposition averaged
;  over the MDN sites in the lon0/lat0 region
;   
;  INPUT:
;    version --> the version number for run
;    resolution --> the resolution of model output: '4x5' or '05x0666'

;  KEYWORDS:
;    LAT0     --> vector containing min/max of latitudes to be plotted. Default
;              is LAT0=[20,57]
;    LON0     --> min/max of longitudes to be plotted.
;    SELECT_YEAR --> range of MDN years to be used for observation. Example: [2004,2008]
;    YEAR     --> Same as SELECT_YEAR, but for model runs
;    rgm      --> only plot that from gaseous phase
;    hgp      --> only plot wet deposition from particulate phase
;    REGION   --> if /region is set then region boundaries are plotted on top of the map and 
;               the seasonal precipitation is plotted for each season.
;    REGRID   --> if /regrid is set then the MDN data will be gridded onto the GEOS-Chem
;                grid based on the resolution of the simulation.
;    SITES    --> plot seasonal cycle at individual sites.
;    VERSION2 --> Name of second model version to plot (this will be plotted with a red line on the seasonal plots).
;    MONTHLY  --> Plot maps of monthly wet deposition
;    SEASONAL -->
;    SNOW correlation
; 
;  EXAMPLES:
;  IDL> plot_mdn,'01-12-11.1','4x5'
;  IDL> plot_mdn,'01-12-11.1','05x0666',year=[2007,2008],select_year=[2007,2008]
; 
; 
;  Created by L. Jaeglé, October 2009.
;  Updated by Y. Zhang, Jan, 2011
;  Updated by L. Jaegle, Feb 2011 (plot wet deposition by region, grid mdn data onto model grid)
;                               - plot individual sites (/site)
;                               - plot maps of monthly (/monthly) and seasonal (/seasonal) wet deposition
; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ;  Default long/lat range for plots
 if not keyword_set(lon0) then lon0=[-140,-40]
 if not keyword_set(lat0) then lat0=[10,70]
 if not keyword_set(snow) then snow=1
 if not keyword_set(select_year) then select_year=[2004,2010]
 if not keyword_set(year) then year=[2004,2010]
 if not keyword_set(mettype) then mettype='GEOS5'
 if not keyword_set(metresolution) then metresolution='05x0666'
 
if  resolution eq '4x5' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=4 )
if  resolution eq '05x0666' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=[2./3.,1./2.] )
GridInfo   = CTM_GRID( ModelInfo )
dlon=gridinfo.di/2. ;  half of di
dlat=gridinfo.dj/2. ;  half of dj
lon_all=gridinfo.xmid
lat_all=gridinfo.ymid

 ;  define plotting limit, taking into account the grid size
 ind0=where(lat_all ge lat0(0))
 ind1=where(lon_all ge lon0(0))
 ind2=where(lat_all ge lat0(1))
 ind3=where(lon_all ge lon0(1))
 lat0=[lat_all(ind0(0))-dlat,lat_all(ind2(0))+dlat]
 lon0=[lon_all(ind1(0))-dlon,lon_all(ind3(0))+dlon]

 if not keyword_set(cven) then cven=1.0
 ; if not keyword_set(lon0) then lon0=[-70,-40]
 ; if not keyword_set(lat0) then lat0=[20,36]

; number of years
myears=max(year)-min(year)+1L

 ;  number of years selected for MDN data
 nyears=max(select_year)-min(select_year)+1

 ;  Strings corresponding to selected years for MDN data and years for model simulations
 ;  will be used in the plotting (jaegle 1/31/2011)
 if n_elements(year) eq 1 then year_str=strcompress(year,/remo) else $
 	year_str=strcompress(min(year),/remo)+'-'+strcompress(max(year),/remo)
 if min(select_year) eq max(select_year) then select_year_str=strcompress(select_year(0),/remo) else $
 	select_year_str=strcompress(select_year(0),/remo)+'-'+strcompress(select_year(1),/remo)

 ;  Select minimum number of points to be used in selecting MDN sites.
 ;  In a full data record in a given year, there are 52 weekly measurements.
 ;  So a 90% temporal requirement would mean mdn_data.npts> 47 points.
 ;     a 80% temporal requirement would mean mdn_data.npts> 41 points.
 ;  Multiply by the number of years.
;  min_npts=41    ; *nyears ;  80% coverage
;  min_npts=35 ;  70% coverage
;   min_npts=20  ; 50% coverage ?
;   min_npts=1  ; 25% coverage ?
   min_npts=39  ; 75% coverage
;   min_npts=0
;   min_npts=10    ; no filter

 smethod='strict'

 ;  Read in the MDN data 
; store the data in some files
 filename = 'mdn_data_' + string(select_year,format='(I4,"-",I4)') + '_' + resolution + '_'  + $
                  string(snow, format='(I1)') + '_' + mettype +  '_' + metresolution + '.sav' 
 status = file_info(filename)
 if status.exists then begin
   print, 'file exist, now reading ...'
   restore, filename=filename 
 endif else begin
   print, 'file not exist, now generating ...'
   read_mdn,mdn_data,select_year=select_year, snow=snow, mettype=mettype, metresolution=metresolution ,/correct
   save, mdn_data, filename=filename
 endelse

; read in original MDN data
 if keyword_set(sn0w) then read_mdn,mdn_data0,select_year=select_year, snow=0, mettype=mettype, metresolution=metresolution, /correct

 ;  Select sites in plotting domain and a minimum number of points available:
 ;  this number will depend on the number of years selected. Typically, for
 ;  a full data record in a given year, there are 52 weekly measurements.

; ind=where(mdn_data.lat ge lat0(0) and mdn_data.lat le lat0(1) and $
;           mdn_data.lon ge lon0(0) and mdn_data.lon le lon0(1) and $
;	   min(mdn_data.npts, dimension=1) ge min_npts)  ; and $

if smethod eq 'strict'  then begin
      ind=where(mdn_data.lat ge lat0(0) and mdn_data.lat le lat0(1) and $
                 mdn_data.lon ge lon0(0) and mdn_data.lon le lon0(1) and $
	             min(mdn_data.npts, dimension=1) ge min_npts ) 
endif else begin
      if smethod eq 'loose'  then begin
          ind=where(mdn_data.lat ge lat0(0) and mdn_data.lat le lat0(1) and $
                       mdn_data.lon ge lon0(0) and mdn_data.lon le lon0(1) )
      endif else begin
          if smethod eq 'hybrid'  then begin
           ind=where(mdn_data.lat ge lat0(0) and mdn_data.lat le lat0(1) and $
                           mdn_data.lon ge lon0(0) and mdn_data.lon le lon0(1) and $
          	             total(mdn_data.npts,1) ge min_npts * nyears * 52 )
          endif else begin
             print, 'no matching criteria for selecting site! stopping ...'
             stop
          endelse
      endelse
endelse

 print, string(n_elements(ind),format='(I3)') + ' sites are selected'

;	   mdn_data.annual_dep gt 0.  and $
;	   min(mdn_data.season_dep, dimension=1) gt 0. and $
;       min(mdn_data.monthly_dep, dimension=1) gt 0. )
	   
 ; ind=where(mdn_data.lat ge lat0(0)-dlat and mdn_data.lat le lat0(1)+dlat and $
 ;           mdn_data.lon ge lon0(0)-dlon and mdn_data.lon le lon0(1)+dlon and $
; 	   mdn_data.npts ge min_npts)
	   ; mdn_data.npts ge 47);  1 year, 90% temporal coverage
	   ; mdn_data.npts ge 41);  1 year, 80% temporal coverage
 mdn_data=mdn_data(ind)
 if keyword_set(sn0w) then  mdn_data0=mdn_data0(ind)

 mdn_data_all=mdn_data

 ;  convert MDN annual deposition data from ng/m2/day to ug/m2/year
 mdn_wetdep=mdn_data.annual_dep*1e-3*365.
 nsites=n_elements(mdn_data.id)

; store nsites
 nsites_all = nsites

 mdn_monthly=mean(mdn_data.monthly_dep,2,/nan)
 mdn_monthly_npts=total(mdn_data.monthly_npts,2,/nan)

 ;  extract monthly mean model values for each individual site 
 mdn_wetdep_monthly_site=mdn_data.monthly_dep 
 mdn_wetdep_annual_site=mdn_data.annual_dep*1e-3*365
 mdn_wetdep_annual_site_std=mdn_data.annual_dep_std*1e-3*365

 ;  Then read in the model output
 ;file_std='./ctm.bpch.'+version
 file_std=file

 ;  option to read a second model
; if keyword_set(version2) then file_std2='/home/disk/p/yanxuz/4-geos/0-hg/Run.'+resolution2+'/dataout/ctm.bpch.'+version2
 if keyword_set(version2) then file_std2='./' + version2
 if keyword_set(version3) then file_std3='/home/disk/p/yanxuz/4-geos/0-hg/Run.'+resolution3+'/dataout/ctm.bpch.'+version3
 if keyword_set(version4) then file_std4='/home/disk/p/yanxuz/4-geos/0-hg/Run.'+resolution4+'/dataout/ctm.bpch.'+version4
 

 if strcmp(resolution,'4x5') then res=4 else res=[2./3.,1./2.]
 if keyword_set(version2) then begin
    if strcmp(resolution2,'4x5') then res2=4 else res2=[2./3.,1./2.]
 endif
 if keyword_set(version3) then begin
    if strcmp(resolution3,'4x5') then res3=4 else res3=[2./3.,1./2.]
 endif
 if keyword_set(version4) then begin
    if strcmp(resolution4,'4x5') then res4=4 else res4=[2./3.,1./2.]
 endif
 
 year_num=1.
 for iyear=min(year),max(year) do begin
   print, 'Reading model file ',file_std, ' for year ',iyear

   read_run_wdep,file=file_std,conc_all,emis_all,drydep_tmp,$
            wetdep_tmp,longitude,latitude,convert,data_model,prof_all,$
	    gridinfo,lat0=[lat0(0)+dlat,lat0(1)-dlat],lon0=[lon0(0)+dlon,lon0(1)-dlon],year=iyear
   lat_all=gridinfo.ymid 
   lon_all=gridinfo.xmid

   if keyword_set(version2) then  begin
      ; special treat for version2, because only 2008 2009 data is available
;      if iyear mod 2 eq 1 then iiyear=2008 else iiyear=2009
      read_run_wdep,file=file_std2,conc_all2,emis_all2,drydep_tmp2,$
            wetdep_tmp2,longitude2,latitude2,convert2,data_model2,prof_all2,$
	    gridinfo2,lat0=[lat0(0)+dlat,lat0(1)-dlat],lon0=[lon0(0)+dlon,lon0(1)-dlon],year=iyear
      lat_all2=gridinfo2.ymid 
      lon_all2=gridinfo2.xmid
   endif
     if keyword_set(version3) then  begin
      read_run_wdep,file=file_std3,conc_all3,emis_all3,drydep_tmp3,$
            wetdep_tmp3,longitude3,latitude3,convert3,data_model3,prof_all3,$
	    gridinfo3,lat0=[lat0(0)+dlat,lat0(1)-dlat],lon0=[lon0(0)+dlon,lon0(1)-dlon],year=iyear
      lat_all3=gridinfo3.ymid 
      lon_all3=gridinfo3.xmid
   endif
     if keyword_set(version4) then  begin
      read_run_wdep,file=file_std4,conc_all4,emis_all4,drydep_tmp4,$
            wetdep_tmp4,longitude4,latitude4,convert4,data_model4,prof_all4,$
	    gridinfo4,lat0=[lat0(0)+dlat,lat0(1)-dlat],lon0=[lon0(0)+dlon,lon0(1)-dlon],year=iyear
      lat_all4=gridinfo4.ymid 
      lon_all4=gridinfo4.xmid
   endif
   
   if iyear eq min(year) then begin
     wetdep_all=wetdep_tmp
     drydep_all=drydep_tmp
     if keyword_set(version2) then begin
     	wetdep_all2=wetdep_tmp2
     	drydep_all2=drydep_tmp2
     endif
     if keyword_set(version3) then begin
     	wetdep_all3=wetdep_tmp3
     	drydep_all3=drydep_tmp3
     endif
     if keyword_set(version4) then begin
     	wetdep_all4=wetdep_tmp4
     	drydep_all4=drydep_tmp4
     endif

     ; anomaly
    if keyword_set(anomaly) then begin
      di=size(wetdep_all.rgm)
      model_anomaly=fltarr(di[1],di[2],myears,12)     ; 12 monthes
      for imonth=0,11 do begin
         model_anomaly[*,*,0,imonth] = (wetdep_tmp.rgm + wetdep_tmp.hgp)[*,*,imonth]		 		 
      endfor
    endif
	 
   endif else begin
   
   year_num=year_num+1.
   
   wetdep_all.rgm=wetdep_all.rgm+wetdep_tmp.rgm
   wetdep_all.hgp=wetdep_all.hgp+wetdep_tmp.hgp
   wetdep_all.vrgm=wetdep_all.vrgm+wetdep_tmp.vrgm
   wetdep_all.vhgp=wetdep_all.vhgp+wetdep_tmp.vhgp
   wetdep_all.vls=wetdep_all.vls+wetdep_tmp.vls
   wetdep_all.vcv=wetdep_all.vcv+wetdep_tmp.vcv

   drydep_all.total=drydep_all.total+drydep_tmp.total
   
   if keyword_set(version2) then begin
     wetdep_all2.rgm=wetdep_all2.rgm+wetdep_tmp2.rgm
     wetdep_all2.hgp=wetdep_all2.hgp+wetdep_tmp2.hgp
     wetdep_all2.vrgm=wetdep_all2.vrgm+wetdep_tmp2.vrgm
     wetdep_all2.vhgp=wetdep_all2.vhgp+wetdep_tmp2.vhgp
     wetdep_all2.vls=wetdep_all2.vls+wetdep_tmp2.vls
     wetdep_all2.vcv=wetdep_all2.vcv+wetdep_tmp2.vcv
     drydep_all2.total=drydep_all2.total+drydep_tmp2.total
   endif
   if keyword_set(version3) then begin
     wetdep_all3.rgm=wetdep_all3.rgm+wetdep_tmp3.rgm
     wetdep_all3.hgp=wetdep_all3.hgp+wetdep_tmp3.hgp
     wetdep_all3.vrgm=wetdep_all3.vrgm+wetdep_tmp3.vrgm
     wetdep_all3.vhgp=wetdep_all3.vhgp+wetdep_tmp3.vhgp
     wetdep_all3.vls=wetdep_all3.vls+wetdep_tmp3.vls
     wetdep_all3.vcv=wetdep_all3.vcv+wetdep_tmp3.vcv
     drydep_all3.total=drydep_all3.total+drydep_tmp3.total
   endif
   if keyword_set(version4) then begin
     wetdep_all4.rgm=wetdep_all4.rgm+wetdep_tmp4.rgm
     wetdep_all4.hgp=wetdep_all4.hgp+wetdep_tmp4.hgp
     wetdep_all4.vrgm=wetdep_all4.vrgm+wetdep_tmp4.vrgm
     wetdep_all4.vhgp=wetdep_all4.vhgp+wetdep_tmp4.vhgp
     wetdep_all4.vls=wetdep_all4.vls+wetdep_tmp4.vls
     wetdep_all4.vcv=wetdep_all4.vcv+wetdep_tmp4.vcv
     drydep_all4.total=drydep_all4.total+drydep_tmp4.total
   endif

     ; anomaly
    if keyword_set(anomaly) then begin
     for imonth=0,11 do begin
       model_anomaly[*,*,iyear-min(year),imonth] = (wetdep_tmp.rgm + wetdep_tmp.hgp)[*,*,imonth]
     endfor
    endif
      
   endelse

endfor

if keyword_set(rgm) and not keyword_set(hgp) then wetdep_all.hgp = wetdep_all.hgp * 0.
if not keyword_set(rgm) and keyword_set(hgp) then wetdep_all.rgm = wetdep_all.rgm * 0.

wetdep_all.rgm=wetdep_all.rgm/year_num
wetdep_all.hgp=wetdep_all.hgp/year_num
wetdep_all.vrgm=wetdep_all.vrgm/year_num
wetdep_all.vhgp=wetdep_all.vhgp/year_num
wetdep_all.vls=wetdep_all.vls/year_num
wetdep_all.vcv=wetdep_all.vcv/year_num
drydep_all.total=drydep_all.total/year_num
model_wetdep =total(wetdep_all.rgm  +  wetdep_all.hgp  ,3)*convert
model_wetdep_plot=model_wetdep
model_drydep =total(drydep_all.total  ,3)*convert

if keyword_set(version2) then begin
  wetdep_all2.rgm=wetdep_all2.rgm/year_num
  wetdep_all2.hgp=wetdep_all2.hgp/year_num
  wetdep_all2.vrgm=wetdep_all2.vrgm/year_num
  wetdep_all2.vhgp=wetdep_all2.vhgp/year_num
  wetdep_all2.vls=wetdep_all2.vls/year_num
  wetdep_all2.vcv=wetdep_all2.vcv/year_num
  drydep_all2.total=drydep_all2.total/year_num
  model_wetdep2 =total(wetdep_all2.rgm  +  wetdep_all2.hgp  ,3)*convert2
  model_wetdep_plot2=model_wetdep2
  model_drydep2 =total(drydep_all2.total  ,3)*convert2
endif
if keyword_set(version3) then begin
  wetdep_all3.rgm=wetdep_all3.rgm/year_num
  wetdep_all3.hgp=wetdep_all3.hgp/year_num
  wetdep_all3.vrgm=wetdep_all3.vrgm/year_num
  wetdep_all3.vhgp=wetdep_all3.vhgp/year_num
  wetdep_all3.vls=wetdep_all3.vls/year_num
  wetdep_all3.vcv=wetdep_all3.vcv/year_num
   drydep_all3.total=drydep_all3.total/year_num
  model_wetdep3 =total(wetdep_all3.rgm  +  wetdep_all3.hgp  ,3)*convert3
  model_wetdep_plot3=model_wetdep3
  model_drydep3 =total(drydep_all3.total  ,3)*convert3
endif
if keyword_set(version4) then begin
  wetdep_all4.rgm=wetdep_all4.rgm/year_num
  wetdep_all4.hgp=wetdep_all4.hgp/year_num
  wetdep_all4.vrgm=wetdep_all4.vrgm/year_num
  wetdep_all4.vhgp=wetdep_all4.vhgp/year_num
  wetdep_all4.vls=wetdep_all4.vls/year_num
  wetdep_all4.vcv=wetdep_all4.vcv/year_num
  drydep_all4.total=drydep_all4.total/year_num
  model_wetdep4 =total(wetdep_all4.rgm  +  wetdep_all4.hgp  ,3)*convert4
  model_wetdep_plot4=model_wetdep4
  model_drydep4 =total(drydep_all4.total  ,3)*convert4
endif

; calculate the anomaly
 if keyword_set(anomaly) then begin
    for iyear=0,myears-1 do begin
      for imonth=0,11 do begin
        ; model_anomaly[*,*,iyear,imonth]=(model_anomaly[*,*,iyear,imonth]-(wetdep_all.rgm+wetdep_all.hgp)[*,*,imonth])*convert    ; absolute anomaly, ug/m2/month
         model_anomaly[*,*,iyear,imonth]=model_anomaly[*,*,iyear,imonth]*convert    ; absolute anomaly, ug/m2/month
      endfor
    endfor
 endif

if strcmp(resolution,'4x5') then begin
  hlat = 2.0
  hlon = 2.5
endif else begin
  hlat = 0.25
  hlon = 1.0/3.0
endelse

if keyword_set(version2) then begin
if strcmp(resolution2,'4x5') then begin
  hlat2 = 2.0
  hlon2 = 2.5
endif else begin
  hlat2 = 0.25
  hlon2 = 1.0/3.0
endelse
endif
if keyword_set(version3) then begin
if strcmp(resolution3,'4x5') then begin
  hlat3 = 2.0
  hlon3 = 2.5
endif else begin
  hlat3 = 0.25
  hlon3 = 1.0/3.0
endelse
endif
if keyword_set(version4) then begin
if strcmp(resolution4,'4x5') then begin
  hlat4 = 2.0
  hlon4 = 2.5
endif else begin
  hlat4 = 0.25
  hlon4 = 1.0/3.0
endelse
endif


;  convert wet deposition from Mg/month to ng/m2/day (annual)

 model_wetdep_monthly=fltarr(12)

 ndays=julday(2+indgen(12),1,2001)-julday(1+indgen(12),1,2001)
 for k=0,12-1 do begin
      model_wetdep_monthly(k)=mean((wetdep_all(k).rgm+ $
      wetdep_all(k).hgp)*convert)*1e3/ndays(k)
 endfor

model_wetdep_monthly_site=fltarr(12,nsites)
model_wetdep_seasonal_site=fltarr(4,nsites)
model_wetdep_annual_site=fltarr(nsites)
model_drydep_annual_site=fltarr(nsites)

if keyword_set(anomaly) then begin
  model_wetdep_monthly_site_anomaly=fltarr(12,myears,nsites)
  model_wetdep_seasonal_site_anomaly=fltarr(4,myears,nsites)
  model_wetdep_annual_site_anomaly=fltarr(myears,nsites)
endif

if keyword_set(version2) then begin
   model_wetdep_monthly_site2=model_wetdep_monthly_site
   model_wetdep_annual_site2=model_wetdep_annual_site
   model_drydep_annual_site2=model_wetdep_annual_site
   model_wetdep_seasonal_site2=model_wetdep_seasonal_site
endif
if keyword_set(version3) then begin
   model_wetdep_monthly_site3=model_wetdep_monthly_site
   model_wetdep_annual_site3=model_wetdep_annual_site
   model_drydep_annual_site3=model_wetdep_annual_site
   model_wetdep_seasonal_site3=model_wetdep_seasonal_site
endif
if keyword_set(version4) then begin
   model_wetdep_monthly_site4=model_wetdep_monthly_site
   model_wetdep_annual_site4=model_wetdep_annual_site
   model_drydep_annual_site4=model_wetdep_annual_site
   model_wetdep_seasonal_site4=model_wetdep_seasonal_site
endif

season=replicate({month:fltarr(3),name:''},4)
season(0).name='DJF' & season(0).month=[12,1,2]-1
season(1).name='MAM' & season(1).month=[3,4,5]-1
season(2).name='JJA' & season(2).month=[6,7,8]-1
season(3).name='SON' & season(3).month=[9,10,11]-1

for k=0,nsites-1 do begin
     ; find location of MDN site
     ctm_index,ctm_type('GEOS5_47L',res=res),i,j,$
	        center=[mdn_data(k).lat,mdn_data(k).lon],/non_interactive
     indlat=where( latitude eq lat_all(j-1) )
     indlon=where( longitude eq lon_all(i-1) )
    
     if keyword_set(debug) then begin
       print,k,mdn_data(k).lat,mdn_data(k).lon
       print,k,latitude(indlat),longitude(indlon)
     endif

     ;  ng/m2/day
     if indlat[0] ne -1 and indlon[0] ne -1 then begin
     model_wetdep_monthly_site(*,k)=(wetdep_all.rgm(indlon(0),indlat(0))+wetdep_all.hgp(indlon(0),indlat(0)))*convert(indlon(0),indlat(0))*1e3/ndays      ; convert: to ug/m2

     if keyword_set(anomaly) then begin
       for iyear=0,myears-1 do begin
          for imonth=0,11 do begin
                 model_wetdep_monthly_site_anomaly[imonth,iyear,k]=model_anomaly[indlon(0),indlat(0),iyear,imonth]    ; ug/m2
          endfor
       endfor
     endif

     model_wetdep_annual_site(k)=total(wetdep_all.rgm(indlon(0),indlat(0))+wetdep_all.hgp(indlon(0),indlat(0)))*convert(indlon(0),indlat(0))
     if keyword_set(anomaly) then begin
         for iyear=0,myears-1 do begin
              model_wetdep_annual_site_anomaly[iyear,k]=total(model_anomaly[indlon(0),indlat(0),iyear, * ])    ; ug/m2
         endfor
     endif

     model_drydep_annual_site(k)=total(drydep_all.total(indlon(0),indlat(0)))*convert(indlon(0),indlat(0))

     ;  calculate seasonal mean  (ug/m2/season)
     for l=0,3 do begin
              model_wetdep_seasonal_site(l,k)=mean(wetdep_all(season(l).month).rgm(indlon(0),indlat(0))+$
                                          wetdep_all(season(l).month).hgp(indlon(0),indlat(0)))*convert(indlon(0),indlat(0))*12./4.
              if keyword_set(anomaly) then begin
                   for iyear=0,myears-1 do begin
                         model_wetdep_seasonal_site_anomaly[l,iyear,k]=total(model_anomaly[indlon(0),indlat(0),iyear, season(l).month ])   ; ug/m2
                   endfor
              endif
     endfor

	 endif
	 
     if keyword_set(version2) then begin
         ctm_index,ctm_type('GEOS5_47L',res=res2),i,j,$
    	        center=[mdn_data(k).lat,mdn_data(k).lon],/non_interactive
         indlat2=where( latitude2 eq lat_all2(j-1) )
         indlon2=where( longitude2 eq lon_all2(i-1) )
         if keyword_set(debug) then begin
           print,k,mdn_data(k).lat,mdn_data(k).lon
           print,k,latitude2(indlat2),longitude2(indlon2)
         endif
         ;  ng/m2/day
         if indlat2[0] ne -1 and indlon2[0] ne -1 then begin
         model_wetdep_monthly_site2(*,k)=(wetdep_all2.rgm(indlon2(0),indlat2(0))+wetdep_all2.hgp(indlon2(0),indlat2(0)))*convert2(indlon2(0),indlat2(0))*1e3/ndays
         model_wetdep_annual_site2(k)=total(wetdep_all2.rgm(indlon2(0),indlat2(0))+wetdep_all2.hgp(indlon2(0),indlat2(0)))*convert2(indlon2(0),indlat2(0))
         model_drydep_annual_site2(k)=total(drydep_all2.total(indlon2(0),indlat2(0)))*convert2(indlon2(0),indlat2(0))
    
         ;  calculate seasonal mean  (ug/m2/season)
         for l=0,3 do model_wetdep_seasonal_site2(l,k)=mean(wetdep_all2(season(l).month).rgm(indlon2(0),indlat2(0))+$
                                              wetdep_all2(season(l).month).hgp(indlon2(0),indlat2(0)))*convert2(indlon2(0),indlat2(0))*12./4.
	     endif
      endif
	  
     if keyword_set(version3) then begin
         ctm_index,ctm_type('GEOS5_47L',res=res3),i,j,$
    	        center=[mdn_data(k).lat,mdn_data(k).lon],/non_interactive
         indlat3=where( latitude3 eq lat_all3(j-1) )
         indlon3=where( longitude3 eq lon_all3(i-1) )
         if keyword_set(debug) then begin
           print,k,mdn_data(k).lat,mdn_data(k).lon
           print,k,latitude3(indlat3),longitude3(indlon3)
         endif
         ;  ng/m2/day
         if indlat3[0] ne -1 and indlon3[0] ne -1 then begin
         model_wetdep_monthly_site3(*,k)=(wetdep_all3.rgm(indlon3(0),indlat3(0))+wetdep_all3.hgp(indlon3(0),indlat3(0)))*convert3(indlon3(0),indlat3(0))*1e3/ndays
         model_wetdep_annual_site3(k)=total(wetdep_all3.rgm(indlon3(0),indlat3(0))+wetdep_all3.hgp(indlon3(0),indlat3(0)))*convert3(indlon3(0),indlat3(0))
         model_drydep_annual_site3(k)=total(drydep_all3.total(indlon3(0),indlat3(0)))*convert3(indlon3(0),indlat3(0))
    
         ;  calculate seasonal mean  (ug/m2/season)
         for l=0,3 do model_wetdep_seasonal_site3(l,k)=mean(wetdep_all3(season(l).month).rgm(indlon3(0),indlat3(0))+$
                                              wetdep_all3(season(l).month).hgp(indlon3(0),indlat3(0)))*convert3(indlon3(0),indlat3(0))*12./4.
	     endif
      endif
	  
	   if keyword_set(version4) then begin
         ctm_index,ctm_type('GEOS5_47L',res=res4),i,j,$
    	        center=[mdn_data(k).lat,mdn_data(k).lon],/non_interactive
         indlat4=where( latitude4 eq lat_all4(j-1) )
         indlon4=where( longitude4 eq lon_all4(i-1) )
         if keyword_set(debug) then begin
           print,k,mdn_data(k).lat,mdn_data(k).lon
           print,k,latitude4(indlat4),longitude4(indlon4)
         endif
         ;  ng/m2/day
         if indlat4[0] ne -1 and indlon4[0] ne -1 then begin
         model_wetdep_monthly_site4(*,k)=(wetdep_all4.rgm(indlon4(0),indlat4(0))+wetdep_all4.hgp(indlon4(0),indlat4(0)))*convert4(indlon4(0),indlat4(0))*1e3/ndays
         model_wetdep_annual_site4(k)=total(wetdep_all4.rgm(indlon4(0),indlat4(0))+wetdep_all4.hgp(indlon4(0),indlat4(0)))*convert4(indlon4(0),indlat4(0))
         model_drydep_annual_site4(k)=total(drydep_all4.total(indlon4(0),indlat4(0)))*convert4(indlon4(0),indlat4(0))
    
         ;  calculate seasonal mean  (ug/m2/season)
         for l=0,3 do model_wetdep_seasonal_site4(l,k)=mean(wetdep_all4(season(l).month).rgm(indlon4(0),indlat4(0))+$
                                              wetdep_all4(season(l).month).hgp(indlon4(0),indlat4(0)))*convert4(indlon4(0),indlat4(0))*12./4.
	     endif
      endif
	  
 endfor

; snow_str='_snow'+string(snow,format='(I1)')
; file_ps='mdn_'+version+'_'+year_str+'_'+resolution+snow_str+'.ps'
; open_device,/ps,/color,/portrait,file=file_ps
; !P.font=0
; Device, /Helvetica, /Isolatin1
 micro=textoidl('\mu')
 c_level=[0,0.5,1,2,3,4,5,6, 8,10,12,14,16,18,20,22,24]
 multipanel,col=1,row=2,omargin=0.08,margin=[0.07,0.09]

 MyCt, /Verbose, /WhGrYlRd,ncolors=20
 ; print,model_wetdep_plot(where(abs(longitude+110.83) le hlon),where(abs(latitude-51.42) le hlat))
 ; limit=[latitude[0]-hlat,longitude[0]-hlon,max(latitude)+hlat,max(longitude)+hlon]
 
 limit=[22,  -130, 55, -53]   ; us
; limit=[25, -125, 50, -64]
; limit=[37,-85, 43, -75] ; ohio
 ctm_overlay,model_wetdep_plot,longitude,latitude,$
     /sample,mdn_wetdep,mdn_data.lon,mdn_data.lat,$
	/continents,thick=0.5, /cbar,$
	div=COLORBAR_NDIV(maxdiv=8),$
	mindata=0.,maxdata=18.0,$
	/usa,limit=limit,$
	t_symbol=1,symsize=1.1,title='Annual Mean Mercury Wet Deposition Flux'+$
       '!cMDN data for ' + select_year_str + '('+year_str+')',$
        cbunit=micro+'g m!u-2!n y!u-1!n',/iso,csfac=1.2


 ;  overplot region boundaries
 ; note: the last region could only be 'ALL' or 'USA'
 if keyword_set(region) then begin
      ; define some regions to average the MDN data
      ; region_all=['PW', 'W', 'C', 'MW', 'NE', 'OH', 'SE', 'Gulf', 'Florida', 'USA']
;       region_all=['PW', 'W', 'C', 'MW', 'NE', 'OH', 'SE', 'Gulf', 'Florida']
;       region_all=['PW', 'W', 'C', 'MW', 'NE', 'OH', 'SE']
;      region_all=['PW', 'W', 'MW', 'NE', 'ORV', 'SE', 'ALL']
        region_all=['MW', 'NE', 'ORV', 'SE']

      ; This is for second paper 
;       region_all=['NORTHEAST', 'MIDWEST', 'SOUTHEAST', 'WEST']

      ; region_all=['MW','NE','OH','SE','USA']
      ; region_all=['SE','Gulf','Florida']
      ; region_all=['Gulf','Florida']
      nregion=n_elements(region_all)
      limit_all=replicate({x:fltarr(5),y:fltarr(5)},nregion)
      for k=0,nregion-1 do begin
         mdn_regions,region_all(k),y,x
         limit_all(k).x=[x(0),x(1),x(1),x(0),x(0)]
	     limit_all(k).y=[y(0),y(0),y(1),y(1),y(0)]
      endfor
      if keyword_set(bound) then begin
      trackd=fltarr(5)
      for k=0,nregion-1 do ctm_overlay,model_wetdep_plot,longitude,latitude,$
         trackd,limit_all(k).x,limit_all(k).y,$
	     t_color=1,t_thick=3,/overplot,$
         limit=limit,/iso
       endif
      ; k=nregion-1 ;  Gulf Coast
      ; ctm_overlay,model_wetdep_plot,longitude,latitude,$
      ;    trackd,limit_all(k).x,limit_all(k).y,$
      ; 	 t_color=2,t_thick=3,/overplot,$
      ;         limit=limit,/iso
 endif

;  also plot seasonal cycle (averaged over entire region selected)
 multipanel,col=2,row=2,omargin=0.03,margin=[0.035,0.08],/noerase
 multipanel,/advance,/noerase,pos=pos
 multipanel,/advance,/noerase,pos=pos
;  multipanel,/advance,/noerase,pos=pos 

 ; yrange=[0,max([max(mean(model_wetdep_monthly_site,2,/nan))*1.2,max(mdn_monthly)*1.2,60])]
 yrange=[0,60]
 plot,1+indgen(12),mdn_monthly,col=1,pos=pos,$
 	ytitle='Monthly deposition (ng m!u-2!nd!u-1!n)',thick=5,yrange=yrange,$
 	xrange=[0.4,12.6],/xst,title='Seasonal cycle of mercury deposition',$
	xtitle='Month',charsize=0.8

 ; for k=0,nsites-1 do oplot,1+indgen(12),mdn_data(k).monthly_dep,col=1,thick=1,line=1
 ; oplot,1+indgen(12),mdn_monthly,col=1,thick=5
 xyouts,1+indgen(12),!y.crange(1)*0.05,string(mdn_monthly_npts,format='(I4)'),col=1,al=0.5,charsize=0.7
 oplot,1+indgen(12),model_wetdep_monthly,col=2,thick=5
 oplot,1+indgen(12),mean(model_wetdep_monthly_site,2,/nan),col=!myct.blue,thick=5
 
 legend,label=['MDN sites '+  string(mean(mdn_monthly), format= '(F5.2)') + ' ' +  'ngm!u-2!nd!u-1!n' ,$
                     'Model all domain '  + string(mean(model_wetdep_monthly), format= '(F5.2)') + ' ' + 'ngm!u-2!nd!u-1!n'  ,$
					 'Model @ MDN sites ' + string(mean(mean(model_wetdep_monthly_site,2,/nan)), format= '(F5.2)') + ' ' + 'nm!u-2!nd!u-1!n'  $
					 ],lcol=[1,2,!myct.blue],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=4,$
 	halign=0.01,valign=0.95,charsize=0.8


 ;  bin data by latitude
;  multipanel,/advance,pos=pos
;  plot,mdn_data.lat,mdn_wetdep,col=1,psym=sym(1),pos=pos,yrange=[0,max([model_wetdep_annual_site,mdn_wetdep])]*1.2,$
;  	xtitle='Latitude (!uo!nN)',ytitle='Annual deposition ('+micro+'g m!u-2!ny!u-1!n)',$
; 	title='Latitudinal variation of!c annual deposition'
;  oplot,latitude,mean(model_wetdep_plot,1),col=2,thick=4
;  oplot,mdn_data.lat,model_wetdep_annual_site,col=!myct.blue,psym=sym(1)
;  legend,label=['MDN sites','Model all domain','Model @ MDN sites'],lcol=[1,2,!myct.blue],$
;  	col=[1,2,!myct.blue],line=[-1,0,-1],symbol=[1,-1,1],boxcolor=-1,thick=4,$
;  	halign=0.1,valign=0.95,charsize=1.

 ;  now do a scatter plot of modeled and observed deposition
 ; multipanel,col=2,row=2,margin=[0.05,0.09],pos=pos,/noerase

;  plot the points in different regions with different colors
 multipanel,/advance,pos=pos
 plot,mdn_wetdep,model_wetdep_annual_site,col=1,xrange=[0,25],yrange=[0,25],$
 	xtitle='MDN deposition ('+micro+'g m!u-2!ny!u-1!n)',$
	ytitle='GEOS-Chem deposition ('+micro+'g m!u-2!ny!u-1!n)',$
	title='Annual deposition comparison',psym=sym(1),pos=pos,/xst,/yst,charsize=0.8,/iso
 
; for ipoint=0, n_elements(mdn_wetdep)-1 do begin
;    if model_wetdep_annual_site[ipoint]/mdn_wetdep[ipoint] gt 1.25  then $
;        xyouts, mdn_wetdep[ipoint], model_wetdep_annual_site[ipoint], /data, $
;               mdn_data[ipoint].id + string(mdn_wetdep[ipoint], format='(f5.2)') + ' ' +$
;                string(model_wetdep_annual_site[ipoint], format='(f5.2)') ,  $
;               col=1,charsize=0.6
; endfor

if keyword_set(region) then begin

  legend_lable=[]
  legend_lcol=[]

; override the points in specified regions with specified colors
ind_other=indgen(n_elements(mdn_wetdep))
for iregion=0,nregion-1 do begin
   ind=where(mdn_data.lat ge limit_all(iregion).y[0] and mdn_data.lat le limit_all(iregion).y[2]  and $
                     mdn_data.lon ge limit_all(iregion).x[0] and mdn_data.lon le limit_all(iregion).x[2] )
   for indind=0, n_elements(ind)-1 do begin
         ind_other[where(ind_other eq ind[indind])]=-1
   endfor
   oplot,mdn_wetdep[ind],model_wetdep_annual_site[ind],col=2+iregion,psym=sym(1)
   legend_lable=[legend_lable, region_all[iregion] + ' '  + string(correlate(mdn_wetdep[ind],model_wetdep_annual_site[ind]), $
                          format= '("r = ",f4.2)' )    + ' ' + $
 		                   string(  mean( (model_wetdep_annual_site[ind]-mdn_wetdep[ind])/mdn_wetdep[ind] )*100.  , format='(F6.2,"%")') ]
   legend_lcol=[legend_lcol,2+iregion]
endfor 

if max(ind_other) gt 0 then begin
ind_other=ind_other(where(ind_other ge 0))

oplot,mdn_wetdep[ind_other],model_wetdep_annual_site[ind_other],col=2+nregion,psym=sym(1)
   legend_lable=[legend_lable, 'OT' + ' '  + string(correlate(mdn_wetdep[ind_other],model_wetdep_annual_site[ind_other]), $
                          format= '("r = ",f4.2)' )]
  ; string(n_elements(ind_other),format='(I2)'  )
   legend_lcol=[legend_lcol,2+nregion]
 legend,label=legend_lable,color=legend_lcol,line=-1,symbol=1+intarr(nregion+1),boxcolor=-1,thick=4,$
 	halign=-0.06,valign=0.98,charsize=0.7
endif else begin
 legend,label=legend_lable,color=legend_lcol,line=-1,symbol=1+intarr(nregion),boxcolor=-1,thick=4,$
 	halign=-0.06,valign=0.98,charsize=0.7
endelse

endif

; if keyword_set(version2) then oplot, mdn_wetdep,model_wetdep_annual_site2,col=2,psym=sym(1)
; if keyword_set(version3) then oplot, mdn_wetdep,model_wetdep_annual_site3,col=3,psym=sym(1)
; if keyword_set(version4) then oplot, mdn_wetdep,model_wetdep_annual_site4,col=4,psym=sym(1)
 
 ;  find AK98 - Kodiak
; qak=where(mdn_data.id eq '"AK98"')
; if qak(0) gt -1 then begin
; 	plots,mdn_wetdep(qak),model_wetdep_annual_site(qak),col=2,psym=sym(1),symsize=1.5
; 	arrow,mdn_wetdep(qak)+7,model_wetdep_annual_site(qak),$
;	      mdn_wetdep(qak)+1,model_wetdep_annual_site(qak),col=2,hthick=4,thick=4,/data
;        xyouts,mdn_wetdep(qak)+7.2,model_wetdep_annual_site(qak)-0.5,'Kodiak Island',/data,col=2
; endif

 rcorr=correlate(mdn_wetdep,model_wetdep_annual_site)
 org_corr,mdn_wetdep,model_wetdep_annual_site,rcorr,$
 	n_elements(mdn_wetdep),slope,intercept
 ; oplot,!x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
 oplot,!x.crange,!x.crange,col=1
 oplot,!x.crange,!x.crange*0.75,col=1,line=2
 oplot,!x.crange,!x.crange*1.25,col=1,line=2
 
  xyouts,!x.crange(1)*0.35,!y.crange(1)*0.05,$
          /data,col=1, 'Model = ' + string(intercept,slope,format= '(F5.2," + ", F4.2, " MDN")' ) +$
          '!cr = '+$
          string(rcorr,format='(F5.2)')+ ';  n = '+$
          string(n_elements(mdn_wetdep), format= '(I2)')+  $
;  	     '!cModel/Obs = ' + string(mean(model_wetdep_annual_site)/mean(mdn_wetdep), format= '(F5.2)')   + $ 
  		 '; Mean Bias = ' + strtrim(string(  mean( (model_wetdep_annual_site-mdn_wetdep)/mdn_wetdep )*100.  , format='(F5.1)'))  + '%'      $
; 		 '!cMean Abs. Bias (%) = ' +  string(   mean(  abs(model_wetdep_annual_site-mdn_wetdep) /mdn_wetdep) *100. , format='(F5.2)' ) + $  
;  		 '!cRMS (' + micro + 'gm!u-2!ny!u-1!n) = ' + string(   sqrt(mean((model_wetdep_annual_site-mdn_wetdep)^2))   , format='(F7.2)' ) + $ 
;          '!cMean Bias (' + micro + 'gm!u-2!ny!u-1!n) =' + string(   mean( model_wetdep_annual_site-mdn_wetdep)   , format='(F7.2)' )  $
  		   ,charsize=0.7

; plot vertical distribution of wet deposition flux
if keyword_set(vertical) then begin
; first plot an annual map
; the whole model grid
;   indlat=where( latitude ge limit[0] and latitude le limit[2])
;   indlon=where( longitude ge limit[1] and longitude le limit[3] )
geos_z_47=[0.071,0.201,0.332,0.466,0.601,0.737,0.875,1.016,1.157,1.301,1.447,1.594,1.769,1.999,2.259,$
2.527,2.801,3.084,3.448,3.904,4.382,4.886,5.419,5.985,6.591,7.241,7.947,8.848,9.938,11.021,12.086,$
13.134,14.17,15.198,16.222,17.243,18.727,20.836,23.02,25.307,28.654,34.024,40.166,47.135,54.834,63.054,$
72.18]

geos_dz_47=[12.189, 8.466,8.138,7.502,6.712,5.858,5.12,2.356,2.236,2.149,2.082,1.02,1.022,1.025,1.032,1.041,1.055,1.073,1.094,$
1.087,0.737,0.677,0.627,0.585,0.549,0.518,0.491,0.467,0.445,0.287,0.278,0.271,0.264,0.257,0.202,0.149,0.146,0.145,0.143,0.141,$
0.139,0.137,0.136,0.134,0.133,0.131,0.129]
geos_dz_47=reverse(geos_dz_47)

   model_vwetdep =total(wetdep_all.vrgm  +  wetdep_all.vhgp  ,4)         ; Mg/yr
   model_vwetdep_plot=model_vwetdep
   model_vwetdep_plot = total(model_vwetdep_plot,1)
   model_vwetdep_plot = total(model_vwetdep_plot,1)

 multipanel, col=1, row=2, pos=pos,omargin=0.03,margin=[0.3,0.04]
   plot, model_vwetdep_plot/geos_dz_47, geos_z_47, xtitle='Wet deposition, Mg/km yr', ytitle='Altitude, km', $
          title='vetical distribution of wet deposition', pos=pos,col=1,yrange=[0,15], xrange=[0, 120],charsize=0.8

   model_vwetdep =total(wetdep_all.vls, 4 )         ; Mg/yr
   model_vwetdep_plot=model_vwetdep
   model_vwetdep_plot = total(model_vwetdep_plot,1)
   model_vwetdep_plot = total(model_vwetdep_plot,1)

   oplot, model_vwetdep_plot/geos_dz_47, geos_z_47,col=!myct.red

   model_vwetdep =total(wetdep_all.vrgm  ,4)         ; Mg/yr
   model_vwetdep_plot=model_vwetdep
   model_vwetdep_plot = total(model_vwetdep_plot,1)
   model_vwetdep_plot = total(model_vwetdep_plot,1)

   oplot, model_vwetdep_plot/geos_dz_47, geos_z_47,col=!myct.green

   legend,label=['Total', 'Large scale', 'RGM'],lcolor=[1, !myct.red, !myct.green],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=4,$
 	halign=0.97,valign=0.95,charsize=0.7

 ; plot for each season
 if keyword_set(seasonal) then begin
    multipanel, col=2, row=2, pos=pos,omargin=0.03,margin=[0.05,0.04]
    for iseason=0,3 do begin
           model_vwetdep =total(wetdep_all(season(iseason).month).vrgm  +  wetdep_all(season(iseason).month).vhgp  ,4)       ; Mg/season
           model_vwetdep_plot=model_vwetdep
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           plot, model_vwetdep_plot/geos_dz_47, geos_z_47, xtitle='Wet deposition, Mg/km season', ytitle='Altitude, km', $
                   title='vetical distribution of !cwet deposition by season ' + season(iseason).name, pos=pos,col=1,yrange=[0,15],xrange=[0, 40],charsize=0.8

           model_vwetdep =total(wetdep_all(season(iseason).month).vls ,4)       ; Mg/season
           model_vwetdep_plot=model_vwetdep
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           oplot, model_vwetdep_plot/geos_dz_47, geos_z_47, col=!myct.red

           model_vwetdep =total(wetdep_all(season(iseason).month).vrgm ,4)       ; Mg/season
           model_vwetdep_plot=model_vwetdep
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           oplot, model_vwetdep_plot/geos_dz_47, geos_z_47, col=!myct.green

           legend,label=['Total', 'Large scale', 'RGM'],lcolor=[1, !myct.red, !myct.green],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=4,$
 	       halign=0.97,valign=0.95,charsize=0.7

           multipanel,/advance, pos=pos
    endfor
 endif

 if keyword_set(region) then begin
    multipanel, col=2, row=2, pos=pos,omargin=0.03,margin=[0.05,0.04]
    for iregion=0,nregion-1 do begin
           indlat=where( latitude ge limit_all(iregion).y[0] and latitude le limit_all(iregion).y[2])
           indlon=where( longitude ge limit_all(iregion).x[0] and longitude le limit_all(iregion).y[2] )

           model_vwetdep =total(wetdep_all.vrgm  +  wetdep_all.vhgp  ,4)         ; Mg/yr
           model_vwetdep_plot=model_vwetdep
           model_vwetdep_plot = model_vwetdep_plot(*,indlat,*)  ; Mg/yr
           model_vwetdep_plot = model_vwetdep_plot(indlon,*,*)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           plot, model_vwetdep_plot/geos_dz_47, geos_z_47, xtitle='Wet deposition, Mg/km yr', ytitle='Altitude, km', $
                   title='vetical distribution of !cwet deposition by region ' + region_all(iregion), pos=pos,col=1,yrange=[0,15],xrange=[0, 30],charsize=0.8

           model_vwetdep =total(wetdep_all.vls  ,4)         ; Mg/yr
           model_vwetdep_plot=model_vwetdep
           model_vwetdep_plot = model_vwetdep_plot(*,indlat,*)  ; Mg/yr
           model_vwetdep_plot = model_vwetdep_plot(indlon,*,*)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           oplot, model_vwetdep_plot/geos_dz_47, geos_z_47, col=!myct.red

           model_vwetdep =total(wetdep_all.vrgm  ,4)         ; Mg/yr
           model_vwetdep_plot=model_vwetdep
           model_vwetdep_plot = model_vwetdep_plot(*,indlat,*)  ; Mg/yr
           model_vwetdep_plot = model_vwetdep_plot(indlon,*,*)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           model_vwetdep_plot = total(model_vwetdep_plot,1)
           oplot, model_vwetdep_plot/geos_dz_47, geos_z_47, col=!myct.green

           legend,label=['Total', 'Large scale', 'RGM'],lcolor=[1, !myct.red, !myct.green],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=4,$
 	       halign=0.97,valign=0.95,charsize=0.7

           multipanel,/advance, pos=pos
    endfor
 endif
endif
; end of plot vertical distribution

 ;  define array that will contain the seasonal total deposition for each region
 if keyword_set(region) then begin
 	mdn_region_season=fltarr(nregion,4)
 	mdn_region_season_std=fltarr(nregion,4)
 	mod_region_season=fltarr(nregion,4)
 	mdn_region_annual=fltarr(nregion)
 	mdn_region_annual_std=fltarr(nregion)
 	mod_region_annual=fltarr(nregion)
 	mod_region_annual_std=fltarr(nregion)
 	mod_region_annual_dry=fltarr(nregion)

	if keyword_set(version2) then begin
 	  mod_region_season2=fltarr(nregion,4)
 	  mod_region_annual2=fltarr(nregion)
 	  mod_region_annual_std2=fltarr(nregion)
 	  mod_region_annual_dry2=fltarr(nregion)
	endif
	if keyword_set(version3) then begin
 	  mod_region_season3=fltarr(nregion,4)
 	  mod_region_annual3=fltarr(nregion)
 	  mod_region_annual_std3=fltarr(nregion)
 	  mod_region_annual_dry3=fltarr(nregion)
	endif
	if keyword_set(version4) then begin
 	  mod_region_season4=fltarr(nregion,4)
 	  mod_region_annual4=fltarr(nregion)
 	  mod_region_annual_std4=fltarr(nregion)
 	  mod_region_annual_dry4=fltarr(nregion)
	endif
	
	endif

;  plot monthly wet deposition maps
if keyword_set(monthly) then begin
 multipanel,col=3,row=4,margin=[0.001,0.001]
 imonth=['Jan','Fen','Mar','Apr','May','June','July','August','September','October','November','December']
 cbar=fltarr(12)
 cbar[10]=1
 cbunit=strarr(12)
 cbunit[10]= 'ng m!u-2!n d!u-1!n'
 for k=0,11 do $
 ctm_overlay,(wetdep_all(k).rgm  +  wetdep_all(k).hgp)*convert*1e3/ndays(k),longitude,latitude,$
       mdn_data.monthly_dep(k),mdn_data.lon,mdn_data.lat,$
	/sample,/continents,thick=1, $
	div=COLORBAR_NDIV(maxdiv=10),$
	mindata=0.,maxdata=100.0,$
	/usa,limit=limit,cbar=cbar(k),$
	t_symbol=1,symsize=0.5,title=imonth(k)+ ' Mercury Wet Deposition Flux'+$
      '!cMDN data for '+select_year_str+  '('+year_str+ ')',$
        cbunit=cbunit[k] , /iso,csfac=0.5,margin=0.01,/nogx,/nogy
endif ;  end /monthly

;  plot seasonal maps of wet deposition (Note: need to check that the calculation is 
;  done correctly...)
if keyword_set(seasonal) then begin
 multipanel,col=2,row=2,margin=[0.01,0.01]
 for k=0,3 do $
 ctm_overlay,mean(wetdep_all(season(k).month).rgm + wetdep_all(season(k).month).hgp,3)*convert*12/4.,longitude,latitude,$
        mdn_data.season_dep(k)*1e-3*365./4.,mdn_data.lon,mdn_data.lat,$
	/sample,/continents,/cbar, thick=1, $
	div=COLORBAR_NDIV(maxdiv=10),$
	mindata=0.,maxdata=6,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,title=season(k).name+' Mercury Wet Deposition Flux'+$
       '!cMDN data for '+select_year_str,$
        cbunit=micro+'g m!u-2!n season!u-1!n',/iso,csfac=0.8,/nogx,/nogy,margin=0.01

;multipanel,col=3,row=4,margin=[0.001,0.001]
 model_data_sites=fltarr(n_elements(mdn_data.season_dep(0)),4)
  for k=0,3 do begin
;  maxdata=5
;tvmap,mean(wetdep_all(season(k).month).rgm + wetdep_all(season(k).month).hgp,3)*convert*12/4.,longitude,latitude,$
;	/sample,/continents, thick=1, $
;	mindata=0,maxdata=maxdata,$
;	/usa,limit=limit,$
;	t_symbol=1,symsize=0.9, $
;        /iso,csfac=0.6,margin=0.01
;
;; plot the model values at MDN sites
 nsites_tmp = n_elements( mdn_data.season_dep(k) )
 for isite_tmp=0, nsites_tmp-1 do begin
     ctm_index,ctm_type('GEOS5_47L',res=res),i,j,$
	        center=[mdn_data[isite_tmp].lat,mdn_data[isite_tmp].lon],/non_interactive
     indlat=where( latitude eq lat_all(j-1) )
     indlon=where( longitude eq lon_all(i-1) )
     model_data_sites[isite_tmp,k]=(mean(wetdep_all(season(k).month).rgm + wetdep_all(season(k).month).hgp,3)*convert*12/4.)[indlon,indlat]
 endfor
; ctm_overlay,0.*mean(wetdep_all(season(k).month).rgm + wetdep_all(season(k).month).hgp,3)*convert*12/4.,longitude,latitude,$
;        model_data_sites[*,k],mdn_data.lon,mdn_data.lat,$
;	/sample,/continents,  $
;       /cbar, div=COLORBAR_NDIV(maxdiv=10),$
;        thick=1, $
;	mindata=0,maxdata=maxdata,$
;	/usa,limit=limit,$
;	t_symbol=1,symsize=0.9, $; title=season(k).name+' Observation',$
;        cbunit=micro+'g m!u-2!n season!u-1!n',$
;       /iso,csfac=0.6,margin=0.01, /nogy
;
; ctm_overlay,0.*mean(wetdep_all(season(k).month).rgm + wetdep_all(season(k).month).hgp,3)*convert*12/4.,longitude,latitude,$
;        mdn_data.season_dep(k)*1e-3*365./4.,mdn_data.lon,mdn_data.lat,$
;	/sample,/continents,  $
;     ;  /cbar, div=COLORBAR_NDIV(maxdiv=10),$
;        thick=1, $
;	mindata=0,maxdata=maxdata,$
;	/usa,limit=limit,$
;	t_symbol=1,symsize=0.9, $; title=season(k).name+' Observation',$
; ;       cbunit=micro+'g m!u-2!n season!u-1!n',$
;       /iso,csfac=0.6,margin=0.01, /nogy
;
  endfor

;    legend,label=[string(year[0],format='(I4)')], lcol=[1], line=[-1],symbol=[-1],boxcolor=-1,smsize=.8, $
;                           halign=-0.03,valign=-0.63,charsize=1.

; plot a series of scatter plot for each season
multipanel,col=2,row=2,margin=[0.1,0.1],pos=pos
seasonstring=['winter','spring','summer','autumn']
for k=0,3 do begin
 xdata=mdn_data.season_dep(k)*1e-3*365./4.
 ydata=model_data_sites[*,k]
 plot,xdata,ydata,col=1,xrange=[0,1.1*max([xdata,ydata])],yrange=[0,1.1*max([xdata,ydata])],$
 	xtitle= 'MDN deposition (' + micro+ 'g m!u-2!nseason!u-1!n)', $
	ytitle= 'GEOS-Chem deposition (' +micro+ 'g m!u-2!nseason!u-1!n)', $
	title = seasonstring[k] + ' deposition comparison',psym=sym(1),pos=pos,/xst,/yst,charsize=0.8, /iso, xminor=10, yminor=10

; override the points in specified regions with specified colors
if keyword_set(region) then begin

  legend_lable=[]
  legend_lcol=[]
ind_other=indgen(n_elements(xdata))
for iregion=0,nregion-1 do begin
   ind=where(mdn_data.lat ge limit_all(iregion).y[0] and mdn_data.lat le limit_all(iregion).y[2]  and $
                     mdn_data.lon ge limit_all(iregion).x[0] and mdn_data.lon le limit_all(iregion).x[2] )
   for indind=0, n_elements(ind)-1 do begin
         ind_other[where(ind_other eq ind[indind])]=-1
   endfor
   oplot,xdata[ind],ydata[ind],col=2+iregion,psym=sym(1)
   legend_lable=[legend_lable, region_all[iregion] + '('  + string(correlate(xdata[ind],ydata[ind]), $
                          format= '(f4.2)' )  + ')']
;                          + ' ' + string(n_elements(ind),format='(I2)'  ) + ' ' + $
;                          string(  100.0*mean((ydata[ind] - xdata[ind]) / xdata[ind]) , format='(f6.2,"%")' ) ]
   legend_lcol=[legend_lcol,2+iregion]
endfor 

if max(ind_other) gt 0 then begin
ind_other=ind_other(where(ind_other ge 0))

oplot,xdata[ind_other],ydata[ind_other],col=2+nregion,psym=sym(1)
   legend_lable=[legend_lable, 'OT' + '('  + string(correlate(xdata[ind_other],ydata[ind_other]), $
                          format= '(f4.2)' )  + ')' ]
;                         + ' ' + string(n_elements(ind_other),format='(I2)'  ) ]
   legend_lcol=[legend_lcol,2+nregion]
 legend,label=legend_lable,color=legend_lcol,line=-1,symbol=1+intarr(nregion+1),boxcolor=-1,thick=4,$
 	halign=-0.1,valign=0.73,charsize=0.75
endif else begin
 legend,label=legend_lable,color=legend_lcol,line=-1,symbol=1+intarr(nregion),boxcolor=-1,thick=4,$
 	halign=-0.1,valign=0.73,charsize=0.75
endelse

endif

 rcorr=correlate(xdata,ydata)
 org_corr,xdata,ydata,rcorr,$
 	n_elements(xdata),slope,intercept
; oplot,!x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
 oplot,!x.crange,!x.crange,col=1
 oplot,!x.crange,!x.crange*0.75,col=1,line=2
 oplot,!x.crange,!x.crange*1.25,col=1,line=2
 
 xyouts,!x.crange(1)/20.,!y.crange(1)*0.92,$
         /data,col=1, 'Model = ' + string(intercept,slope,format= '(F4.2," + ", F4.2, " MDN")' ) +$
         '!cr= '+$
         string(rcorr,format='(F5.2)')+ ';  n= '+$
         string(n_elements(xdata), format= '(I2)')  + $
; 	     '!cModel/Obs = ' + string(mean(ydata)/mean(xdata), format= '(F5.2)')   + $ 
		 '!cMean Bias = ' + strtrim(string(  mean( (ydata-xdata)/xdata )*100.  , format='(F6.1)')) + '%'        $
;		 '!cMean Abs. Bias (%) = ' +  string(   mean(  abs(ydata-xdata) /xdata) *100. , format='(F5.2)' ) + $  
;		 '!cRMS (' + micro + 'gm!u-2!ny!u-1!n) = ' + string(   sqrt(mean((ydata-xdata)^2))   , format='(F7.2)' ) + $ 
;         '!cMean Bias (' + micro + 'gm!u-2!ny!u-1!n) =' + string(   mean( ydata-xdata)   , format='(F7.2)' )  $
		 , charsize=0.8

 multipanel,/advance,pos=pos
endfor
endif

 ;  plot the seasonal cycle of monthly deposition for each separate region
;  if keyword_set(region) then maxdata = findgen(n_elements(region_all))*0.0+80.0 $
;     else maxdata=[80, 50, 80, 120]
if keyword_set(region) then $
if nregion eq 4 then maxdata=[80, 80, 80, 120] else maxdata = findgen(n_elements(region_all))*0.0+120.0

yrange=[0,120]

if keyword_set(region) then begin
 multipanel,col=4,row=5,margin=[0.02,0.03],pos=pos
 !p.charsize=1.5

 nsites_regionall = intarr(nregion)

; plot the multi-year average data
 for k=0,nregion-1 do begin
    mdn_data=mdn_data_all ;  reset data to the whole domain
;    ind=where(mdn_data.lat ge limit_all(k).y(0)-dlat and mdn_data.lat lt limit_all(k).y(2)+dlat and $
;              mdn_data.lon ge limit_all(k).x(0)-dlon and mdn_data.lon lt limit_all(k).x(1)+dlon and $
;   	      mdn_data.npts ge min_npts,nsites)
    ind=where(mdn_data.lat ge limit_all(k).y(0) and mdn_data.lat lt limit_all(k).y(2) and $
              mdn_data.lon ge limit_all(k).x(0) and mdn_data.lon lt limit_all(k).x(2) and $
   	          min(mdn_data.npts,dimension=1) ge min_npts, nsites)
    nsites_regionall(k) = nsites

    if dimension(mdn_data(ind).monthly_dep) eq 2 then begin
          mdn_monthly=mean(mdn_data(ind).monthly_dep,2,/nan)
          if keyword_set(sn0w) then mdn_monthly0=mean(mdn_data0(ind).monthly_dep,2,/nan) 
    endif else begin
          mdn_monthly=mdn_data(ind).monthly_dep
          if keyword_set(sn0w) then mdn_monthly0=mdn_data0(ind).monthly_dep
    endelse

    mdn_monthly_std=mdn_monthly*0.
    ;  calculate seasonal mean wet deposition mdn_data.season_dep is in ng/m2/day --> convert
    ;  to ug/m2/season (where season=total number of days in a season - assume equally weighted months)
    if n_elements(ind) gt 1 then begin
       mdn_region_season(k,*)=mean(mdn_data(ind).season_dep,2,/nan)*1e-3*365./4.
       mdn_region_season_std(k,*)=mean(mdn_data(ind).season_dep_std,2,/nan)*1e-3*365./4.
       mod_region_season(k,*)=mean(model_wetdep_seasonal_site(*,ind),2,/nan)
    endif else begin
       mdn_region_season(k,*)=mdn_data(ind).season_dep*1e-3*365./4.
       mdn_region_season_std(k,*)=mdn_data(ind).season_dep_std*1e-3*365./4.
       mod_region_season(k,*)=model_wetdep_seasonal_site(*,ind)
    endelse
    mdn_region_annual(k)=mean(mdn_wetdep_annual_site(ind),/nan)
    mdn_region_annual_std(k)=mean(mdn_wetdep_annual_site_std(ind),/nan)
    mod_region_annual(k)=mean(model_wetdep_annual_site(ind),/nan)
    mod_region_annual_std(k)=stddev(model_wetdep_annual_site(ind),/nan)
    mod_region_annual_dry(k)=mean(model_drydep_annual_site(ind),/nan)
    if keyword_set(version2) then begin
       mod_region_season2(k,*)=mean(model_wetdep_seasonal_site2(*,ind),2,/nan)
       mod_region_annual2(k)=mean(model_wetdep_annual_site2(ind),/nan)
       mod_region_annual_std2(k)=stddev(model_wetdep_annual_site2(ind),/nan)
       mod_region_annual_dry2(k)=mean(model_drydep_annual_site2(ind),/nan)
    endif
    if keyword_set(version3) then begin
       mod_region_season3(k,*)=mean(model_wetdep_seasonal_site3(*,ind),2,/nan)
       mod_region_annual3(k)=mean(model_wetdep_annual_site3(ind),/nan)
       mod_region_annual_std3(k)=stddev(model_wetdep_annual_site3(ind),/nan)
       mod_region_annual_dry3(k)=mean(model_drydep_annual_site3(ind),/nan)
    endif
    if keyword_set(version4) then begin
       mod_region_season4(k,*)=mean(model_wetdep_seasonal_site4(*,ind),2,/nan)
       mod_region_annual4(k)=mean(model_wetdep_annual_site4(ind),/nan)
       mod_region_annual_std4(k)=stddev(model_wetdep_annual_site4(ind),/nan)
       mod_region_annual_dry4(k)=mean(model_drydep_annual_site4(ind),/nan)
    endif
	
    for i=0,12-1 do mdn_monthly_std(i)=stddev(mdn_data(ind).monthly_dep(i),/nan)
    ; yrange=[0,120]
    if region_all(k) eq 'Florida' then yrange=[0,150]

    ; debug
;	if k eq 0 then print,  'MW MDN ', mdn_monthly
;    if k eq 0 then print, 'MW model', mean(model_wetdep_monthly_site(*,ind),2,/nan)

    plot,1+indgen(12),mdn_monthly,col=1,pos=pos,$
 	ytitle='Monthly dep. (ng m!u-2!nd!u-1!n)',thick=5,yrange=[0,maxdata[k]],$
 	xrange=[0.4,12.6],/xst,title='Region: '+region_all(k),$
	xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
        xtickv=indgen(12)+1,xminor=1,charsize=1.1
    polyfill,[1+indgen(12),reverse(1+indgen(12)),1],$
       [mdn_monthly-mdn_monthly_std,reverse(mdn_monthly+mdn_monthly_std),mdn_monthly(0)-mdn_monthly_std(0)],col=!myct.lightgray
    ; oplot,1+indgen(12),mdn_monthly0,col=!myct.orange,thick=5
    oplot,1+indgen(12),mdn_monthly,col=1,thick=5
    xyouts, 1, !y.crange(1)*0.70,string(nsites,format='(I4," sites")'),col=1,al=0,charsize=0.6
    if n_elements(ind) gt 1 then oplot,1+indgen(12),mean(model_wetdep_monthly_site(*,ind),2,/nan),col=!myct.red,thick=5 $
                                       else oplot,1+indgen(12),model_wetdep_monthly_site(*,ind),col=!myct.red,thick=5
    if keyword_set(version2) then  oplot,1+indgen(12),mean(model_wetdep_monthly_site2(*,ind),2,/nan),col=!myct.green,thick=5
	if keyword_set(version3) then  oplot,1+indgen(12),mean(model_wetdep_monthly_site3(*,ind),2,/nan),col=!myct.green,thick=5,linestyle=2
    if keyword_set(version4) then  oplot,1+indgen(12),mean(model_wetdep_monthly_site4(*,ind),2,/nan),col=!myct.green,thick=5,linestyle=3

	if resolution eq '4x5' then res_str='4' else res_str='05'
    if n_elements(ind) gt 1 then $
	label_accumulate=['MDN' + $
      strtrim(string(mean( 1e-3*365.25 * mdn_monthly,/nan),stddev(1e-3*365.25 *mdn_monthly,/nan),format= '(F5.2,"+/-",F5.2)') )   + ' ' + micro + 'g m!u-2!ny!u-1!n'           , $
;                                  'MDN (original)'  + $
;      strtrim(string(mean(mdn_monthly0,/nan),stddev(mdn_monthly0,/nan),format= '(F4.1,"+/-",F4.1)') )   + ' ng m!u-2!nd!u-1!n'           , $
                                 'NEW simulation'  + $
      strtrim( string(mean(mean(1e-3*365.25 *model_wetdep_monthly_site(*,ind),2,/nan),/nan),stddev(mean(1e-3*365.25 *model_wetdep_monthly_site(*,ind),2,/nan),/nan),format= '(F5.2,"+/-",F5.2)' )  ) ]  $
    else $
	label_accumulate=['MDN' + $
      strtrim(string(mean( 1e-3*365.25 * mdn_monthly,/nan),stddev(1e-3*365.25 *mdn_monthly,/nan),format= '(F5.2,"+/-",F5.2)') )   + ' ' + micro + 'g m!u-2!ny!u-1!n'           , $
;                                  'MDN (original)'  + $
;      strtrim(string(mean(mdn_monthly0,/nan),stddev(mdn_monthly0,/nan),format= '(F4.1,"+/-",F4.1)') )   + ' ng m!u-2!nd!u-1!n'           , $
                                 'NEW simulation'  + $
      strtrim( string(mean(1e-3*365.25 *model_wetdep_monthly_site(*,ind),/nan),stddev(1e-3*365.25 *model_wetdep_monthly_site(*,ind),/nan),format= '(F5.2,"+/-",F5.2)' )  ) ] 

	lcol_accumulate=[1,!myct.red]
	line_accumulate=[0,0]
	symbol_accumulate=[-1,-1]
	
	if keyword_set(version2) then begin
;	   if resolution2 eq '4x5' then res_str='4' else res_str='05'
	   label_accumulate=[label_accumulate, 'OLD simulation'  + $
           ' '+ strtrim(string(mean(mean(1e-3*365.25 *model_wetdep_monthly_site2(*,ind),2,/nan),/nan),stddev(mean(1e-3*365.25 *model_wetdep_monthly_site2(*,ind),2,/nan),/nan),format= '(F5.2,"+/-",F5.2)' )  )       ]
	   lcol_accumulate=[lcol_accumulate,!myct.green]
	   line_accumulate=[line_accumulate,0]
	   symbol_accumulate=[symbol_accumulate,-1]
 	  if keyword_set(version3) then begin
	   if resolution3 eq '4x5' then res_str='4' else res_str='05'
	   label_accumulate=[label_accumulate, 'GC'+res_str+  ' w/o inplume red. ' + strtrim( string(mean(mean(model_wetdep_monthly_site3(*,ind),2,/nan),/nan),stddev(mean(model_wetdep_monthly_site3(*,ind),2,/nan),/nan),format= '(F4.1,"+/-",F4.1)' )  )      ]
	      lcol_accumulate=[lcol_accumulate,!myct.green]
	      line_accumulate=[line_accumulate,2]
	      symbol_accumulate=[symbol_accumulate,-1]
		  	  if keyword_set(version4) then begin
			      if resolution4 eq '4x5' then res_str='4' else res_str='05'
	              label_accumulate=[label_accumulate, 'GC'+res_str+  ' high inplume red. ' + strtrim( string(mean(mean(model_wetdep_monthly_site4(*,ind),2,/nan),/nan),stddev(mean(model_wetdep_monthly_site4(*,ind),2,/nan),/nan),format= '(F4.1,"+/-",F4.1)' )  )      ]
	              lcol_accumulate=[lcol_accumulate,!myct.green]
	              line_accumulate=[line_accumulate,3]
	              symbol_accumulate=[symbol_accumulate,-1]
	          endif
	  endif
   endif
	
   legend,label=label_accumulate,lcol=lcol_accumulate,line=line_accumulate,symbol=symbol_accumulate,boxcolor=-1,thick=5,$
 	         halign=0.01,valign=0.95,charsize=0.6
   
	multipanel,/advance,pos=pos
 endfor
 
; ==========================================================================================================================
; plot the precipitation data for each region
if keyword_set(region) then begin
yrange=[0,30]
 !p.charsize=1.5
 mdn_precipitation_monthly_std=fltarr(12)

 for k=0,nregion-1 do begin
    ind=where(mdn_data.lat ge limit_all(k).y(0) and mdn_data.lat lt limit_all(k).y(2) and $
              mdn_data.lon ge limit_all(k).x(0) and mdn_data.lon lt limit_all(k).x(2) and $
   	      min(mdn_data.npts,dimension=1) ge min_npts,nsites_region)

    if dimension(mdn_data(ind).monthly_precipitation) eq 2 then begin
          mdn_precipitation_monthly=mean(mdn_data(ind).monthly_precipitation,2,/nan)
          for imonth=0,11 do begin
               mdn_precipitation_monthly_std[imonth]=stddev(mdn_data(ind).monthly_precipitation[imonth])
          endfor
          mdn_precipitation_monthly_gctt=mean(mdn_data(ind).monthly_precipitation_gctt,2,/nan)
          mdn_precipitation_monthly_gccv=mean(mdn_data(ind).monthly_precipitation_gccv,2,/nan)
    endif else begin
          mdn_precipitation_monthly=mdn_data(ind).monthly_precipitation
    endelse
	; transfer mm/day to cm/month

    ; debug
	if k eq 0 then print,  'MW MDN ', 3*mdn_precipitation_monthly
    if k eq 0 then print, 'MW model prec ', 0.1*mdn_precipitation_monthly_gctt

    plot,1+indgen(12),3*mdn_precipitation_monthly,col=1,pos=pos,$
 	ytitle='Monthly prec. (cm month!u-1!n)',thick=5,yrange=yrange,$
 	xrange=[0.4,12.6],/xst,title='Region: '+region_all(k),$
	xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
        xtickv=indgen(12)+1,xminor=1,charsize=1.1

    polyfill,[1+indgen(12),reverse(1+indgen(12)),1],$
       3*[mdn_precipitation_monthly-mdn_precipitation_monthly_std,reverse(mdn_precipitation_monthly+mdn_precipitation_monthly_std),$
           mdn_precipitation_monthly[0]-mdn_precipitation_monthly_std[0]],col=!myct.lightgray

     oplot,1+indgen(12),3*mdn_precipitation_monthly,col=1,thick=5    
     oplot,1+indgen(12),0.1*mdn_precipitation_monthly_gctt,col=!myct.blue,thick=5
     oplot,1+indgen(12),0.1*mdn_precipitation_monthly_gccv,col=!myct.purple,thick=5

     legend,label=['MDN'+string(mean(3*mdn_precipitation_monthly),format='(f5.2)')+ ' '+string(stddev(3*mdn_precipitation_monthly),format='(f5.2)'), $
                         mettype + ' total',+string(mean(0.1*mdn_precipitation_monthly_gctt),format='(f5.2)')+ ' '+string(stddev(0.1*mdn_precipitation_monthly_gctt),format='(f5.2)'), $
                         mettype + ' convective'],lcol=[!myct.black,!myct.blue,!myct.purple],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=5,$
 	         halign=0.01,valign=0.95,charsize=0.6
;   if k eq 0 then legend,label=['MDN', 'GEOS-5 total', 'GEOS-5 convective'],$
;             lcol=[!myct.black,!myct.blue,!myct.purple],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=5,$
; 	         halign=0.01,valign=0.95,charsize=0.6

	multipanel,/advance,pos=pos
endfor

endif
; ===========================================================================================================================

  ;  now plot seaonal means by region
   multipanel,col=1,row=2,margin=[0.15,0.05],pos=pos
   !p.thick=1
   !x.thick=1 & !y.thick=1
   ; sub=indgen(4) & nsub=4 ;  for first 4 resgions
   sub=indgen(nregion) & nsub=nregion ;  for all regions
   ; sub=indgen(2) & nsub=2 ;  SE, Florida, and Gulf

   width=1./5. 
   yrange=[0,8]
   if nsub eq 2 then yrange=[0,13]
   !p.charsize=2
   !p.charsize=1.7
 ;  bargraph,mdn_region_season(sub,0),barcolor=!myct.gray67, barXoffset=-1.5*width, $
 ;  	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),yrange=yrange,$
;	ytitle= 'Deposition (' + micro+ 'g m!u-2!n season!u-1!n)',/xst,/yst, charsize=0.8
 ;  bargraph,mdn_region_season(sub,1),barcolor=!myct.LIGHTgreen,BARWIDTH=width,pos=pos,/overplot,BarXOffset=-0.5*width,/no_label
 ;  bargraph,mdn_region_season(sub,2),barcolor=!myct.red,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width*0.5,/no_label
 ;  bargraph,mdn_region_season(sub,3),barcolor=!myct.blue,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width*1.5,/no_label
 ;  legend,label=['Winter','Spring','Summer','Fall'],lcol=[!myct.gray67,!myct.LIGHTgreen,!myct.red,!myct.blue],$
 ;  	line=[0,0,0,0],symbol=[-1,-1,-1,-1],boxcolor=-1,thick=14,$
 ;	halign=0.01,valign=0.95,charsize=1.
   
 ;  multipanel,/advance,pos=pos
   
   bargraph,mod_region_season(sub,0),barcolor=!myct.gray67, barXoffset=-1.5*width, $
   	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),yrange=yrange,$
	ytitle= 'Deposition (' + micro + 'g m!u-2!n season!u-1!n)',/yst, charsize=0.8
   oplot,indgen(nsub)+1-1.5*width,mdn_region_season(sub,0),col=1,psym=sym(1),symsize=2.
   errorbar,indgen(nsub)+1-1.5*width,mdn_region_season(sub,0),mdn_region_season_std(sub,0),col=1,thick=4
   oplot,indgen(nsub)+1-1.5*width,mdn_region_season(sub,0),col=!myct.gray67,psym=sym(1),symsize=1.5

   bargraph,mod_region_season(sub,1),barcolor=!myct.LIGHTgreen,BARWIDTH=width,pos=pos,/overplot,BarXOffset=-0.5*width,/no_label
   oplot,indgen(nsub)+1-0.5*width,mdn_region_season(sub,1),col=1,psym=sym(1),symsize=2.
   errorbar,indgen(nsub)+1-0.5*width,mdn_region_season(sub,1),mdn_region_season_std(sub,1),col=1,thick=4
   oplot,indgen(nsub)+1-0.5*width,mdn_region_season(sub,1),col=!myct.lightgreen,psym=sym(1),symsize=1.5

   bargraph,mod_region_season(sub,2),barcolor=!myct.red,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width*0.5,/no_label
   oplot,indgen(nsub)+1+width*0.5,mdn_region_season(sub,2),col=1,psym=sym(1),symsize=2.
   errorbar,indgen(nsub)+1+width*0.5,mdn_region_season(sub,2),mdn_region_season_std(sub,2),col=1,thick=4
   oplot,indgen(nsub)+1+width*0.5,mdn_region_season(sub,2),col=!myct.red,psym=sym(1),symsize=1.5

   bargraph,mod_region_season(sub,3),barcolor=!myct.blue,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width*1.5,/no_label
   oplot,indgen(nsub)+1+width*1.5,mdn_region_season(sub,3),col=1,psym=sym(1),symsize=2.
   errorbar,indgen(nsub)+1+width*1.5,mdn_region_season(sub,3),mdn_region_season_std(sub,3),col=1,thick=4
   oplot,indgen(nsub)+1+width*1.5,mdn_region_season(sub,3),col=!myct.blue,psym=sym(1),symsize=1.5

   legend,label=['Winter','Spring','Summer','Fall!cCircle: MDN!cBar: Model'],$
      lcol=[!myct.gray67,!myct.LIGHTgreen,!myct.red,!myct.blue],$
   	line=[0,0,0,0],symbol=[-1,-1,-1,-1],boxcolor=-1,thick=14,$
 	halign=0.01,valign=0.95,charsize=1.

  ; if keyword_set(version2) then begin
  ;   multipanel,/advance,pos=pos
  ;   multipanel,/advance,pos=pos
  ;   us_mod_region_season=mod_region_season-mod_region_season2 ;  difference between std run and run without US emissions
  ; bargraph,mod_region_season(sub,0),barcolor=!myct.gray67,$
  ; 	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),yrange=yrange,$
	;ytitle= 'Deposition (' + micro + 'g m!u-2!n season!u-1!n)',/yst, charsize=0.8
  ; bargraph,us_mod_region_season(sub,0),barcolor=!myct.green,$
  ; 	BARWIDTH=width,pos=pos,$
	;/yst,/overplot,barlabel=us_mod_region_season(sub,0)/mod_region_season(sub,0)*100.,l_format='(I2,"%")',$
	;barcharsize=1.2
  ; oplot,indgen(nsub)+1,mdn_region_season(sub,0),col=1,psym=sym(1),symsize=2.
  ; errorbar,indgen(nsub)+1,mdn_region_season(sub,0),mdn_region_season_std(sub,0),col=1,thick=4
  ; oplot,indgen(nsub)+1,mdn_region_season(sub,0),col=!myct.gray67,psym=sym(1),symsize=1.5
  ; bargraph,mod_region_season(sub,1),barcolor=!myct.LIGHTgreen,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width,/no_label
  ; bargraph,us_mod_region_season(sub,1),barcolor=!myct.green,$
  ; 	BARWIDTH=width,pos=pos,/yst,/overplot,barlabel=us_mod_region_season(sub,1)/mod_region_season(sub,1)*100.,l_format='(I2,"%")',$
	;barcharsize=1.2,BarXOffset=width
  ; oplot,indgen(nsub)+1+width,mdn_region_season(sub,1),col=1,psym=sym(1),symsize=2.
  ; errorbar,indgen(nsub)+1+width,mdn_region_season(sub,1),mdn_region_season_std(sub,1),col=1,thick=4
  ; oplot,indgen(nsub)+1+width,mdn_region_season(sub,1),col=!myct.lightgreen,psym=sym(1),symsize=1.5
  ; bargraph,mod_region_season(sub,2),barcolor=!myct.red,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width*2.,/no_label
  ; bargraph,us_mod_region_season(sub,2),barcolor=!myct.green,$
  ; 	BARWIDTH=width,pos=pos,/yst,/overplot,barlabel=us_mod_region_season(sub,2)/mod_region_season(sub,2)*100.,l_format='(I2,"%")',$
	;barcharsize=1.2,BarXOffset=width*2.
  ; oplot,indgen(nsub)+1+width*2,mdn_region_season(sub,2),col=1,psym=sym(1),symsize=2.
  ; errorbar,indgen(nsub)+1+width*2,mdn_region_season(sub,2),mdn_region_season_std(sub,2),col=1,thick=4
  ; oplot,indgen(nsub)+1+width*2,mdn_region_season(sub,2),col=!myct.red,psym=sym(1),symsize=1.5
  ; bargraph,mod_region_season(sub,3),barcolor=!myct.blue,BARWIDTH=width,pos=pos,/overplot,BarXOffset=width*3.,/no_label
  ; bargraph,us_mod_region_season(sub,3),barcolor=!myct.green,$
  ; 	BARWIDTH=width,pos=pos,/yst,/overplot,barlabel=us_mod_region_season(sub,3)/mod_region_season(sub,3)*100.,l_format='(I2,"%")',$
	;barcharsize=1.2,BarXOffset=width*3.
  ; oplot,indgen(nsub)+1+width*3,mdn_region_season(sub,3),col=1,psym=sym(1),symsize=2.
  ; errorbar,indgen(nsub)+1+width*3,mdn_region_season(sub,3),mdn_region_season_std(sub,3),col=1,thick=4
  ; oplot,indgen(nsub)+1+width*3,mdn_region_season(sub,3),col=!myct.blue,psym=sym(1),symsize=1.5
  ; endif

   ;  plot annual mean
   multipanel,/advance,pos=pos
   width=0.4
   yrange=[0,20]
   if nsub eq 2 then yrange=[0,25]
   bargraph,mdn_region_annual(sub),barcolor=!myct.red,$
   	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),$
	ytitle= 'Wet deposition (' + micro+ 'g m!u-2!n year!u-1!n)',yrange=yrange,BarXOffset=-0.2, charsize=0.8
   errorbar,indgen(nsub)+1-0.2,mdn_region_annual(sub),mdn_region_annual_std(sub),col=1,thick=8
   ; errorbar,indgen(nsub)+1+0.05,mdn_region_annual(sub),mdn_region_annual_std(sub),col=!myct.red,thick=4
   ; oplot,indgen(nsub)+1+0.05,mdn_region_annual(sub),col=!myct.red,psym=sym(1),symsize=2
   bargraph,mod_region_annual(sub),barcolor=!myct.gray67,$
   	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),$
	ytitle='Wet deposition ('+micro+ 'g m!u-2!n year!u-1!n)',yrange=yrange,BarXOffset=0.2,/overplot
   errorbar,indgen(nsub)+1+0.2,mod_region_annual(sub),mod_region_annual_std(sub),col=1,thick=8
   legend,label=['MDN observations','GEOS-Chem model'],lcol=[!myct.red,!myct.gray67],$
   	line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=14,$
 	halign=0.01,valign=0.95,charsize=1.
		
   ;  plot annual mean (US vs global contributions)
 ;  if keyword_set(version2) then begin
 ;    us_mod_region_annual=mod_region_annual-mod_region_annual2 ;  difference between std run and run without US emissions
 ;    global_mod_region_annual=mod_region_annual2
 ;  multipanel,/advance,pos=pos
 ;  width=0.3
 ;  yrange=[0.,30.]
 ;  bargraph,mdn_region_annual(sub),barcolor=!myct.red,$
 ;  	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),$
;	ytitle='Wet/Dry deposition ('+micro+ 'g m!u-2!n year!u-1!n)',yrange=yrange,BarXOffset=-width,$
;	charsize=1.0	
 ;  errorbar,indgen(nsub)+1-width,mdn_region_annual(sub),mdn_region_annual_std(sub),col=1,thick=8
 ;  bargraph,mod_region_annual(sub),barcolor=!myct.gray67,$
 ;  	BARWIDTH=width,pos=pos,/no_label,$
;	BarXOffset=0.,/overplot
 ;  bargraph,us_mod_region_annual(sub),barcolor=!myct.blue,$
 ;  	BARWIDTH=width,pos=pos,xlabel=region_all(sub),$
;	BarXOffset=0.,/overplot,$
;	barlabel=us_mod_region_annual(sub)/mod_region_annual(sub)*100.,l_format='(f4.1,"%")',$
;	barcharsize=0.6
 ;  axis,xaxis=1,col=1,xcharsize=0.0001
 ;  ; errorbar,indgen(nsub)+1+0.,mod_region_annual(sub),mod_region_annual_std(sub),col=1,thick=8
 ;
 ;  ; plot dry deposition in the same plot
 ;  us_mod_region_annual_dry=mod_region_annual_dry-mod_region_annual_dry2 ;  difference between std run and run without US emissions
 ;  bargraph,mod_region_annual_dry(sub),barcolor=!myct.gray,$
 ;  	BARWIDTH=width,BarXOffset=width,pos=pos,/no_label,xlabel=region_all(sub),$
;	/overplot
 ;  bargraph,us_mod_region_annual_dry(sub),barcolor=!myct.purple,$
 ;  	BARWIDTH=width,BarXOffset=width,pos=pos,xlabel=region_all(sub),/overplot,$
;	barlabel=us_mod_region_annual_dry(sub)/mod_region_annual_dry(sub)*100.,l_format='(f4.1,"%")',$
;	barcharsize=0.6
 ;  axis,xaxis=1,col=1,xcharsize=0.0001
 ;  
 ;  legend,label=['Wet deposition, MDN  observations','Wet deposition, US contribution','Wet deposition, Global contribution',  $
 ;                                                                            'Dry deposition, US contribution', 'Dry deposition, Global contribution'],$
 ;  lcol=[!myct.red,!myct.blue,!myct.gray67,!myct.gray,!myct.purple],$
 ;  	line=[0,0,0,0,0],symbol=[-1,-1,-1,-1,-1],boxcolor=-1,thick=14,$
 ;	halign=0.01,valign=0.95,charsize=0.8
 ;
 ;  ;  also plot dry deposition (only model)
 ;   multipanel,/advance,pos=pos
 ;  us_mod_region_annual_dry=mod_region_annual_dry-mod_region_annual_dry2 ;  difference between std run and run without US emissions
 ;  bargraph,mod_region_annual_dry(sub),barcolor=!myct.gray,$
 ;  	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),$
;	ytitle= 'Dry Deposition (' + micro+ 'g m!u-2!n year!u-1!n)',yrange=[0,30]
 ;  bargraph,us_mod_region_annual_dry(sub),barcolor=!myct.blue,$
 ;  	BARWIDTH=width,pos=pos,xlabel=region_all(sub),/overplot,$
;	barlabel=us_mod_region_annual_dry(sub)/mod_region_annual_dry(sub)*100.,l_format='(f4.1,"%")',$
;	barcharsize=1.2
 ;  axis,xaxis=1,col=1,xcharsize=0.0001
 ;
 ;  ;  and total deposition
 ;  multipanel,/advance,pos=pos
 ;  us_mod_region_annual_total = mod_region_annual+mod_region_annual_dry-$
 ;                     (mod_region_annual2+mod_region_annual_dry2) ;  difference between std run and run without US emissions
 ;  bargraph,mod_region_annual(sub)+ mod_region_annual_dry(sub),barcolor=!myct.gray,$
 ;  	BARWIDTH=width,pos=pos,/no_label,xlabel=region_all(sub),$
;	ytitle='Total Deposition ('+micro+ 'g m!u-2!n year!u-1!n)',yrange=[0,45],/yst
 ;  bargraph,us_mod_region_annual_total(sub),barcolor=!myct.blue,$
 ;  	BARWIDTH=width,pos=pos,xlabel=region_all(sub),/overplot,$
;	barlabel=us_mod_region_annual_total(sub)/(mod_region_annual(sub)+mod_region_annual_dry(sub))*100.,l_format='(f4.1,"%")',$
;	barcharsize=1.2
 ;  axis,xaxis=1,col=1,xcharsize=0.0001
 ;
 ;
 ;  endif

   !p.charsize=1

endif ;  end of plotting seasonal cycle by region

  if keyword_set(sites) and keyword_set(region) then begin
   !p.charsize=1.5
   for k=0,nregion-1 do begin
      mdn_data=mdn_data_all ;  reset data to the whole domain
      ind=where(mdn_data.lat ge limit_all(k).y(0)-dlat and mdn_data.lat lt limit_all(k).y(2)+dlat and $
              mdn_data.lon ge limit_all(k).x(0)-dlon and mdn_data.lon lt limit_all(k).x(2)+dlon and $
   	      min(mdn_data.npts,dimension=1) ge min_npts,nsites)
      mdn_monthly=mean(mdn_data(ind).monthly_dep,2,/nan)
      mdn_monthly_std=mdn_monthly*0.
      multipanel,col=1,row=2,margin=[0.02,0.02],pos=pos
      map_set,/usa,/continents,$
      	limit=[limit_all(k).y(0)-dlat,limit_all(k).x(0)-dlon,limit_all(k).y(2)+dlat,limit_all(k).x(2)+dlon],col=1,pos=pos,/iso
      oplot,mdn_data(ind).lon,mdn_data(ind).lat,col=2,psym=sym(1),symsize=1
      xyouts,mdn_data(ind).lon,mdn_data(ind).lat,mdn_data(ind).id,/data,col=1,charsize=1.
      multipanel,col=3,row=5,margin=[0.04,0.04],pos=pos
      yrange=[0,max(mdn_data(ind).monthly_dep,/nan)]
      if region_all(k) eq 'Gulf' then yrange=[0,150]
      if region_all(k) eq 'Florida' then yrange=[0,200]
      for i=0,nsites-1 do begin
        ; yrange=[0,max(mdn_data(ind(i)).monthly_dep,/nan)>max(model_wetdep_monthly_site(*,ind(i)))]
        plot,1+indgen(12),mdn_data(ind(i)).monthly_dep,col=1,pos=pos,$
 	  ytitle='Monthly deposition (ng m!u-2!nd!u-1!n)',thick=5,yrange=yrange,$
 	  xrange=[0.4,12.6],/xst,title='Region: '+region_all(k)+' site:'+mdn_data(ind(i)).id,$
	  xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
          xtickv=indgen(12)+1,xminor=1
        oplot,1+indgen(12),mdn_data(ind(i)).monthly_dep,col=1,thick=5
	errorbar,1+indgen(12),mdn_data(ind(i)).monthly_dep,mdn_data(ind(i)).monthly_dep_std,col=1,thick=5
        xyouts,1.1,!y.crange(1)*0.9,mdn_data(ind(i)).name,col=1,al=0,charsize=0.5
        oplot,1+indgen(12),model_wetdep_monthly_site(*,ind(i)),col=!myct.red,thick=5
        if keyword_set(version2) then oplot,1+indgen(12),model_wetdep_monthly_site2(*,ind(i)),col=!myct.green,thick=5
        if keyword_set(version3) then oplot,1+indgen(12),model_wetdep_monthly_site3(*,ind(i)),col=!myct.green,thick=5,linestyle=2
		 if keyword_set(version4) then oplot,1+indgen(12),model_wetdep_monthly_site4(*,ind(i)),col=!myct.green,thick=5,linestyle=3
		 
        multipanel,/advance,pos=pos
      endfor
   endfor
   !p.charsize=1.
  endif

 ;  plot map of annual wet deposition regridded onto the GEOS-Chem grid
 if keyword_set(regrid) then begin
    multipanel,col=2,row=3,mar=0.01
    cbpos=[0.1,0.05,0.9,0.1]
    indlon=where(mdn_data_regrid.lon ge lon0(0) and mdn_data_regrid.lon le lon0(1))
    indlat=where(mdn_data_regrid.lat ge lat0(0) and mdn_data_regrid.lat le lat0(1))
    mdn_plot=mdn_data_regrid.mdn*1e-3*365.
    mdn_plot=mdn_plot(indlon,*) & mdn_plot=mdn_plot(*,indlat)

    ;  plot MDN wet deposition mapped onto the GEOS-Chem grid
    qind=where(finite(mdn_plot) eq 0)
    mdn_plot(qind)=-9999.
    ctm_overlay,mdn_plot,mdn_data_regrid.lon(indlon),mdn_data_regrid.lat(indlat),$
         mdn_wetdep,mdn_data.lon,mdn_data.lat,$
         /sample,/continents,/usa,thick=1, $
         mindata=0.,maxdata=24.0,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!cMDN data for '+string(select_year,format='(I4,"-",I4)'),$
         cbunit=micro+'g m!u-2!n y!u-1!n',/iso,csfac=0.8,limit=limit,/cbar,$
	     min_valid=0.,botoutofrange=!myct.mediumgray,tria=1,$    ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    tvmap,model_wetdep_plot,longitude,latitude,$
         /sample,/continents,/usa,thick=1, $
         mindata=0.,maxdata=24.0,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!cGEOS-Chem '+year_str,$
         cbunit=micro+'g m!u-2!n y!u-1!n',/iso,csfac=0.8,limit=limit,/cbar,$  ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    ;  set model grid points where no data is available to NaN
    model_wetdep_plot_sample=model_wetdep_plot
    model_wetdep_plot_sample(qind)=-9999.

    ;  plot model wet deposition only where there is MDN data
    tvmap,model_wetdep_plot_sample,longitude,latitude,$
         /sample,/continents,/usa,thick=1, $
         mindata=0.,maxdata=24.0,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!cGEOS-Chem '+year_str,$
         cbunit=micro+'g m!u-2!n y!u-1!n',/iso,csfac=0.8,limit=limit,/cbar,$
	     min_valid=0.,botoutofrange=!myct.mediumgray,tria=1,  $  ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    Myct,/buwhwhrd,ncolors=10
    ;  plot the difference MDN-model
    wetdep_diff=model_wetdep_plot_sample-mdn_plot
    wetdep_diff(qind)=-9999.
    tvmap,wetdep_diff,longitude,latitude,$
         /sample,/continents,/usa,thick=1,$
         mindata=-15,maxdata=15.0,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!cGEOS-Chem - MDN',$
         cbunit=micro+'g m!u-2!n y!u-1!n',/iso,csfac=0.8,limit=limit,/cbar,$
	     min_valid=-15.,botoutofrange=!myct.mediumgray,tria=1,$   ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    ;  plot the percent difference
    qind=where(mdn_plot ne -9999)
    wetdep_diff_perc=mdn_plot*0.-9999.
    wetdep_diff_perc(qind)=-(mdn_plot(qind)-model_wetdep_plot_sample(qind))/mdn_plot(qind)*100. ;  %
    tvmap,wetdep_diff_perc,longitude,latitude,$
         /sample,/continents,/usa,thick=1,$
         mindata=-100,maxdata=100.0,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!c(GEOS-Chem - MDN) / MDN',$
         cbunit='%',/iso,csfac=0.8,limit=limit,/cbar,$
	 min_valid=-100.,botoutofrange=!myct.mediumgray,tria=1,$    ; cbpos=cbpos,$
	 div=COLORBAR_NDIV(maxdiv=8)

    multipanel,/noerase,margin=[0.07,0.03,0.02,0.03],pos=pos 
    q=where(mdn_plot ne -9999)
    !p.charsize=1.5
    !x.charsize=!p.charsize & !y.charsize=!p.charsize
    plot,mdn_plot(q),model_wetdep_plot_sample(q),col=1,xrange=[0,30],yrange=[0,30],$
 	xtitle= 'MDN deposition ('+micro+'g m!u-2!ny!u-1!n)', thick=1, $
	ytitle='GEOS-Chem deposition ('+micro+'g m!u-2!ny!u-1!n)',$
	title='Annual deposition !ccomparison',psym=sym(1),pos=pos,/xst,/yst
    rcorr=correlate(mdn_plot(q),model_wetdep_plot_sample(q))
    org_corr,mdn_plot(q),model_wetdep_plot_sample(q),rcorr,$
 	n_elements(mdn_plot(q)),slope,intercept
    oplot,!x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
    oplot,!x.crange,!x.crange,col=1
    xyouts,!x.crange(1)/20.,!y.crange(1)*0.92,$
         /data,col=1, 'Model = '+string(intercept,slope,format='(F6.3,"+ MDN x",F5.2)')+$
         '!c(r= '+$
         string(rcorr,format='(F5.2)')+ ';  npts= '+$
         string(n_elements(mdn_plot(q)),format='(I6)')+ ')'+$
	     '!cModel/Obs='+string(mean(model_wetdep_plot_sample(q)/mdn_plot(q)),format='(F5.2)'),charsize=0.9

     MyCt, /Verbose, /WhGrYlRd,ncolors=20
    !p.charsize=1.
    !x.charsize=!p.charsize & !y.charsize=!p.charsize
 endif

;===========================================================================================
; plot the anomaly
;===========================================================================================
if keyword_set(anomaly) then begin
 myct,/BuYlRd, ncolors=20
; myct, /diff, ncolors=10

  multipanel,col=2,row=4,margin=[0.02,0.02]
  model_anomaly_annual=reform(mean(model_anomaly,4))
  model_anomaly_annual_clim=mean(model_anomaly_annual,3)
  cbar=intarr(7)  & cbar[6]=1
  cbunit=strarr(7) & cbunit(6)='%'
  for iyear=0,myears-1 do begin
     tvmap, (model_anomaly_annual[*,*,iyear]-model_anomaly_annual_clim)/model_anomaly_annual_clim*100,longitude,latitude,$
	/sample,/continents, thick=1, $
	mindata=-50.,maxdata=50,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,  title='Modeled annual anomaly ' + string(iyear+min(select_year),format='(I4)'), $
        /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01, $
        cbar=cbar(iyear), 	div=COLORBAR_NDIV(maxdiv=10),  $
       cbunit=cbunit(iyear)
  endfor

if keyword_set(seasonal) then begin
 for iseason=0,3 do begin
  multipanel,col=2,row=4,margin=[0.02,0.02]
   model_anomaly_season_clim =   mean(reform(mean(model_anomaly[*,*,*,season[iseason].month],4)),3)
   for iyear=0,myears-1 do begin
      tvmap, (mean(model_anomaly[*,*,iyear,season[iseason].month],4)/model_anomaly_season_clim-1.0)*100,longitude,latitude,$
 	/sample,/continents, thick=1, $
 	mindata=-50,maxdata=50,$
 	/usa,limit=limit,$
 	t_symbol=1,symsize=0.9, title='Modeled '+season[iseason].name+' anomaly' + string(iyear+min(select_year),format='(I4)'), $
         /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01,$
         cbar=cbar(iyear), 	div=COLORBAR_NDIV(maxdiv=10),  $
         cbunit=cbunit(iyear)
   endfor
 endfor
endif

; model  value at MDN sites
multipanel,col=2,row=4,margin=[0.02,0.02]
model_wetdep_annual_site_anomaly_clim=mean(model_wetdep_annual_site_anomaly,1)

for iyear=0,myears-1 do begin
percent=reform((model_wetdep_annual_site_anomaly[iyear,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 )
ctm_overlay, 0.0*wetdep_all(0).rgm-50.01,longitude,latitude,$
       [(percent), -50] , [(mdn_data.lon), -123] , [(mdn_data.lat), 27] ,$
	/sample,/continents,  $
      cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	min_valid=-50.,maxdata=50,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,  title='Modeled annual wetd anomaly !cat MDN sites ' + string(iyear+min(year),format='(I4)'), $
        cbunit=cbunit(iyear),$
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
endfor
percent= ( reform((model_wetdep_annual_site_anomaly[4,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 ) + $
                 reform((model_wetdep_annual_site_anomaly[5,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 ) + $
                 reform((model_wetdep_annual_site_anomaly[6,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 ) ) / 3.0 - $
               ( reform((model_wetdep_annual_site_anomaly[1,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 ) + $
                 reform((model_wetdep_annual_site_anomaly[2,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 ) + $
                 reform((model_wetdep_annual_site_anomaly[3,*]-model_wetdep_annual_site_anomaly_clim[*])/model_wetdep_annual_site_anomaly_clim[*]*100 ) ) / 3.0
ctm_overlay, 0.0*wetdep_all(0).rgm-50.01,longitude,latitude,$
       [(percent), -50] , [(mdn_data.lon), -123] , [(mdn_data.lat), 27] ,$
	/sample,/continents,  $
     div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	min_valid=-50.,maxdata=50,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,  title='Modeled annual wetd (2008-2010 - 2005-2007) !cat MDN sites ' + string(iyear+min(year),format='(I4)'), $        
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
 
if keyword_set(seasonal) then begin
 for iseason=0,3 do begin
 multipanel,col=2,row=4,margin=[0.02,0.02]
 for iyear=0,myears-1 do begin
 percent=reform((model_wetdep_seasonal_site_anomaly[iseason,iyear,*]/model_wetdep_seasonal_site[iseason,*]-1.0)*100)
  ctm_overlay,0.*wetdep_all(0).rgm-50.01,longitude,latitude,$
       [(percent), -50] , [(mdn_data.lon), -123] , [(mdn_data.lat), 27] ,$
 	/sample,/continents,  $
        cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
         thick=1, $
	min_valid=-50.,maxdata=50,$
 	/usa,limit=limit,$
 	t_symbol=1,symsize=0.9, title='Modeled '+season(iseason).name+ ' wetd anomaly !cat MDN sites ' + string(iyear+min(year),format='(I4)'), $
         cbunit=cbunit(iyear),$
        /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
 endfor
 endfor
endif
 
; MDN observations
multipanel,col=2,row=4,margin=[0.02,0.02]
for iyear=0,nyears-1 do begin
 percent=reform(mdn_data.annual_dep_anomaly[iyear]*100)
 ctm_overlay,0.*wetdep_all(0).rgm-50.01,longitude,latitude,$
       [(percent), -50] , [(mdn_data.lon), -123] , [(mdn_data.lat), 27] ,$    ; 1e-3*365.25* is not needed to transfer ng/m2/d to ug/m2/y
	/sample,/continents,  $
     cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	min_valid=-50.,maxdata=50,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,  title='Observed annual wetd anomaly !cat MDN sites ' + string(iyear+min(select_year),format='(I4)'), $
        cbunit=cbunit(iyear),$
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
endfor
 percent= ( reform(mdn_data.annual_dep_anomaly[4]*100) + $
                  reform(mdn_data.annual_dep_anomaly[5]*100) + $
                  reform(mdn_data.annual_dep_anomaly[6]*100) ) / 3.0  - $
                ( reform(mdn_data.annual_dep_anomaly[1]*100) + $
                  reform(mdn_data.annual_dep_anomaly[2]*100) + $
                  reform(mdn_data.annual_dep_anomaly[3]*100) ) / 3.0
 ctm_overlay,0.*wetdep_all(0).rgm-50.01,longitude,latitude,$
       [(percent), -50] , [(mdn_data.lon), -123] , [(mdn_data.lat), 27] ,$
	/sample,/continents,  $
    div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	min_valid=-50.,maxdata=50,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,  title='Observed annual wetd  (2008-2010 - 2005-2007) !cat MDN sites ' + string(iyear+min(select_year),format='(I4)'), $
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01


if keyword_set(seasonal) then begin
 for iseason=0,3 do begin
 multipanel,col=2,row=4,margin=[0.02,0.02]
 for iyear=0,myears-1 do begin
 percent=reform(mdn_data.season_dep_anomaly[iseason,iyear]*100)
  ctm_overlay,0.*wetdep_all(0).rgm-50.01,longitude,latitude,$
        [(percent), -50] , [(mdn_data.lon), -123] , [(mdn_data.lat), 27] ,$                   ; 1e-3*365.25*
 	/sample,/continents,  $
        cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
         thick=1, $
	min_valid=-50.,maxdata=50,$
 	/usa,limit=limit,$
 	t_symbol=1,symsize=0.9, title='Observed '+season[iseason].name+ ' wetd anomaly !cat MDN sites ' + string(iyear+min(select_year),format='(I4)'), $
         cbunit=cbunit(iyear)  ,$
        /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
 endfor
 endfor
endif

; plot the anomaly of model bias
; restore the nsites
nsites=nsites_all

model_annual_bias_sites=fltarr(myears,nsites)
model_seasonal_bias_sites=fltarr(4,myears,nsites)
for iyear=0,myears-1 do begin
   for isite=0,nsites-1 do begin
       model_annual_bias_sites[iyear,isite]=model_wetdep_annual_site_anomaly[iyear,isite]- mdn_data[isite].annual_dep*(1.0+mdn_data[isite].annual_dep_anomaly[iyear])*1e-3*365.25    ; need to transfer ng/m2/d to ug/m2/y
       model_seasonal_bias_sites[*,iyear,isite]=model_wetdep_seasonal_site[*,isite]-mdn_data[isite].season_dep_anomaly[*,iyear]*1e-3*365.25
   endfor
endfor

;  model bias
 cbunit(6)=micro+'g m!u-2!n'
 multipanel,col=2,row=4,margin=[0.001,0.001]
 for iyear=0,nyears-1 do begin
  ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
         model_annual_bias_sites[iyear,*],mdn_data.lon,mdn_data.lat,$
 	/sample,/continents,  $
       cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
         thick=1, $
 	mindata=-4.,maxdata=4,$
 	/usa,limit=limit,$
 	t_symbol=1,symsize=0.9, title='Annual model bias !cat MDN sites ' + string(iyear+min(select_year),format='(I4)'), $
        cbunit=cbunit(iyear),$
        /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
 endfor

; slope of model bias, for annual data
slope_mdn=0.*model_annual_bias_sites(0,*)
slope_mdn_err=slope_mdn
rcorr=slope_mdn
chisquare=slope_mdn

for isite=0, n_elements(slope_mdn)-1 do begin
    xdata=findgen(nyears)
    ydata=-model_annual_bias_sites[*,isite]
    r=correlate(xdata,ydata)
;    org_corr,xdata,ydata,r,nyears,slope,intercept,slope_err     
    result=poly_fit(xdata,ydata,1,sigma=sigma,chisq=chisq,yfit=yfit)
    rcorr[isite]=r
;    slope_mdn[isite]=slope
;    slope_mdn_err[isite]=slope_err
    slope_mdn[isite]=result[1]
    slope_mdn_err[isite]=sigma[1]
    chisquare[isite]=chisq
        
endfor
; chisq_threshold=10000.0
 chisq_threshold=12.59  ; alpha=5%, df = 7-1
; chisq_threshold=10.64 ; 10%
; chisq_threshold=16.81 ; 1%
; chisq_threshold=18.54 ; 0.5%

ind_significance=where(chisquare le chisq_threshold)

multipanel,col=2,row=2,margin=[0.001,0.001]
 ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
        rcorr[ind_significance],mdn_data[ind_significance].lon,mdn_data[ind_significance].lat,$
	/sample,/continents,  $
      /cbar, div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	mindata=-1.,maxdata=1,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9, title='Correlation coefficient between model !cresidue and year at MDN sites', $
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01

 ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
        100.0*slope_mdn[ind_significance]/(0.3625*mdn_data[ind_significance].annual_dep),mdn_data[ind_significance].lon,mdn_data[ind_significance].lat,$
	/sample,/continents,  $
      /cbar, div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	mindata=-10,maxdata=10,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9, title='Slope of model residue against year at MDN sites', cbunit='%', $
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01

if keyword_set(version2) then begin
 ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
;   assuming decrease by 2.5% per year
        ; 0.3625 transfer ng/m2/day to ug/m2/year
        100.0*(slope_mdn[ind_significance] + model_wetdep_annual_site2[ind_significance]*0.025) / $
          (0.36525*mdn_data[ind_significance].annual_dep), mdn_data[ind_significance].lon,mdn_data[ind_significance].lat,$
	/sample,/continents,  $
      /cbar, div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
       mindata=-10,maxdata=10,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9, title='Slope of model residue against year at MDN sites !c subtract the global trend', cbunit='%', $
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
endif

myct,/WhGrYlRd 
 ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
        chisquare[ind_significance],mdn_data[ind_significance].lon,mdn_data[ind_significance].lat,$
	/sample,/continents,  $
      /cbar, div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	mindata=min(chisquare[ind_significance]),maxdata=max(chisquare[ind_significance]),$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9, title='Chi-square goodness of fit', $
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01 

ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
        100.0*slope_mdn_err[ind_significance]/(0.36525*mdn_data[ind_significance].annual_dep),mdn_data[ind_significance].lon,mdn_data[ind_significance].lat,$
	/sample,/continents,  $
      /cbar, div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
;	mindata=min(slope_mdn_err[ind_significance]),maxdata=max(slope_mdn_err[ind_significance]),$
    mindata=0, maxdata=5, $
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9, title='Standard errror of slope of model !cresidue against year at MDN sites', cbunit='%', $
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
 myct,/BuWhRd 

; fraction of model bias anomaly
model_annual_bias_sites_anomaly=model_annual_bias_sites[*,*]*0.0
model_seasonal_bias_sites_anomaly=model_seasonal_bias_sites[*,*,*]*0.0
 for isite=0,nsites-1 do begin
      model_annual_bias_sites_anomaly[*,isite]=(model_annual_bias_sites[*,isite]-mean(model_annual_bias_sites[*,isite]))/mean(model_annual_bias_sites[*,isite])*100
      for iseason=0,3 do model_seasonal_bias_sites_anomaly[iseason,*,isite]=(model_seasonal_bias_sites[iseason,*,isite]-mean(model_seasonal_bias_sites[iseason,*,isite]))/ $
                                                            mean(model_seasonal_bias_sites[iseason,*,isite])*100
 endfor

cbunit(6)='%'
multipanel,col=2,row=4,margin=[0.02,0.02]
for iyear=0,nyears-1 do begin
 ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
        model_annual_bias_sites_anomaly[iyear,*],mdn_data.lon,mdn_data.lat,$
	/sample,/continents,  $
       cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
        thick=1, $
	mindata=-100,maxdata=100,$
	/usa,limit=limit,$
	t_symbol=1,symsize=0.9,  title='Annual model bias anomaly !cat MDN sites ' + string(iyear+min(select_year),format='(I4)'), $
        cbunit=cbunit(iyear),$
       /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
endfor

if keyword_set(seasonal) then begin
 for iseason=0,3 do begin
 multipanel,col=2,row=4,margin=[0.02,0.02]
 for iyear=0,myears-1 do begin
  ctm_overlay,0.*wetdep_all(0).rgm,longitude,latitude,$
        4.0*model_seasonal_bias_sites_anomaly[iseason,iyear,*] ,mdn_data.lon,mdn_data.lat,$
 	/sample,/continents,  $
        cbar=cbar(iyear), div=COLORBAR_NDIV(maxdiv=10),$
         thick=1, $
 	mindata=-50.,maxdata=50,$
 	/usa,limit=limit,$
 	t_symbol=1,symsize=0.9, title=season[iseason].name+ ' model bias anomaly !cat MDN sites ' + string(iyear+min(select_year),format='(I4)'), $
         cbunit=cbunit(iyear),$
        /iso,csfac=0.8,/NOGXLABELS ,/NOGYLABELS,margin=0.01
 endfor
 endfor
endif

myct,/WhGrYlRd 
endif

; ====================================================================================
; end of plot anomaly
; ====================================================================================

; plot interannual variability of different regions
if keyword_set(interannual) and keyword_set(region) then begin

; modify the year range, which is different from the standard years
; select_year=[1996,2010]  ; year range for MDN data
; year=[1996,2010]            ; year range for GC model

; tyear=[1996,2010]          ; a larger year range covers both the observation and model year range
  tyear = year

   nsyear=max(select_year)-min(select_year)+1
   ntyear=max(tyear)-min(tyear) +1
   nmonth=ntyear*12
   mdn_region_interannual_site=fltarr(nregion, nmonth)   &  mdn_region_interannual_site[*]=!values.f_nan
   mdn_region_interannual_site_std=fltarr(nregion, nmonth) & mdn_region_interannual_site_std[*]=!values.f_nan
   mod_region_interannual=fltarr(nregion, nmonth) & mod_region_interannual[*]=!values.f_nan
 
    for iyear=min(tyear),max(tyear) do begin
      
      ; mdn_data
      if iyear ge min(select_year) and iyear le max(select_year) then begin
            print, 'reading MDN data for year'+string(iyear)
            read_mdn,file,mdn_data_tmp,mdn_data_regrid_tmp,mdn_data_regrid_monthly_tmp,select_year=iyear,res=resolution,min_npts=min_npts , snow=snow, mettype=mettype, metresolution=metresolution
   
            for iregion=0,nregion-1 do begin
                ind=where(mdn_data_tmp.lat ge limit_all(iregion).y(0) and mdn_data_tmp.lat lt limit_all(iregion).y(2) and $
                    mdn_data_tmp.lon ge limit_all(iregion).x(0) and mdn_data_tmp.lon lt limit_all(iregion).x(2) and $
                    min(mdn_data_tmp.npts,dimension=1) ge min_npts)     
                if n_elements(ind) eq 1 then $
                mdn_region_interannual_site(iregion, 12*(iyear-min(tyear) ) : 12*(iyear-min(tyear) )+11 )=mdn_data_tmp(ind).monthly_dep $
                else mdn_region_interannual_site(iregion, 12*(iyear-min(tyear) ) : 12*(iyear-min(tyear) )+11 )=mean(mdn_data_tmp(ind).monthly_dep,2,/nan)
            
                for imonth=1,12 do begin
                      if n_elements(ind) gt 1 then $
                               mdn_region_interannual_site_std( iregion,  (iyear-min(tyear))*12+imonth-1)=stddev(mdn_data_tmp(ind).monthly_dep(imonth-1),/nan)  $
                               else mdn_region_interannual_site_std( iregion,  (iyear-min(tyear))*12+imonth-1)=0
                      if not finite(mdn_region_interannual_site_std( iregion,  (iyear-min(tyear))*12+imonth-1)) then mdn_region_interannual_site_std( iregion,  (iyear-min(tyear))*12+imonth-1) = 0.
                endfor
            
            endfor
       endif

       ; mod_data
       if iyear ge min(year) and iyear le max(year) then begin
               print, 'reading GEOS-Chem wet deposition flux data for year'+string(iyear)
               read_run_wdep,file=file_std,dummy1,dummy2,dummy5,$
                       wetdep_all_tmp,longitude,latitude,convert,dummy3,dummy4,$
    	                gridinfo,lat0=[lat0(0)+dlat,lat0(1)-dlat],lon0=[lon0(0)+dlon,lon0(1)-dlon],year=iyear
                 
                 nsites=n_elements(mdn_data_tmp)
                 model_wetdep_monthly_site=fltarr(12,nsites)+!values.f_nan
                 for isite=0,nsites-1 do begin
                      ; find location of MDN site
                      ctm_index,ctm_type('GEOS5_47L',res=res),i,j,$
                 	                  center=[mdn_data_tmp(isite).lat,mdn_data_tmp(isite).lon],/non_interactive
                      indlat=where( latitude eq lat_all(j-1) )
                      indlon=where( longitude eq lon_all(i-1) )
                     
                      ;  ng/m2/day
                      if indlat[0] ne -1 and indlon[0] ne -1 then $
                      model_wetdep_monthly_site(*,isite)=(wetdep_all_tmp.rgm(indlon(0),indlat(0),*)+wetdep_all_tmp.hgp(indlon(0),indlat(0),* ))*convert(indlon(0),indlat(0))*1e3/ndays
                 endfor
           
                 for iregion=0,nregion-1 do begin
                     ind=where(mdn_data_tmp.lat ge limit_all(iregion).y(0) and mdn_data_tmp.lat lt limit_all(iregion).y(2) and $
                         mdn_data_tmp.lon ge limit_all(iregion).x(0) and mdn_data_tmp.lon lt limit_all(iregion).x(2) and $
              	          min(mdn_data_tmp.npts,dimension=1) ge min_npts) 
                     if n_elements(ind) gt 1 then $
                       mod_region_interannual(iregion, 12*(iyear-min(tyear) ) : 12*(iyear-min(tyear) )+11 )=mean(model_wetdep_monthly_site(*,ind),2,/nan) $
                     else $
                       mod_region_interannual(iregion, 12*(iyear-min(tyear) ) : 12*(iyear-min(tyear) )+11 )=model_wetdep_monthly_site(*,ind)

                  endfor
       endif
   
   endfor

; fix NAN with model values
; ind = where( finite(mdn_region_interannual_site) eq 0)
; mdn_region_interannual_site(ind) = mod_region_interannual(ind)

; plot the interannual variation
multipanel,col=1,row=8,margin=[0.03,0.005],pos=pos
 !p.charsize=1.5

; which year to start?
start_year=1996

xrange=[-1,  nmonth - 12 * (start_year - min(select_year)) ]
for iregion=0,nregion-1 do begin
       yrange=[0, 4]
       if region_all(iregion) eq 'SE' or region_all(iregion) eq 'Gulf'  then yrange=[0,6.]
       if region_all(iregion) eq 'Florida' then yrange=[0,7.]
       if region_all(iregion) eq 'WEST' then yrange=[-0.5, 2]
       if region_all(iregion) eq 'NORTHEAST' then yrange=[0, 2]

       if iregion eq nregion-1 then begin
                  xtickname=string( indgen(ntyear-(start_year-min(select_year)))+start_year, format='(6x,I4)' ) 
                  xtitle='Year'
       endif else begin
                  xtickname=replicate(' ', ntyear-(start_year-min(select_year)))
                  xtitle=' '
       endelse
       if iregion eq 1 then begin
                   ytitle='Wet deposition, ' + micro + 'g m!u-2!nmonth!u-1!n'
       endif else begin
                   ytitle=' '  
       endelse

       plot, indgen( nmonth - 12 * (start_year - min(select_year)) ) , 0.03*mdn_region_interannual_site(iregion, 12*(start_year-min(select_year)):nmonth-1), col=1, pos=pos,$
 	        ytitle=ytitle, thick=1,yrange=yrange,$
        	/xst, /yst, xrange=xrange, xtitle=xtitle , $
            xticks=ntyear-(start_year-min(select_year))-1, xtickname=xtickname, $
		    xtickv=indgen(ntyear-(start_year-min(select_year)))*12+5.5, xminor=12,yminor=2, xticklen=0.05, yticklen=0.01

       print,region_all(iregion), transpose(0.03*mdn_region_interannual_site(iregion,12*(start_year-min(select_year)):nmonth-1))
       polyfill,[ indgen( (nsyear-(start_year-min(select_year)))*12 ), reverse(indgen( (nsyear-(start_year-min(select_year)))*12 )),0], $
                0.03*[ reform(mdn_region_interannual_site(iregion,12*(start_year-min(select_year)):nmonth-1) - mdn_region_interannual_site_std(iregion,12*(start_year-min(select_year)):nmonth-1)), $
                reverse(reform(mdn_region_interannual_site(iregion,12*(start_year-min(select_year)):nmonth-1)+mdn_region_interannual_site_std(iregion,12*(start_year-min(select_year)):nmonth-1))), $
                 mdn_region_interannual_site[iregion,12*(start_year-min(select_year))]-mdn_region_interannual_site_std[iregion,12*(start_year-min(select_year))]]  $
                ,col=!myct.lightgray
 
       oplot, indgen( nmonth - 12 * (start_year - min(select_year)) ), 0.03*mdn_region_interannual_site(iregion,12*(start_year-min(select_year)):nmonth-1),col=1,thick=1

       oplot, indgen( nmonth - 12 * (start_year - min(select_year)) ), 0.03*reform(mod_region_interannual(iregion,12*(start_year-min(select_year)):nmonth-1)) ,col=!myct.red,thick=1
	   if iregion eq 0 then legend,label=['MDN', 'GEOS-Chem model '],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=1,$
 	             halign=0.01,valign=0.95,charsize=.75
       xyouts, !x.crange(1)*0.6, !y.crange(1)*0.85, $
         /data,col=1, region_all(iregion) + ' (' +  string(nsites_regionall(iregion),format='(I2," sites")') + ')', charsize=0.7

       multipanel,/advance,pos=pos
;       multipanel,/advance,pos=pos

endfor

; large array to hold the time series
meta_ydata = fltarr(nmonth, nregion, 2)  ; month,  region,  direct observation or subtract model

; directly regression the observation, also filter out the seasonal cycle
multipanel,col=1,row=8,margin=[0.045,0.005,0.015,0.005],pos=pos
 !p.charsize=1.5

xrange=[-1,  nmonth ]

for iregion=0,nregion-1 do begin
       yrange=[-1,1]
        xdata=indgen( nmonth )
        ydata= 0.03*mdn_region_interannual_site(iregion,*)

;       del the seasonal cycle
        ydata_ref=reform(ydata,12,ntyear)
        ydata_mean=mean(ydata_ref,2,/nan)
        for imonth=0,11 do begin
            ydata_ref[imonth,*]=(ydata_ref[imonth,*]-ydata_mean(imonth))
        endfor
        ydata=reform(ydata_ref,nmonth)

        ind_finite=indgen(nmonth)
        ydata=ydata[ind_finite]
        rcorr=correlate(xdata,ydata)
        result=poly_fit(xdata,ydata,1,sigma=sigma,yfit=yfit)
        slope=result[1]
        intercept=result[0]
        t=slope / ( sqrt(total((ydata-yfit)^2)/(nmonth-2.0))/sqrt(total((xdata-mean(xdata))^2))  )       

       if iregion eq nregion-1 then begin
                  xtickname=string( indgen(ntyear)+min(year), format='(I4)' )
                  xtitle='Year'
       endif else begin
                  xtickname=replicate(' ', ntyear)
                  xtitle=' '
       endelse
  
       ; store data
       meta_ydata(*,iregion,0)=ydata

       plot, xdata , ydata , col=1, pos=pos, $
 	        ytitle='Hg wet deposition !c(' + micro + 'g m!u-2!nmonth!u-1!n)',thick=1,yrange=yrange,$
        	/xst,  xrange=xrange,xtitle=xtitle , $
            xticks=ntyear-1, xtickname=xtickname, $
		    xtickv=indgen(ntyear)*12+5.5, xminor=6
 
       oplot,!x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
       xyouts,!x.crange(1)/20.,!y.crange(1)*0.8,$
         /data,col=1, 'MDN = ' +         string(144*intercept,format='(F8.5)') + '+ month x' + $
        string(144*slope,format='(F8.5)') + '+/-' + '(' + $
        string(144*sigma[1],format='(F8.5)') + ')' +  $
        string(rcorr,format='(" r = ",F5.2)') +  $
        string(t,format='(" t = ",F5.2)') + $ 
        string(t_cvf(0.025,84-2), format='(" t_cutoff = ",F5.2)'), charsize=0.7

       xyouts, !x.crange(1)*0.85, !y.crange(1)*0.70,$
         /data,col=1, region_all(iregion) + ' (' +  string(nsites_regionall(iregion),format='(I2)') + ')', charsize=0.7
       
       multipanel,/advance,pos=pos
       
endfor

; plot the model bias
multipanel,col=1,row=8,margin=[0.045,0.005,0.015,0.005],pos=pos
 !p.charsize=1.5

xrange=[-1,  nmonth]

for iregion=0,nregion-1 do begin
       yrange=[-1,1]
;       if region_all(iregion) eq 'SE' or region_all(iregion) eq 'Gulf'  then yrange=[0,6.]
;       if region_all(iregion) eq 'Florida' then yrange=[0,7.]
        xdata=indgen( nmonth )
        ydata=-(0.03*reform(mod_region_interannual(iregion,*)) - 0.03*mdn_region_interannual_site(iregion,*))

;       del the seasonal cycle
        ydata_ref=reform(ydata,12,ntyear)
        ydata_mean=mean(ydata_ref,2,/nan)
        for imonth=0,11 do begin
            ydata_ref[imonth,*]=(ydata_ref[imonth,*]-ydata_mean(imonth))
        endfor
        ydata=reform(ydata_ref,nmonth)

        ind_finite=indgen(nmonth)
        ydata=ydata[ind_finite]
;        rcorr=correlate(xdata,ydata)
;        org_corr,xdata[ind_finite],ydata[ind_finite],rcorr,nmonth-48,slope,intercept
        result=poly_fit(xdata,ydata,1,sigma=sigma,yfit=yfit)
        slope=result[1]
        intercept=result[0]         
        t=slope / ( sqrt(total((ydata-yfit)^2)/(nmonth-2.0))/sqrt(total((xdata-mean(xdata))^2))  )     

;       use mann-kendall test
;      result=mkstrend(ydata)
;      slope=result.slope
;      intercept=result.intercept
;      p=result.p

       if iregion eq nregion-1 then begin
                  xtickname=string( indgen(ntyear)+min(year), format='(I4)' )
                  xtitle='Year'
       endif else begin
                  xtickname=replicate(' ', ntyear)
                  xtitle=' '
       endelse

       ; store data
       meta_ydata(*,iregion,1)=ydata

       plot, xdata , ydata , col=1, pos=pos, $
 	        ytitle='Hg wet deposition !c(' + micro + 'g m!u-2!nm!u-1!n)',thick=1,yrange=yrange,$
        	/xst, xrange=xrange,xtitle=xtitle , $
            xticks=ntyear-1,xtickname=xtickname, $
		    xtickv=indgen(ntyear)*12+5.5, xminor=6
 
       oplot,!x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
;       xyouts,!x.crange(1)/20.,!y.crange(1)*0.8,$
;         /data,col=1, 'Model residue = ' + string(144*intercept,format='(F8.5)') + '+ month x' + $
;        string(144*slope,format='(F8.5)') + '+/-' + '(' + $
;        string(144*sigma[1],format='(F8.5)') + ')' +  $
;        string(rcorr,format='(" r = ",F5.2)') + $
;        string(t,format='(" t = ",F5.2)') + $ 
;        string(t_cvf(0.025,84-2), format='(" t_cutoff = ",F5.2)'), charsize=0.7
       xyouts,!x.crange(1)/60.,!y.crange(1)*0.70,$
         /data,col=1,  $
        'slope: '+ string(144*slope,format='(F6.3)') + ' ' + $
        '+/-' + string(144*sigma[1],format='(F5.3)') + ' ' + $
         micro + 'g m!u-2!n yr!u-2!n' + '(' +  string(100* 144*slope /mean(12*0.03*mdn_region_interannual_site(iregion,*),/nan),  $
                     format="(f5.2)"  )+ '% yr!u-1!n)!c' + $ 
        string(t,format='(" t = ",F5.2)') + ', ' + $ 
        string(t_cvf(0.025,nmonth-2), format='(" t!lcutoff!n = ",F5.2)'), $
;        string(p,format='(" p = ",F5.3)'),   $ 
         charsize=0.9

       xyouts, !x.crange(1)*0.85, !y.crange(1)*0.70,$
         /data,col=1, region_all(iregion) + ' (' +  string(nsites_regionall(iregion),format='(I2)') + ')', charsize=0.7
      
       multipanel,/advance,pos=pos
       
endfor

; plot the two time series together
multipanel,col=1,row=8,margin=[0.03,0.005],pos=pos
!p.charsize=1.5
xdata=indgen(nmonth)

xrange=[-1,  nmonth ]
for iregion=0,nregion-1 do begin
       yrange=[-1, 1]
 if region_all(iregion) eq 'NORTHEAST' then yrange=[-0.5, 0.5]
 if region_all(iregion) eq 'WEST' then yrange=[-0.5, 0.5]

       if iregion eq nregion-1 then begin
                  xtickname=string( indgen(ntyear)+min(year), format='(I4)' ) 
                  xtitle='Year'
       endif else begin
                  xtickname=replicate(' ', ntyear)
                  xtitle=' '
       endelse
       if iregion eq 1 then begin
                   ytitle='Hg wet deposition anomaly, ' + micro + 'g m!u-2!nmonth!u-1!n'
       endif else begin
                   ytitle=' '  
       endelse

       plot, indgen( nmonth ) , meta_ydata(*,iregion,0), col=1, pos=pos,$
 	        ytitle=ytitle, thick=1,yrange=yrange,$
        	/xst, /yst, xrange=xrange, xtitle=xtitle , $
            xticks=ntyear-1, xtickname=xtickname, $
		    xtickv=indgen(ntyear)*12+5.5, xminor=12,yminor=2, xticklen=0.05, yticklen=0.01

        result=poly_fit(xdata,meta_ydata(*,iregion,0),1,sigma=sigma,yfit=yfit)
        slope1=result[1]
        intercept1=result[0]         
        t1=slope1 / ( sqrt(total((meta_ydata(*,iregion,0)-yfit)^2)/(nmonth-2.0))/sqrt(total((xdata-mean(xdata))^2))  )   
        p1=1-t_pdf(abs(t1),nmonth-2)  

       oplot, indgen( nmonth ), meta_ydata(*,iregion,1) ,col=!myct.red,thick=1

        result=poly_fit(xdata,meta_ydata(*,iregion,1),1,sigma=sigma,yfit=yfit)
        slope2=result[1]
        intercept2=result[0]         
        t2=slope2 / ( sqrt(total((meta_ydata(*,iregion,1)-yfit)^2)/(nmonth-2.0))/sqrt(total((xdata-mean(xdata))^2))  )     
        p2=1-t_pdf(abs(t2),nmonth-2)

	   if iregion eq 0 then legend,label=['MDN', 'MDN - GC'],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=1,$
 	             halign=0.01,valign=0.95,charsize=.75
       xyouts, !x.crange(0)+(!x.crange(1)-!x.crange(0))*0.6, !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.85, $
         /data,col=1, region_all(iregion) + ' (' +  string(nsites_regionall(iregion),format='(I2," sites")') + ')', charsize=0.7

       xyouts,!x.crange(1)/60.,  !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.20, /data,col=1,  $
        string(100* 144 * slope1/mean(12*0.03*mdn_region_interannual_site(iregion,*),/nan), format='(F5.2)') +  $
         '% yr!u-1!n' + string(p1,format='(" p = ",F5.3)'),   $ 
          charsize=0.7

       xyouts,!x.crange(1)/60.,  !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.10, /data,col=!myct.red,  $
        string(100* 144 * slope2/mean(12*0.03*mdn_region_interannual_site(iregion,*),/nan), format='(F5.2)') +  $
         '% yr!u-1!n' + string(p2,format='(" p = ",F5.3)'),   $ 
          charsize=0.7

       multipanel,/advance,pos=pos
 ;      multipanel,/advance,pos=pos

endfor

; plot the interannual variation, by monthly anomaly
print, 'monthly values'
multipanel,col=1,row=6,margin=[0.03,0.03],pos=pos
 !p.charsize=1.5

xrange=[-1,  nmonth ]

for iregion=0,nregion-1 do begin

       ; calculate the anomaly
       data_monthly_mean=fltarr(12)
       data_monthly_stdev=fltarr(12)
       data_nyear=n_elements(mdn_region_interannual_site(iregion,*))/12
       data_monthly_anomaly=fltarr(12*data_nyear)       
      for imonth=0,11 do begin
            data_monthly_mean[imonth]=mean(0.03*mdn_region_interannual_site(iregion,imonth+indgen(data_nyear)*12))
            data_monthly_stdev[imonth]=stddev(0.03*mdn_region_interannual_site(iregion,imonth+indgen(data_nyear)*12))
            for iyear=0,data_nyear-1 do begin
                 data_monthly_anomaly(imonth+12*iyear)=(0.03*mdn_region_interannual_site(iregion,imonth+12*iyear)-data_monthly_mean[imonth]) $
                                                                   / data_monthly_stdev[imonth]
            endfor
       endfor

       plot, indgen( nmonth ) , data_monthly_anomaly, col=1, pos=pos,$
 	        ytitle='Monthly deposition !cstandardized anomaly',thick=1,yrange=[-3,3],$
        	/xst, title='Region: ' + region_all(iregion), xrange=xrange,xtitle='Year' , $
            xticks=ntyear-1,xtickname=string( indgen(ntyear)+min(tyear), format='(I4)' ), $
		    xtickv=indgen(ntyear)*12+5.5, xminor=6
        print, region_all(iregion) + ' obs '
        print, data_monthly_anomaly

       ; calculate the anomaly
       data_monthly_anomaly=fltarr(12*data_nyear)       
      for imonth=0,11 do begin
            for iyear=0,data_nyear-1 do begin
                 data_monthly_anomaly(imonth+12*iyear)=(0.03*mod_region_interannual(iregion,imonth+12*iyear)-data_monthly_mean[imonth]) $
                                                                   / data_monthly_stdev[imonth]
            endfor
       endfor

       oplot, indgen( nmonth )   , data_monthly_anomaly ,col=!myct.red,thick=1
        print, region_all(iregion) + ' mod '
        print, data_monthly_anomaly
       oplot, indgen( nmonth )   , 0. * indgen( nmonth ), col=!myct.black,thick=1,line=2
	   legend,label=['MDN', 'GC ' + resolution],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=1,$
 	             halign=0.01,valign=0.95,charsize=.75
       multipanel,/advance,pos=pos

endfor

; plot summer only
multipanel,col=1,row=6,margin=[0.03,0.03],pos=pos
for iregion=0,nregion-1 do begin
 
       xrange=[-1,nmonth/2]       

       ; calculate the anomaly
       data_monthly_mean=fltarr(12)
       data_monthly_stdev=fltarr(12)
       data_nyear=n_elements(mdn_region_interannual_site(iregion,*))/12
       data_monthly_anomaly=fltarr(12*data_nyear)       
      for imonth=0,11 do begin
            data_monthly_mean[imonth]=mean(0.03*mdn_region_interannual_site(iregion,imonth+indgen(data_nyear)*12))
            data_monthly_stdev[imonth]=stddev(0.03*mdn_region_interannual_site(iregion,imonth+indgen(data_nyear)*12))
            for iyear=0,data_nyear-1 do begin
                 data_monthly_anomaly(imonth+12*iyear)=(0.03*mdn_region_interannual_site(iregion,imonth+12*iyear)-data_monthly_mean[imonth]) $
                                                                   / data_monthly_stdev[imonth]
            endfor
       endfor

       month_range=[5, 6, 7, 8, 9, 10]
       ind_range=month_range
       for iyear=1,data_nyear-1 do begin
         ind_range=[ind_range,12*iyear+month_range]
       endfor

       plot, indgen( nmonth/2 ) , data_monthly_anomaly[ind_range], col=1, pos=pos,$
 	        ytitle='Monthly deposition !cstandardized anomaly',thick=1,yrange=[-3,3],$
        	/xst, title='Region: ' + region_all(iregion) + ' for summer and autumn', xrange=xrange,xtitle='Year' , $
            xticks=ntyear-1,xtickname=string( indgen(ntyear)+min(tyear), format='(I4)' ), $
		    xtickv=indgen(ntyear)*6+2.5, xminor=6


       ; calculate the anomaly
       data_monthly_anomaly=fltarr(12*data_nyear)       
      for imonth=0,11 do begin
            for iyear=0,data_nyear-1 do begin
                 data_monthly_anomaly(imonth+12*iyear)=(0.03*mod_region_interannual(iregion,imonth+12*iyear)-data_monthly_mean[imonth]) $
                                                                   / data_monthly_stdev[imonth]
            endfor
       endfor

       oplot, indgen( nmonth/2 )   , data_monthly_anomaly[ind_range] ,col=!myct.red,thick=1

       oplot, indgen( nmonth )   , 0. * indgen( nmonth ), col=!myct.black,thick=1,line=2
	   legend,label=['MDN', 'GC ' + resolution],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=1,$
 	             halign=0.01,valign=0.95,charsize=.75
       multipanel,/advance,pos=pos

endfor

; plot winter only
multipanel,col=1,row=6,margin=[0.03,0.03],pos=pos
for iregion=0,nregion-1 do begin
 
       xrange=[-1,nmonth/2]       

       ; calculate the anomaly
       data_monthly_mean=fltarr(12)
       data_monthly_stdev=fltarr(12)
       data_nyear=n_elements(mdn_region_interannual_site(iregion,*))/12
       data_monthly_anomaly=fltarr(12*data_nyear)       
      for imonth=0,11 do begin
            data_monthly_mean[imonth]=mean(0.03*mdn_region_interannual_site(iregion,imonth+indgen(data_nyear)*12))
            data_monthly_stdev[imonth]=stddev(0.03*mdn_region_interannual_site(iregion,imonth+indgen(data_nyear)*12))
            for iyear=0,data_nyear-1 do begin
                 data_monthly_anomaly(imonth+12*iyear)=(0.03*mdn_region_interannual_site(iregion,imonth+12*iyear)-data_monthly_mean[imonth]) $
                                                                   / data_monthly_stdev[imonth]
            endfor
       endfor

       month_range=[0, 1, 2, 3, 4, 11]
       ind_range=month_range
       for iyear=1,data_nyear-1 do begin
         ind_range=[ind_range,12*iyear+month_range]
       endfor

       plot, indgen( nmonth/2 ) , data_monthly_anomaly[ind_range], col=1, pos=pos,$
 	        ytitle='Monthly deposition !cstandardized anomaly',thick=1,yrange=[-3,3],$
        	/xst, title='Region: ' + region_all(iregion) + ' for winter and spring', xrange=xrange,xtitle='Year' , $
            xticks=ntyear-1,xtickname=string( indgen(ntyear)+min(tyear), format='(I4)' ), $
		    xtickv=indgen(ntyear)*6+2.5, xminor=6


       ; calculate the anomaly
       data_monthly_anomaly=fltarr(12*data_nyear)       
      for imonth=0,11 do begin
            for iyear=0,data_nyear-1 do begin
                 data_monthly_anomaly(imonth+12*iyear)=(0.03*mod_region_interannual(iregion,imonth+12*iyear)-data_monthly_mean[imonth]) $
                                                                   / data_monthly_stdev[imonth]
            endfor
       endfor

       oplot, indgen( nmonth/2 )   , data_monthly_anomaly[ind_range] ,col=!myct.red,thick=1

       oplot, indgen( nmonth )   , 0. * indgen( nmonth ), col=!myct.black,thick=1,line=2
	   legend,label=['MDN', 'GC ' + resolution],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=1,$
 	             halign=0.01,valign=0.95,charsize=.75
       multipanel,/advance,pos=pos

endfor

; plot each year seperately
multipanel,col=2,row=5,margin=[0.03,0.03],pos=pos
 !p.charsize=1.5
for iregion=0,nregion-1 do begin
       yrange=[0,3.5]
       if region_all(iregion) eq 'MW' then yrange=[0,2.5]
       if region_all(iregion) eq 'NE' then yrange=[0,2]
       if region_all(iregion) eq 'OH' then yrange=[0,2.5]     
       if region_all(iregion) eq 'SE' or region_all(iregion) eq 'Gulf'  then yrange=[0,4.]
       if region_all(iregion) eq 'Florida' then yrange=[0,7.]
     for iyear=0,nmonth/12-1 do begin
       if iyear eq 0 then plot, indgen( 12 ) , 0.03*mdn_region_interannual_site(iregion,iyear*12:iyear*12+11), col=iyear+1, pos=pos,$
 	        ytitle='Observed Monthly !cdeposition (' + micro + 'g m!u-2!nm!u-1!n)',thick=5,yrange=yrange,/yst,$
        	/xst, title='Region: ' + region_all(iregion), xrange=[0,11],xtitle='Month'  else $
             oplot, indgen( 12 ) , 0.03*mdn_region_interannual_site(iregion,iyear*12:iyear*12+11), col=iyear+1,thick=5
        legend,label=[string(indgen(nmonth/12)+select_year[0])],lcol=indgen(nmonth/12)+1,line=0*indgen(nmonth/12),symbol=0*indgen(nmonth/12)-1,boxcolor=-1,thick=5,$
 	             halign=0.01,valign=0.95,charsize=.5,nlines=4
     endfor
     multipanel,/advance,pos=pos
endfor

multipanel,col=2,row=5,margin=[0.03,0.03],pos=pos
 !p.charsize=1.5
for iregion=0,nregion-1 do begin
       yrange=[0,3.5]
       if region_all(iregion) eq 'MW' then yrange=[0,2.5]
       if region_all(iregion) eq 'NE' then yrange=[0,2]
       if region_all(iregion) eq 'OH' then yrange=[0,2.5]     
       if region_all(iregion) eq 'SE' or region_all(iregion) eq 'Gulf'  then yrange=[0,4.]
       if region_all(iregion) eq 'Florida' then yrange=[0,7.]
     for iyear=0,nmonth/12-1 do begin
       if iyear eq 0 then plot, indgen( 12 ) , 0.03*reform(mod_region_interannual(iregion,iyear*12:iyear*12+11)), col=iyear+1, pos=pos,$
 	        ytitle='Modeled Monthly !cdeposition (' + micro + 'g m!u-2!nm!u-1!n)',thick=5,yrange=yrange,$
        	/xst, title='Region: ' + region_all(iregion), xrange=[0,11],xtitle='Month' else $
           oplot, indgen( 12 ) , 0.03*reform(mod_region_interannual(iregion,iyear*12:iyear*12+11)),col=iyear+1,thick=5
        legend,label=[string(4+indgen(nmonth/12-4)+select_year[0])],lcol=indgen(nmonth/12-4)+5,line=0*indgen(nmonth/12-4),symbol=0*indgen(nmonth/12-4)-1,boxcolor=-1,thick=5,$
 	             halign=0.01,valign=0.95,charsize=.5
     endfor
     multipanel,/advance,pos=pos
endfor

; plot the interannual variation, integrated to years
print, 'integrated to years'
multipanel,col=2,row=5,margin=[0.03,0.03],pos=pos
 !p.charsize=1.5

for iregion=0,nregion-1 do begin
       xrange=[1995, 2011]

       ydata=0.36525*shrink_data(mdn_region_interannual_site(iregion,*),12)
       print,region_all(iregion)
       print,ydata
       ydata=ydata-mean(ydata)
       yrange=[-4,4]
;       yrange=[0,15.]
;       if region_all(iregion) eq 'MW' then yrange=[5,15]
;       if region_all(iregion) eq 'NE' then yrange=[5,15]
;       if region_all(iregion) eq 'OH' then yrange=[5,20]     
;       if region_all(iregion) eq 'SE'  then yrange=[5,25]
;       if region_all(iregion) eq 'SE' or region_all(iregion) eq 'Gulf'  then yrange=[0,30]
;       if region_all(iregion) eq 'Florida' then yrange=[0,30]

       plot, indgen( ntyear ) + min(select_year), ydata,   col=1, pos=pos,$
 	        ytitle='Anomaly for annual!cdeposition (' + micro + 'g m!u-2!ny!u-1!n)' , yrange=yrange, thick=5, $
        	/xst, title='Region: ' + region_all(iregion), xrange=xrange,xtitle='Year', xminor=2 , yminor=2

        xdata=indgen(ntyear)+min(year)
        ydata=ydata[indgen(ntyear)+4]
        result=poly_fit(xdata,ydata,1,yfit=yfit,sigma=sigma)
        slope=result[1]
        intercept=result[0]
        oplot,[min(year),2010],intercept+[min(year),2010]*slope,col=1,line=1,thick=4
        t=slope / ( sqrt(total((ydata-yfit)^2)/(nmonth-2.0))/sqrt(total((xdata-mean(xdata))^2))  )

         xyouts,1999.5, 3.4,/data,col=1, 'MDN = ' + $
        string(intercept,format='(F7.1)') + '+ year x' + $
        string(slope,format='(F5.2)') + '!c' + $
        'Standard error of the slope: '+string(sigma[1],format='(F8.5)') + '!cr = ' + $
        string(rcorr,format='(F5.2)') + ' t = ' + $
        string(t,format='(F6.2)') + ' t-cutoff (0.05) = ' + $
        string(t_cvf(0.05,ntyear-2),format='(F5.2)'), charsize=0.7

;       print, 'MDN', 0.36525*shrink_data(mdn_region_interannual_site(iregion,*),12)
;       polyfill,[ indgen( nsyear*12 ), reverse(indgen( nsyear*12 )),0], $
;                [ reform(mdn_region_interannual_site(iregion,0:nsyear*12-1) - mdn_region_interannual_site_std(iregion,0:nsyear*12-1)), $
;                reverse(reform(mdn_region_interannual_site(iregion,0:nsyear*12-1)+mdn_region_interannual_site_std(iregion,0:nsyear*12-1))), $
;                 mdn_region_interannual_site[iregion,0]-mdn_region_interannual_site_std[iregion,0]]  $
;                ,col=!myct.lightgray
 
;       oplot,indgen( nmonth ),mdn_region_interannual_site(iregion,*),col=1,thick=1
        ydata=0.36525*shrink_data(reform(mod_region_interannual(iregion,*)),12)
       print,region_all(iregion), ' mod'
       print,ydata
        ydata=ydata-mean(ydata,/nan)
       oplot, indgen( ntyear ) + min(select_year)  , ydata ,col=!myct.red,thick=5
       ; print, 'GC', 0.36525*shrink_data(reform(mod_region_interannual(iregion,*)),12)
	   legend,label=['MDN', 'GEOS-Chem ' + resolution ],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=5,$
 	             halign=0.8,valign=0.05,charsize=.75
       multipanel,/advance,pos=pos
endfor

; plot the interannual variation, integrated to years, anomaly of (model - observation)
multipanel,col=2,row=5,margin=[0.03,0.03],pos=pos
 !p.charsize=1.5

for iregion=0,nregion-1 do begin
       yrange=[-3,3]

       xrange=[1995, 2011]
       bias_data=-(0.36525*shrink_data(reform(mod_region_interannual(iregion,*)),12) - 0.36525*shrink_data(mdn_region_interannual_site(iregion,*),12))
       bias_data_mean=mean(bias_data,/nan)
       bias_data=bias_data-bias_data_mean

       plot, indgen( ntyear ) + min(select_year), bias_data ,   col=1, pos=pos,$
 	        ytitle='Anomaly of model !c residue (' + micro + 'g m!u-2!ny!u-1!n)' , yrange=yrange, thick=5, $
        	/xst, title='Region: ' + region_all(iregion), xrange=xrange,xtitle='Year', xminor=1 , yminor=1

        rcorr=correlate(indgen(ntyear)+min(year),bias_data[indgen(ntyear)])
;        org_corr,indgen(7)+min(year),bias_data[indgen(7)+4],rcorr,7,slope,intercept
        xdata=indgen(ntyear)+min(year)
        ydata=bias_data[indgen(ntyear)]
        result=poly_fit(xdata,ydata,1,yfit=yfit,sigma=sigma)
        slope=result[1]
        intercept=result[0]
        t=slope / ( sqrt(total((ydata-yfit)^2)/(nmonth-2.0))/sqrt(total((xdata-mean(xdata))^2))  )
        oplot,!x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
         xyouts,2003.3,3,/data,col=1, 'MDN = ' + $
        string(intercept,format='(F7.1)') + '+ year x' + $
        string(slope,format='(F8.5)') + '!c' + $
        'Standard error of the slope: '+string(sigma[1],format='(F8.5)') + '!cr = ' + $
        string(rcorr,format='(F5.2)') + ' t = ' + $
        string(t,format='(F6.2)') + ' t-cutoff (0.05) = ' + $
        string(t_cvf(0.05,ntyear-2),format='(F5.2)'), charsize=0.7

       multipanel,/advance,pos=pos

endfor


; plot the interannual variation, integrated to seasons
if keyword_set(seasonal) then begin
for iseason=0,3 do begin
multipanel,col=2,row=5,margin=[0.03,0.03],pos=pos
 !p.charsize=1.5
for iregion=0,nregion-1 do begin
       ; yrange=[0,15.]
       ; if region_all(iregion) eq 'MW' then yrange=[0,5]
       ; if region_all(iregion) eq 'NE' then yrange=[0,5]
       ; if region_all(iregion) eq 'OH' then yrange=[0,7]     
       ; if region_all(iregion) eq 'SE'  then yrange=[0,8]

       xrange=[1995, 2011]
       if region_all(iregion) eq 'SE' or region_all(iregion) eq 'Gulf'  then yrange=[0,30]
       if region_all(iregion) eq 'Florida' then yrange=[0,30]
       
       mdn_season_data=fltarr(ntyear)+!values.f_nan
       mod_season_data=fltarr(ntyear)+!values.f_nan
       for iyear=0,ntyear-1 do begin
           mdn_season_data[iyear] = 0.36525*mean(mdn_region_interannual_site(iregion, iyear*12+season[iseason].month))/4.0
           if iyear+tyear[0] ge year[0] then mod_season_data[iyear] = 0.36525*mean(mod_region_interannual(iregion, iyear*12+season[iseason].month))/4.0
       endfor
       yrange=[min([mdn_season_data,mod_season_data]), max([mdn_season_data,mod_season_data])]

       plot, indgen( ntyear ) + min(select_year), mdn_season_data,   col=1, pos=pos,$
 	        ytitle='Seasonal deposition!c (' + micro + 'g m!u-2!nseason!u-1!n)' , yrange=yrange, thick=5, $
        	/xst, title='Season: '+season(iseason).name+ ' Region: ' + region_all(iregion), xrange=xrange,xtitle='Year', xminor=2 , yminor=2
;       print, 'MDN', 0.36525*shrink_data(mdn_region_interannual_site(iregion,*),12)
;       polyfill,[ indgen( nsyear*12 ), reverse(indgen( nsyear*12 )),0], $
;                [ reform(mdn_region_interannual_site(iregion,0:nsyear*12-1) - mdn_region_interannual_site_std(iregion,0:nsyear*12-1)), $
;                reverse(reform(mdn_region_interannual_site(iregion,0:nsyear*12-1)+mdn_region_interannual_site_std(iregion,0:nsyear*12-1))), $
;                 mdn_region_interannual_site[iregion,0]-mdn_region_interannual_site_std[iregion,0]]  $
;                ,col=!myct.lightgray
 
;       oplot,indgen( nmonth ),mdn_region_interannual_site(iregion,*),col=1,thick=1

       oplot, indgen( ntyear ) + min(select_year)  , mod_season_data , col=!myct.red,thick=5
;       print, 'GC', 0.36525*shrink_data(reform(mod_region_interannual(iregion,*)),12)
	   legend,label=['MDN', 'GEOS-Chem ' + resolution ],lcol=[1,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=5,$
 	             halign=0.01,valign=0.95,charsize=.75
       multipanel,/advance,pos=pos

endfor

endfor
endif

endif

; close_device
 print, '================================================'
print, ' Congratulations! Program terminate normally.'
print, '================================================'

  if keyword_Set(debug) then stop
;stop
end
