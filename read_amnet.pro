pro read_amnet,site,data, plot=plot, version,select_year=select_year,year=year,read=read, $
               model=model,ts=ts,enable_4x5=enable_4x5,enable_05x0666=enable_05x0666, $
			   version2=version2,new=new
			   
; select_year: year range for observation
; year: year range for model

; only read model year 2005 - 2007 !

; ====================================================================================
; read AMNET site observations
; ====================================================================================

 if not keyword_set(select_year) then select_year=[2000,2010]

  dir='./AMNet_data/'
  filename=dir+site+'.txt'
  if (not file_test(filename)) then begin
    print,site,' file not exist'
    return
  endif

  ; first find out number of lines in file
  nlines=file_lines(filename)
  
  ; Open and read the filen
  OPEN_FILE,filename,ilun
  data_str=strarr(nlines)
  readf,ilun,data_str
  free_lun,ilun
  nvar=17
  array=strarr(nlines,nvar)

  for k=0,nlines-1 do begin
    str=strsplit(data_str(k),",",count=count,/extract)
    if count gt 0 then array(k,0:count-1)=str
  endfor

  ftmp=fltarr(nlines-1) & ftmp(*)=!values.f_nan
  itmp=intarr(nlines-1) & itmp(*)=!values.f_nan
  stmp=strarr(nlines-1)
  data={siteid:'',year:itmp,month:itmp,day:itmp,julday:ftmp,gem:ftmp,gom:ftmp,pbm:ftmp, $
	 gem_monthly:fltarr(12),gem_monthly_std:fltarr(12),gem_num:fltarr(12),$
	 gom_monthly:fltarr(12),gom_monthly_std:fltarr(12),gom_num:fltarr(12),$
	 pbm_monthly:fltarr(12),pbm_monthly_std:fltarr(12),pbm_num:fltarr(12), $
     mod_gem_4x5: fltarr(12), mod_gom_4x5: fltarr(12), mod_pbm_4x5: fltarr(12), $
     mod_gem_4x52: fltarr(12), mod_gom_4x52: fltarr(12), mod_pbm_4x52: fltarr(12), $
     mod_gem_05x0666: fltarr(12), mod_gom_05x0666: fltarr(12), mod_pbm_05x0666: fltarr(12), $
     mod_gem_05x06662: fltarr(12), mod_gom_05x06662: fltarr(12), mod_pbm_05x06662: fltarr(12) $
        }

  ; read in data table
  data.siteid=array(1,0)
  for k=0,nlines-2 do begin
     if (array(k+1,2) ne 'NULL') then readdate,array(k+1,2),year1,month1,day1,hour1,minute1 else continue
     if (array(k+1,3) ne 'NULL') then readdate,array(k+1,3),year2,month2,day2,hour2,minute2 else continue

	 if year2 lt min(select_year) or year1 gt max(select_year) then continue
	 
     julday1=julday(month1,day1,year1,hour1,minute1,0.)
     julday2=julday(month2,day2,year2,hour2,minute2,0.)

     data.year(k) = year1
     data.month(k)= month1
     data.day(k)=day1
     data.julday(k)=0.5*(julday1+julday2)

     if (array(k+1,6) eq 'V') then begin
          data.pbm(k)=array(k+1,5)
          data.pbm_monthly(month1-1)=data.pbm_monthly(month1-1)+data.pbm(k)
          data.pbm_num(month1-1)=data.pbm_num(month1-1)+1.
     endif

     if (array(k+1,8) eq 'V') then begin
          data.gom(k)=array(k+1,7)
          data.gom_monthly(month1-1)=data.gom_monthly(month1-1)+data.gom(k)
          data.gom_num(month1-1)=data.gom_num(month1-1)+1.
     endif

     if (fix(array(k+1,10)) gt 0) then begin
          data.gem(k)=array(k+1,9)
          data.gem_monthly(month1-1)=data.gem_monthly(month1-1)+data.gem(k)
          data.gem_num(month1-1)=data.gem_num(month1-1)+1.
     endif          
  endfor

     ; calculate monthly average data
     FOR j=0,11 DO BEGIN
          data.gem_monthly(j)=data.gem_monthly(j)/data.gem_num(j)
          data.gom_monthly(j)=data.gom_monthly(j)/data.gom_num(j)
          data.pbm_monthly(j)=data.pbm_monthly(j)/data.pbm_num(j)
		  
          ind=where(data.month eq j+1)
          if (ind(0) ne -1) then begin
            data.gem_monthly_std(j)=stddev(data.gem(ind),/nan)    ; ng/m3
            data.gom_monthly_std(j)=stddev(data.gom(ind),/nan)    ; pg/m3
            data.pbm_monthly_std(j)=stddev(data.pbm(ind),/nan)    ; pg/m3
          endif
     endfor

; ====================================================================================
 ; read model data
; ====================================================================================
	 
 if keyword_set(model) then begin

 ; only when time series is plotted
 if keyword_set(ts) then begin
 ; get the daily concentration of the observation
 ; needed only when model data is also plotted
 ind_year=where(data.year ge min(year) and data.year le max(year))

 if isleap(year) then ndays=366 else ndays=365

 daily_gem=fltarr(ndays)
 daily_gom=fltarr(ndays)
 daily_pbm=fltarr(ndays)

 n_gem=fltarr(ndays)
 n_gom=fltarr(ndays)
 n_pbm=fltarr(ndays)

 julday_1stday=julday(1,1,year,0,0,0)
 for k=0,nlines-2 do begin
   if data.year[k] ne year then continue
   n_day=fix(data.julday[k]-julday_1stday)

   if data.gem[k] gt 0. then begin
     daily_gem[n_day]=daily_gem[n_day]+data.gem[k]
     n_gem[n_day]=n_gem[n_day]+1.
   endif
   if data.gom[k] gt 0. then begin
     daily_gom[n_day]=daily_gom[n_day]+data.gom[k]
     n_gom[n_day]=n_gom[n_day]+1.
   endif
   if data.pbm[k] gt 0. then begin
     daily_pbm[n_day]=daily_pbm[n_day]+data.pbm[k]
     n_pbm[n_day]=n_pbm[n_day]+1.
   endif
 endfor

 daily_gem=daily_gem/n_gem
 daily_gom=daily_gom/n_gom
 daily_pbm=daily_pbm/n_pbm
endif



  ; read monthly mean data
    conc_convert=200.6/(8.31*273.2)*1D2 ; convert from pptv to ng/m3 STP
    read_amnet_siteinfo,siteinfo, mode='c'
    nsites=n_elements(siteinfo)

	found=0
    for k=0,nsites-1 do begin
      ; print,site,siteinfo(k).siteid
      if (site eq siteinfo(k).siteid)then begin
         lat=siteinfo(k).lat
         lon=siteinfo(k).lon
		 ele=siteinfo(k).ele
		 sitename=siteinfo(k).sitename
		 found=1
         ; print,lat,lon
         break
      endif
    endfor
	if found eq 0 then return

  if keyword_set(enable_4x5) then begin

    mod_gem_4x5=fltarr(12)    ; ng/3
    mod_gom_4x5=fltarr(12)    ; pg/m3
    mod_pbm_4x5=fltarr(12)    ; pg/m3
	
	if keyword_set( version2 ) then begin
            mod_gem_4x52=fltarr(12)   ; ng/m3
            mod_gom_4x52=fltarr(12)   ; pg/m3
            mod_pbm_4x52=fltarr(12)   ; pg/m3	
	endif

   lev=alt2lev(lat,lon,ele,resolution='4x5')
	
	; read the data
    file='~/4-geos/0-hg/Run.4x5/dataout/ctm.bpch.'+version

	if keyword_set( version2 ) then begin
           file2='~/4-geos/0-hg/Run.4x5/dataout/ctm.bpch.'+version2
	endif

    for iyear=min(year),max(year) do begin	
;    for iyear=2005,2007 do begin	

       for month=1,12 do begin
         success=ctm_get_datablock(tmp1, 'ij-avg-$',  file=file,tracer=1, $
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1))) 
  
         if (iyear eq min(year)) then mod_gem_4x5(month-1)=tmp1*conc_convert else $
		                                       mod_gem_4x5(month-1)=mod_gem_4x5(month-1)+tmp1*conc_convert

        ; old simulation results
        if not keyword_set(new) then begin		 
         success=ctm_get_datablock(tmp2, 'ij-avg-$',file=file,tracer=2,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	     if (iyear eq min(year)) then mod_gom_4x5(month-1)=500.0*tmp2*conc_convert else $
                                                mod_gom_4x5(month-1)=mod_gom_4x5(month-1)+500.0*tmp2*conc_convert
												
   	     success=ctm_get_datablock(tmp3, 'ij-avg-$',file=file,tracer=3,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	      if (iyear eq min(year)) then mod_pbm_4x5(month-1)=  (0.5* tmp2 +  tmp3) * conc_convert*1000. else $
		                                         mod_pbm_4x5(month-1)=mod_pbm_4x5(month-1) + (0.5* tmp2 +  tmp3) * conc_convert*1000.
        endif else begin
         success=ctm_get_datablock(tmp2, 'pl-hg2-$',file=file,tracer=12,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	     if (iyear eq min(year)) then mod_gom_4x5(month-1)=1000.0*tmp2*conc_convert else $
                                                mod_gom_4x5(month-1)=mod_gom_4x5(month-1)+1000.0*tmp2*conc_convert
												
   	     success=ctm_get_datablock(tmp3, 'pl-hg2-$',file=file,tracer=11,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	      if (iyear eq min(year)) then mod_pbm_4x5(month-1)=  tmp3 * conc_convert*1000. else $
		                                         mod_pbm_4x5(month-1)=mod_pbm_4x5(month-1) + tmp3 * conc_convert*1000.
        endelse

       endfor
     endfor
	 
	 mod_gem_4x5=mod_gem_4x5/ (  max(year) - min(year) +1.0)     ; ng/m3
	 mod_gom_4x5=mod_gom_4x5/ (  max(year) - min(year) +1.0)     ; ng/m3
	 mod_pbm_4x5=mod_pbm_4x5/ (  max(year) - min(year) +1.0)     ; ng/m3

     data.mod_gem_4x5=mod_gem_4x5
     data.mod_gom_4x5=mod_gom_4x5
     data.mod_pbm_4x5=mod_pbm_4x5

	if keyword_set( version2 ) then begin
       	for iyear=min(year),max(year) do begin	

                for month=1,12 do begin
                  success=ctm_get_datablock(tmp1, 'ij-avg-$',file=file2,tracer=1,$
                                            lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
            								
                   if (iyear eq min(year)) then mod_gem_4x52(month-1)=tmp1*conc_convert else $
       		                                       mod_gem_4x52(month-1)=mod_gem_4x52(month-1)+tmp1*conc_convert
            	  
                  
        ; old simulation results
        if not keyword_set(new) then begin		 
         success=ctm_get_datablock(tmp2, 'ij-avg-$',file=file2,tracer=2,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	     if (iyear eq min(year)) then mod_gom_4x5(month-1)=500.0*tmp2*conc_convert else $
                                                mod_gom_4x5(month-1)=mod_gom_4x5(month-1)+500.0*tmp2*conc_convert
												
   	     success=ctm_get_datablock(tmp3, 'ij-avg-$',file=file2,tracer=3,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	      if (iyear eq min(year)) then mod_pbm_4x5(month-1)=  (0.5* tmp2 +  tmp3) * conc_convert*1000. else $
		                                         mod_pbm_4x5(month-1)=mod_pbm_4x5(month-1) + (0.5* tmp2 +  tmp3) * conc_convert*1000.
        endif else begin
         success=ctm_get_datablock(tmp2, 'pl-hg2-$',file=file2,tracer=12,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	     if (iyear eq min(year)) then mod_gom_4x52(month-1)=1000.0*tmp2*conc_convert else $
                                                mod_gom_4x52(month-1)=mod_gom_4x52(month-1)+1000.0*tmp2*conc_convert
												
   	     success=ctm_get_datablock(tmp3, 'pl-hg2-$',file=file2,tracer=11,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	      if (iyear eq min(year)) then mod_pbm_4x52(month-1)=  tmp3 * conc_convert*1000. else $
		                                         mod_pbm_4x52(month-1)=mod_pbm_4x52(month-1) + tmp3 * conc_convert*1000.
        endelse

                endfor
       	endfor
       
         mod_gem_4x52=mod_gem_4x52/ (  max(year) - min(year) +1.0)
       	 mod_gom_4x52=mod_gom_4x52/ (  max(year) - min(year) +1.0)
       	 mod_pbm_4x52=mod_pbm_4x52/ (  max(year) - min(year) +1.0)

         data.mod_gem_4x52=mod_gem_4x52
         data.mod_gom_4x52=mod_gom_4x52
         data.mod_pbm_4x52=mod_pbm_4x52

	endif

	 
	endif

  if keyword_set(enable_05x0666) then begin
;    conc_convert=200.6/(8.31*273.2)*1D2 ; convert from pptv to ng/m3 STP
 ;   read_amnet_siteinfo,siteinfo
 ;   nsites=n_elements(siteinfo)

    mod_gem_05x0666=fltarr(12)   ; ng/m3
    mod_gom_05x0666=fltarr(12)   ; pg/m3
    mod_pbm_05x0666=fltarr(12)   ; pg/m3
	
	if keyword_set( version2 ) then begin
            mod_gem_05x06662=fltarr(12)   ; ng/m3
            mod_gom_05x06662=fltarr(12)   ; pg/m3
            mod_pbm_05x06662=fltarr(12)   ; pg/m3	
	endif
    print , 'good'

    lev=alt2lev(lat,lon,ele,resolution='05x0666')
	
    file='~/4-geos/0-hg/Run.05x0666/dataout/ctm.bpch.'+version
	
	if keyword_set( version2 ) then begin
           file2='~/4-geos/0-hg/Run.05x0666/dataout/ctm.bpch.'+version2
	endif
	
	for iyear=min(year),max(year) do begin	
;	for iyear=2005,2007 do begin	
    
         for month=1,12 do begin
           success=ctm_get_datablock(tmp1, 'ij-avg-$',file=file,tracer=1,$
                                     lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
     								
            if (iyear eq min(year)) then mod_gem_05x0666(month-1)=tmp1*conc_convert else $
		                                       mod_gem_05x0666(month-1)=mod_gem_05x0666(month-1)+tmp1*conc_convert
     	  
        ; old simulation results
        if not keyword_set(new) then begin		 
           success=ctm_get_datablock(tmp2, 'ij-avg-$',file=file,tracer=2,$
                                     lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
     								
    	     if (iyear eq min(year)) then mod_gom_05x0666(month-1)=500.0*tmp2*conc_convert else $
                                                mod_gom_05x0666(month-1)=mod_gom_05x0666(month-1)+500.0*tmp2*conc_convert
     	  
           success=ctm_get_datablock(tmp3, 'ij-avg-$',file=file,tracer=3,$
                                     lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
     								
            if (iyear eq min(year)) then mod_pbm_05x0666(month-1)=(0.5* tmp2 +  tmp3) * conc_convert*1000. else $
		                                         mod_pbm_05x0666(month-1)=mod_pbm_05x0666(month-1) + (0.5* tmp2 +  tmp3) * conc_convert*1000.
        endif else begin
         success=ctm_get_datablock(tmp2, 'pl-hg2-$',file=file,tracer=12,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	     if (iyear eq min(year)) then mod_gom_05x0666(month-1)=1000.0*tmp2*conc_convert else $
                                                mod_gom_05x0666(month-1)=mod_gom_05x0666(month-1)+1000.0*tmp2*conc_convert
												
   	     success=ctm_get_datablock(tmp3, 'pl-hg2-$',file=file,tracer=11,$
                                   lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
   
   	      if (iyear eq min(year)) then mod_pbm_05x0666(month-1)=  tmp3 * conc_convert*1000. else $
		                                         mod_pbm_05x0666(month-1)=mod_pbm_05x0666(month-1) + tmp3 * conc_convert*1000.
        endelse

         endfor
	endfor

 	 mod_gem_05x0666=mod_gem_05x0666/ (  max(year) - min(year) +1.0)
	 mod_gom_05x0666=mod_gom_05x0666/ (  max(year) - min(year) +1.0)
	 mod_pbm_05x0666=mod_pbm_05x0666/ (  max(year) - min(year) +1.0)

     data.mod_gem_05x0666=mod_gem_05x0666
     data.mod_gom_05x0666=mod_gom_05x0666
     data.mod_pbm_05x0666=mod_pbm_05x0666


	if keyword_set( version2 ) then begin
       	for iyear=min(year),max(year) do begin	

                for month=1,12 do begin
                  success=ctm_get_datablock(tmp1, 'ij-avg-$',file=file2,tracer=1,$
                                            lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
            								
                   if (iyear eq min(year)) then mod_gem_05x06662(month-1)=tmp1*conc_convert else $
       		                                       mod_gem_05x06662(month-1)=mod_gem_05x06662(month-1)+tmp1*conc_convert
            	  
                     ; old simulation results
                     if not keyword_set(new) then begin		 
                        success=ctm_get_datablock(tmp2, 'ij-avg-$',file=file2,tracer=2,$
                                                  lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
                  								
                 	     if (iyear eq min(year)) then mod_gom_05x06662(month-1)=500.0*tmp2*conc_convert else $
                                                             mod_gom_05x06662(month-1)=mod_gom_05x06662(month-1)+500.0*tmp2*conc_convert
                  	  
                        success=ctm_get_datablock(tmp3, 'ij-avg-$',file=file2,tracer=3,$
                                                  lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
                  								
                         if (iyear eq min(year)) then mod_pbm_05x06662(month-1)=(0.5* tmp2 +  tmp3) * conc_convert*1000. else $
             		                                         mod_pbm_05x06662(month-1)=mod_pbm_05x06662(month-1) + (0.5* tmp2 +  tmp3) * conc_convert*1000.
                     endif else begin
                      success=ctm_get_datablock(tmp2, 'pl-hg2-$',file=file2,tracer=12,$
                                                lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
                
                	     if (iyear eq min(year)) then mod_gom_05x06662(month-1)=1000.0*tmp2*conc_convert else $
                                                             mod_gom_05x06662(month-1)=mod_gom_05x06662(month-1)+1000.0*tmp2*conc_convert
             												
                	     success=ctm_get_datablock(tmp3, 'pl-hg2-$',file=file2,tracer=11,$
                                                lat=lat,lon=lon,lev=lev,tau0=nymd2tau(ymd2date(iyear,month,1)))
                
                	      if (iyear eq min(year)) then mod_pbm_05x06662(month-1)=  tmp3 * conc_convert*1000. else $
             		                                         mod_pbm_05x06662(month-1)=mod_pbm_05x06662(month-1) + tmp3 * conc_convert*1000.
                     endelse

                endfor
       	endfor
       
         mod_gem_05x06662=mod_gem_05x06662/ (  max(year) - min(year) +1.0)
       	 mod_gom_05x06662=mod_gom_05x06662/ (  max(year) - min(year) +1.0)
       	 mod_pbm_05x06662=mod_pbm_05x06662/ (  max(year) - min(year) +1.0)

         data.mod_gem_05x06662=mod_gem_05x06662
         data.mod_gom_05x06662=mod_gom_05x06662
         data.mod_pbm_05x06662=mod_pbm_05x06662

	endif
	 
  endif

   ; read time series data
   if keyword_set(ts)  then begin

   if keyword_set(read) then begin
   mod_julday=fltarr(ndays)
   mod_gem_ts_4x5=fltarr(ndays)
   mod_gom_ts_4x5=fltarr(ndays)
   mod_pbm_ts_4x5=fltarr(ndays)
   mod_gem_ts_05x0666=fltarr(ndays)
   mod_gom_ts_05x0666=fltarr(ndays)
   mod_pbm_ts_05x0666=fltarr(ndays)

   idata=0
   for iday=0,ndays-1 do begin
     itau=nymd2tau(ymd2date(year,1,1))+24.*iday
     dummy=tau2yymmdd(itau)
     if dummy.day eq 1 then ctm_cleanup
     mod_julday[iday]=julday(dummy.month,dummy.day,dummy.year,12,0.,0.)
     if keyword_set(enable_4x5) then begin
       mod_gem_ts_4x5[iday]=get_conc(version,'4x5',lat,lon,dummy.year,dummy.month,dummy.day,1)
       mod_gom_ts_4x5[iday]=0.5*get_conc(version,'4x5',lat,lon,dummy.year,dummy.month,dummy.day,2)
       mod_pbm_ts_4x5[iday]=mod_gom_ts_4x5[iday]+get_conc(version,'4x5',lat,lon,dummy.year,dummy.month,dummy.day,3)
     endif

     if keyword_set(enable_05x0666) then begin
       mod_gem_ts_05x0666[iday]=get_conc(version,'05x0666',lat,lon,dummy.year,dummy.month,dummy.day,1)
       mod_gom_ts_05x0666[iday]=0.5*get_conc(version,'05x0666',lat,lon,dummy.year,dummy.month,dummy.day,2)
       mod_pbm_ts_05x0666[iday]=mod_gom_ts_05x0666[iday]+get_conc(version,'05x0666',lat,lon,dummy.year,dummy.month,dummy.day,3)
     endif
   endfor
   
     ; save data for next time running
     save,mod_julday,mod_gem_ts_4x5,mod_gom_ts_4x5,mod_pbm_ts_4x5, $
        mod_gem_ts_05x0666,mod_gom_ts_05x0666,mod_pbm_ts_05x0666, $
        filename='TS_'+version+'_4x5_'+string(year,format='(I4)')+'_'+$
                 site+'.sav'
   endif else begin
      ; restore the saved ts file
      restore,filename='TS_'+version+'_4x5_'+string(year,format='(I4)')+'_'+site+'.sav'

   endelse


   endif

 endif

  if keyword_set(plot) then begin

     open_device,/ps,/color,/portrait,file='./sites/AMNET_TGM+TOM_'+site+'_'+version+'_'+string(min(year),format='(I4)')+'-'+string(max(year),format='(I4)')+ '.ps'
     !P.font=0
     Device, /Helvetica, /Isolatin1
 ;    multipanel,col=2,row=3,margin=[0.05,0.05],pos=pos
     !p.charsize=0.8
     !P.MULTI  = [0 ,2 ,2]

     ; dummy=label_date(date_format='%M 1!C%Y')
        ; hg0
     ;	plot,julday_1stday+indgen(ndays),daily_gem,col=1,$
     ;		psym=-4,symsize=0.3,xtickformat='label_date',$
     ;  	ytitle='GEM (ng/m!u3!n)',pos=pos,/xst,xtickinterval=30.,$
     ;      yrange=[0.8,2.8]
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_4x5) then oplot,mod_julday,conc_convert*mod_gem_ts_4x5,color=!myct.red,thick=1
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_05x0666) then oplot,mod_julday,conc_convert*mod_gem_ts_05x0666,color=!myct.green,thick=1
     ;
     ;  xyouts,0.5,0.90,/normal,data.siteid,col=!myct.black,al=0.5,charsize=1.1
     ;  multipanel,/advance,pos=pos
     ;  
     ;  ; gom
     ;	plot,julday_1stday+indgen(ndays),daily_gom,col=1,$
     ;		psym=-4,symsize=.3,xtickformat='label_date',$
     ;  	ytitle='GOM (pg/m!u3!n)',pos=pos,/xst,xtickinterval=30.,$
     ;      yrange=[0.,100.]
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_4x5) then oplot,mod_julday,1000.*conc_convert*mod_gom_ts_4x5,color=!myct.red,thick=1
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_05x0666) then oplot,mod_julday,1000.*conc_convert*mod_gom_ts_05x0666,color=!myct.green,thick=1
     ;
     ;  multipanel,/advance,pos=pos
     ;
     ;  ; pbm
     ;	plot,julday_1stday+indgen(ndays),daily_pbm,col=1,$
     ;		psym=-4,symsize=.3,xtickformat='label_date',$
     ;  	ytitle='PBM (pg/m!u3!n)',pos=pos,/xst,xtickinterval=30.,$
     ;      yrange=[0.,100.]
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_4x5) then  oplot,mod_julday,1000.*conc_convert*mod_pbm_ts_4x5,color=!myct.red,thick=1
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_05x0666) then oplot,mod_julday,1000.*conc_convert*mod_pbm_ts_05x0666,color=!myct.green,thick=1
     ;
     ;  multipanel,/advance,pos=pos

        ; tgm
     ;	plot,julday_1stday+indgen(ndays),daily_gem+daily_gom/1000.,col=1,$
     ;	     psym=-4,symsize=.3,xtickformat='label_date', $
     ;   	 ytitle='TGM (ng/m!u3!n)',pos=pos,/xst,xtickinterval=30., $
     ;       yrange=[0.8,2.8]
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_4x5) then  oplot,mod_julday,conc_convert*(mod_gem_ts_4x5+mod_pbm_ts_4x5),color=!myct.red,thick=1
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_05x0666) then oplot,mod_julday,conc_convert*(mod_gem_ts_05x0666+mod_pbm_ts_05x0666),color=!myct.green,thick=1
     ;  
     ;  multipanel,/advance,pos=pos
     ;
     ;  ; TOM = RGM+PBM
     ;	plot,julday_1stday+indgen(ndays),daily_gom+daily_pbm,col=1,$
     ;	     psym=-4,symsize=.3,xtickformat='label_date', $
     ;   	 ytitle='TOM (pg/m!u3!n)',pos=pos,/xst,xtickinterval=30., $
     ;       yrange=[0.,140.]
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_4x5) then  oplot,mod_julday,1000.*conc_convert*(mod_gom_ts_4x5+mod_pbm_ts_4x5),color=!myct.red,thick=1
     ;  if keyword_set(ts) and keyword_set(model) and keyword_set(enable_05x0666) then oplot,mod_julday,1000.*conc_convert*(mod_gom_ts_05x0666+mod_pbm_ts_05x0666),color=!myct.green,thick=1
     ;  
     ;  multipanel,/advance,pos=pos
        
      ;  multipanel,col=3,row=7,omargin=0.1,margin=[0.03,0.02],pos=pos,/noerase
      ;
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos
      ;  multipanel,/advance,pos=pos

     ;	plot,indgen(12)+1,data.gem_monthly,$
     ;		symsize=1.0,color=!myct.black,thick=4,$
     ;   	ytitle='GEM (ng/m!u3!n)',pos=pos,yrange=[1,2]
     ;   if keyword_set(model) then begin
     ;       if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gem_4x5,symsize=1.0,color=!myct.red,thick=4
	;		if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gem_05x0666,symsize=1.0,color=!myct.green,thick=4
     ;   endif        
     ;   multipanel,/advance,pos=pos
     ;
     ;	plot,indgen(12)+1,data.gom_monthly,$
     ;		symsize=1.0,color=!myct.black,thick=4,$
     ;   	ytitle='GOM (pg/m!u3!n)',pos=pos,yrange=[0,40]
     ;   if keyword_set(model) then begin
     ;       if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gom_4x5,symsize=1.0,color=!myct.red,thick=4
     ;       if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gom_05x0666,symsize=1.0,color=!myct.green,thick=4
     ;   endif        
     ;   multipanel,/advance,pos=pos
     ;
     ;	plot,indgen(12)+1,data.pbm_monthly,$
     ;		symsize=1.0,color=!myct.black,thick=4,$
     ;   	ytitle='PBM (pg/m!u3!n)',pos=pos,yrange=[0,40]
     ;   if keyword_set(model) then begin
     ;       if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_pbm_4x5,symsize=1.0,color=!myct.red,thick=4
     ;       if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_pbm_05x0666,symsize=1.0,color=!myct.green,thick=4
     ;   endif
     ;   multipanel,/advance,pos=pos
	 
	   ; tom
;   	    tom_monthly=(data.pbm_monthly+data.gom_monthly)
;		tom_monthly_std=sqrt(data.gom_monthly_std^2+data.pbm_monthly_std^2)
;		
		title=data.siteid+ ' ' + sitename + ' ' + $
                   string(lat,format='(F5.2)')+'!uo!nN;'+$
	             string(lon,format='(F7.2)')+'!uo!nE'
;	             
;    	plot,indgen(12)+1,tom_monthly,$
;    		symsize=1.0,color=!myct.black,thick=4,$
;       	    ytitle='TOM (pg/m!u3!n)',pos=pos,yrange=[0,200.],$
;		    xrange=[0.4,12.6],/xst,title= title, $
;		    xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
;		    xtickv=indgen(12)+1,xminor=1
;
;       errorbar,indgen(12)+1,tom_monthly, tom_monthly_std,col=1,thick=4
;
;       oplot,indgen(12)+1,tom_monthly,col=1,thick=4	 
;		if keyword_set(model) then begin
;            if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gom_4x5+mod_pbm_4x5,color=!myct.red,thick=4
;            if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gom_05x0666+mod_pbm_05x0666,color=!myct.green,thick=4
;        endif
;		legend,label=['AMNET'+string(mean(tom_monthly,/nan),stddev(tom_monthly,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ),  $
;		                    'GEOS-Chem 4x5'  + string(mean(mod_gom_4x5+mod_pbm_4x5,/nan),stddev(mod_gom_4x5+mod_pbm_4x5,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ), $
;							'GEOS-Chem 05x0666'  + string(mean(mod_gom_05x0666+mod_pbm_05x0666,/nan),stddev(mod_gom_05x0666+mod_pbm_05x0666,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ) ], $
;		lcol=[1,!myct.red,!myct.green],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=4,$
; 	                           halign=0.01,valign=0.95,charsize=0.8
;	 
;        multipanel,/advance,pos=pos

;       oplot,indgen(12)+1,oneyear_gem,col=1,thick=5	 
		
       ; polyfill,[1+indgen(12),reverse(1+indgen(12)),1],$
        ;     [data.gem_monthly-data.gem_monthly_std,reverse(data.gem_monthly+data.gem_monthly_std),data.gem_monthly[0]-data.gem_monthly_std[0]],col=!myct.lightblack
		
		
;      plot,indgen(12)+1,oneyear_gem,$
;     		symsize=1.0,color=!myct.black,thick=4,$
;        	ytitle='GEM (ng/m!u3!n)',pos=pos,yrange=[min(oneyear_gem-oneyear_gem_std),max(oneyear_gem+oneyear_gem_std)],xrange=[0.4,12.6],/xst,title=data.siteid+'_'+string(year,format='(I4)'),$
;			xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
;			xtickv=indgen(12)+1,xminor=1
;		
;        polyfill,[1+indgen(12),reverse(1+indgen(12)),1],$
;             [oneyear_gem-oneyear_gem_std,reverse(oneyear_gem+oneyear_gem_std),oneyear_gem(0)-oneyear_gem_std(0)],col=!myct.lightblack


		; tgm
	    tgm_monthly=data.gem_monthly+data.gom_monthly/1000.   ; ng/m3
		tgm_monthly_std=sqrt(data.gem_monthly_std^2+(data.gom_monthly_std/1000.)^2)     ; ng/m3

     	plot,indgen(12)+1,tgm_monthly,$
     		symsize=1.0,color=!myct.black,thick=4,$
        	ytitle='TGM (ng/m!u3!n)',pos=[0.3500050  ,   0.83   ,  0.750005  ,   0.95000],yrange=[0.5,3.0],  $
			            xrange=[0.4,12.6],/xst,title=title, $
			            xtitle=' ',xticks=11,xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '],$
			            xtickv=indgen(12)+1,xminor=1
		
        ; polyfill,[1+indgen(12),reverse(1+indgen(12)),1],$
         ;    [tgm_monthly-tgm_monthly_std,reverse(tgm_monthly+tgm_monthly_std),tgm_monthly[0]-tgm_monthly_std[0]],col=!myct.lightblack
        errorbar,indgen(12)+1,tgm_monthly, tgm_monthly_std,col=!myct.black,thick=1

       xyouts,9.7, 2.7,string(min(select_year),format='(I4)')+'-'+string(max(select_year),format='(I4)'),charsize=0.8,/data,col=1

;        oneyear_tgm=oneyear_gem+oneyear_gom/1000.
;		oneyear_tgm_std=sqrt(oneyear_gem_std^2+oneyear_gom_std^2/1.e6)
;        plot,indgen(12)+1,oneyear_tgm,$
;     		symsize=1.0,color=!myct.black,thick=4,$
;        	ytitle='TGM (ng/m!u3!n)',pos=pos,yrange=[min(oneyear_tgm-oneyear_tgm_std),max(oneyear_tgm+oneyear_tgm_std)],xrange=[0.4,12.6],/xst,title=data.siteid+'_'+string(year,format='(I4)'),$
;			xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
;			xtickv=indgen(12)+1,xminor=1
;		
;        polyfill,[1+indgen(12),reverse(1+indgen(12)),1],$
;             [oneyear_tgm-oneyear_tgm_std,reverse(oneyear_tgm+oneyear_tgm_std),oneyear_tgm(0)-oneyear_tgm_std(0)],col=!myct.lightblack

		oplot,indgen(12)+1,tgm_monthly,col=1,thick=4
;		oplot,indgen(12)+1,oneyear_tgm,col=1,thick=5	 
		if keyword_set(model) then begin
            if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gem_4x5+mod_gom_4x5/1000.,color=!myct.red,thick=4
            if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gem_05x0666+mod_gom_05x0666/1000.,color=!myct.red,thick=4
			if keyword_set(version2) then begin
                 if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gem_4x52+mod_gom_4x52/1000.,color=!myct.green,thick=4
                 if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gem_05x06662+mod_gom_05x06662/1000.,color=!myct.green,thick=4
            endif
        endif
		
		label_accumulate=['AMNET '+string(mean(tgm_monthly,/nan),stddev(tgm_monthly,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ) ]
		lcol_accumulate=1
        line_accumulate=0
        symbol_accumulate=-1

        if keyword_set(enable_4x5) then begin
           label_accumulate=[label_accumulate, 'GC4 '+ string(mean(mod_gem_4x5+mod_gom_4x5/1000.,/nan), $
                                     stddev(mod_gem_4x5+mod_gom_4x5/1000.,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' )]
           lcol_accumulate=[lcol_accumulate,!myct.red]
           line_accumulate=[line_accumulate,0]
           symbol_accumulate=[symbol_accumulate,-1]
        endif
        
		if keyword_set(enable_05x0666) then begin
              	label_accumulate=[label_accumulate, 'STD simulation' + string(mean(mod_gem_05x0666+mod_gom_05x0666/1000.,/nan),$
                                     stddev(mod_gem_05x0666+mod_gom_05x0666/1000.,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ) ]
           lcol_accumulate=[lcol_accumulate,!myct.red]
           line_accumulate=[line_accumulate,0]
           symbol_accumulate=[symbol_accumulate,-1]
        endif

		if keyword_set (version2) then begin
           if keyword_set(enable_05x0666) then $
               label_accumulate = [		label_accumulate , $
			   version2 + string(mean(mod_gem_05x06662+mod_gom_05x06662/1000.,/nan),stddev(mod_gem_05x06662+mod_gom_05x06662/1000.,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ) ]
           if keyword_set(enable_4x5) then $
               label_accumulate = [		label_accumulate , $
			   version2 + string(mean(mod_gem_4x52+mod_gom_4x52/1000.,/nan),stddev(mod_gem_4x52+mod_gom_4x52/1000.,/nan),format= '(F4.2," +/- ",F4.2," ng/m!u3!n")' ) ]

			   lcol_accumulate = [lcol_accumulate, !myct.green ]
			   line_accumulate = [ line_accumulate, 0]
			   symbol_accumulate = [symbol_accumulate, -1]
		endif

;		legend,label=label_accumulate, $
;							lcol=lcol_accumulate,line=line_accumulate,symbol=symbol_accumulate,boxcolor=-1,thick=4,$
; 	                           halign=0.01,valign=0.95,charsize=.8
	
;        multipanel,/advance,pos=pos

     ;	plot,indgen(12)+1,data.pbm_monthly+data.gom_monthly,$
     ;		symsize=1.0,color=!myct.black,thick=4,$
     ;  	ytitle='TOM (pg/m!u3!n)',pos=pos,yrange=[0,80]
     ;  if keyword_set(model) then begin
     ;      if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_pbm_4x5+mod_gom_4x5,symsize=1.0,color=!myct.red,thick=4
     ;      if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_pbm_05x0666+mod_gom_05x0666,symsize=1.0,color=!myct.green,thick=4
     ;  endif
     ;  multipanel,/advance,pos=pos
     ;
     ;	plot,indgen(12)+1,data.pbm_monthly/(data.pbm_monthly+data.gom_monthly),$
     ;		symsize=1.0,color=!myct.black,thick=4,$
     ;  	ytitle='Part. Frac.',pos=pos,yrange=[0,1]
        
     ; close_Device
 ;    multipanel,/advance,pos=pos
  ;   multipanel,/advance,pos=pos

	 ; open_device,/ps,/color,/portrait,file=' AMNET_TOM_'+site+'_'+version+'_'+string(min(year),format='(I4)')+'-'+string(max(year),format='(I4)')+ '.ps'
     ; !P.font=0
     ; Device, /Helvetica, /Isolatin1
     ; multipanel,col=2,row=3,margin=[0.05,0.05],pos=pos
     ; !p.charsize=1.5


 ; ===============================================================================================
  	    gom_monthly=data.gom_monthly
		gom_monthly_std=sqrt(data.gom_monthly_std^2)	
	                                                    
    	plot,indgen(12)+1,gom_monthly,$
    		symsize=1.0,color=!myct.black,thick=4,$
       	    ytitle='RGM (pg/m!u3!n)',pos=[0.3500050 ,    0.705   ,  0.750005   ,  0.825],yrange=[0,40.],$
		    xrange=[0.4,12.6],/xst, $
            xtitle=' ',xticks=11,xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '],$
            xtickv=indgen(12)+1,xminor=1
                                                          
       errorbar,indgen(12)+1,gom_monthly, gom_monthly_std,col=!myct.black,thick=1

       oplot,indgen(12)+1,gom_monthly,col=1,thick=4	 
		if keyword_set(model) then begin
            if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gom_4x5,color=!myct.red,thick=4
            if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gom_05x0666,color=!myct.red,thick=4
			if keyword_set(version2) then begin
                    if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gom_05x06662,color=!myct.green,thick=4
                    if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gom_4x52,color=!myct.green,thick=4
            endif
        endif

		label_accumulate='AMNET '+string(mean(gom_monthly,/nan),stddev(gom_monthly,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' )
		lcol_accumulate=1
        line_accumulate=0
        symbol_accumulate=-1

        if keyword_set(enable_4x5) then begin
           label_accumulate=[label_accumulate,  'GC4 '  + string(mean(mod_gom_4x5,/nan),$
                                          stddev(mod_gom_4x5,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' )]
           lcol_accumulate=[lcol_accumulate,!myct.red]
           line_accumulate=[line_accumulate,0]
           symbol_accumulate=[symbol_accumulate,-1]
        endif
        
		if keyword_set(enable_05x0666) then begin
              	label_accumulate=[label_accumulate, 'STD simulation' + string(mean(mod_gom_05x0666,/nan),$
                                           stddev(mod_gom_05x0666,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' )]
           lcol_accumulate=[lcol_accumulate,!myct.red]
           line_accumulate=[line_accumulate,0]
           symbol_accumulate=[symbol_accumulate,-1]
        endif

		if keyword_set (version2) then begin
              if keyword_set(enable_05x0666) then $
               label_accumulate = [		label_accumulate , $
			  version2 + string(mean(mod_gom_05x06662,/nan),stddev(mod_gom_05x06662,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ) ]
              if keyword_set(enable_4x5) then $
               label_accumulate = [		label_accumulate , $
			  version2 + string(mean(mod_gom_4x52,/nan),stddev(mod_gom_4x52,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ) ]

			   lcol_accumulate = [lcol_accumulate, !myct.green ]
			   line_accumulate = [ line_accumulate, 0]
			   symbol_accumulate = [symbol_accumulate, -1]
		endif

;		legend,label=label_accumulate, $
;							lcol=lcol_accumulate,line=line_accumulate,symbol=symbol_accumulate,boxcolor=-1,thick=4,$
; 	                           halign=0.01,valign=0.95,charsize=.8
         
 ; ===============================================================================================
  	    pbm_monthly=data.pbm_monthly
		pbm_monthly_std=sqrt(data.pbm_monthly_std^2)	
	             
    	plot,indgen(12)+1,pbm_monthly,$
    		symsize=1.0,color=!myct.black,thick=4,$
       	    ytitle='PBM (pg/m!u3!n)',pos= [0.3500050  ,   0.58   ,  0.750005  ,   0.70000],yrange=[0,40.],$
		    xrange=[0.4,12.6],/xst, $   ; title= title, $
		    xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
		    xtickv=indgen(12)+1,xminor=1

       errorbar,indgen(12)+1,pbm_monthly, pbm_monthly_std,col=!myct.black,thick=1

       oplot,indgen(12)+1,pbm_monthly,col=1,thick=4	 
		if keyword_set(model) then begin
            if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_pbm_4x5,color=!myct.red,thick=4
            if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_pbm_05x0666,color=!myct.red,thick=4
			if keyword_set(version2) then begin
                    if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_pbm_05x06662,color=!myct.green,thick=4
                    if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_pbm_4x52,color=!myct.green,thick=4
            endif
        endif

		label_accumulate='AMNET '+string(mean(pbm_monthly,/nan),stddev(pbm_monthly,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' )
		lcol_accumulate=1
        line_accumulate=0
        symbol_accumulate=-1

        if keyword_set(enable_4x5) then begin
           label_accumulate=[label_accumulate,  'GC4 '  + string(mean(mod_pbm_4x5,/nan),$
                                          stddev(mod_pbm_4x5,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' )]
           lcol_accumulate=[lcol_accumulate,!myct.red]
           line_accumulate=[line_accumulate,0]
           symbol_accumulate=[symbol_accumulate,-1]
        endif
        
		if keyword_set(enable_05x0666) then begin
              	label_accumulate=[label_accumulate, 'STD simulation' + string(mean(mod_pbm_05x0666,/nan),$
                                           stddev(mod_pbm_05x0666,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' )]
           lcol_accumulate=[lcol_accumulate,!myct.red]
           line_accumulate=[line_accumulate,0]
           symbol_accumulate=[symbol_accumulate,-1]
        endif

		if keyword_set (version2) then begin
              if keyword_set(enable_05x0666) then $
               label_accumulate = [		label_accumulate , $
			  version2 + string(mean(mod_pbm_05x06662,/nan),stddev(mod_pbm_05x06662,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ) ]
              if keyword_set(enable_4x5) then $
               label_accumulate = [		label_accumulate , $
			  version2 + string(mean(mod_pbm_4x52,/nan),stddev(mod_pbm_4x52,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ) ]

			   lcol_accumulate = [lcol_accumulate, !myct.green ]
			   line_accumulate = [ line_accumulate, 0]
			   symbol_accumulate = [symbol_accumulate, -1]
		endif

;		legend,label=label_accumulate, $
;							lcol=lcol_accumulate,line=line_accumulate,symbol=symbol_accumulate,boxcolor=-1,thick=4,$
; 	                           halign=0.01,valign=0.95,charsize=.8

 ;    print,pos
     close_Device
;	 spawn,'convert -density 500 psfile.ps file.png' 
;	  file_move, 'file.png',' AMNET_TOM_'+site+'_'+version+'_'+string(min(year),format='(I4)')+'-'+string(max(year),format='(I4)')+ '.png'
	 
	 
;	  open_device,/ps,/color,/portrait,file=' AMNET_GOM_'+site+'_'+version+'_'+string(min(year),format='(I4)')+'-'+string(max(year),format='(I4)')+ '.ps'
;    !P.font=0
;    Device, /Helvetica, /Isolatin1
;    multipanel,col=2,row=3,margin=[0.05,0.05],pos=pos
;    !p.charsize=1.5
;
; 	    gom_monthly=data.gom_monthly
;		gom_monthly_std=data.gom_monthly_std
;		
;	             
;   	plot,indgen(12)+1,gom_monthly,$
;   		symsize=1.0,color=!myct.black,thick=4,$
;      	    ytitle='GOM (pg/m!u3!n)',pos=pos,yrange=[0,100.],$
;		    xrange=[0.4,12.6],/xst,title= title, $
;		    xtitle='Month',xticks=11,xtickname=['J','F','M','A','M','J','J','A','S','O','N','D'],$
;		    xtickv=indgen(12)+1,xminor=1
;
;      errorbar,indgen(12)+1,gom_monthly, gom_monthly_std,col=!myct.black,thick=1
;
;      oplot,indgen(12)+1,gom_monthly,col=1,thick=4	 
;		if keyword_set(model) then begin
;           if keyword_set(enable_4x5) then oplot,indgen(12)+1,mod_gom_4x5,color=!myct.red,thick=4
;           if keyword_set(enable_05x0666) then oplot,indgen(12)+1,mod_gom_05x0666,color=!myct.green,thick=4
;			if keyword_set(version2) then oplot,indgen(12)+1,mod_gom_05x06662,color=!myct.green,thick=4,linestyle=2
;       endif
		
;		label_accumulate=['AMNET '+string(mean(gom_monthly,/nan),stddev(gom_monthly,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ),  $
;		                    'GEOS-Chem 4x5 '  + string(mean(mod_gom_4x5,/nan),stddev(mod_gom_4x5,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ), $
;							'GEOS-Chem 05x0666 '  + string(mean(mod_gom_05x0666,/nan),stddev(mod_gom_05x0666,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ) ]
;		lcol_accumulate=[1,!myct.red,!myct.green]
;        line_accumulate=[0,0,0]
;        symbol_accumulate=	[-1,-1,-1]	
;
;		if keyword_set (version2) then begin
;               label_accumulate = [		label_accumulate , $
;			   'GC05 w/o inplume red. ' + string(mean(mod_gom_05x06662,/nan),stddev(mod_gom_05x06662,/nan),format= '(F5.1," +/- ",F5.1," pg/m!u3!n")' ) ]
;			   lcol_accumulate = [lcol_accumulate, !myct.green ]
;			   line_accumulate = [ line_accumulate, 2]
;			   symbol_accumulate = [symbol_accumulate, -1]
;		endif
;		legend,label=label_accumulate, $
;							lcol=lcol_accumulate,line=line_accumulate,symbol=symbol_accumulate,boxcolor=-1,thick=4,$
; 	                           halign=0.01,valign=0.95,charsize=.8

 ;    close_Device
     ; file_move,'idl.ps','psfile.ps',/overwrite    ;this will rename the file.png into any other name you wish
    ; spawn,'convert -density 500 psfile.ps file.png'            ; this will transfer the ps file into a png file named file.png
	; spawn,'convert -extract 3900x1450-70+140 file.png out.png'
    ; if keyword_set(model) then begin
    ;   file_move,'out.png','AMNET_'+site+'_'+version+'_'+string(min(year),format='(I4)')+'-'+string(max(year),format='(I4)')+ '.png',/overwrite    ; this will rename the file.png into any other name you wish
    ; endif else begin 
    ;   file_move,'out.png','AMNET_'+site+'.png',/overwrite    ; this will rename the file.png into any other name you wish
    ; endelse
  
  endif

end

