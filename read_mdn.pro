pro read_mdn, data, select_year=select_year, snow=snow, mettype=mettype, metresolution=metresolution, correct=correct

; +++++++++++++++++++++++++++++++++++++++++++++
;
; Program reading in the MDN data from ~/ALASKA/data/MDN/mdnic.txt
; file containing all MDF data (downloaded in Oct 2009) or from an individual site.
; INPUTS:
;            FILE --> name of MDN file to be read. If not given, then the default
;                     is to read the file containing all the site data:
;                     '~/ALASKA/data/MDN/mdnic.txt'. Another possibility is 
;                     to read in individual files from  '~/ALASKA/data/MDN/site/'
; OUTPUT:
;            DATA --> structure containing the MDN data
;            DATA_REGRID --> structure containing gridded annual MDN data with resolutions determined
;                            by RES keyword.
;            DATA_REGRID_MONTHLY --> gridded monthly MDN data
; KEYWORDS:
;            SELECT_YEAR --> 2 element array containing start and end year of data to 
;                            be read in. If not given than all the data is read in.
;            Rgres     --> if set, then the program will map the MDN data onto a regular grid of
;                        resolution RES. RES='4x5' or RES='05x0666'. The result will be in data_regrid.
;            MIN_NPTS  --> Minimum number of points to consider for regridding. Only sites with 
;                          model than MIN_NPTS valid points (typically 80% temporal coverage) will 
;                          be used in the regridding. If MIN_NPTS is not passed, then all the data
;                          from all the sites will be used.
;             outlier: only influence the running monthly mean
;            type: GEOS5 or MERRA met fields
;
; Example:
; IDL> read_mdn,file,mdn_data,select_year=[2005,2005],res='4x5',/plot
; Created: Lyatt Jaegle October 6 2009.
; Updated: Lyatt Jaegle January 30, 2011: 
;                   - replaced mean with mean_nan for the monthly and annual mean deposition calculation
;                   - added gridding option
;                   - also calculate seasonal means: Winter (DJF), Spring (MAM), Summer (JJA), Fall (SON)
; ++++++++++++++++++++++++++++++++++++++++++++++
; number of years to read
nyears=max(select_year)- min(select_year)+1L

if metresolution eq '4x5' then res = [5.0, 4.0] else res = [2./3., 1./2.]

; read snow correlation ratio
restore, './precipitation/precipitation_data_2D_7years_' + mettype + '_' + metresolution +  '.sav'   ; 151x121x84x3,   monthly mean through 2004 to 2010, 
precipitation=mod_data

if mettype eq 'MERRAL' then begin
       precipitation_clim=(precipitation[*,*,0:11,*] + precipitation[*,*,12:23,*] + precipitation[*,*,24:35,*] + $
                            precipitation[*,*,36:47,*] + precipitation[*,*,48:59,*] + precipitation[*,*,60:71,*] + $
                            precipitation[*,*,72:83,*] + precipitation[*,*,84:95,*] + precipitation[*,*,96:107,*] + $
                            precipitation[*,*,108:119,*] + precipitation[*,*,120:131,*] + precipitation[*,*,132:143,*] + $
                            precipitation[*,*,144:155,*] + precipitation[*,*,156:167,*] + precipitation[*,*,168:179,*] $
                             ) / 15.0   ; climatological mean
       precipitation_start = 1996
       precipitation_end = 2010
endif else begin
       precipitation_clim=(precipitation[*,*,0:11,*] + precipitation[*,*,12:23,*] + precipitation[*,*,24:35,*] + $
                            precipitation[*,*,36:47,*] + precipitation[*,*,48:59,*] + precipitation[*,*,60:71,*] + $
                            precipitation[*,*,72:83,*] ) / 7.0   ; climatological mean
       precipitation_start = 2004
       precipitation_end = 2010
endelse

  ; first read in site info
  read_mdn_siteinfo,sites

  file='./MDN_data/mdnic.txt'

  ; first find out number of lines in file
  nlines=file_lines(file)

  ; Open and read the filen
  OPEN_FILE,file,ilun
  data_str=strarr(nlines)
  readf,ilun,data_str
  free_lun,ilun
  nvar=13
  array=strarr(nlines,nvar)

  for k=0L,nlines-1L do begin
      str=strsplit(data_str(k),",",count=count,/extract)
      if count gt 0 then array(k,0:count-1)=str
  endfor
  header=strlowcase(reform(array(0,*)))
  data=array(1:*,*)
  nlines=nlines-1

  ; replace missing values (-9.00 with NaN)
  q=where(data eq '-9.00' or data eq '-7.00')
  data(q)=!values.f_nan

  ; test to see if this is a new version of the MDN files or an
  ; older version. The difference is in the way "dateon" and "dateoff"
  ; are specified (with or without " ").
  test=strpos(array(1,1),'"')
  if test(0) eq -1 then begin
    print, 'New version of the file '
    new=1
  endif else begin
    print, 'Old version of the file '
    new=0
  endelse

  ; create structure containing data for the sites
  name=data(*,0)
  
  ; make it obselete, yanxu, 3/31/11
  ; if keyword_set(new) then name='"'+name+'"'
  
  uniq_sites=name(uniq(name,sort(name)))

  nsites=n_elements(uniq_sites)
  tmp=fltarr(1000)+!values.f_nan
  stmp=strarr(1000)
  mtmp=fltarr(12)+!values.f_nan
  nmtmp=fltarr((max(select_year)-min(select_year)+1)*12)+!values.f_nan
  setmp=fltarr(4)+!values.f_nan
  sites_all=replicate({id:'',name:'',lat:0.,lon:0.,$
                       start_day:stmp,start_month:stmp,start_year:stmp,start_hour:tmp,start_jday:tmp,$
                       end_day:stmp,end_month:stmp,end_year:stmp,end_hour:tmp,end_jday:tmp,$
                       day:stmp,month:stmp,year:stmp,hour:tmp,jday:tmp,$
                       hgdep:tmp,hgconc:tmp,rgppt:tmp,svol:tmp,subppt:tmp,$
		               qr:stmp,$
		               daily_dep:tmp,$
                       daily_precipitation:tmp,$
                       monthly_precipitation:mtmp,$
                       monthly_precipitation_gctt:mtmp,$
                       monthly_precipitation_gccv:mtmp,$
		               monthly_dep:mtmp,monthly_dep_std:mtmp,monthly_dep_median:mtmp,monthly_npts:fltarr(12),monthly_dep_anomaly:fltarr(12,nyears)+!values.f_nan,$
		               season_dep:setmp,season_dep_std:setmp,season_npts:fltarr(4),season_dep_anomaly:fltarr(4,nyears)+!values.f_nan,$
		               annual_dep:!values.f_nan,annual_dep_std:!values.f_nan,annual_dep_anomaly:fltarr(nyears)+!values.f_nan ,$
                       running_monthly_dep:nmtmp, running_monthly_prec:nmtmp, $    ; ng/m2/day,  mm/day
                       running_monthly_prec_gctt:nmtmp, running_monthly_prec_gccv:nmtmp, $   ; mm/month
		               npts: intarr(max(select_year)- min(select_year)+1L)  },nsites)
			   
  ; read in date/time
  dateon=data(*,where(header eq '"dateon"' or header eq 'dateon' )) 
  dateoff=data(*,where(header eq '"dateoff"' or header eq 'dateoff' ))
  if not keyword_set(new) then begin
   start_month=strmid(dateon,1,2) & start_day=strmid(dateon,4,2) & start_year=strmid(dateon,7,4)
   start_hour=float(strmid(dateon,12,2))+float(strmid(dateon,15,2))/60.
   end_month=strmid(dateoff,1,2) & end_day=strmid(dateoff,4,2) & end_year=strmid(dateoff,7,4)
   end_hour=float(strmid(dateoff,12,2))+float(strmid(dateoff,15,2))/60.
  endif else begin
   start_month=strmid(dateon,0,2) & start_day=strmid(dateon,3,2) & start_year=strmid(dateon,6,4)
   start_hour=float(strmid(dateon,11,2))+float(strmid(dateon,14,2))/60.
   end_month=strmid(dateoff,0,2) & end_day=strmid(dateoff,3,2) & end_year=strmid(dateoff,6,4)
   end_hour=float(strmid(dateoff,11,2))+float(strmid(dateoff,14,2))/60.
  endelse
  start_jday=julday(start_month,start_day,start_year,0)-julday(1,1,start_year,0)+1+start_hour/24.
  end_jday=julday(end_month,end_day,end_year,0)-julday(1,1,end_year,0)+1+end_hour/24.

  ; jday=(start_jday+end_jday)/2.
  del_day=end_jday-start_jday
  q=where(del_day lt 0)
  del_day(q)=del_day(q)+(julday(12,31,start_year(q))-julday(1,1,start_year(q))+1)
  
  jday=start_jday+del_day/2.  
  caldat,jday+julday(1,1,start_year)-1,month,day,year,hour,minute

  ; read in quality of data ("A" = fully qualified with no problems
  ;                          "B" = valid data with minor problems, used for summary statistics
  ;                          "C" = invalid data, not used for summary statistics
  ; Initially try to use only A and B data
  QR=data(*,where(header eq '"qr"' or header eq 'qr' ))

  for k=0,nsites-1 do begin
    sites_all(k).id=uniq_sites(k)
    ind=where(sites.id eq uniq_sites(k))

    ; skip when site not found
    if ind eq -1 then continue

    sites_all(k).name=sites(ind).name
    sites_all(k).lat=float(sites(ind).lat)
    sites_all(k).lon=float(sites(ind).lon)

    ; Try using data with A or B quality.
    ind0=where(name eq uniq_sites(k),count0)

    if keyword_set(new) then $
    	ind=where(name eq uniq_sites(k) and (qr eq 'A' or qr eq 'B') and $
                  year ge min(select_year) and year le max(select_year)  ,count)  $
                        else $
	    ind=where(name eq uniq_sites(k) and (qr eq '"A"' or qr eq '"B"') and $
                  year ge min(select_year) and year le max(select_year)  ,count)
    ; ind=where(name eq uniq_sites(k) and $
    ;          year ge select_year(0) and year le select_year(1)  ,count)

	for i_select_year=min(select_year), max(select_year) do begin
	  if keyword_set(new) then begin
    	dummy=where(name eq uniq_sites(k) and (qr eq 'A' or qr eq 'B') and $
                  year eq i_select_year  , tmp_count   ) 
		sites_all(k).npts( i_select_year-min(select_year)) = tmp_count
      endif else begin
	    dummy=where(name eq uniq_sites(k) and (qr eq '"A"' or qr eq '"B"') and $
                  year eq i_select_year  , tmp_count      )
		sites_all(k).npts( i_select_year-min(select_year)) = tmp_count
	  endelse
	endfor
	
    ; sites_all(k).npts=count
    print, 'Site ',uniq_sites(k),count0,count
    ; HGDEP = Total mercury deposition, ng/m2 . The product of SUBPPT and HGCONC. 
    if count eq 0 then goto,skip_this
    hg=data(*,where(header eq '"hgdep"'  or header eq 'hgdep' )) & hg(where(hg lt 0.0))=!values.f_nan  & sites_all(k).hgdep(0:count-1)=hg(ind)
    hg=data(*,where(header eq '"hgconc"' or header eq 'hgconc' )) & hg(where(hg lt 0.0))=!values.f_nan & sites_all(k).hgconc(0:count-1)=hg(ind)
    hg=data(*,where(header eq '"svol"' or header eq 'svol'  )) & hg(where(hg lt 0.0))=!values.f_nan & sites_all(k).svol(0:count-1)=hg(ind)
    hg=data(*,where(header eq '"rgppt"' or header eq 'rgppt' ))& hg(where(hg lt 0.0))=!values.f_nan & sites_all(k).rgppt(0:count-1)=hg(ind)   ; precipitation in mm
    hg=data(*,where(header eq '"subppt"'  or header eq 'subppt' ))& hg(where(hg lt 0.0))=!values.f_nan & sites_all(k).subppt(0:count-1)=hg(ind)  ; precipitation in mm

    ; note for the data
    ;  SITE CODE: 2-letter state or province designator followed by a two digit number.
    ;  START DATE: (mm/dd/yyyy hh:mm), GMT
    ;  END DATE: (mm/dd/yyyy hh:mm), GMT
    ;  RGPPT: Precipitation amount as measured by the rain gage in millimeters.
    ;              Trace amounts are indicated by -7.00 and missing amount by -9.00.
    ;  SVOL: Sample Volume, ml. Missing amounts are indicated by a -9.00.
    ;  SUBPPT: Rain gage precipitation amount, if available, in mm. If the rain gage
    ;                value (RGPPT) is missing, the precipitation amount in mm is calculated from the
    ;                net sample volume caught in the sample bottle. A value of 0.127 is inserted for
    ;                Trace sample types. Missing amounts are indicated by a -9.00
    ;  HGCONC: Total mercury concentration reported by the lab in ng/L. Missing
    ;                  amounts are indicated by a -9.00
    ;  HGDEP: Total mercury deposition, ng/m2 . The product of SUBPPT and
    ;               HGCONC. Missing amounts are indicated by a -9.00

    sites_all(k).start_day(0:count-1)=start_day(ind)
    sites_all(k).start_month(0:count-1)=start_month(ind)
    sites_all(k).start_year(0:count-1)=start_year(ind)
    sites_all(k).start_hour(0:count-1)=start_hour(ind)
    sites_all(k).start_jday(0:count-1)=start_jday(ind)
    sites_all(k).end_day(0:count-1)=end_day(ind)
    sites_all(k).end_month(0:count-1)=end_month(ind)
    sites_all(k).end_year(0:count-1)=end_year(ind)
    sites_all(k).end_hour(0:count-1)=end_hour(ind)
    sites_all(k).end_jday(0:count-1)=end_jday(ind)
    sites_all(k).day(0:count-1)=end_day(ind)
    sites_all(k).month(0:count-1)=month(ind)
    sites_all(k).year(0:count-1)=year(ind)
    sites_all(k).hour(0:count-1)=hour(ind)
    sites_all(k).jday(0:count-1)=jday(ind)
    sites_all(k).qr(0:count-1)=qr(ind)
	
    ; calculate daily deposition in ng/m2/day as HGDEP/(jday_end-jday_start)
    ; sites_all(k).daily_dep(0:count-1)=sites_all(k).hgdep(0:count-1)/(end_jday(ind)-start_jday(ind))
    ctm_index,ctm_type('GEOS5_47L',res=res),i1,j1,$
        center=[sites_all(k).lat,sites_all(k).lon],/non_interactive
    
    if metresolution eq '05x0666' then begin
     i0 = 60
     j0 = 200
    endif else begin
;     i0 = 8
;     j0 = 25
     i0 = 0
     j0 = 0 
    endelse

if snow eq 1 and sites_all(k).lat ge 10. and sites_all(k).lat le 70. and sites_all(k).lon ge -140. and sites_all(k).lon le -40.  then begin
;    accu_snowfrac=[]
    for i_point=0,count-1 do begin
          if year(ind(i_point)) ge precipitation_start and year(ind(i_point)) le precipitation_end then begin
               i_month = 12* (year(ind(i_point)) - precipitation_start) + month(ind(i_point)) - 1
               snow_fraction = precipitation( -i0+i1-1, -j0+j1-1,i_month,2  )/ precipitation( -i0+i1-1, -j0+j1-1,i_month,0  )
          endif else begin
               i_month = month(ind(i_point)) - 1 
               snow_fraction = precipitation_clim( -i0+i1-1, -j0+j1-1,i_month,2  )/ precipitation_clim( -i0+i1-1, -j0+j1-1,i_month,0  )
          endelse

          ; assume only capture half of the snow by MDN sites
 ;         accu_snowfrac=[accu_snowfrac,snow_fraction]
;          snow_correction_factor = snow_fraction * 2.0 + 1.0 - snow_fraction
          snow_correction_factor = snow_fraction * 1.13 + 1.0 - snow_fraction
	      sites_all(k).daily_dep(i_point)=sites_all(k).hgdep(i_point)/del_day(ind(i_point))*snow_correction_factor
	      sites_all(k).rgppt(i_point)=sites_all(k).rgppt(i_point)*snow_correction_factor
	      sites_all(k).subppt(i_point)=sites_all(k).subppt(i_point)*snow_correction_factor		  
    endfor
;    print, mean(accu_snowfrac)
endif 

if snow eq 2 and sites_all(k).lat ge 10. and sites_all(k).lat le 70. and sites_all(k).lon ge -140. and sites_all(k).lon le -40.  then begin
    year_last=0L & month_last=0L   
    year_current=0L & month_current=0L   
    
    for i_point=0,count-1 do begin
        year_last=year_current & month_last=month_current
        year_current=year(ind(i_point)) & month_current=month(ind(i_point))

          if year(ind(i_point)) ge precipitation_start and year(ind(i_point)) le precipitation_end then begin
               i_month = 12* (year(ind(i_point)) - precipitation_start) + month(ind(i_point)) - 1
               precipitation_gc=precipitation( -i0+i1-1, -j0+j1-1,i_month,0  )
               snow_fraction = precipitation( -i0+i1-1, -j0+j1-1,i_month,2  )/ precipitation( -i0+i1-1, -j0+j1-1,i_month,0  )
          endif else begin
               i_month = month(ind(i_point)) - 1 
               precipitation_gc=precipitation_clim( -i0+i1-1, -j0+j1-1,i_month,0  )
               snow_fraction = precipitation_clim( -i0+i1-1, -j0+j1-1,i_month,2  )/ precipitation_clim( -i0+i1-1, -j0+j1-1,i_month,0  )
          endelse

          if year_last ne year_current or month_current ne month_last then begin
               qind = where(sites_all(k).month eq month(ind(i_point)) and sites_all(k).year eq year(ind(i_point)), points)
               ;precipitation_mdn=total(sites_all(k).rgppt(qind))
			   precipitation_mdn=total(sites_all(k).subppt(qind))
          endif
 
          ; assume only capture half of the snow by MDN sites
          ; snow_correction_factor = snow_fraction * 2.0 + 1.0 - snow_fraction
          snow_correction_factor=snow_fraction/(1.-(1-precipitation_mdn/precipitation_gc)/snow_fraction)+1.0-snow_fraction
          if snow_correction_factor lt 0. or snow_correction_factor gt 5. then snow_correction_factor=precipitation_gc/precipitation_mdn
	      sites_all(k).daily_dep(i_point)=sites_all(k).hgdep(i_point)/del_day(ind(i_point))*snow_correction_factor
	      sites_all(k).rgppt(i_point)=sites_all(k).rgppt(i_point)*snow_correction_factor
		  sites_all(k).subppt(i_point)=sites_all(k).subppt(i_point)*snow_correction_factor
    endfor
endif 

if snow eq 0 then begin
	sites_all(k).daily_dep(0:count-1)=sites_all(k).hgdep(0:count-1)/del_day(ind)
endif
;sites_all(k).daily_precipitation(0:count-1)=sites_all(k).rgppt(0:count-1)/del_day(ind)
sites_all(k).daily_precipitation(0:count-1)=sites_all(k).subppt(0:count-1)/del_day(ind)

for i_month=1,12 do begin     
        month_index = 12*indgen(nyears) + 12* (min(select_year) - precipitation_start) + i_month - 1
        if snow gt 0 then begin
            sites_all(k).monthly_precipitation_gctt(i_month-1)=mean(precipitation( -i0+i1-1, -j0+j1-1, month_index  ,0  ))
            sites_all(k).monthly_precipitation_gccv(i_month-1)=mean(precipitation( -i0+i1-1, -j0+j1-1, month_index  ,1  ))
        endif
endfor

sites_all(k).running_monthly_prec_gctt=precipitation(-i0+i1-1, -j0+j1-1,  indgen(12*nyears) + 12* (min(select_year) - precipitation_start)   ,0  )
sites_all(k).running_monthly_prec_gccv=precipitation(-i0+i1-1, -j0+j1-1, indgen(12*nyears) + 12* (min(select_year) - precipitation_start)   ,1  )
	 
    ; calculate monthly mean deposition in ng/m2/day
    monthly_dep=fltarr(12,nyears)+!values.f_nan
    monthly_prec=fltarr(12,nyears)+!values.f_nan
    season_dep=fltarr(4,nyears)+!values.f_nan
    annual_dep=fltarr(nyears)+!values.f_nan

    for i=0,12-1 do begin
        qind = where(sites_all(k).month eq i+1 and sites_all(k).daily_dep ge 0.0, points)
	; modified to have at least 2 points for a monthly mean (jaegle 1/5/2010)
	if points gt 1 then begin
	        ; replace mean with mean_nan --> if only NaN values are available, then
		; the result of mean_nan will be NaN instead of 0 (jaegle, 1/30/11)
        ; correct it by the number of points, 4.35 is the average number of weeks per month, yanxu
        if keyword_set(correct) then time_correct_factor=min([points / (4.35*nyears),1.0])  $
                                            else  time_correct_factor=1.0   ; not correct
		sites_all(k).monthly_dep(i)=mean_nan(sites_all(k).daily_dep(qind),/nan) * time_correct_factor
		sites_all(k).monthly_precipitation(i)=mean_nan(sites_all(k).daily_precipitation(qind),/nan) * time_correct_factor
		sites_all(k).monthly_dep_median(i)=median(sites_all(k).daily_dep(qind))
		sites_all(k).monthly_dep_std(i)=stddev(sites_all(k).daily_dep(qind),/nan)
		sites_all(k).monthly_npts(i)=points

        for iyear=0,nyears-1 do begin           
              qind_iyear=where(sites_all(k).month eq i+1 and sites_all(k).year eq iyear+select_year[0] and sites_all(k).daily_dep ge 0.0, points_iyear)
              if points_iyear ge 1 then begin
                 if keyword_set(correct) then time_correct_factor=min([points_iyear / 4.35, 1.0])  $ ; also need to correct
                                            else  time_correct_factor=1.0   ; not correct
                       monthly_dep(i,iyear)=mean_nan(sites_all(k).daily_dep(qind_iyear),/nan) * time_correct_factor                   ; ng/m2/day
                       monthly_prec(i,iyear)=mean_nan(sites_all(k).daily_precipitation(qind_iyear),/nan) * time_correct_factor  ; mm/day
              endif
        endfor   ; iyear
        
        ; correct by the multiple year mean
        if (mean(monthly_dep(i,*), /nan) gt 0.0) then time_correct_factor = sites_all(k).monthly_dep(i) / mean(monthly_dep(i,*), /nan) else time_correct_factor=0.0
        monthly_dep(i,*) = monthly_dep(i, *) * time_correct_factor
        if (mean(monthly_prec(i,*), /nan) gt 0.0) then  time_correct_factor = sites_all(k).monthly_precipitation(i) / mean(monthly_prec(i,*), /nan)  else time_correct_factor=0.0
        monthly_prec(i,*) = monthly_prec(i, *) * time_correct_factor
 
	endif

    endfor ; end loop over months

    sites_all(k).running_monthly_dep=reform(monthly_dep,12*nyears)       ; ng/m2/day
    sites_all(k).running_monthly_prec=reform(monthly_prec,12*nyears)     ; mm/day

    ; calculate seasonal means in ng/m2/day
    season=fltarr(4,3)
    season(0,*)=[12,1,2]-1
    season(1,*)=[3,4,5]-1
    season(2,*)=[6,7,8]-1
    season(3,*)=[9,10,11]-1
    season_string=['DJF','MAM','JJA','SON']

    ; calculate the anomaly
    for iyear=0,nyears-1 do begin
          annual_dep[iyear]=mean(monthly_dep[*,iyear])
          for iseason=0,3 do begin
              season_dep[iseason,iyear]=mean(monthly_dep[season[iseason,*],iyear])
          endfor
    endfor

    clim_annual_dep=mean(annual_dep)
    if  nyears gt 1 then clim_season_dep=mean(season_dep,2) else clim_season_dep=season_dep
    if  nyears gt 1 then clim_month_dep=mean(monthly_dep,2) else clim_month_dep=monthly_dep

    for iyear=0,nyears-1 do begin
          sites_all(k).annual_dep_anomaly[iyear]=(annual_dep[iyear]-clim_annual_dep)/clim_annual_dep
          for iseason=0,3 do begin
              sites_all(k).season_dep_anomaly[iseason,iyear]=(season_dep[iseason,iyear]-clim_season_dep[iseason])/clim_season_dep[iseason]
          endfor
          for imonth=0,11 do begin
              sites_all(k).monthly_dep_anomaly[imonth,iyear]=(monthly_dep[imonth,iyear]-clim_month_dep[imonth])/clim_month_dep[imonth]
          endfor
    endfor

    ; clim mean
    sites_all(k).season_dep= mean_nan( sites_all(k).monthly_dep(season), 2, /nan)
    sites_all(k).season_npts= total( sites_all(k).monthly_npts(season), 2, /nan)
    for i=0,3 do begin
            sites_all(k).season_dep_std(i)= stddev( sites_all(k).monthly_dep(season(i,*)), /nan)
            for iyear=0,nyears-1 do begin
                     sites_all(k).season_dep_anomaly(i,iyear)=mean_nan( sites_all(k).monthly_dep_anomaly(season(i,*),iyear), /nan)
            endfor
    endfor

    ; calculate annual mean deposition in ng/m2/day
    ; as above, replace mean with mean_nan  (jaegle, 1/30/11)
    sites_all(k).annual_dep=mean_nan(sites_all(k).monthly_dep,/nan)
    sites_all(k).annual_dep_std=stddev(sites_all(k).monthly_dep,/nan)
    for iyear=0,nyears-1 do begin
             sites_all(k).annual_dep_anomaly(iyear)=mean_nan( sites_all(k).monthly_dep_anomaly(*,iyear), /nan)
    endfor
    skip_this:
  endfor ; end loop over sites

  data=sites_all

end
