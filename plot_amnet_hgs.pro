pro plot_amnet_hgs, file, resolution, year,hg2frac=hg2frac,title=title

; plot gom from amnet and camnet observations, comparison with model
; spatial distribution, scatter plot, monthly mean

lat0=[9.75,70.25]
lon0=[-140.-1./3.,-40.+1./3.]

limit=[24,-130,55,-60]

if strcmp(resolution,'4x5') then begin
  hlat = 2.0
  hlon = 2.5
endif else begin
  hlat = 0.25
  hlon = 1.0/3.0
endelse

; read model output
conc_convert=200.6/(8.31*273.2)*1D2 ; convert from pptv to ng/m3 STP
year_num=1.

; file='~/4-geos/0-hg/Run.'+resolution+'/dataout/ctm.bpch.'+version

for iyear=min(year),max(year) do begin

if keyword_set(hg2frac) then $
read_run,file=file,conc_all,emis_all,drydep_all,$
            wetdep_all,longitude,latitude,convert,data_model,prof_all,$
	        gridinfo,lat0=lat0,lon0=lon0,year=iyear,/hg2frac $
else $
read_run,file=file,conc_all,emis_all,drydep_all,$
            wetdep_all,longitude,latitude,convert,data_model,prof_all,$
	        gridinfo,lat0=lat0,lon0=lon0,year=iyear


if iyear eq min(year) then begin 
 model_gem=conc_all.hg0*conc_convert
 model_newrgm=conc_all.newrgm*conc_convert
 model_newrpm=conc_all.newrpm*conc_convert
endif else begin
 year_num=year_num+1.
 model_gem=model_gem+conc_all.hg0*conc_convert
 model_newrgm=model_newrgm+conc_all.newrgm*conc_convert
 model_newrpm=model_newrpm+conc_all.newrpm*conc_convert
endelse

endfor

 model_gem=reform(model_gem[*,*,0,*])/year_num
 model_newrgm=reform(model_newrgm[*,*,0,*])/year_num
 model_newrpm=reform(model_newrpm[*,*,0,*])/year_num

; using pg/m3 as the unit
 model_newrgm=model_newrgm*1000.
 model_newrpm=model_newrpm*1000.

; read amnet observation
read_amnet_siteinfo,siteinfo

nsites=n_elements(siteinfo)
ndata=fltarr(12)+nsites

obs_monthly_gem=fltarr(12)
obs_annual_gem=fltarr(nsites)  & obs_annual_gem(*)=!values.f_nan
mod_monthly_gem=fltarr(12)
mod_annual_gem=fltarr(nsites) 

obs_monthly_rgm=fltarr(12)
obs_annual_rgm=fltarr(nsites) & obs_annual_rgm(*)=!values.f_nan
mod_monthly_newrgm=fltarr(12)
mod_annual_newrgm=fltarr(nsites)

obs_monthly_rpm=fltarr(12)
obs_annual_rpm=fltarr(nsites)  & obs_annual_rpm(*)=!values.f_nan
mod_monthly_newrpm=fltarr(12)
mod_annual_newrpm=fltarr(nsites)

site_lat=fltarr(nsites)
site_lon=fltarr(nsites)

; handle the amnet data

for k=0,nsites-1 do begin
  print, 'now reading... ',siteinfo(k).siteid
  read_amnet,siteinfo(k).siteid,data

  site_lat(k)=siteinfo(k).lat
  site_lon(k)=siteinfo(k).lon
  obs_annual_gem(k)=mean(data.gem_monthly,/nan)
  obs_annual_rgm(k)=mean(data.gom_monthly,/nan)
  obs_annual_rpm(k)=mean(data.pbm_monthly,/nan)
  
  for imonth=0,11 do begin
    if data.gom_monthly[imonth] ne data.gom_monthly[imonth] then begin
      ndata[imonth]=ndata[imonth]-1
      continue
    endif
    obs_monthly_gem[imonth]=obs_monthly_gem[imonth]+data.gem_monthly[imonth]
    obs_monthly_rgm[imonth]=obs_monthly_rgm[imonth]+data.gom_monthly[imonth]
    obs_monthly_rpm[imonth]=obs_monthly_rpm[imonth]+data.pbm_monthly[imonth]
  endfor

  i_ind=where(abs(longitude-site_lon[k]) le hlon )
  j_ind=where(abs(latitude-site_lat[k]) le hlat )
  if (i_ind lt 0 or j_ind lt 0) then begin
    ndata=ndata-1
    continue
  endif

  mod_annual_gem(k)=total(model_gem(i_ind,j_ind,*))/12.
  mod_annual_newrgm(k)=total(model_newrgm(i_ind,j_ind,*))/12.
  mod_annual_newrpm(k)=total(model_newrpm(i_ind,j_ind,*))/12.

  for imonth=0,11 do begin
    if data.gom_monthly[imonth] ne data.gom_monthly[imonth] then continue
    mod_monthly_gem[imonth]=mod_monthly_gem[imonth]+model_gem[i_ind,j_ind,imonth]
    mod_monthly_newrgm[imonth]=mod_monthly_newrgm[imonth]+model_newrgm[i_ind,j_ind,imonth]
    mod_monthly_newrpm[imonth]=mod_monthly_newrpm[imonth]+model_newrpm[i_ind,j_ind,imonth]
  endfor
endfor

; average the data

mod_monthly_gem=mod_monthly_gem/ndata
obs_monthly_gem=obs_monthly_gem/ndata
model_gem_annual = total(model_gem,3)/12.

mod_monthly_newrgm=mod_monthly_newrgm/ndata
obs_monthly_rgm=obs_monthly_rgm/ndata
model_newrgm_annual = total(model_newrgm,3)/12.

mod_monthly_newrpm=mod_monthly_newrpm/ndata
obs_monthly_rpm=obs_monthly_rpm/ndata
model_newrpm_annual = total(model_newrpm,3)/12.

; ploting
;open_device,/ps,/color,/portrait,file='AMNet_HGS_'+version+'_'+string(min(year),format='(I4)')+'_'+string(max(year),format='(I4)')+'_'+resolution+'.ps'
;!P.font=0
;Device, /Helvetica, /Isolatin1
 
multipanel,col=1,row=2, omargin=0.0
ctm_overlay, model_gem_annual, longitude, latitude,$
	   /sample,/continents,/cbar,DIVISIONS=5, $
	   mindata=1.0,maxdata=2.0,$
	   /usa,limit=limit,  $
       cbunit='GEM, ng m!u-3!n',$   
       obs_annual_gem,site_lon,site_lat,$
       t_symbol=1,symsize=2.6,/iso,csfac=1.0,title=title

ctm_overlay, model_newrgm_annual, longitude, latitude,$
	   /sample,/continents,/cbar,DIVISIONS=5, $
	   mindata=0.,maxdata=40,$
	   /usa,limit=limit,  $
       cbunit='RGM, pg m!u-3!n',$   
       obs_annual_rgm,site_lon,site_lat,$
       t_symbol=1,symsize=2.6,/iso,csfac=1.0

multipanel,col=1,row=2, omargin=0.0
ctm_overlay, model_newrpm_annual, longitude, latitude,$
	   /sample,/continents,/cbar,DIVISIONS=5, $
	   mindata=0,maxdata=20,$
	   /usa,limit=limit,  $
       cbunit='PBM, pg m!u-3!n',$   
       obs_annual_rpm,site_lon,site_lat,$
       t_symbol=1,symsize=2.6,/iso,csfac=1.0

multipanel,col=2,row=2, pos=pos, margin=[0.05,0.05]
; hg0
 plot, obs_annual_gem, mod_annual_gem,$
  	ytitle='GEOS-Chem Model Hg0, ngm!u-3!n',yrange=[1.3,2.4],$
  	xrange=[1.3,2.4],xtitle='AMNET Hg0, ngm!u-3!n',$
     color=!myct.black,psym=sym(1), pos=pos,/iso

 rcorr=correlate(obs_annual_gem,mod_annual_gem)
 org_corr,obs_annual_gem,mod_annual_gem,rcorr,$
  	n_elements(obs_annual_gem),slope,intercept
 oplot, !x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
 oplot, !x.crange,!x.crange,col=1
 xyouts,!x.crange(0)+(!x.crange(1)-!x.crange(0))/20.,$
        !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.92,$
          /data,col=1, 'Model = ' + string(intercept,slope,format= '(F6.3," +  AMNET x", F5.2)' ) + $
          '!c(r= ' + $
          string(rcorr,format='(F5.2)')+ '; npts= '+$
          string(n_elements(obs_annual_gem),format= '(I6)' )+ ')' +$
 	      '!c Model/Obs=' + string(mean(mod_annual_gem/obs_annual_gem),format='(F5.2)') + '!c' + $
          string(mean(mod_annual_gem),format='(F4.2)') + '+/-' + $
          string(stddev(mod_annual_gem),format='(F4.2)') + 'ng/m3!c' + $
          string(mean(obs_annual_gem),format='(F4.2)') + '+/-' + $
          string(stddev(obs_annual_gem),format='(F4.2)') + 'ng/m3!c', $
          charsize=0.9
multipanel, /advance, pos=pos

; newrgm
 plot, obs_annual_rgm, mod_annual_newrgm,$
  	ytitle='GEOS-Chem Model RGM, pgm!u-3!n',yrange=[0,40],$
  	xrange=[0,40],xtitle='AMNET RGM, pgm!u-3!n',$
     color=!myct.black,psym=sym(1), pos=pos,/iso

 rcorr=correlate(obs_annual_rgm,mod_annual_newrgm)
 org_corr,obs_annual_rgm,mod_annual_newrgm,rcorr,$
  	n_elements(obs_annual_rgm),slope,intercept
 oplot, !x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
 oplot, !x.crange,!x.crange,col=1
 xyouts,!x.crange(0)+(!x.crange(1)-!x.crange(0))/20.,$
        !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.92,$
          /data,col=1, 'Model = ' + string(intercept,slope,format= '(F6.3," +  AMNET x", F5.2)' ) + $
          '!c(r= ' + $
          string(rcorr,format='(F5.2)')+ '; npts= '+$
          string(n_elements(obs_annual_rgm),format= '(I6)' )+ ')' +$
 	      '!c Model/Obs=' + string(mean(mod_annual_newrgm/obs_annual_rgm),format='(F5.2)') + '!c' + $
          string(mean(mod_annual_newrgm),format='(F5.2)') + '+/-' + $
          string(stddev(mod_annual_newrgm),format='(F5.2)') + 'pg/m3!c' + $
          string(mean(obs_annual_rgm),format='(F5.2)') + '+/-' + $
          string(stddev(obs_annual_rgm),format='(F5.2)') + 'pg/m3!c', $
          charsize=0.9
multipanel, /advance, pos=pos

; newrpm
 plot, obs_annual_rpm, mod_annual_newrpm,$
  	ytitle='GEOS-Chem Model RPM, pgm!u-3!n',yrange=[0,20],$
  	xrange=[0,20],xtitle='AMNET RPM, pgm!u-3!n',$
     color=!myct.black,psym=sym(1), pos=pos,/iso

 rcorr=correlate(obs_annual_rpm,mod_annual_newrpm)
 org_corr,obs_annual_rpm,mod_annual_newrpm,rcorr,$
  	n_elements(obs_annual_rpm),slope,intercept
 oplot, !x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
 oplot, !x.crange,!x.crange,col=1
 xyouts,!x.crange(0)+(!x.crange(1)-!x.crange(0))/20.,$
        !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.92,$
          /data,col=1, 'Model = ' + string(intercept,slope,format= '(F6.3," +  AMNET x", F5.2)' ) + $
          '!c(r= ' + $
          string(rcorr,format='(F5.2)')+ '; npts= '+$
          string(n_elements(obs_annual_rpm),format= '(I6)' )+ ')' +$
 	      '!c Model/Obs=' + string(mean(mod_annual_newrpm/obs_annual_rpm),format='(F5.2)') + '!c' + $
          string(mean(mod_annual_newrpm),format='(F5.2)') + '+/-' + $
          string(stddev(mod_annual_newrpm),format='(F5.2)') + 'pg/m3!c' + $
          string(mean(obs_annual_rpm),format='(F5.2)') + '+/-' + $
          string(stddev(obs_annual_rpm),format='(F5.2)') + 'pg/m3!c', $
          charsize=0.9
;multipanel, /advance, pos=pos

;close_device
end