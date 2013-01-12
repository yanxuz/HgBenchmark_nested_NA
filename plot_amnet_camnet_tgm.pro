pro plot_amnet_camnet_tgm,file,resolution,year, regrid=regrid,hg2frac=hg2frac,title=title

; plot tgm from amnet and camnet observations, comparison with model
; spatial distribution, scatter plot, monthly mean

; lat0=[26,58]
; lon0=[-125,-60]

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
 model_gom=conc_all.newrgm*conc_convert
 model_pbm=conc_all.newrpm*conc_convert
endif else begin
 year_num=year_num+1.
 model_gem=model_gem+conc_all.hg0*conc_convert
 model_gom=model_gom+conc_all.newrgm*conc_convert
 model_pbm=model_pbm+conc_all.newrpm*conc_convert
endelse

endfor

 model_gem=reform(model_gem[*,*,0,*])/year_num
 model_gom=reform(model_gom[*,*,0,*])/year_num
 model_pbm=reform(model_pbm[*,*,0,*])/year_num
 model_tgm=model_gem+model_gom

; read amnet observation
read_amnet_siteinfo,siteinfo
; read camnet data
; read_camnet_allsite_data,version,resolution,tmpsite,year=year,data_camnet
restore,'./CAMNet_data/data_camnet.sav'
nsites=n_elements(siteinfo)+n_elements(data_camnet)

ndata=fltarr(12)+nsites

obs_monthly_tgm=fltarr(12)
mod_monthly_tgm=fltarr(12)
obs_annual_tgm=fltarr(nsites)
mod_annual_tgm=fltarr(nsites)
site_name=strarr(nsites)
site_lat=fltarr(nsites)
site_lon=fltarr(nsites)

; handle the amnet data
for k=0,n_elements(siteinfo)-1 do begin
  print, 'now reading... ' , siteinfo(k).siteid
  read_amnet,siteinfo(k).siteid,data
  
  site_lat(k)=siteinfo(k).lat
  site_lon(k)=siteinfo(k).lon
  obs_annual_tgm(k)=mean(data.gem_monthly,/nan)+mean(data.gom_monthly,/nan)/1000.0
  site_name(k)=siteinfo(k).siteid
  for imonth=0,11 do begin
    if data.gem_monthly[imonth]+data.gom_monthly[imonth]/1000.0 ne data.gem_monthly[imonth]+data.gom_monthly[imonth]/1000.0 then begin
      ndata[imonth]=ndata[imonth]-1   
      continue
    endif
    obs_monthly_tgm[imonth]=obs_monthly_tgm[imonth]+data.gem_monthly[imonth]+data.gom_monthly[imonth]/1000.0
  endfor

  i_ind=where(abs(longitude-site_lon[k]) le hlon )
  j_ind=where(abs(latitude-site_lat[k]) le hlat )
  if (i_ind lt 0 or j_ind lt 0) then begin
    ndata=ndata-1
    continue
  endif
  mod_annual_tgm(k)=total(model_tgm(i_ind,j_ind,*))/12.
  for imonth=0,11 do begin
    if data.gem_monthly[imonth]+data.gom_monthly[imonth]/1000.0 ne data.gem_monthly[imonth]+data.gom_monthly[imonth]/1000.0 then continue
    mod_monthly_tgm[imonth]=mod_monthly_tgm[imonth]+model_tgm[i_ind,j_ind,imonth]
  endfor
endfor

; handle the camnet data
for k=n_elements(siteinfo),nsites-1 do begin

  site_lat(k)=data_camnet[k-n_elements(siteinfo)].lat
  site_lon(k)=data_camnet[k-n_elements(siteinfo)].lon
  obs_annual_tgm(k)=mean(data_camnet[k-n_elements(siteinfo)].tgm,/nan)
  site_name(k)=data_camnet[k-n_elements(siteinfo)].sitename
  for imonth=0,11 do begin
    obs_monthly_tgm[imonth]=obs_monthly_tgm[imonth]+data_camnet[k-n_elements(siteinfo)].tgm[imonth]
  endfor

  i_ind=where(abs(longitude-site_lon[k]) le hlon )
  j_ind=where(abs(latitude-site_lat[k]) le hlat )
  if (i_ind lt 0 or j_ind lt 0) then begin
    ndata=ndata-1
    continue
  endif
  mod_annual_tgm(k)=total(model_tgm(i_ind,j_ind,*))/12.
  for imonth=0,11 do begin
    mod_monthly_tgm[imonth]=mod_monthly_tgm[imonth]+model_tgm[i_ind,j_ind,imonth]
  endfor
endfor

; average the data
mod_monthly_tgm=mod_monthly_tgm/ndata
obs_monthly_tgm=obs_monthly_tgm/ndata
; print,ndata

; background_level=3.0
background_distance=100.0  ; km
background_emission=2000.0
; background_emission=200000.0

obs_annual_tgm_amnet=obs_annual_tgm[0:n_elements(siteinfo)-1]
site_lon_amnet=site_lon[0:n_elements(siteinfo)-1]
site_lat_amnet=site_lat[0:n_elements(siteinfo)-1]
; amnet_background=where( obs_annual_tgm_amnet lt background_level )
amnet_background=[]
for ipoint = 0, n_elements(site_lon_amnet)-1 do begin
     if  vicinity_emission(site_lat_amnet[ipoint],site_lon_amnet[ipoint],background_distance)  lt background_emission  and $
          obs_annual_tgm_amnet(ipoint) lt 2.0 $
               then amnet_background=[amnet_background, ipoint]
     print,vicinity_emission(site_lat_amnet[ipoint],site_lon_amnet[ipoint],background_distance)
endfor

obs_annual_tgm_camnet=obs_annual_tgm[n_elements(siteinfo):nsites-1]
site_lon_camnet=site_lon[n_elements(siteinfo):nsites-1]
site_lat_camnet=site_lat[n_elements(siteinfo):nsites-1]
; camnet_background=where( obs_annual_tgm_camnet lt background_level )
camnet_background=[]
for ipoint = 0, n_elements(site_lon_camnet)-1 do begin
     if  vicinity_emission(site_lat_camnet[ipoint],site_lon_camnet[ipoint],background_distance)  lt background_emission  then camnet_background=[camnet_background, ipoint]
     print,vicinity_emission(site_lat_camnet[ipoint],site_lon_camnet[ipoint],background_distance)
endfor

model_tgm_annual=total(  model_tgm,3  )/12.

; ploting
;open_device,/ps,/color,/portrait,file='psfile.ps'
;!P.font=0
;Device, /Helvetica, /Isolatin1
conc_range=[1.0,2.0]

 multipanel,col=1,row=2, omargin=0.0
ctm_overlay,model_tgm_annual,longitude,latitude,$
	   /sample,/continents,/cbar,DIVISIONS=6, $
	   mindata=conc_range[0],maxdata=conc_range[1], $
	   /usa, /iso, csfac=1.0,  $
       limit=limit, $
       ; position=[0.1,0.5,0.9,0.8],$
       cbunit='TGM, ng m!u-3!n',$   
       obs_annual_tgm_amnet[amnet_background],site_lon_amnet[amnet_background],site_lat_amnet[amnet_background],$
;       obs_annual_tgm,site_lon,site_lat,$
       t_symbol=1,symsize=2.0,title=title
;       title='Comparison of TGM with AMNET-CAMNET sites!c'+version+'_'+resolution+' during '+$
;       string(year[0],format='(I4)')+'-'+string(year[1],format='(I4)')

ctm_overlay,model_tgm_annual,longitude,latitude,$
	   /sample,/continents, $
	   mindata=conc_range[0],maxdata=conc_range[1], $
	   /usa,limit=limit,$
       position=[0.1,0.5,0.9,0.8],$
       obs_annual_tgm_camnet[camnet_background],site_lon_camnet[camnet_background],site_lat_camnet[camnet_background],$
       t_symbol=4,symsize=2.0,/overplot,/conus

legend,label=['AMNET', 'CAMNET'],lcol=[1,1],line=[-1,-1],symbol=[1,4],boxcolor=-1,smsize=1.2,$
                           halign=0.05,valign=1.3,charsize=1.
					   
; determine the range of plot
; ind_background=where(obs_annual_tgm lt background_level)
 ind_background=[amnet_background,camnet_background+n_elements(siteinfo)]
 obs_annual_tgm=obs_annual_tgm(ind_background)
 mod_annual_tgm=mod_annual_tgm(ind_background)


 xmin=min(obs_annual_tgm)
 xmax=max(obs_annual_tgm)
 ymin=min(mod_annual_tgm)
 ymax=max(mod_annual_tgm)
 zmin=min([xmin,ymin])
 zmax=max([xmax,ymax])
; 
 zmin=zmin-0.1*(zmax-zmin)
 zmax=zmax+0.09*(zmax-zmin)
; 
;
 plot,obs_annual_tgm,mod_annual_tgm,$
  	ytitle='GEOS-Chem Model TGM, ngm!u-3!n',yrange=[1.2,2.0],$
  	xrange=[1.2,2.0],xtitle='AMNET-CAMNET TGM, ngm!u-3!n',$
     color=!myct.black,psym=sym(1),$
	 pos=[0.15,0.1,0.75,0.5],$
;     pos=pos, $
	 /noerase
; 
; 
 rcorr=correlate(obs_annual_tgm,mod_annual_tgm)
 org_corr,obs_annual_tgm,mod_annual_tgm,rcorr,$
  	n_elements(obs_annual_tgm),slope,intercept
 oplot, !x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
 oplot, !x.crange,!x.crange,col=1
 
 xyouts, !x.crange(0)+(!x.crange(1)-!x.crange(0))/20.,$
        !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.92,$
          /data, col=1,  'Model = ' + string(intercept,slope,format= '(F6.3,"+ (C)AMNET x",F5.2)'   ) + $
          '!c(r= ' + string(rcorr,format= '(F5.2)') + '; npts= '+$
          string(n_elements(obs_annual_tgm),format= '(I6)')+ ')'  + $
 	     '!c Model/Obs = ' + string(mean(mod_annual_tgm)/mean(obs_annual_tgm), format= '(F5.2)')   + $ 
		 '!c Mean Bias (%) = ' + string(  mean( (mod_annual_tgm-obs_annual_tgm)/obs_annual_tgm )*100.  , format='(F5.2)')    +    $
		 '!c Mean Abs. Bias (%) = ' +  string(   mean(  abs(mod_annual_tgm-obs_annual_tgm) /obs_annual_tgm) *100. , format='(F5.2)' ) + $  
		 '!c RMS (ng/m!u3!n) = ' + string(   sqrt(mean((mod_annual_tgm-obs_annual_tgm)^2))   , format='(F7.2)' ) + $ 
         '!c Mean Bias (ng/m!u3!n) = ' + string(   mean( mod_annual_tgm-obs_annual_tgm)   , format='(F7.2)' )  $
		 , charsize=1.0
		 
 xyouts, !x.crange(0)+(!x.crange(1)-!x.crange(0))*0.5, $
        !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.05, $
		/data,col=1, 'GC : '+ string(mean(mod_annual_tgm), stddev(mod_annual_tgm), format= '(F5.2,"+/-",F5.2)') + ' ng/m!u3!n' + $
        '!cObs: ' + string(mean(obs_annual_tgm), stddev(obs_annual_tgm), format= '(F5.2,"+/-",F5.2)') + ' ng/m!u3!n' $
		, charsize=1.0

;  plot,indgen(12)+1,obs_monthly_tgm,color=!myct.black, $
;      ytitle='Monthly Mean!cTGM concentration, ngm!u-3!n',$
;      yrange=[min([mod_monthly_tgm,obs_monthly_tgm]),max([mod_monthly_tgm,obs_monthly_tgm])],$
;      xtitle='Month', $
; ;	 pos=[0.6,0.02,0.95,0.27],$
;      pos=pos, $
; 	 thick=4,/noerase
; 
;  oplot,indgen(12)+1,mod_monthly_tgm,color=!myct.red,thick=4
;  legend,label=['Obs','Model'],lcol=[!myct.black,!myct.red],line=[0,0],symbol=[-1,-1],boxcolor=-1,thick=4,$
;  	halign=0.01,valign=0.95,charsize=1.
;  xyouts,indgen(12)+1,!y.crange(0)+(!y.crange(1)-!y.crange(0))*0.05,string(ndata,format='(I2)'),col=1,al=0.5,charsize=0.7

 ;  plot map of annual TGM regridded onto the GEOS-Chem grid
 if keyword_set(regrid) then begin
         if  resolution eq '4x5' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=4 )
         if  resolution eq '05x0666' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=[2./3.,1./2.] )
         GridInfo   = CTM_GRID( ModelInfo )

    ; regrid annual mean wet deposition
    data_regrid=compute_grid_mdn(  [site_lon_amnet[amnet_background], site_lon_camnet[camnet_background] ], $
                                                    [site_lat_amnet[amnet_background],  site_lat_camnet[camnet_background] ], $
                                                    [obs_annual_tgm_amnet[amnet_background],obs_annual_tgm_camnet[camnet_background] ], $
                                                     res=resolution,             gridinfo=gridinfo)

    multipanel,col=2,row=3,mar=0.01
    cbpos=[0.1,0.05,0.9,0.1]

    ;  plot MDN wet deposition mapped onto the GEOS-Chem grid
    qind=where(finite(data_regrid.avg_data) eq 0)
    data_regrid.avg_data(qind)=-9999.
    ctm_overlay,data_regrid.avg_data,data_regrid.lon,data_regrid.lat,$
          [obs_annual_tgm_amnet[amnet_background],obs_annual_tgm_camnet[camnet_background] ], $
          [site_lon_amnet[amnet_background], site_lon_camnet[camnet_background] ], $
          [site_lat_amnet[amnet_background],  site_lat_camnet[camnet_background] ], $
         /sample,/continents,/usa,thick=1, $
         mindata=1.2,maxdata=1.8,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!c(C)AMNet data for '+string(year,format='(I4,"-",I4)'),$
         cbunit= 'ng m!u-3!n',/iso,csfac=0.8,limit=limit,/cbar,$
	     min_valid=1.2,botoutofrange=!myct.mediumgray,tria=1,$    ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8),/conus

    tvmap, model_tgm_annual, longitude, latitude,$
         /sample,/continents,/usa,thick=1, $
         mindata=1.2, maxdata=1.8,limit=limit, /iso, $
         t_symbol=1,t_color=1,symsize=.3, /cbar, cbunit='ng m!u-3!n', $
         title='!c!cGEOS-Chem '+string(year,format='(I4,"-",I4)'),$
	     min_valid=1.2,botoutofrange=!myct.mediumgray,tria=1,$    ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    indlon=where(data_regrid.lon ge lon0(0) and data_regrid.lon le lon0(1))
    indlat=where(data_regrid.lat ge lat0(0) and data_regrid.lat le lat0(1))
    obs_tgm_regrid=data_regrid.avg_data
    obs_tgm_regrid=obs_tgm_regrid(indlon,*) & obs_tgm_regrid=obs_tgm_regrid(*,indlat)

    ;  set model grid points where no data is available to NaN
    qind=where( obs_tgm_regrid lt 0. )
    model_tgm_annual_sample=model_tgm_annual
    model_tgm_annual_sample(qind)=-9999.
 
    ;  plot model wet deposition only where there is AMNET data
     tvmap, model_tgm_annual_sample, longitude, latitude,$
         /sample,/continents,/usa,thick=1, $
         mindata=1.2, maxdata=1.8, limit=limit, /iso, /cbar, cbunit='ng m!u-3!n', $
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!cGEOS-Chem '+string(year,format='(I4,"-",I4)'),$
	     min_valid=1.2,botoutofrange=!myct.mediumgray,tria=1,$    ; cbpos=cbpos,$
	    div=COLORBAR_NDIV(maxdiv=8)

    Myct,/buwhwhrd,ncolors=10
    ;  plot the difference MDN-model
    model_tgm_annual_diff=model_tgm_annual_sample-obs_tgm_regrid
    model_tgm_annual_diff(qind)=-9999.
    tvmap, model_tgm_annual_diff, longitude, latitude,$
         /sample,/continents,/usa,thick=1,$
         mindata=-0.3,maxdata=0.3,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!cGEOS-Chem - (C)AMNet',$
         cbunit= 'ng m!u-3!n',/iso,csfac=0.8,limit=limit,/cbar,$
	     min_valid=-0.3,botoutofrange=!myct.mediumgray,tria=1,$   ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    ;  plot the percent difference
    qind=where( obs_tgm_regrid gt 0. )
    model_tgm_annual_diff_perc=obs_tgm_regrid*0.-9999.
    model_tgm_annual_diff_perc(qind)=-(obs_tgm_regrid(qind)-model_tgm_annual_sample(qind))/obs_tgm_regrid(qind)*100. ;  %    
    tvmap,model_tgm_annual_diff_perc,longitude,latitude,$
         /sample,/continents,/usa,thick=1,$
         mindata=-10,maxdata=10.0,$
         t_symbol=1,t_color=1,symsize=.3,$
         title='!c!c(GEOS-Chem - (C)AMNet) / (C)AMNet',$
         cbunit='%',/iso,csfac=0.8,limit=limit,/cbar,$
	     min_valid=-10,botoutofrange=!myct.mediumgray,tria=1,$    ; cbpos=cbpos,$
	     div=COLORBAR_NDIV(maxdiv=8)

    multipanel,/noerase,margin=[0.07,0.03,0.02,0.03],pos=pos 
    q=where(obs_tgm_regrid gt 0.)
    !p.charsize=1.5
    !x.charsize=!p.charsize & !y.charsize=!p.charsize
    plot, obs_tgm_regrid(q), model_tgm_annual_sample(q),col=1,xrange=[1.2,2.0],yrange=[1.2, 2.0],$
 	xtitle= '(C)AMNet ('+ 'ng m!u-3!n)', thick=1, $
	ytitle='GEOS-Chem TGM ('+ 'ng m!u-3!n)', $
	title='Annual Mean TGM',psym=sym(1),pos=pos,/xst,/yst
    rcorr=correlate(obs_tgm_regrid(q),model_tgm_annual_sample(q))
    org_corr,obs_tgm_regrid(q),model_tgm_annual_sample(q),rcorr,$
 	n_elements(obs_tgm_regrid(q)),slope,intercept
    oplot, !x.crange,intercept+!x.crange*slope,col=1,line=1,thick=4
    oplot, !x.crange,!x.crange,col=1
    xyouts, !x.crange(0)+(!x.crange(1)-!x.crange(0))/20., !y.crange(0)+(!y.crange(1)-!y.crange(0))*0.93,$
         /data,col=1, 'Model = '+string(intercept,slope,format='(F6.3,"+ (C)AMNet x",F5.2)')+$
         '!c(r= '+$
         string(rcorr,format='(F5.2)')+ ';  npts= '+$
         string(n_elements(obs_tgm_regrid(q)),format='(I6)')+ ')'+$
	     '!cModel/Obs='+string(mean(model_tgm_annual_sample(q)/obs_tgm_regrid(q)),format='(F5.2)'),charsize=0.9

     MyCt, /Verbose, /WhGrYlRd,ncolors=20
    !p.charsize=1.
    !x.charsize=!p.charsize & !y.charsize=!p.charsize
 endif



;close_device

; file_move,'idl.ps','psfile.ps',/overwrite    ;this will rename the file.png into any other name you wish
; spawn,'convert -density 500 psfile.ps file.png'            ;this will transfer the ps file into a png file named file.png
; file_move,'psfile.ps','CAMNet_AMNet_TGM_'+version+'_'+string(min(year),format='(I4)')+'_'+string(max(year),format='(I4)')+'_'+resolution+'.ps',/overwrite    ; this will rename the file.png into any other name you wish

; for i=0,n_elements(site_name)-1 do begin
;   print,site_name[i],obs_annual_tgm[i],mod_annual_tgm[i]
; endfor

;stop
end