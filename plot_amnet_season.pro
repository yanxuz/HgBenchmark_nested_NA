pro plot_amnet_season,conc_all_bm, conc_all_rf
; plot seasonal cycle of speciated Hg concentrations

 multipanel,row=4, col=3, margin=[0.03,0.05], pos=pos
 x=indgen(12)
 monthstring=['J','F','M','A','M','J','J','A','S','O','N','D']

sites=['CA48', 'MD08', 'MD98', 'MD99', 'MS12', 'MS99', 'NH06', 'NJ05', 'NJ30', 'NJ32', 'NJ54', 'NS01', 'NY06', 'NY20', 'NY43', 'NY95', 'OH02', 'UT96', 'UT97', 'VT99', 'WV99']
lat=[36.81, 39.7053, 39.028, 39.028, 30.4294, 30.4294, 43.11, 39.4, 40.4728, 40.7876, 40.6414, 44.4328, 40.868, 43.9731, 43.1463, 43.1463, 39.3078, 41.0467, 40.7118, 44.5283, 39.0636]
lon=[-121.78, -79.0122, -76.8171, -76.8171, -88.4277, -88.4277, -70.95, -74.37, -74.4226, -74.6763, -74.2084, -65.2056, -73.8782, -74.2231, -77.5481, -77.5481, -82.1182, -112.0248, -111.9609, -72.8684, -79.4222]
nsite=n_elements(sites)

for isite=0,nsite-1 do begin
 print,'new reading ... ' + sites(isite)
 read_amnet,sites(isite),data
 y=[data.gem_monthly+data.gom_monthly/1000.]
 y_std=sqrt(data.gem_monthly_std^2+(data.gom_monthly_std/1000.)^2)
 plot, pos=pos, x , y, symsize=1.0,color=!myct.black,thick=4,$
        	yrange=[0.5,3.0],  xrange=[-0.5, 11.5], $
			            xtitle='Month', xticks=11,xtickname=monthstring, $
			            xtickv=indgen(12), xminor=1,ytitle='Hg0 concentration, ng/m3'
 ctm_index,ctm_type('GEOS5_47L',res=[2./3.,1./2.]),i,j,center=[lat(isite),lon(isite)],/non_interactive
 y=conc_all_bm.hg0(i-61,j-201)
 oplot, thick=4, x, y, color=!myct.red
 y=conc_all_rf.hg0(i-61,j-201)
 oplot, thick=4, x, y, color=!myct.green
 legend,label=['Observation', 'New model', 'Old model'],lcol=[1, 2, 3],line=[0,0,0],symbol=[-1,-1,-1],boxcolor=-1,thick=4,$
 	                 halign=0.01,valign=0.9,charsize=0.6
 xyouts,/data,9,2.7, sites(isite),col=1,charsize=0.8
 if isite eq 0 then xyouts,/data,21,3.5, 'Seasonal cycle of speciated Hg concentration at AMNet sites',col=1,charsize=1.2, align=0.5


multipanel,/advance, pos=pos

 y=[data.gom_monthly]
 y_std=[data.gom_monthly_std]
 plot, pos=pos, x , y, symsize=1.0,color=!myct.black,thick=4,$
        	yrange=[0,60],  xrange=[-0.5, 11.5], $
			            xtitle='Month', xticks=11,xtickname=monthstring, $
			            xtickv=indgen(12), xminor=1,ytitle='RGM concentration, pg/m3'
 y=1000.0*conc_all_bm.newrgm(i-61,j-201)
 oplot, thick=4, x, y, color=!myct.red
 y=1000.0*conc_all_rf.newrgm(i-61,j-201)
 oplot, thick=4, x, y, color=!myct.green
multipanel,/advance, pos=pos

 y=[data.pbm_monthly]
 y_std=[data.pbm_monthly_std]
 plot, pos=pos, x , y, symsize=1.0,color=!myct.black,thick=4,$
        	yrange=[0,60],  xrange=[-0.5, 11.5], $
			            xtitle='Month', xticks=11,xtickname=monthstring, $
			            xtickv=indgen(12), xminor=1,ytitle='PBM concentration, pg/m3'
 y=1000.0*conc_all_bm.newrpm(i-61,j-201)
 oplot, thick=4, x, y, color=!myct.red
 y=1000.0*conc_all_rf.newrpm(i-61,j-201)
 oplot, thick=4, x, y, color=!myct.green
multipanel,/advance, pos=pos

endfor

end