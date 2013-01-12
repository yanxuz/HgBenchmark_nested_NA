;
;-----------------------------------------------------------------------

pro mercury_benchmark_nested_na, filename=filename, reference=reference, $
       psfilename=psfilename


   ; conversion factors
   pptv_ngm3 = 8.93d0
   ngm3_ppqv = 112d0

   ; reference
   if not keyword_set(reference) then reference='default.05x0666.NA.2009'

   ; print ps file
   if not keyword_set(psfilename) then psfilename='benchmark_nested_na.ps'
   open_device,/ps,/color,/portrait,file=psfilename
   !p.font=0
   device, /helvetica, /isolatin1
   multipanel,col=2,row=3,margin=[0.05,0.05],pos=pos
   !p.charsize=1.5

   ; read model data
   read_run_all,file=filename,$
	conc_all_bm,emis_all_bm,drydep_all_bm,wetdep_all_bm,chem_all_bm,$
	longitude,latitude,convert,prof_all_bm,gridinfo,transp_all_bm

   read_run_all,file=reference,$
	conc_all_rf,emis_all_rf,drydep_all_rf,wetdep_all_rf,chem_all_rf,$
	longitude,latitude,convert,prof_all_rf,gridinfo,transp_all_rf,/hg2frac

   ; tgm
   plot_amnet_camnet_tgm,filename,'05x0666',2009,title='new model TGM'
   plot_amnet_camnet_tgm,reference,'05x0666',2009,/hg2frac,title='old model TGM'

   ; speciated hg concentration
   plot_amnet_hgs,filename,'05x0666',2009,title='new model speciated Hg'
   plot_amnet_hgs,reference,'05x0666',2009,/hg2frac,title='old model speciated Hg'

   ; speciated hg season cycle
   plot_amnet_season,conc_all_bm,conc_all_rf

   ; mdn
   plot_mdn,filename,'05x0666',sel=2009,year=2009,mett='GEOS5',metr='05x0666',$
            version2=reference,resolution2='05x0666',/season,/region
   ;plot_mdn,reference,'05x0666',sel=2009,year=2009,mett='GEOS5',metr='05x0666'

   ; compare between models
   ; concentrations
   plotdiff,mean(conc_all_bm.hg0,3),mean(conc_all_rf.hg0,3),longitude,latitude,$
                  mindata=0.3,maxdata=2.0,mindiff=-0.5,maxdiff=0.5,minperc=-50,maxperc=50,$
                  unit='ng/m3',title='Hg0 concentration at surface',/showmean
   plotdiff,1000.0*mean(conc_all_bm.newrgm,3),1000.0*mean(conc_all_rf.newrgm,3),longitude,latitude,$
                  mindata=0,maxdata=40,mindiff=-20,maxdiff=20,minperc=-50,maxperc=50,$
                  unit='pg/m3',title='RGM concentration at surface',/showmean
   plotdiff,1000.0*mean(conc_all_bm.newrpm,3),1000.0*mean(conc_all_rf.newrpm,3),longitude,latitude,$
                  mindata=0,maxdata=20,mindiff=-20,maxdiff=20,minperc=-50,maxperc=50,$
                  unit='pg/m3',title='RPM concentration at surface',/showmean

   ; emissions
   plotdiff,total(emis_all_bm.anthr,3),total(emis_all_rf.anthr,3),longitude,latitude,$
           mindata=0,maxdata=1.0,mindiff=-0.5,maxdiff=0.5,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Anthropogenic Hg0 emission'
   plotdiff,total(emis_all_bm.anthr_rgm,3),total(emis_all_rf.anthr_rgm,3),longitude,latitude,$
           mindata=0,maxdata=0.1,mindiff=-0.05,maxdiff=0.05,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Anthropogenic Hg2 emission'
   plotdiff,total(emis_all_bm.anthr_hgp,3),total(emis_all_rf.anthr_hgp,3),longitude,latitude,$
           mindata=0,maxdata=0.1,mindiff=-0.05,maxdiff=0.05,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Anthropogenic HgP emission'
   plotdiff,total(emis_all_bm.ocean,3),total(emis_all_rf.ocean,3),longitude,latitude,$
           mindata=0,maxdata=0.1,mindiff=-0.02,maxdiff=0.02,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Oceanic Hg0 emission'
   plotdiff,total(emis_all_bm.land,3),total(emis_all_rf.land,3),longitude,latitude,$
           mindata=0,maxdata=0.05,mindiff=-0.01,maxdiff=0.01,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Soil Hg0 emission (reflective deposition)'
   plotdiff,total(emis_all_bm.land_natural,3),total(emis_all_rf.land_natural,3),longitude,latitude,$
           mindata=0,maxdata=0.05,mindiff=-0.01,maxdiff=0.01,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Geogenic Hg0 emission'
   plotdiff,total(emis_all_bm.bb,3),total(emis_all_rf.bb,3),longitude,latitude,$
           mindata=0,maxdata=0.05,mindiff=-0.01,maxdiff=0.01,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Biomass burning Hg0 emission'
   plotdiff,total(emis_all_bm.soil,3),total(emis_all_rf.soil,3),longitude,latitude,$
           mindata=0,maxdata=0.05,mindiff=-0.01,maxdiff=0.01,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Soil Hg0 emission'

   ; dry depositions
   plotdiff,total(drydep_all_bm.hg0,3),total(drydep_all_rf.hg0,3),longitude,latitude,$
           mindata=0,maxdata=0.1,mindiff=-0.02,maxdiff=0.02,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Hg0 dry deposition to land'
   plotdiff,total(drydep_all_bm.rgm,3)+total(drydep_all_bm.hgp,3),$
            total(drydep_all_rf.rgm,3)+total(drydep_all_rf.hgp,3),$
           longitude,latitude,$
           mindata=0,maxdata=0.1,mindiff=-0.05,maxdiff=0.05,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='HgII dry deposition'
   plotdiff,total(drydep_all_bm.rgm_ss,3),total(drydep_all_rf.rgm_ss,3),longitude,latitude,$
           mindata=0,maxdata=0.05,mindiff=-0.01,maxdiff=0.01,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='HgII seasalt loss'

   ; wet deposition
   plotdiff,total(wetdep_all_bm.rgm,3)+total(wetdep_all_bm.hgp,3),$
            total(wetdep_all_rf.rgm,3)+total(wetdep_all_rf.hgp,3),$
           longitude,latitude,$
           mindata=0,maxdata=0.1,mindiff=-0.02,maxdiff=0.02,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='HgII wet deposition'

   ; chemistry
   plotdiff,total(total(chem_all_bm.net_ox,4),3),$
            total(total(chem_all_rf.net_ox,4),3),$
           longitude,latitude,$
           mindata=0,maxdata=0.05,mindiff=-0.01,maxdiff=0.01,minperc=-50,maxperc=50,$
           unit='Mg/yr',title='Hg0 nex oxidation'

   geos_z=[0.071,0.201,0.332,0.466,0.601,0.737,0.875,1.016,1.157,1.301,1.447,$
      1.594,1.769,1.999,2.259,2.527,2.801,3.084,3.448,3.904,4.382,4.886,$
      5.419,5.985,6.591,7.241,7.947,8.848,9.938,11.021,12.086,13.134,14.17,$
      15.198,16.222,17.243,18.727,20.836,23.02,25.307,28.654,34.024,40.166,$
      47.135,54.834,63.054,72.18]
   volume=ctm_boxsize(ctm_grid(ctm_type('GEOS5_47L',res=[2.0/3.0,0.5])),/m3,/volume)
   volume=volume(60:210,*,*)
   volume=volume(*,200:320,*)

   plotdiff,1e15*mean(total(chem_all_bm.net_ox,4)/volume,1),$
            1e15*mean(total(chem_all_rf.net_ox,4)/volume,1),$
           latitude,geos_z,$
           mindata=0,maxdata=2,mindiff=-0.2,maxdiff=0.2,minperc=-50,maxperc=50,$
           unit='ng/m3yr',title='Zonal mean Hg0 nex oxidation',/zonal,/showmean

   plotdiff,mean(total(chem_all_bm.br,4),1)/1e8,$
            mean(total(chem_all_rf.br,4),1)/1e8,$
           latitude,geos_z,$
           mindata=0,maxdata=1,mindiff=-1,maxdiff=1,minperc=-100,maxperc=100,$
           unit='10!u8!n molec/cm3',title='Zonal mean Br concentration',/zonal,/showmean

   close_device
stop
end
