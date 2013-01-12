pro read_amnet_siteinfo,data_all,plot=plot, mode=mode
 
 ; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ; This procedure reads in header info for all the AMNET files
 ; and prints out info about sites, number of years for each site (if /debug is set)
 ; Also, plots map with location of sites on top (if /plot is set).
 ; ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 dir='./AMNet_data/'

if keyword_set(mode) then begin
 if mode eq 's' then filename=dir+'siteinfo_select.txt'
 if mode eq 'c' then  filename=dir+'siteinfo_complete.txt'
endif else begin
 ; default case is select, which is selected for 1st paper
 filename=dir+'siteinfo_select.txt'
endelse
 nlines=file_lines(filename)
 nsites=nlines-1
 open_file,filename,ilun
 
 data_str=strarr(nlines)
 readf,ilun,data_str
 free_lun,ilun
 nvar=5
 array=strarr(nlines,nvar)

 for k=0,nlines-1 do begin
   str=strsplit(data_str(k),",",count=count,/extract)
   if count gt 0 then array(k,0:count-1)=str
 endfor

 ; find begining of tables
 struc={siteid:'',sitename:'',lat:0.,lon:0.,ele:0.}
 data_all=replicate(struc,nsites)

 for k=0,nsites-1 do begin
   data_all(k).siteid=array(k+1,0)
   data_all(k).sitename=array(k+1,1)
   data_all(k).lat=array(k+1,2)
   data_all(k).lon=array(k+1,3)
   data_all(k).ele=array(k+1,4)
 endfor
 
 ; plot the sites location
 if keyword_Set(plot) then begin
  open_device,/ps,/color
  ;multipanel,col=2,row=1,margin=[0.02,0.01],pos=pos
  map_set,50,-100,/ortho,/iso,col=1,/continents,pos=pos,limit=[20,-140,85,-60],title='AMNET sites'
  map_continents,/countries,/continents,col=1
  oplot,data_all.lon,data_all.lat,psym=sym(1),col=2,symsize=1.5
  close_Device
 endif


end
