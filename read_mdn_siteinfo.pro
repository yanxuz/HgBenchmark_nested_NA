pro read_mdn_siteinfo,sites,plot=plot

;+++++++++++++++++++++++++++++++++++++++++++++
; program reading in the MDN site info from ~/ALASKA/data/MDN/mdnsites.txt
;++++++++++++++++++++++++++++++++++++++++++++++

  file='./MDN_data/mdnsites.txt'

  ; first find out number of lines in file
  nlines=file_lines(file)

  ; Open and read the filen
  OPEN_FILE,file,ilun
  data_str=strarr(nlines)
  readf,ilun,data_str
  free_lun,ilun
  nvar=15
  array=strarr(nlines,nvar)

  for k=0L,nlines-1L do begin
      str=strsplit(data_str(k),",",count=count,/extract,/preserve)
      if count gt 0 then array(k,0:count-1)=str
  endfor
  header=reform(array(0,*))
  data=array(1:*,*)
  nsites=nlines-1

  ; create structure containing data for all the sites
  sites=replicate({name:'',id:'',status:'',lat:'',lon:'',alt:'',$
  		   start_day:'',start_month:'',start_year:'',$
  		   end_day:'',end_month:'',end_year:''},nlines-1)
  sites.name=data(*,1)
  sites.id=data(*,0)

  sites.status=data(*,where(header eq '"status"'  or header eq 'status'  ))
  sites.lat=data(*,where(header eq '"latitude"'  or header eq 'latitude'     ))
  sites.lon=data(*,where(header eq '"longitude"'  or header eq 'longitude' ))
  sites.alt=data(*,where(header eq '"elevation"'  or header eq 'elevation' ))
  start_date=data(*,where(header eq '"startdate"' or header eq 'startdate' ))
  end_date=data(*,where(header eq '"stopdate"' or header eq 'stopdate' ))
  
  for k=0,nsites-1 do begin
      str=strsplit(start_date(k),"/",count=count,/extract)
      if count eq 3 then begin
      	sites(k).start_day=str(1)
      	sites(k).start_month=str(0)
      	sites(k).start_year=str(2)
      endif
      str=strsplit(end_date(k),"/",count=count,/extract)
      if count eq 3 then begin
      	sites(k).end_day=str(1)
      	sites(k).end_month=str(0)
      	sites(k).end_year=str(2)
      endif
  endfor


  if keyword_Set(plot) then begin
    ; only plot sites north of 40N
    ind=where(sites.lat ge 40)
    sites=sites(ind)
    open_device,/ps,/color,file='MDN_sites.ps'
    multipanel,col=1,row=1,margin=[0.02,0.01],pos=pos
    map_set,50,-100,/ortho,/iso,col=1,pos=pos,limit=[40,-140,85,-50],title='MDN sites'
    map_continents,/countries,col=1,/continents,/hires
    oplot,sites.lon,sites.lat,psym=sym(3),col=2,symsize=0.5
    oplot,sites.lon,sites.lat,psym=sym(2),col=2,symsize=0.5
    for k=0,n_elements(sites.name)-1 do begin
      xyouts,sites(k).lon,sites(k).lat,sites(k).id,col=2,/data,charsize=0.8,al=(k mod 2)
      ;xyouts,0.7,0.99-0.03*k,$
    ;	string(sites(k).id,sites(k).lat,sites(k).lon,sites(k).start_year,sites(k).end_year,$
     ;   format='(A4,1x,F5.2,"N",F8.2,"E",1x,A4,"-",A4)'),col=1,/normal,charsize=0.7,al=0.5
    endfor
    close_Device
  endif

end
