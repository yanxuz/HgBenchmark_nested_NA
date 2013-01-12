function compute_grid_mdn, Longitude, Latitude, Data, Lon, Lat, $
                      Bin_Data=Bin_Data, Bin_Ct=Bin_Ct, $
		      Var_Data=Var_Data,RES=RES,$
		      gridinfo=gridinfo,Verbose=Verbose

; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;  This function will grid the MDN data onto a regular grid
;  
;  INPUTS:
;       LATITUDE     --> Latitudes of input data points 
;       LONGITUDE    --> Longitudes of input data points
;       DATA         --> Values of input data points
;       
;  KEYWORDS
;       RES         --> '4x5' or '05x0666' resolution of grid
; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   

   if not keyword_set(RES) then RES='4x5' ;  default resolution

   ;  determine grid size
   if  RES eq '4x5' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=4 ) 
   if  RES eq '2x25' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=2 ) 
   if  RES eq '05x0666' then ModelInfo  = CTM_TYPE( 'GEOS5_47L', res=[2./3.,1./2.] ) 

   GridInfo   = CTM_GRID( ModelInfo )
   Lon=gridinfo.xmid & Lat=gridinfo.ymid
   Lon_edge=gridinfo.xedge & Lat_edge=gridinfo.yedge
   Lat_res=gridinfo.dj & Lon_res=gridinfo.di
   nLat=n_elements(Lat)
   nLon=n_elements(Lon)

   ;  Define arrays
   Avg_Data = DblArr( GridInfo.IMX, GridInfo.JMX ) + !values.f_nan
   Bin_Data = DblArr( GridInfo.IMX, GridInfo.JMX ) 
   Bin_Ct   = DblArr( GridInfo.IMX, GridInfo.JMX )

   ;  Loop over sites
   For N = 0,n_elements(Longitude)-1 do begin
    if finite(Data(N) ) then begin
      ctm_index,Modelinfo,i,j,$
              center=[Latitude(N),Longitude(N)],/non_interactive
      if RES eq '05x0666' then begin
        Bin_Ct[I-1,J-1]=Bin_Ct[I-1,J-1]+1L
        Bin_Data[I-1,J-1]=Bin_Data[I-1,J-1]+Data(N)
      endif else begin
        ;  check whether grids around should be included (within 0.1*grid)
        xdis  = abs( Longitude(N)-Lon_edge(I-1) )
        xdis1 = abs( Longitude(N)-Lon_edge(I) )
        ydis  = abs( Latitude(N)-Lat_edge(J-1) )
        ydis1 = abs( Latitude(N)-Lat_edge(J) )
        points=1.
        iuse=0 & juse=0
        if mean(ydis)  lt 0.1*Lat_res then juse=[juse,-1]
        if mean(ydis1) lt 0.1*Lat_res then juse=[juse,+1]
        if mean(xdis)  lt 0.1*Lon_res then iuse=[iuse,-1]
        if mean(xdis1) lt 0.1*Lon_res then iuse=[iuse,+1]
        counts=n_elements(iuse)+n_elements(juse)-1.
        if n_elements(iuse) eq 2 and n_elements(juse) eq 2 then counts=counts+1.
        for k=0,n_elements(iuse)-1 do Bin_Ct[I-1+iuse(k),J-1+juse]=Bin_Ct[I-1+iuse(k),J-1+juse]+1./float(counts)
        for k=0,n_elements(iuse)-1 do Bin_Data[I-1+iuse(k),J-1+juse]=Bin_Data[I-1+iuse(k),J-1+juse]+Data(N)/float(counts)
      endelse

      if keyword_set(verbose) then begin
        print,Latitude(N),Longitude(N),Lat(j-1),Lon(i-1),Bin_Ct[I-1,J-1],Bin_Data[I-1,J-1]
      endif
    endif
   Endfor

  ;  Compute the average -- exclude missing values (negatives)
  for J=0L, GridInfo.JMX-1L do begin
  for I=0L, GridInfo.IMX-1L do begin
       if ( Bin_Ct[I,J] gt 0L ) then begin
                 Avg_Data[I,J]= Bin_Data[I,J] / Float( Bin_Ct[I,J] )
       endif 
  endfor
  endfor

  data_regrid={lat:lat,lon:lon,avg_data:avg_data,npts:bin_ct}

  ;  Return the average data to the calling program
  return, data_regrid
end
