function vicinity_emission,lat,lon,distance
; return the Hg0, Hg2, and HgP emission near a given location lat, lon
; emission: the emission amount
; distance: a circle near the point within the given distance
; created by yanxu, 09/13/2011

; emission file names
filename0 = './emission/NEI2005_NPRI2005_NewSpeciation_Hg0.geos5.05x0666.2005'
filename2 = './emission/NEI2005_NPRI2005_NewSpeciation_Hg2.geos5.05x0666.2005'
filenamep = './emission/NEI2005_NPRI2005_NewSpeciation_HgP.geos5.05x0666.2005'

; initiate emission
emission = {hg0:0.0,hg2:0.0,hgp:0.0}

; get emission data
diagname='hg-srce'
tracer=1
lev=1
infile=filename0
success=ctm_get_datablock(data0,diagname,                    $
			 x=longitude, y=latitude,                             $
             tracer=tracer,                         $
             filename=infile,                       $
             tau0=nymd2tau(ymd2date(2005,1,1)),     $
             lev=lev,lat=[10.0, 70.0],lon=[-140.333, -39.667]               $
             )
if not success then stop

diagname='hg-srce'
tracer=6
lev=1
infile=filename2
success=ctm_get_datablock(data2,diagname,                    $
			 x=longitude, y=latitude,                             $
             tracer=tracer,                         $
             filename=infile,                       $
             tau0=nymd2tau(ymd2date(2005,1,1)),     $
             lev=lev,lat=[10.0, 70.0],lon=[-140.333, -39.667]               $
             )
if not success then stop

diagname='hg-srce'
tracer=9
lev=1
infile=filenamep
success=ctm_get_datablock(datap,diagname,                    $
			 x=longitude, y=latitude,                             $
             tracer=tracer,                         $
             filename=infile,                       $
             tau0=nymd2tau(ymd2date(2005,1,1)),     $
             lev=lev,lat=[10.0, 70.0],lon=[-140.333, -39.667]               $
             )
if not success then stop

; add the emission together within the given distance
pi = 4.0 * atan(1.0)
deg2rad = pi/180.0
; range = 2*fix(distance/50.0)

for i=0,150 do begin
  for j=0,120 do begin
;            if abs(latitude[j] - lat) gt range or abs(longitude[i] - lon) gt range then continue
            d = 6372.8*acos(sin(latitude[j]*deg2rad)*sin(lat*deg2rad)+cos(latitude[j]*deg2rad)*cos(lat*deg2rad)*cos((lon-longitude[i])*deg2rad))  ; km
            if ( d le distance) then begin
                emission.hg0=emission.hg0+data0[i,j]
                emission.hg2=emission.hg2+data2[i,j]
                emission.hgp=emission.hgp+datap[i,j]
            endif
  endfor
endfor

thg=emission.hg0+emission.hg2+emission.hgp
return,thg
end