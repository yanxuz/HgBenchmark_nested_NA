function alt2lev,lat,lon,alt,resolution=resolution
; determine the level from the surface elevation

; input
; lat : latitude
; lon: longitude
; alt: altitude, meters

zedge=[123.,254.,387.,521.,657.,795.,934.,1075.,1218.,1363.,1510.,1659.,1860.,2118.,2382.,2654.,2932.,3219.,3665.,4132.,4623.,5142.,$
         5692.,6277.,6905.,7582.,8320.,9409.]

if not keyword_set(resolution) then resolution= '4x5'

	; read surface elevation data
    filename='./phis/Phis.geos5.' + resolution

    success0=ctm_get_datablock(phis,'dao-flds',tracer=15,filename=filename,$
                    tau0=131472.0,lev=1,lat=float(lat),lon=float(lon))
	alt_effect=float(alt)-phis       ; m

	if alt_effect lt zedge[0] then begin
	    print,alt,phis
	    return,1
	endif
	
    for i=0,26 do begin
        if alt_effect ge zedge[i] and alt_effect lt zedge[i+1] then break
    endfor
    lev=i+2

	; print,alt,phis,alt_effect,lev
    return,lev
end