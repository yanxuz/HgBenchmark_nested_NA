pro mdn_regions,region,lat0,lon0,$
	latitude=latitude,longitude=longitude,$
	ind_lat,ind_lon

; +++++++++++++++++++++++++++++++++++++++++++++++++
; This procedure defines MDN regions and
; returns indices for lon/lat for each region.
; Regions: Midwest (MW), Norteast (NE), Ohio River Valley (OH)
;          Southeast (SE), Gulf Coast (Gulf)
; +++++++++++++++++++++++++++++++++++++++++++++++++
case region of
    'W':begin
	   ; West
	   lat0=[30,50]
	   lon0=[-118,-104]
	   end
	'C':begin
       ; Central
       lat0=[24,54]
       lon0=[-104,-95]
	   end
    'PW':begin
       ;pacific west
       lat0=[32,50]
       lon0=[-125,-118]
	   end
	'MW':begin
	    ; Midwest
		lat0=[42,50]
		lon0=[-95,-80]
	     end
	'NE':begin
		lat0=[42,47]
		lon0=[-80,-55]
	     end
	'ORV':begin
		lat0=[35,42] 
		lon0=[-95,-65]
	     end
	'SE':begin
		lat0=[24,35]
		lon0=[-95,-76]
	     end
	'Gulf':begin
		lat0=[24.5,32.5]
		lon0=[-95,-85]
		;lat0=[24.5,33]
		;lon0=[-95,-78]
	     end
	'Florida':begin
		lat0=[24.5,30]
		lon0=[-85,-78]
             end
     'USA': begin
		lat0=[24.5,48]
		lon0=[-125,-60]
	  end
      'ALL':begin
        lat0=[10,70]
        lon0=[-140,-40]
      end
; now we define some new regions for the second paper
      'WEST':begin
         lat0=[24,49]
         lon0=[-125,-100]
      end
      'MIDWEST':begin
         lat0=[36,49]
         lon0=[-100,-82.5]
      end
      'NORTHEAST':begin
         lat0=[36,49]
         lon0=[-82.5,-60]
      end
      'SOUTHEAST':begin
         lat0=[24,36]
         lon0=[-100,-75]
      end
      'PENNSYLVANIA':begin
         lat0=[39.2,42.5]
         lon0=[-81.3,-73.7]
      end
	  
	 else: begin
	 	print,'Region not found'
		stop
    		end
endcase
 if keyword_set(latitude) then $
 	ind_lat=where(latitude ge lat0(0) and latitude le lat0(1))
 if keyword_set(longitude) then $
 	ind_lon=where(longitude ge lon0(0) and longitude le lon0(1))
end
