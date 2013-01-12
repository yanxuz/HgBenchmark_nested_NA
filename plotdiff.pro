pro plotdiff,data1,data2,longitude,latitude,$
              mindata=mindata,maxdata=maxdata,$
              mindiff=mindiff,maxdiff=maxdiff,$
              minperc=minperc,maxperc=maxperc,$
              unit=unit,title=title,showmean=showmean,zonal=zonal

multipanel,col=2,row=2

if not keyword_set(showmean) then begin
total1=total(data1,/nan)
total2=total(data2,/nan)
endif else begin
total1=mean(data1,/nan)
total2=mean(data2,/nan)
endelse
diff=total1-total2
perc=diff/total2*100.0

if keyword_set(zonal) then begin
 myct,/WhGrYlRd 
tvplot,data1,longitude,latitude,/sample,/cbar,division=4,cbunit=unit,mindata=mindata,maxdata=maxdata,title='new benchmark '+string(total1,format='(f7.3)'),yrange=[0,30]
tvplot,data2,longitude,latitude,/sample,/cbar,division=4,cbunit=unit,mindata=mindata,maxdata=maxdata,title='old benchmark '+string(total2,format='(f7.3)'),yrange=[0,30]
 myct,/BuWhRd 
tvplot,(data1-data2),longitude,latitude,/sample,/cbar,division=4,cbunit=unit,mindata=mindiff,maxdata=maxdiff,title='new - old '+string(diff,format='(f7.3)'),yrange=[0,30]
tvplot,((data1-data2)/data2),longitude,latitude,/sample,/cbar,division=4,cbunit=unit,mindata=minperc,maxdata=maxperc,title='100(new - old)/old '+string(perc,format='(f7.3)')+'%',yrange=[0,30]
 myct,/WhGrYlRd 
xyouts,10.0,9.5,title,col=1,align=0.5,charsize=1.2

endif else begin

 myct,/WhGrYlRd 
tvmap,data1,longitude,latitude,mindata=mindata,maxdata=maxdata,/cbar,div=5,/usa,title='new benchmark '+string(total1,format='(f7.3)'),/sample,/coast,cbunit=unit,/iso
tvmap,data2,longitude,latitude,mindata=mindata,maxdata=maxdata,/cbar,div=5,/usa,title='old benchmark '+string(total2,format='(f7.3)'),/sample,/coast,cbunit=unit,/iso
 myct,/BuWhRd 
tvmap,data1-data2,longitude,latitude,mindata=mindiff,maxdata=maxdiff,/cbar,div=5,/usa,title='new - old '+string(diff,format='(f7.3)'),/sample,/coast,cbunit=unit,/iso
temp=100.0*(data1-data2)/data2
ind=where(finite(temp) eq 0)
temp(ind)=0.0
tvmap,temp,longitude,latitude,mindata=minperc,maxdata=maxperc,/cbar,div=5,/usa,title='100(new - old)/old '+string(perc,format='(f7.3)')+'%',/sample,/coast,cbunit='%',/iso
 myct,/WhGrYlRd 
xyouts,10.0,9.0,title,col=1,align=0.5,charsize=1.2

endelse

end