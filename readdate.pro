pro readdate,str,year,month,day,hour,minute

len=strlen(str)

pos=strpos(str,'/')
month=fix(strmid(str,0,pos))
len=len-pos-1
str=strmid(str,pos+1,len)

pos=strpos(str,'/')
day=fix(strmid(str,0,pos))
len=len-pos-1
str=strmid(str,pos+1,len)

pos=4
year=fix(strmid(str,0,pos))
len=len-pos-1
str=strmid(str,pos+1,len)

pos=strpos(str,':')
hour=fix(strmid(str,0,pos))
len=len-pos-1
str=strmid(str,pos+1,len)

minute=fix(strmid(str,0,len))
end