'sdfopen http://nomads.ncep.noaa.gov:80/dods/gfs_0p50/gfs20191119/gfs_0p50_00z'
*'set t 1 17'
'define pres = 'pressfc
'set sdfwrite pres.nc'
'sdfwrite 'pres
'define z1000 = 'hgtprs'(lev=1000)'
'set sdfwrite z1000.nc'
'sdfwrite 'z1000
'define z500 = 'hgtprs'(lev=500)'
'set sdfwrite z500.nc'
'sdfwrite 'z500
'define u250 = 'ugrdprs'(lev=250)'
'set sdfwrite u.nc'
'sdfwrite 'u250
'define v250 = 'vgrdprs'(lev=250)'
'set sdfwrite v.nc'
'sdfwrite 'v250
'define sfcz = 'hgtsfc
'set sdfwrite sfcz.nc'
'sdfwrite 'sfcz
'define t2m = 'tmp2m
'set sdfwrite t2m.nc'
'sdfwrite 't2m
