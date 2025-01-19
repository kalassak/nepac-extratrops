import sys, math
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import minimum_filter, maximum_filter
from mpl_toolkits.basemap import Basemap, addcyclic

print "this script has even started"

def floodfill(a, x, y, min_pres, lows, fillx, filly, origx, origy, traversed, invalid):
	if (origx, origy) in invalid:
		#print "closing all recursive instances, already found to be invalid"
		return 0
	if calc_dist(x, y, origx, origy) > 2500:
		#is 2500 km too close?
		print "too far from center, assume low not deep enough"
		invalid.append((origx, origy))
		return 0
	if x != origx and y != origy:
		if (x, y) in lows:
			print "this low is not 8 hPa deep! %s %s within p_min+8" % (x, y)
			invalid.append((origx, origy))
			return 0
	if a[y][x] < min_pres+8 and (x, y) not in traversed:
		traversed.append((x, y))
		if x > 0:
			floodfill(a,x-1,y, min_pres, lows, fillx, filly, origx, origy, traversed, invalid)
		if x < len(a[y]) - 1:
			floodfill(a,x+1,y, min_pres, lows, fillx, filly, origx, origy, traversed, invalid)
		if y > 0:
			floodfill(a,x,y-1, min_pres, lows, fillx, filly, origx, origy, traversed, invalid)
		if y < len(a) - 1:
			floodfill(a,x,y+1, min_pres, lows, fillx, filly, origx, origy, traversed, invalid)
		fillx.append(x)
		filly.append(y)

def calc_dist(lat1,lon1,lat2,lon2):
	R = 6371.0088
	lat1 = np.radians([lat1])
	lon1 = np.radians([lon1])
	lat2 = np.radians([lat2])
	lon2 = np.radians([lon2])

	dlat = lat2 - lat1
	dlon = lon2 - lon1
	a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2) **2
	c = 2 * np.arctan2(a**0.5, (1-a)**0.5)
	d = R * c
	return d

def fmt_lon(n):
	if n > 180.0:
		s = "%.1fW" % math.fabs(n-360)
	else:
		s = "%.1fE" % n
	return s

def fmt_lat(n):
	if n > 0.0:
		s = "%.1fN" % n
	else:
		s = "%.1fS" % math.fabs(n)
	return s

PLOT_DEPTH_FILL = False
YYYYMMDD = sys.argv[1]
HH = sys.argv[2]

nc = netCDF4.Dataset("pres.nc")
nc2 = netCDF4.Dataset("z1000.nc")
nc3 = netCDF4.Dataset("z500.nc")
nc4 = netCDF4.Dataset("u.nc")
nc5 = netCDF4.Dataset("v.nc")
nc6 = netCDF4.Dataset("sfcz.nc")
nc7 = netCDF4.Dataset("t2m.nc")
nc8 = netCDF4.Dataset("u10m.nc")
nc9 = netCDF4.Dataset("v10m.nc")

lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]

psfc = nc.variables['pres'][:]
hgt1000 = nc2.variables['z1000'][:]
hgt500 = nc3.variables['z500'][:]
u = nc4.variables['u250'][:]
v = nc5.variables['v250'][:]
sfcz = nc6.variables['sfcz'][:]
t2m = nc7.variables['t2m'][:]
u10m = nc8.variables['u10m'][:]
v10m = nc9.variables['v10m'][:]

g = 9.81
R_d = 287

thicks = hgt500 - hgt1000
wind250 = (u**2 + v**2)**0.5*1.94384449 #m/s to kt
wind10m = (u10m**2 + v10m**2)**0.5*1.94384449 #m/s to kt
mslp = psfc * np.exp(sfcz*g/R_d/t2m)/100

print np.amax(thicks)
print np.amin(thicks)

WIND_SAMPLE_RATE = 6

nc.close()
nc2.close()
nc3.close()
nc4.close()
nc5.close()
nc6.close()
nc7.close()
nc8.close()
nc9.close()

# find low centers
data_ext = minimum_filter(mslp, 50, mode='nearest')
mxy, mxx = np.where(data_ext == mslp)

# get lat/lon of all low centers
low_lons = lons[mxx]
low_lats = lats[mxy]

nepac_low_lats = []
nepac_low_lons = []
nepac_xs = []
nepac_ys = []
lows = [] #for low coordinate pairs
for lon, lat, x, y in zip(low_lons, low_lats, mxx, mxy):
	lows.append((x, y))
	if lon >= 140 and lon <= 250 and lat >= 30 and lat <= 70:
		nepac_low_lats.append(lat)
		nepac_low_lons.append(lon)
		nepac_xs.append(x)
		nepac_ys.append(y)

#print nepac_low_lats
#print nepac_low_lons

#find wind maxima
data_ext = maximum_filter(wind10m, 50, mode='nearest')
mwy, mwx = np.where(data_ext == wind10m)

wind_lons = lons[mwx]
wind_lats = lats[mwy]

#m = Basemap(projection='cyl', llcrnrlat=20,urcrnrlat=70,llcrnrlon=360-180,urcrnrlon=360-110,resolution='l') #-90 90 -180 180
m = Basemap(width=9000000,height=6000000,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=50,lon_0=-160.)

ucolors = [(1, 1, 1), (0, 1, 1), (0, 0.84, 0.84), (0, 0.74, 0.74), (0, 1, 0.64), (0, 0.9, 0.59), (1, 1, 0), (1, 0.75, 0), (1, 0.5, 0),
	(1, 0, 0), (1, 0, 1), (1, 0.5, 1)]
ulevs = [0, 10, 20, 27, 34, 42, 49, 56, 64, 85, 100, 120]
jetcolors = [(1, 1, 1), (0, 0.93, 1), (0, 0.7, 1), (0, 0.35, 1), (0, 0, 1), (0.3, 0, 1), (0.5, 0, 1), (0.88, 0.27, 1), (1, 0.7, 1)] #(0, 1, 0.82), 
jetlevs = [0, 64, 85, 100, 120, 140, 160, 180, 200, 225, 250]
qcolors = [(1, 1, 1), (0.5, 1, 0.5), (0.3, 0.8, 0.3), (0.1, 0.6, 0.1), (0, 0.2, 0.6), (0, 0.3, 0.8)]
qlevs = [0, 0.5, 1, 2, 4, 8, 16]
thicklevels = np.arange(0, 12000, 60)
preslevs = np.arange(900, 1060, 4)

parallels = np.arange(0, 90, 10)
meridians = np.arange(0, 360, 10)

fig = plt.figure(figsize=(18.6, 10.5))
ax = fig.add_axes((0,0,1,1))

m.drawcoastlines()
m.drawparallels(parallels)
m.drawmeridians(meridians)

# for basemap
lons2, lats2 = np.meshgrid(lons, lats)
xx, yy = m(lons2, lats2)

MXX, MXY = m(lons[mxx], lats[mxy]) #convert pressure mins to map coordinates
MWX, MWY = m(lons[mwx], lats[mwy]) #convert wind maxs

# plot contours
print "plotting thickness"
con_thick = m.contour(xx, yy, thicks, levels=thicklevels, linewidths=1, linestyles='dashed', colors='r')
plt.clabel(con_thick, con_thick.levels, fmt='%d')
print "plotting psfc"
m.contour(xx, yy, mslp, levels=preslevs, colors='k')
print "plotting winds"
m.contourf(xx, yy, wind250, levels=jetlevs, colors=jetcolors, zorder=0)

#plot cyclones
plt.plot(MXX, MXY, 'ro')
#determine valid nepac cyclones
#lows.append((451,305))
valid = []
for nepac_x, nepac_y in zip(nepac_xs, nepac_ys):
	print 'checking depth of %s %s' % (nepac_x, nepac_y)
	fillx = []
	filly = []
	traversed = []
	invalid = []
	floodfill(mslp, nepac_x, nepac_y, mslp[nepac_y,nepac_x], lows, fillx, filly, nepac_x, nepac_y, traversed, invalid)
	if (nepac_x, nepac_y) in invalid:
		print 'this is not nameable cyclone'
		valid.append(['name', lats[nepac_y], lons[nepac_x], 0, mslp[nepac_y,nepac_x], 'EX', 0., 0.])
		plotopts = 'mo'
	else:
		#check if the storm is associated with gales
		#determine max 10m wind speed
		last_wind = 0
		storm_wind = 0
		maxwlo = 0.
		maxwla = 0.
		for wlo, wla, wx, wy in zip(wind_lons, wind_lats, mwx, mwy):
			dist = calc_dist(wla, wlo, lats[nepac_y], lons[nepac_x]-360)
			wind = wind10m[wy,wx]
			if dist < 1000 and wind > last_wind: #find highest wind maximum within 1000 km of center
				storm_wind = wind
				maxwlo = wlo
				maxwla = wla
				last_wind = wind

		#if the storm is associated with gales
		if storm_wind >= 34:
			print 'this is nameable'
			valid.append(['name', lats[nepac_y], lons[nepac_x], storm_wind, mslp[nepac_y,nepac_x], 'WS', maxwla, maxwlo])
			plotopts = 'go'
		else:
			print 'not nameable, no gales'
			valid.append(['name', lats[nepac_y], lons[nepac_x], storm_wind, mslp[nepac_y,nepac_x], 'EX', maxwla, maxwlo])
			plotopts = 'co'
	if PLOT_DEPTH_FILL == True:
		FILLX, FILLY = m(lons[fillx], lats[filly])
		plt.plot(FILLX, FILLY, plotopts)
#read in most recent data - name them via continuity
f = open('current', 'r')
k = f.readlines()
f.close()
laststorms = []
for line in k:
	l = line.strip('\n').split(',')
	laststorms.append((l[0], float(l[1]), float(l[2]), int(l[3]), int(l[4])))
#output new data
f = open('current', 'w')
f3 = open('all', 'a')
f4 = open('wind', 'a')
for storm in valid:
	#match storms
	i = 0
	named = False
	for laststorm in laststorms:
		dist = calc_dist(laststorm[1], laststorm[2], storm[1], storm[2]-360) #convert current valid lon to negative (west)
		if dist < 648: #within 648 km (moving at less than ~30 m/s)
			storm[0] = laststorm[0] #assume continuity of this storm and keep name
			named = True
			laststorms.pop(i) #remove from list, storms cannot split
			break
		i += 1

	#name the storm if it didn't exist before
	if named == False and storm[5] == 'EX':
		pass #we're still tracking this low but it doesn't meet windstorm criteria
	elif named == False and storm[2] >= 180:
		f2 = open('names_ne', 'r')
		names_ne = f2.readlines()
		f2.close()
		storm[0] = names_ne.pop(0).strip('\n')
		f2 = open('names_ne', 'w')
		f2.writelines(names_ne)
		f2.close()
	elif named == False:
		f2 = open('names_nw', 'r')
		names_nw = f2.readlines()
		f2.close()
		storm[0] = names_nw.pop(0).strip('\n')
		f2 = open('names_nw', 'w')
		f2.writelines(names_nw)
		f2.close()

	#write data to current file
	if storm[0] != 'name': #if actually named
		s = '%s,%s,%s,%d,%d\n' % (storm[0], storm[1], storm[2]-360, storm[3], storm[4])
		f.writelines(s)
		s2 = 'BB,%s,%s%s,,,%s,%s,%d,%d,%s\n' % (storm[0], YYYYMMDD, HH, fmt_lat(float(storm[1])), fmt_lon(float(storm[2])), storm[3], storm[4], storm[5])
		f3.writelines(s2)
		s3 = '%s,%s%s,%s,%s\n' % (storm[0], YYYYMMDD, HH, storm[6], storm[7]-360)
		f4.writelines(s3)

	#plot the storm name
	if storm[5] == 'WS': #if actually a windstorm
		STORMX, STORMY = m(storm[2], storm[1]) #lon, lat
		print STORMX
		if STORMX > 0:
			an = '%s\n%d kts, %d hPa' % (storm[0], storm[3], storm[4])
			plt.annotate(an, xy=(STORMX, STORMY), xytext=(STORMX+50000, STORMY+50000))
f.close()
f3.close()

plt.colorbar(fraction=0.025, pad=0.01)
print 'saving figure...'
plt.savefig("/home/kalassak/nepac/1000-500_thickness_mslp_250_wind_%s%s.png" % (YYYYMMDD, HH), bbox_inches='tight', pad_inches=0, dpi=100)
print 'figure saved!'
plt.close()

#SECOND PLOT --- plot mslp & 10m winds
fig = plt.figure(figsize=(18.6, 10.5))
ax = fig.add_axes((0,0,1,1))

m.drawcoastlines(linewidth=1.2)
m.drawparallels(parallels)
m.drawmeridians(meridians)

# plot contours
print "plotting psfc"
con_mslp = m.contour(xx, yy, mslp, levels=preslevs, linewidths=0.75, colors='k')
plt.clabel(con_mslp, con_mslp.levels, fmt='%d')
print "plotting sfc winds"
m.contourf(xx, yy, wind10m, levels=ulevs, colors=ucolors, zorder=0)
m.barbs(xx[::WIND_SAMPLE_RATE,::WIND_SAMPLE_RATE], yy[::WIND_SAMPLE_RATE,::WIND_SAMPLE_RATE], u10m[::WIND_SAMPLE_RATE,::WIND_SAMPLE_RATE]*1.94384449, v10m[::WIND_SAMPLE_RATE,::WIND_SAMPLE_RATE]*1.94384449) #make sure to convert m/s to kt

#plot cyclones
plt.plot(MXX, MXY, 'ro')
#plt.plot(MWX, MWY, 'yo')
for storm in valid:
	if storm[5] == 'WS': #if actually a windstorm
		STORMX, STORMY = m(storm[2], storm[1]) #lon, lat
		if STORMX+50000 > 0:	
			an = '%s\n%d kts, %d hPa' % (storm[0], storm[3], storm[4])
			t = plt.text(STORMX+500000, STORMY+500000, an, fontsize=12, zorder=10)
			#t = plt.text(an, xy=(STORMX, STORMY), xytext=(STORMX+50000, STORMY+50000), fontsize=14, zorder=10)
			t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='red'))
	

plt.colorbar(fraction=0.025, pad=0.01)
print 'saving figure 2...'
plt.savefig("/home/kalassak/nepac/mslp_10m_wind_%s%s.png" % (YYYYMMDD, HH), bbox_inches='tight', pad_inches=0, dpi=100)
print 'figure 2 saved!'
plt.close()

print nepac_low_lats
print nepac_low_lons

#m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
#m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
