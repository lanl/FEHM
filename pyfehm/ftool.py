"""Various tools for use with PyFEHM."""

"""
Copyright 2013.
Los Alamos National Security, LLC. 
This material was produced under U.S. Government contract DE-AC52-06NA25396 for 
Los Alamos National Laboratory (LANL), which is operated by Los Alamos National 
Security, LLC for the U.S. Department of Energy. The U.S. Government has rights 
to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES 
ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce 
derivative works, such modified software should be clearly marked, so as not to 
confuse it with the version available from LANL.

Additionally, this library is free software; you can redistribute it and/or modify 
it under the terms of the GNU Lesser General Public License as published by the 
Free Software Foundation; either version 2.1 of the License, or (at your option) 
any later version. Accordingly, this library is distributed in the hope that it 
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General 
Public License for more details.
"""

import numpy as np
import os,math,platform,string,difflib
WINDOWS = platform.system()=='Windows'
if WINDOWS: slash = '\\'
else: slash = '/'

from fdflt import*

dflt = fdflt()

# dictionary of unit name, number to multiply by to yield SI measure
units = dict([
# flow
['gpm',1/15.85],
['tph',1000./3600.],

# pressure
['pa',1.e-6],
['kpa',1.e-3],
['mpa',1.],
['psi',1./145.05],
['bar',0.1],

# volume
['gal',3.78541],
['gallon',3.78541],
['gallons',3.78541],

# length
['ft',0.3048],
['feet',0.3048],
['mile',1609.34],
['miles',1609.34],

# time
['day',3600.*24],
['days',3600.*24],
['year',3600.*24*365.25],
['yr',3600.*24*365.25],
['years',3600.*24*365.25],

# mass
['ton',1000.],
['tons',1000.],
['Mt',1.e9],

# temperature
['F',[5./9,32.]],
['farenheit',[5./9,32.]],
]
)
#-----------------------------------------------------------------------------------------------------
#-------------------------------------- USEFUL TOOLS, FUNCTION ---------------------------------------
#-----------------------------------------------------------------------------------------------------
def UTM_to_latlong(x,y,zone,hemisphere=1):
	'''Return latitude and longitude corresponding to supplied UTM (Universal Transverse Mercator) coordinates.
	
	:param x: Easting inside of zone.
	:type x: fl64
	:param y: Northing inside of zone.
	:type y: fl64
	:param zone: UTM zone (integer between 1 and 60).
	:type zone: int
	:param hemisphere: Integer denoting hemisphere, 1 = northern, -1 = southern.
	:type hemisphere: bool
	
	:returns: Two column list of latitudes and longitudes.
	'''

	x = np.array(x); y = np.array(y); zone = np.array(zone); hemisphere = np.array(hemisphere)

	x=500e3-x
	k0 = 0.9996
	M = y/k0
	e = .081819191
	ep2 = 0.006739497
	a = 6378137.
	mu = M/(a*(1-e**2/4-3*e**4/64-5*e**6/256))
#	print mu
	e1 = (1-np.sqrt(1-e**2))/(1+np.sqrt(1-e**2))
#	print e1
	J1 = 3*e1/2.-27*e1**3/32
#	print J1
	J2 = 21*e1**2/16-55*e1**4/32
#	print J2
	J3 = 151*e1**3/96
#	print J3
	J4 = 1097*e1**4/512
#	print J4
	fp = mu+J1*np.sin(2*mu)+J2*np.sin(4*mu)+J3*np.sin(6*mu)+J4*np.sin(8*mu)
#	print fp
	C1=ep2*np.cos(fp)**2
#	print C1
	T1 = np.tan(fp)**2
#	print T1
	R1 = a*(1-e**2)/np.sqrt((1-e**2*np.sin(fp)**2)**3)
#	print R1
	N1 = a/np.sqrt(1-e**2*np.sin(fp)**2)
#	print N1
	D = x/(N1*k0)
#	print D
	Q1 = N1*np.tan(fp)/R1
#	print Q1
	Q2 = D**2/2
#	print Q2
	Q3 = (5+3*T1+10*C1-4*C1**2-9*ep2)*D**4/24
#	print Q3
	Q4 = (61+90*T1+298*C1+45*T1**2-3*C1**2-252*ep2)*D**6/720
#	print Q4
	lat = fp-Q1*(Q2-Q3+Q4)
#	print lat
	Q5 = D
#	print Q5
	Q6 = (1+2*T1+C1)*D**3/6
#	print Q6
	Q7 = (5-2*C1+28*T1-3*C1**2+8*ep2+24*T1**2)*D**5/120
#	print Q7
	long0 = (zone*6-183)/180.*math.pi
#	print long0
	long = long0-(Q5-Q6+Q7)/np.cos(fp)
#	print long
	
	return [hemisphere*lat/math.pi*180.,long/math.pi*180.]
def latlong_to_UTM(lat,long):
	'''Return UTM easting, northing, zone number and hemisphere corresponding to supplied latitude and longitude.
	
	:param lat: Latitude, in deceimal degrees.
	:type lat: fl64
	:param long: Longitude, in decimal degrees.
	:type long: fl64
	
	:returns: Four column list of UTM easting, northing, zone and hemisphere (1 = northern, -1 = southern).
	'''
	
	k0=0.9996
	e=0.081819191
	ep2=0.006739497
	a = 6378137.
	A0=6367449.146
	B0=16038.42955
	C0=16.83261333
	D0=0.021984404
	E0=0.000312705
	
	UTM_zone=31+np.floor(long/6.)
#	print UTM_zone
	UTM_zone_CM = 6.*UTM_zone-183.
#	print UTM_zone_CM
	
	P=(long-UTM_zone_CM)*math.pi/180.
#	print P
	Q=lat*math.pi/180.
#	print Q
	R=long*math.pi/180.
#	print R
	S=a*(1.-e*e)/(np.sqrt(1.-((e*np.sin(Q))**2))**3)
#	print S
	T=a/(np.sqrt(1-(e*np.sin(Q))**2))
#	print T
	V=A0*Q - B0*np.sin(2*Q) + C0*np.sin(4*Q) - D0*np.sin(6*Q) + E0*np.sin(8*Q)
#	print V
	X=V*k0
#	print X
	Y=T*np.sin(Q)*np.cos(Q)/2
#	print Y
	Z=((T*np.sin(Q)*np.cos(Q)**3)/24)*(5-np.tan(Q)**2+9*ep2*np.cos(Q)**2+4*ep2**2*np.cos(Q)**4)*k0
#	print Z
	AA=T*np.cos(Q)*k0
#	print AA
	AB=(np.cos(Q))**3*(T/6)*(1-np.tan(Q)**2+ep2*np.cos(Q)**2)*k0
#	print AB
	AC=(P**6*T*np.sin(Q)*np.cos(Q)**5/720.)*(61.-58.*np.tan(Q)**2+np.tan(Q)**4+270*ep2*np.cos(Q)**2-330*ep2*np.sin(Q)**2)*k0
#	print AC
	
	y=(X+Y*P*P+Z*P**4)
#	print 
	if y<0: y += 10000000
	
	x=500000+(AA*P+AB*P**3)
	
	if lat<0: hemisphere = -1
	else: hemisphere = 1
	
	return [x,y, int(UTM_zone), hemisphere]
def powspace(x0,x1,N=10,power=1):
	'''Returns a sequence of numbers spaced according to the power law (x1-x0)**(1-power)*linspace(0,(x1-x0),N)**base + x0
	
	:param x0: First number in sequence.
	:type x0: fl64
	:param x1: Last number in sequence.
	:type x1: fl64
	:param N: Total items in sequence.
	:type N: int
	:param power: Index of power law. If negative, spacing order will be reversed from "big-to-small".
	:type power: fl64
	'''
	if power>0:
		return (x1-x0)**(1-power)*np.linspace(0,x1-x0,N)**power+x0
	elif power<0:
		return np.sort(x1-((x1-x0)**(1-abs(power))*np.linspace(0,x1-x0,N)**abs(power)))
def SI(quantity,unit=None):
	'''Returns an SI measurement corresponding to a supplied non standard unit, e.g., '31.6 gpm'. Output reverts to FEHM format, e.g., pressures and stresses quoted in MPa.
	
	:param quanityt: Number with unit. If a float or list of floats provided, then unit must 
	:type quantity: str (fl64)
	:param unit: Flag indicating unit of float supplied in quantity.
	:type unit: str
	
	:returns: Quantity in SI units.
	'''
	
	if isinstance(quantity,list):
		return [SI(q,unit) for q in quantity]
	elif isinstance(quantity,np.ndarray):
		return np.array([SI(q,unit) for q in quantity])
	
	if unit == None:	
		#remove any spaces
		quantity = ''.join(quantity.split())
		#separate number and unit
		unit = ''.join([i for i in quantity if not i.isdigit() and not i=='.'])
		quantity = float(''.join([i for i in quantity if i.isdigit() or i=='.']))
		
	try:
		if unit not in ['F','farenheit']:
			return quantity*units[unit]
		else:
			return (quantity-units[unit][1])*units[unit][0]
	
	except KeyError:
		print unit +' not a valid unit for conversion'
		return	
def pyfehm_print(s):
	if not dflt.silent: print s
#-----------------------------------------------------------------------------------------------------
#------------------------------ FUNCTIONS AND CLASSES FOR INTERNAL USE -------------------------------
#-----------------------------------------------------------------------------------------------------
class fpath(object):
	def __init__(self,filename = None,work_dir = None,parent = None):
		self._filename = filename
		self.absolute_to_file = None		# location where originally read DOES NOT CHANGE
		self.absolute_to_workdir = None	# working directory CAN CHANGE
		self.parent = parent
	def update(self, wd):
		'''called when work_dir is updated'''		
		if wd == None: 
			self.absolute_to_workdir = None
			return
		
		if WINDOWS: wd = wd.replace('/','\\')
		else: wd = wd.replace('\\','/')
		
		absolute = False
		if WINDOWS and wd[1]==':': absolute = True
		if not WINDOWS and wd[0]=='/': absolute = True
		if absolute:
			self.absolute_to_workdir = wd
		else:
			self.absolute_to_workdir = os.getcwd()+slash+wd	
	def _get_filename(self): return self._filename
	def _set_filename(self,value): 
		# ensure path specification consistent with OS
		if WINDOWS: value = value.replace('/','\\')
		else: value = value.replace('\\','/')
		# check if any slashes exist
		if slash in value:
			self._filename = value.split(slash)[-1]
		else: 
			self._filename = value
			self.absolute_to_file = os.getcwd()
			return
		# check if absoulte or relative specification
		path = value.split(slash)[:-1]
		path = string.join(path,slash)
		
		absolute = False
		if WINDOWS and path[1]==':': absolute = True
		if not WINDOWS and path[0]=='/': absolute = True
		
		if absolute:
			self.absolute_to_file = path
		else:
			self.absolute_to_file = os.getcwd()+slash+path
	filename = property(_get_filename, _set_filename) #: (**)
	def _get_full_path(self): 
		return self.absolute_to_file+slash+self.filename
	full_path = property(_get_full_path) #: (**)
def dict_key_check(dict,keys,dict_name):
	'''Return False if dict contains only the supplied keys and no extras.
	'''
	returnFlag = False
	ws = 'Key error in '+dict_name+'.\n'
	for k in dict.keys():
		if k not in keys:
			ws += 'No such key \''+k+'\''
			if len(k)>2:
				matches = difflib.get_close_matches(k,keys)
			else:
				matches = difflib.get_close_matches(k,keys,cutoff = 0.5)
			print k
			print keys
			if len(matches)>0:
				ws+=', did you mean?\n'
				for match in matches:
					ws+='- '+match+'\n'
			else:
				ws+='.\n'
			returnFlag = True
	if returnFlag: print ws
	return returnFlag
def os_path(path):
	if WINDOWS: path = path.replace('/','\\')
	else: path = path.replace('\\','/')
	return path
def _title_string(s,n): 						#prepends headers to sections of FEHM input file
	if not n: return
	ws = '# '
	pad = int(np.floor((n - len(s) - 2)/2))
	for i in range(pad): ws+='-'
	ws+=s
	for i in range(pad): ws+='-'
	ws+='\n'
	return ws
def _zone_ind(indStr): return abs(int(indStr))-(int(indStr)+abs(int(indStr)))/2
def flatten(l):
	'''Takes a nested list and returns a flattened list (upon list conversion of output).
	:param l: Nested list.
	:type l: lst
	'''
	import collections
	for el in l:
		if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
			for sub in flatten(el):
				yield sub
		else:
			yield el
def sub_cubes(cube): 			# returns sub-cubes formed by sub-dividing the given cube into eight equal
	c = 0.5*(cube[0]+cube[1]) 		# cube centre
	
	dp = c - cube[0]
	p1,p8 = cube[0], c
	p2 = p1+np.array([dp[0],0,0])
	p3 = p1+np.array([0,dp[1],0])
	p4 = p1+np.array([dp[0],dp[1],0])
	p5 = p1+np.array([0,0,dp[2]])
	p6 = p1+np.array([dp[0],0,dp[2]])
	p7 = p1+np.array([0,dp[1],dp[2]])
	
	#c0=[cube[0], c]
	#
	#c1=[np.array([c[0],cube[0][1],cube[0][2]]),np.array([cube[1][0],c[1],c[2]])]
	#
	#c2=[np.array([cube[0][0],cube[0][1],c[2]]),np.array([c[0],c[1],cube[1][2]])]
	#
	#c3=[np.array([c[0],cube[0][1],c[2]]),np.array([cube[1][0],c[1],cube[1][2]])]
	#
	#c4=[np.array([cube[0][0],c[1],cube[0][2]]),np.array([c[0],cube[1][1],c[2]])]
	#
	#c5=[np.array([c[0],c[1],cube[0][2]]),np.array([cube[1][0],cube[1][1],c[2]])]
	#
	#c6=[np.array([cube[0][0],c[1],c[2]]),np.array([c[0],cube[1][1],cube[1][2]])]
	#
	#c7=[c, cube[1]]
	
	c0=[p1,p1+dp]
	c1=[p2,p2+dp]
	c2=[p3,p3+dp]
	c3=[p4,p4+dp]
	c4=[p5,p5+dp]
	c5=[p6,p6+dp]
	c6=[p7,p7+dp]
	c7=[p8,p8+dp]
	
	return [c0,c1,c2,c3,c4,c5,c6,c7]
def in_cube(pos,cube):
	"""Tests if the 3-D point lies in an axis-aligned cube, defined as a two-element list of arrays
	[bottom left front, top right back]."""
	return all([cube[0][i]<=pos[i]<cube[1][i] for i in xrange(3)])
def cubes_intersect(cube1,cube2):
	"""Returns True if two cubes intersect."""
	return all([(cube1[1][i]>=cube2[0][i]) and (cube2[1][i]>=cube1[0][i]) for i in xrange(2)])
def save_name(save='',variable='',time=0., node=0): 		# returns file name and extension for saving
	pdf = False
	if save:
		save = save.split('.')
		if len(save)==1: 
			print 'No extension specified, default to .png'
			ext = 'png'
		elif len(save)>2: print 'Too many dots!'; return
		else:
			if save[1] in ['png','eps','pdf']:
				ext = save[1]
			else: print 'Unrecognized extension'; return
		if ext == 'pdf': ext = 'eps'; pdf = True
		save_fname=save[0]+'.'+ext
	else: 
		from glob import glob
		ext = 'png'
		if time:
			varStr = variable+'_time'+str(time)
		elif node:
			varStr = variable+'_node'+str(node)
		else:
			varStr = variable
		files=glob('pyfehm_sliceplot_'+varStr+'_*.png')
		if not files: ind = 1
		else:
			inds = []
			for file in files:
				file = file.split('pyfehm_sliceplot_'+varStr+'_')
				inds.append(int(file[1].split('.png')[0]))
			ind = np.max(inds)+1	
		save_fname='pyfehm_sliceplot_'+varStr+'_'+str(ind)+'.png'
	return ext, save_fname, pdf	
def make_directory(fname):
	fname = fname.split('\\')
	fname = fname[:-1]
	dirname = ''
	for f in fname: dirname += f
	if not os.path.isdir(dirname): os.system('mkdir '+dirname)
	