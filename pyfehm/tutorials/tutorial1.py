# 7.1 PyFEHM tutorial 1
# Cube with fixed pressures/temperatures
# 7.1.1 Getting started

#import os,sys
#sys.path.append('c:\\python\\pyfehm')

from fdata import*
from fpost import*

root = 'tut1'
dat = fdata(work_dir = root) 					# creates an 'empty' input file

# 7.1.2 Grid generation
x = np.linspace(0,10,11)
dat.grid.make(root+'_GRID.inp',x=x,y=x,z=x)
#dat.grid.plot(root+'_GRID.png',color='r',angle=[45,45])

# 7.1.3 Zone creation
zn = fzone(index=1,name='lower') 			# (a) create the zone object, with index, name
zn.rect([-0.,-0.,-0.1],[10.1,10.1,3.1])		# (b) assign a bounding box for a rectangular zone
dat.add(zn)									# (c) add the zone to the input file
dat.zone[1].what 							# print out some information about the zone

dat.new_zone(2,'middle',rect=[[-0.,-0.,3.1],[10.1,10.1,6.1]]) 	# the .new_zone() method implements steps 
dat.new_zone(3,'upper',rect=[[-0.,-0.,6.1],[10.1,10.1,10.1]])	# (a)-(c) above in one step.
dat.zone['upper'].what 											# print out some information about the zone

# 7.1.4 Adding macros
# note, if we had omitted this macro entirely, PyFEHM would have added in a default with these values
rm = fmacro('rock',param=(('density',2500),('specific_heat',1000),('porosity',0.1)))
dat.add(rm)

pm = fmacro('perm',zone=1,param=(('kx',1.e-15),('ky',1.e-15),('kz',1.e-16))) 	# (a) create a permeability macro, assign to zone 1
dat.add(pm) 																	# (b) add perm macro to input file

kmid = 1.e-20
dat.zone[2].permeability = kmid 			# use material property attributes to do steps (a)-(b) in one step
dat.zone['upper'].permeability = [1.e-14,1.e-14,1.e-15] 	# three item list indicates anisotropic

# 7.1.5 Initial and boundary conditions
pres = fmacro('pres',param=(('pressure',5.),('temperature',60.),('saturation',1)))
dat.add(pres)

hflx = fmacro('hflx',zone='ZMAX',param=(('heat_flow',30),('multiplier',1.e10))) 			# (a) create fixed temperature hflx macro
dat.add(hflx)																				# (b) add macro to input file

dat.zone['ZMIN'].fix_temperature(80) 			# does steps (a)-(b) above for ZMIN

flow = fmacro('flow',zone='YMIN',param=(('rate',6.),('energy',-60.),('impedance',1.e6))) 	# (a) create fixed pressure flow macro
dat.add(flow) 																				# (b) add macro to input file

dat.zone['YMAX'].fix_pressure(P=4.,T=60.) 		# does steps (a)-(b) above for YMAX

# 7.1.6 Running the simulation
dat.cont.variables.append(['xyz','temperature','pressure','grid','mat'])
dat.cont.format = 'surf'
#dat.cont.time_interval = 2.

dat.hist.nodelist.append([100,400,800])
#dat.hist.nodelist.append([100])
dat.hist.variables.append(['pressure','temperature'])

dat.tf=100.
#dat.time['max_time_TIMS']=10.

dat.iter['machine_tolerance_TMCH'] = -0.5e-5

dat.files.root = root
dat.run(root+'_INPUT.dat', diagnostic = True) 		# note, because no executable path is specified, PyFEHM retrieves the executable specified in the default path

# 7.1.7 Visualisation
#c = fcontour(dat.work_dir+'\\*.csv')
#c.slice_plot(save='Tslice.png',cbar=True,levels=11,slice=['x',5],variable='T',method='linear',title='temperature / degC',
#xlabel='y / m', ylabel = 'z / m')
#c.slice_plot(save='Pslice.png',cbar=True,levels=np.linspace(4,6,9),slice=['x',5],variable='P',method='linear',title='pressure / MPa',
#xlabel='y / m', ylabel = 'z / m')

