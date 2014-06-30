# 7.1 PyFEHM tutorial 4
# Dynamic simulation monitoring and kill conditions
# 7.1.1 Getting started

from fdata import*
from fpost import*

root = 'tut4'
dat = fdata(work_dir=root) 					# creates an 'empty' input file

# 7.1.2 Grid generation
x = np.linspace(0,10,11)
dat.grid.make(root+'_GRID.inp',x=x,y=x,z=x)
dat.grid.plot(root+'_GRID.png',color='r',angle=[45,45])

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
dat.cont.variables.append(['xyz','temperature','pressure'])

dat.tf=20.
#dat.time['max_time_TIMS']=10.

dat.files.root = root

########################################### DEFINE A FUNCTION ###########################################
def kill_condition(dat):
	# function must have a single input, and fdata object
	# when the function is called within dat.run(), it will be passed dat as an input
	# this function returns True when the simulation is to be halted
	# during function execution, state information at each node is available
	
	# for example, kill the simulation when the node at [7.,7.,7.] has a temperature less than 50degC
	nd = dat.grid.node_nearest_point([7.,7.,9.])
		
	return nd.T<50.

dat.run(root+'_INPUT.dat',until=kill_condition) # pass the function to the argument 'until'
