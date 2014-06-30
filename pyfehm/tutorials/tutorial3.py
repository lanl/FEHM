# 7.3 PyFEHM tutorial 3
# Cube with fixed pressures/temperatures - batch submission
from fdata import*
from fpost import*

# 7.3.2 Setting up the simulations
import multiprocessing

models = range(1,13)
processors = 5

params = [
#[kmid,Tbase], 		# legend
[1.e-20,80], 			# 1 
[1.e-20,120], 			# 2 
[1.e-16,80], 			# 3 
[1.e-16,120], 			# 4 

[1.e-20,80], 			# 1 
[1.e-20,120], 			# 2 
[1.e-16,80], 			# 3 
[1.e-16,120], 			# 4 

[1.e-20,80], 			# 1 
[1.e-20,120], 			# 2 
[1.e-16,80], 			# 3 
[1.e-16,120], 			# 4 
]

def execute(j):		
	simulation(j,params[j-1])

if __name__ == '__main__':	
	p = multiprocessing.Pool(processes = processors)
	p.map(execute,models)

def simulation(j,param):

	kmid,Tbase = param
	
	root = 'tut3_'+str(j)
	dat = fdata(work_dir = root)

	# 7.1.2 Grid generation
	x = np.linspace(0,10,11)
	dat.grid.make(root+'_GRID.inp',x=x,y=x,z=x)
	dat.grid.plot(root+'_GRID.png',color='r',angle=[45,45])

	# 7.1.3 Zone creation
	zn = fzone(index=1,name='lower') 			
	zn.rect([-0.,-0.,-0.1],[10.1,10.1,3.1])		
	dat.add(zn)									
	dat.zone[1].what 							

	dat.new_zone(2,'middle',rect=[[-0.,-0.,3.1],[10.1,10.1,6.1]]) 	
	dat.new_zone(3,'upper',rect=[[-0.,-0.,6.1],[10.1,10.1,10.1]])	
	
	# 7.1.4 Adding macros
	rm = fmacro('rock',param=(('density',2500),('specific_heat',1000),('porosity',0.1)))
	dat.add(rm)

	pm = fmacro('perm',zone=1,param=(('kx',1.e-15),('ky',1.e-15),('kz',1.e-16))) 	
	dat.add(pm) 																	

	dat.zone[2].permeability = kmid 			
	dat.zone['upper'].permeability = [1.e-14,1.e-14,1.e-15] 	

	# 7.1.5 Initial and boundary conditions
	pres = fmacro('pres',param=(('pressure',5.),('temperature',60.),('saturation',1)))
	dat.add(pres)

	hflx = fmacro('hflx',zone='ZMAX',param=(('heat_flow',30),('multiplier',1.e10)))
	dat.add(hflx)

	dat.zone['ZMIN'].fix_temperature(Tbase)

	flow = fmacro('flow',zone='YMIN',param=(('rate',6.),('energy',-60.),('impedance',1.e6)))
	dat.add(flow)

	dat.zone['YMAX'].fix_pressure(P=4.,T=60.)

	# 7.1.6 Running the simulation
	dat.cont.variables.append(['xyz','temperature','pressure'])

	dat.tf=10.

	dat.files.root = root
	dat.run(root+'_INPUT.dat',exe='c:\\users\\264485\\fehm\\source\\src\\fehm.exe',files=['outp'])

	# 7.1.7 Visualisation
	c = fcontour(root+'\\*.csv',latest=True)
	c.slice_plot(save=root+'\\Tslice.png',cbar=True,levels=11,slice=['x',5],variable='T',method='linear',title='temperature / degC',
	xlabel='y / m', ylabel = 'z / m')
	c.slice_plot(save=root+'\\Pslice.png',cbar=True,levels=np.linspace(4,6,9),slice=['x',5],variable='P',method='linear',title='pressure / MPa',
	xlabel='y / m', ylabel = 'z / m')
