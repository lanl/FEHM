
print 'Testing imports'
from fdata import*

# test reading of internode fluxes
ndflx = fnodeflux('run.internode_fluxes.out')

# test pyfehm capabilities
path = fpath()
path.filename = 'test.py'
print path.filename
print path.absolute_to_file
print path.full_path
print ''
path.filename = '/test/test.py'
print path.filename
print path.absolute_to_file
print path.full_path
print ''
path.filename = 'c:\\test\\test.py'
print path.filename
print path.absolute_to_file
print path.full_path
print ''

def test_grid(geo):

	geo = fgrid()
	geo.read('mesh_avs.avs')

	x = np.linspace(0,10,11)
	geo.make('pyfehm_unittest_GRID.inp',x=x,y=x,z=x)
	if not os.path.isfile('pyfehm_unittest_GRID.inp'): print 'where is file?'; return False
	xs = np.sort(np.unique([nd.position[0] for nd in geo.nodelist]))
	ys = np.sort(np.unique([nd.position[1] for nd in geo.nodelist]))
	zs = np.sort(np.unique([nd.position[2] for nd in geo.nodelist]))
	
	if len(geo.nodelist) != 1331: print 'not enough nodes'; return False
	
	for x1,x2 in zip(x,xs):
		if x1 != x2: print 'x node in wrong place'; return False
		
	for x1,x2 in zip(x,ys):
		if x1 != x2: print 'y node in wrong place'; return False
		
	for x1,x2 in zip(x,zs):
		if x1 != x2: print 'z node in wrong place'; return False
		
	geo.read('pyfehm_unittest_GRID.inp',full_connectivity = True)
	
	xL,yL,zL = len(xs), len(ys), len(zs)
	nd6 = []; nd5 = []; nd4 = []; nd3 = []
	for nd in geo.nodelist:
		if len(nd.connected_nodes) == 6: nd6.append(nd)
		if len(nd.connected_nodes) == 5: nd5.append(nd)
		if len(nd.connected_nodes) == 4: nd4.append(nd)
		if len(nd.connected_nodes) == 3: nd3.append(nd)
	
	if len(nd3) != 8: print 'wrong number corner nodes'; return False
	if len(nd4) != ((xL-2)+(yL-2)+(zL-2))*4: print 'wrong number edge nodes'; return False
	if len(nd5) != ((xL-2)*(yL-2)+(yL-2)*(zL-2)+(xL-2)*(zL-2))*2: print 'wrong number face nodes'; return False
	if len(nd6) != (xL-2)*(yL-2)*(zL-2): print 'wrong number internal nodes'; return False
	
	x,y,z = geo.node_nearest_point([5.8,5.8,5.8]).position
	if x!=6. or y!=6. or z!=6.: print 'node nearest point failed'; return False
	geo.what
	geo.plot(save='pyfehm_grid_unit.png',color='r',cutaway = 'middle', angle=[25,45])
	os.system('del pyfehm_grid_unit.png')
	
	return True

def test_input():
	# create from scratch
	dat = fdata()
	x = np.linspace(0,10,11)
	dat.grid.make('pyfehm_unittest_GRID.inp',x=x,y=x,z=x)
	#dat.grid.lagrit_stor(overwrite=True)
	dat.grid.what
	# zone creation
	# 1. rect
	zn = fzone(index=1,name='zone1')
	zn.rect([1,1,1],[2,2,2])
	dat.add(zn)
	dat.zone[1].what
	if len(dat.zone['zone1'].nodelist) != 8: print 'rect zone wrong number of nodes'; return False
	dat.delete(dat.zonelist[-1])
	# 2. nodelist
	zn = fzone(1,nodelist=5)
	dat.add(zn)
	dat.delete(dat.zonelist[-1])
	zn = fzone(1,nodelist=[5,8])
	dat.add(zn)
	dat.delete(dat.zonelist[-1])
	nd1,nd2 = dat.grid.node[12],dat.grid.node[13]
	nd1.what
	dat.add(fzone(1,nodelist=nd1))
	dat.add(fzone(1,nodelist=nd1))
	dat.add(fzone(1,nodelist=nd1),overwrite=True)
	dat.delete(dat.zonelist[-1])
	dat.add(fzone(1,nodelist=[nd1,nd2]))
	dat.delete(dat.zonelist[-1])
	# 3. new_zone()
	dat.new_zone(1,'hi',rect=[[2.,2.,2.],[3.,3.,3.]],permeability = 1.e-14, density = 2500., porosity = 0.1, specific_heat = 800.,
		Pi = 1., Ti = 50., youngs_modulus = 10.e3, poissons_ratio = 1., thermal_expansion = 1., pressure_coupling = 1.)
	if dat.zone[1].permeability != 1.e-14: print 'perm error'; return False
	if dat.zone[1].density != 2500.: print 'density error'; return False
	if dat.zone[1].specific_heat != 800.: print 'spec heat error'; return False
	if dat.zone[1].porosity != 0.1: print 'por error'; return False
	if dat.zone[1].Pi != 1.: print 'Pi error'; return False
	if dat.zone[1].Ti != 50.: print 'Ti error'; return False
	if dat.zone[1].youngs_modulus != 10.e3: print 'youngs mod error'; return False
	if dat.zone[1].poissons_ratio != 1.: print 'poissons ratio error'; return False
	if dat.zone[1].thermal_expansion != 1.: print 'alpha error'; return False
	if dat.zone[1].pressure_coupling != 1.: print 'biot error'; return False
	dat.new_zone(1,'hi',nodelist = 5)
	dat.new_zone(1,'hi',nodelist = 5,overwrite=True)
	dat.new_zone(1,'hi',nodelist = [nd1,nd2],permeability =[1.e-14,2.e-14,1.e-14], density = 2500., Pi = 1., youngs_modulus = 1., thermal_expansion = 1.,
		overwrite=True)
	#dat.delete(dat.zonelist[-1])
	
	# macro creation
	dat.zone[0].permeability = 1.e-14
	dat.add(fmacro('cond',zone=0,param=(('cond_x',1.),('cond_y',1.),('cond_z',1.))))
	dat.add(fmacro('rock',zone=0,param=(('density',2500.),('porosity',.1),('specific_heat',800.))))
	dat.add(fmacro('elastic',zone=0,param=(('youngs_modulus',1.),('poissons_ratio',1.))))
	dat.add(fmacro('biot',zone=0,param=(('pressure_coupling',1.),('thermal_expansion',1.))))
	dat.add(fmacro('flow',zone=0,param=(('rate',1.),('energy',1.),('impedance',1.))))
	dat.add(fmacro('grad',zone=0,param=(('reference_coord',0),('direction',3),
		('variable',1),('reference_value',1),('gradient',-9.81*1e3/1e6))))
	dat.add(fmacro('grad',zone=0,param=(('reference_coord',0),('direction',3),
		('variable',2),('reference_value',1),('gradient',-0.04))))
	dat.add(fmacro('co2frac',param=(('water_rich_sat',1.),('co2_rich_sat',0.),('co2_mass_frac',0.),
		('init_salt_conc',0.),('override_flag',1))))
	dat.add(fmacro('co2pres',param=(('pressure',1),('temperature',1),('phase',1))))
		
	dat.add(fmacro('co2flow',zone=1,param=(('rate',1),('energy',1),('impedance',1),('bc_flag',1))))
	dat.add(fmacro('co2flow',zone=1,param=(('rate',1),('energy',1),('impedance',1),('bc_flag',1))))
	dat.add(fmacro('co2flow',zone='hi',param=(('rate',1),('energy',1),('impedance',1),('bc_flag',1))))
	dat.add(fmacro('co2flow',zone=1,param=(('rate',1),('energy',1),('impedance',1),('bc_flag',1))),overwrite=True)
	
	dat.carb.on(iprtype = 3)
	dat.strs.on()
	dat.nfinv = True
	dat.nobr = True
	dat.adif = 200.
	
	#output stuff
	dat.hist.variables.append(['flow','temperature'])
	dat.hist.format = 'surf'
	dat.hist.nodelist.append(4)
	dat.hist.nodelist.append(nd1)
	
	dat.cont.variables.append(['xyz','pressure'])
	dat.cont.format = 'surf'

	dat.flxo.append((2,5))
	dat.flxo.append((dat.grid.node[8],10))
	dat.flxo.append((dat.grid.node[8],dat.grid.node[15]))
	
	dat.add(fmacro('stressboun',zone=1,subtype='fixed',param=(('direction',1),('value',0))))
	dat.add(fboun(type='ti_linear',zone=1,times=[1,2,3],variable=[['ft']+[2,3,4,],['pw']+[4,5,6]]))	

	# rlp stuff
	rlp = fmodel('rlp',index=17,param=[.05,1,1,0,1,1,0,0,1,1,1,0,1,0])
	dat.add(rlp)
	newrlpm = frlpm_table(group=2,zone=dat.zone[0])
	
	s = np.linspace(0,1,11)
	
	newrlpm.saturation = s
	
	Smin = 0.2
	rlp_w = ((s-Smin)/(1-Smin))**3.1
	rlp_w[np.where(s<Smin)[0]] = 0.
	
	rlp_c = 0.8*((1-s)/(1-Smin))**3.1
	rlp_c[np.where(s<Smin)[0]] = 0.8
	
	newrlpm.phase1 = ['water',rlp_w]
	newrlpm.phase2 = ['co2_liquid',rlp_c]
	newrlpm.capillary = np.zeros((1,len(s)))[0]
	
	dat.add(newrlpm)
	
	# other models
	rlp = fmodel('vcon',index=1,param=[1.,1.,1.])
	dat.add(rlp)
	rlp = fmodel('ppor',index=1,zonelist = 1, param=(('compressibility',1)))
	dat.add(rlp)
	pm = fmodel('permmodel',index=25)
	pm.zonelist=dat.zone[1]
	pm.param['shear_frac_tough'] = 1.
	pm.param['static_frict_coef'] = 1.
	pm.param['dynamic_frict_coef'] = 1.
	pm.param['frac_num'] = 1.
	pm.param['onset_disp'] = 1.
	pm.param['disp_interval'] = 1.
	pm.param['max_perm_change'] = 1.
	pm.param['frac_cohesion'] = 1.
	dat.add(pm)
	
	dat.zone[0].what
	
	print 'test writing'
	dat.write('pyfehm_unittest_INPUT.dat')
	
	print 'test reading'
	dat0 = fdata('pyfehm_unittest_INPUT.dat','pyfehm_unittest_GRID.inp')
	
	# test run
	dat.tf = 1.
	dat.dtn = 1.
	dat.dtx = 1.5
	dat.dti = 1.
	dat.dtmax = 2.
	dat.dtmin = 0.1
	dat.change_timestepping(10.)
	dat.change_timestepping(20.)
	dat.change_timestepping(30.,new_dti = .1)
	dat.change_timestepping(40.,new_dtmax = 10.)
	dat.change_timestepping(50.,new_dtx = 1.1)
	dat.change_timestepping(60.)
	dat.cont.format = 'avsx'
	dat.files.rsto='pyfehm_unittest_OUT.rsto'
	dat.verbose = True
	
	print 'test running'
	#dat.permlist[0].param['Kx'] = 1.e-5
	dat.run(verbose = False)
	
	dat.verbose = True
	
	dat.output_times = [10,20,30]
	dat.output_times = np.linspace(1,10,20)
	
	print 'test incon reading'
	dat.incon.read('pyfehm_unittest_OUT.rsto')
	
	dat0 = fdata('pyfehm_unittest_INPUT.dat','pyfehm_unittest_GRID.inp','pyfehm_unittest_OUT.rsto')
	dat.incon.stressgrad(1.,1.,1.)
	dat.incon.critical_stress()
	
	print 'test work dir change'
	dat.work_dir = 'test\\test'
	dat.tf = 5.
	
	print 'test running again'
	dat.run(autorestart = 2,use_paths = True)
	os.system('rmdir test /s /q')
	dat.run()
	os.system('rmdir test /s /q')
	dat.work_dir = None
	dat.run()
	
	return True
#################### CONSTRUCTORS #########################
print 'Testing grid'
geo = fgrid()
if not test_grid(geo): 
	print 'ERROR: grid constructor'
	
dat = fdata()
dat = fdata(work_dir = 'test/test')
x = np.linspace(0,10,11)
dat.grid.make('pyfehm_unittest_GRID.inp',x=x,y=x,z=x)
dat.grid.lagrit_stor(overwrite=True)

if not os.path.isdir(dat.work_dir): print 'where work dir?'
if not os.path.isfile(dat.work_dir+'\\pyfehm_unittest_GRID.inp'): print 'where grid in work dir?'

print 'Testing input'
os.system('rmdir test /s /q')
if not test_input():
	print 'ERROR: input files constructor'
os.system('del pyfehm_unittest_INPUT*.*')
os.system('del pyfehm_unittest_GRID*.*')
os.system('del pyfehm_unittest_OUT*.*')

###################
print 'No errors!'

