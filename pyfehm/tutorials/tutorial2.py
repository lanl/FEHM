# 7.2 PyFEHM tutorial 2
# Set up a five spot EGS injection/production simulation, quarter domain
# 7.2.1 First steps
from fdata import*
from fpost import*

dat = fdata(work_dir = 'tutorial2')

# 7.2.2 Grid generation
X0,X1 = 0,1.e3
Z0,Z1 = -1.5e3,-0.5e3

injX,injY = 0.,0.
proX,proY = 300.,300.

base = 2
dx = proX/2.
x = dx**(1-base)*np.linspace(0,dx,8)**base
dx2 = X1 - proX
x2 = dx2**(1-base)*np.linspace(0,dx2,10)**base
X = np.sort(list(x) + list(2*dx-x)[:-1] + list(2*dx+x2)[1:])

Za_base = -1.1e3
Za_top = -0.8e3

Z = list(np.linspace(Z0,Za_base,5)) + list(np.linspace(Za_base,Za_top,11))[1:] + list(np.linspace(Za_top,Z1,5))[1:]

dat.grid.make('quarterGrid.inp', x= X, y = X, z = Z)

dat.grid.plot('quarterGrid.png', angle = [45,45], color = 'b', cutaway = [proX,proY,-1000])

# 7.2.3 Zone creation
#zone_res=fzone(index=10,name='reservoir')
#zone_res.rect()
#dat.add(zone_res)

dat.new_zone(10,name='reservoir',rect = [[X0-0.1,X0-0.1,Za_base+0.1], [X1+0.1,X1+0.1,Za_top-0.1]])

dat.zone['reservoir'].plot('reservoirZone.png',angle = [30,30], equal_axes=False,color = 'g')

#zone_con = fzone(index = 20, name='confining_lower')
#zone_con.rect()
#dat.add(zone_con)
dat.new_zone(20,name='confining_lower',rect = [[X0-0.1,X0-0.1,Z0-0.1], [X1+0.1,X1+0.1,Za_base+0.1]])
dat.new_zone(21,name='confining_upper',rect = [[X0-0.1,X0-0.1,Za_top-0.1], [X1+0.1,X1+0.1,Z1+0.1]])


#zone_con = fzone(index = 21, name='confining_upper')
#zone_con.rect()
#dat.add(zone_con)

# 7.2.4 Material property assignment
perm_res, perm_con = 1.e-14, 1.e-16
rho_res, rho_con = 2300., 2500.
phi_res, phi_con = 0.1, 0.01
cond = 2.5
H = 1.e3

# method 1
perm = fmacro('perm')
perm.zone = 'reservoir'
perm.param['kx'] = perm_res
perm.param['ky'] = perm_res
perm.param['kz'] = perm_res
dat.add(perm)

# method 2
perm=fmacro('perm', zone=21, param=(('kx',perm_con), ('ky',perm_con), ('kz',perm_con)))
dat.add(perm)

# method 3
dat.zone['confining_lower'].permeability = perm_con

dat.add(fmacro('rock', zone=dat.zone['reservoir'], param=(('density',rho_con), ('specific_heat',H), ('porosity',phi_res))))
dat.add(fmacro('rock', zone=dat.zone['confining_lower'], param=(('density',rho_con), ('specific_heat',H), ('porosity',phi_con))))

dat.zone['confining_upper'].density = rho_con
dat.zone['confining_upper'].specific_heat = H
dat.zone['confining_upper'].porosity = phi_con

dat.zone[0].conductivity = cond

# 7.2.5 Injectors and producers
#inj_zone = fzone(index=30, name='injection')
#inj_zone.rect()
#dat.add(inj_zone)

dat.new_zone(30,'injection',rect = [[injX-0.1,injY-0.1,-1000.1],[injX+0.1,injY+0.1,-899.9]])

pro_node=dat.grid.node_nearest_point([proX,proY,-950])
dat.new_zone(40,'production',nodelist = pro_node)
#pro_zone=fzone(index=40, name='production', type='nnum', nodelist=[pro_node])
#dat.add(pro_zone)

flow=fmacro('flow', zone = 'production', param=(('rate',6), ('energy',30), ('impedance',1.)))
dat.add(flow)

injRate = 4.
injTemp = 60.

boun = fboun(zone=['injection'],times=[0,365.],variable=[['dsw',-injRate,-injRate],['ft',injTemp,injTemp]])
dat.add(boun)

# 7.2.6 Initial conditions
dat.add(fmacro('grad',zone=0,param=(('reference_coord',0.),('direction',3),
('variable',2),('reference_value',25.),('gradient',-0.06))))

dat.add(fmacro('grad',zone=0,param=(('reference_coord',0.),('direction',3),
('variable',1),('reference_value',0.1),('gradient',-9.81*1e3/1e6))))

#upper_zone = fzone(index=50,name='zmax')
#upper_zone.rect([X0-0.1,X0-0.1,Z1-0.1],[X1+0.1,X1+0.1,Z1+0.1])
#dat.add(upper_zone)

#upper_flow = fmacro('flow',zone='ZMAX',param=(('rate',0.1+Z1*-9.81*1e3/1e6),('energy',-(25.+Z1*-0.06)),('impedance',100.)))
#dat.add(upper_flow)
dat.zone['ZMAX'].fix_pressure(P=0.1+Z1*-9.81*1e3/1e6,T=25.+Z1*-0.06)

# 7.2.7 Setting up stresses
dat.dtn = 1
dat.files.rsto='EGS_INCON.ini'
dat.run('EGS_flow_INPUT.dat')

dat.incon.read('tutorial2\\EGS_INCON.ini')

dat.incon.stressgrad(xgrad=0.6,ygrad=0.8,zgrad=2500*abs(Z1)*9.81/1e6,calculate_vertical=True,vertical_fraction=True)

dat.strs.on()
E_res,E_con = 2e3,2e4
nu_res,nu_con = 0.15,0.35
#dat.add(fmacro('elastic',zone='reservoir',param=(('youngs_modulus',E_res),('poissons_ratio',nu_res))))
#dat.add(fmacro('elastic',zone='confining_lower',param=(('youngs_modulus',E_con),('poissons_ratio',nu_con))))
#dat.add(fmacro('elastic',zone='confining_upper',param=(('youngs_modulus',E_con),('poissons_ratio',nu_con))))
#dat.add(fmacro('biot',param=(('thermal_expansion',3e-5),('pressure_coupling',1.))))

dat.zone['reservoir'].youngs_modulus = E_res
dat.zone['confining_lower'].youngs_modulus = E_con
dat.zone['confining_upper'].youngs_modulus = E_con

dat.zone['reservoir'].poissons_ratio = nu_res
dat.zone['confining_lower'].poissons_ratio = nu_con
dat.zone['confining_upper'].poissons_ratio = nu_con

dat.zone[0].thermal_expansion = 3e-5
dat.zone[0].pressure_coupling = 1.

#zn = fzone(index=60,name='xmin'); zn.rect([X0-0.1,X0-0.1,Z0-0.1],[X0+0.1,X1+0.1,Z1+0.1]); dat.add(zn)
#zn = fzone(index=61,name='ymin'); zn.rect([X0-0.1,X0-0.1,Z0-0.1],[X1+0.1,X0+0.1,Z1+0.1]); dat.add(zn)
#zn = fzone(index=62,name='zmin'); zn.rect([X0-0.1,X0-0.1,Z0-0.1],[X1+0.1,X1+0.1,Z0+0.1]); dat.add(zn)
dat.strs.fem = True
dat.strs.bodyforce=False
dat.add(fmacro('stressboun',zone='XMIN',subtype='fixed',param=(('direction',1),('value',0)),write_one_macro=True))
dat.add(fmacro('stressboun',zone='YMIN',subtype='fixed',param=(('direction',2),('value',0)),write_one_macro=True))
dat.add(fmacro('stressboun',zone='ZMIN',subtype='fixed',param=(('direction',3),('value',0)),write_one_macro=True))

# 7.2.8 Running the model
dat.cont.variables.append(['xyz','pressure','liquid','temperature','stress','displacement','permeability'])
dat.cont.format = 'surf'
dat.cont.timestep_interval = 1000.
dat.cont.time_interval = 365.25/2

dat.hist.variables.append(['temperature','pressure','flow','zfl'])
dat.hist.format = 'surf'
dat.hist.nodelist.append(dat.zone['production'].nodelist[0])
dat.hist.zoneflux.append(dat.zone['injection'])

dat.dtn = 1000
dat.tf = 365.25*10.
dat.dtmax = 365.25

dat.files.root = 'EGS'

dat.nobr = True

dat.run('EGS_stress_INPUT.dat')

# 7.2.9 Post-processing

cont = fcontour(dat.files.root+'*.csv',latest=True)

cont.slice_plot(save='temperature_slice.png',cbar=True, levels = 10,
	slice = ['z',-1000],divisions=[100,100],variable = 'T',xlims = [0,500],ylims=[0,500],title='final temperature distribution')

cont.profile_plot(save='pressure_profile.png',profile = np.array([[0,0,-1000],[400,400,-1000]]),variable = 'P',
ylabel = 'pressure / MPa', title='pressure with distance from injector', color = 'g', marker = 'o--',method = 'linear')

cont.profile_plot(save='stress_profile.png',profile = np.array([[600,600,-500],[600,600,-1500]]),variable = 'strs_xx',
xlabel='depth / m',ylabel = 'pressure / MPa', title='horizontal stress with depth', color = 'k', marker = '-',method = 'linear',elevationPlot = True)

cont.cutaway_plot(variable='T', save='temperature_cutaway.png',xlims=[0,400],
ylims = [0,400],zlims = [-1000,-800],cbar=True,levels=np.linspace(60,90,11),
grid_lines='k:',title='temperature contours / $^o$C')

hist = fhistory(dat.files.root+'_flow_his.dat')

hist.time_plot(variable='flow',save='extraction_plot.png',node=hist.nodes[0],scale_t=1./365.25,
xlabel='time / years',ylabel='flow kg s$^{-1}$')











