"""Various templates for use with PyFEHM."""

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
import scipy as sp
import os,math
from fdata import*
from ftool import*
from matplotlib import pyplot as plt
from fvars import*
from fdflt import*

dflt = fdflt()

#-----------------------------------------------------------------------------------------------------
#------------------------------------- FEHM MODEL TEMPLATES ------------------------------------------
#-----------------------------------------------------------------------------------------------------
class wellbore_model(fdata):
	'''Create a simple 2D radial well bore model. Returns an fdata object corresponding to the model.
	
	User inputs dimensions of the wellbore, steel pipe, casing and reservoir, their material properties and 
	permeability, information about initial reservoir conditions and the injection operation. 
	The model grid is generated on initialisation
	but this can be replaced with a more complex mesh by calling ``read_grid()`` - zones will be automatically 
	reassigned.
	
	The simulation is executed by calling the ``run()`` method.
	
	Output data are visualised by calling ``plot()`` or ``summarise()`` methods.
	
	:param xL: Horizontal dimension of model.
	:type xL: fl64
	:param zL: Vertical dimension of model.
	:type zL: fl64
	:param wellbore_radius: Well-bore radius.
	:type wellbore_radius: fl64
	:param wellbore_xdiv: Horizontal grid divisions in wellbore.
	:type wellbore_xdiv: int
	:param pipe_width: Steel pipe width.
	:type pipe_width: fl64
	:param pipe_xdiv: Horizontal grid divisions in steel pipe.
	:type pipe_xdiv: int
	:param casing_width: Casing width.
	:type casing_width: fl64
	:param casing_xdiv: Horizontal grid divisions in casing.
	:type casing_xdiv: int
	:param reservoir_xdiv: Horizontal grid divisions in reservoir.
	:type reservoir_xdiv: int
	
	:param zdiv: Number of vertical divisions for the grid (not including the feedzone).
	:type zdiv: lst[int]
	:param zprop: Proportioning of vertical dimension for gridding.
	:type zprop: lst[fl64]
	
	:param injection_temperature: Temperature of fluid injected at the wellhead.
	:type injection_temperature: fl64
	:param injection_flow_rate: Flow rate of fluid injected at the wellhead.
	:type injection_flow_rate: fl64
	
	:param initial_temperature: Specifies initial temperature conditions in the reservoir. If positive, initial temperature is interpreted as isotropic. If negative, initial temperature is interpreted as a vertical gradient. If a string, initial temperature corresponds to a text file containing a temperature-depth profile and applies this as the initial reservoir temperatures.
	:type initial_temperature: fl64, str
	
	:param simulation_time: Length of simulation in days.
	:type simulation_time: fl64
	
	:param inputfilename: Name of input file.
	:type inputfilename: str
	:param gridfilename: Name of grid file.
	:type gridfilename: str
	
	:param pipe_density: Density of steel pipe.
	:type pipe_density: fl64
	:param pipe_specific_heat: Specific heat of steel pipe.
	:type pipe_specific_heat: fl64
	:param casing_density: Density of casing.
	:type casing_density: fl64
	:param casing_specific_heat: Specific heat of casing.
	:type casing_specific_heat: fl64
	:param reservoir_density: Density of reservoir.
	:type reservoir_density: fl64
	:param reservoir_specific_heat: Specific heat of reservoir.
	:type reservoir_specific_heat: fl64
	:param reservoir_permeability: Reservoir permeability, specified as either isotropic k0 or anisotropic [kx,ky,kz].
	:type reservoir_permeability: fl64, lst[fl64]
	:param wellbore_permeability: Wellbore permeability, surrogate representing rapid transport down well. Set to high value.
	:type wellbore_permeability: fl64, lst[fl64]
	'''
	def __init__(self,
		xL, zL, wellbore_radius, wellbore_xdiv, pipe_width, pipe_xdiv, casing_width, casing_xdiv, 
		reservoir_xdiv, zdiv, zprop,
		
		injection_temperature, injection_flow_rate,
		
		initial_temperature,
		
		simulation_time,
		
		inputfilename='', gridfilename='', work_dir = None,
		
		pipe_density = 2500., pipe_specific_heat=1.e3, pipe_conductivity = 2., 
		casing_density = 2500., casing_specific_heat=1.e3, casing_conductivity = 0.5, 
		reservoir_density = 2500., reservoir_specific_heat=1.e3, reservoir_conductivity = 1., reservoir_porosity = 0.1, reservoir_permeability = 1.e-15, 
		wellbore_permeability=1.e-4,
		
		surface_pressure=0.1, surface_temperature=25.
		):
		
		# 1. inherit properties and initialise
		super(wellbore_model,self).__init__(filename='',gridfilename='',inconfilename='',sticky_zones=dflt.sticky_zones,associate=dflt.associate,work_dir = None,
			full_connectivity=dflt.full_connectivity,skip=[],keep_unknown=dflt.keep_unknown)
		self.nobr = True
		
		self._injection_temperature = injection_temperature
		self._injection_flow_rate = injection_flow_rate
		
		self._initial_temperature = initial_temperature
		
		self._simulation_time = simulation_time
		
		self._pipe_density = pipe_density
		self._casing_density = casing_density
		self._reservoir_density = reservoir_density
		
		self._pipe_specific_heat = pipe_specific_heat
		self._casing_specific_heat = casing_specific_heat
		self._reservoir_specific_heat = reservoir_specific_heat
		
		self._pipe_conductivity = pipe_conductivity
		self._casing_conductivity = casing_conductivity
		self._reservoir_conductivity = reservoir_conductivity
		
		self._reservoir_porosity = reservoir_porosity
		self._reservoir_permeability = reservoir_permeability
		self._wellbore_permeability = wellbore_permeability
		
		self._surface_pressure = surface_pressure
		self._surface_temperature = surface_temperature
		
		self._remember_inittemp = None
		
		if inputfilename: self._filename = inputfilename
		if gridfilename: self.grid._path.filename = gridfilename
		
		self.work_dir = work_dir
				
		# 2. create a grid
		if xL == 0 or zL == 0: print 'Error: grid dimensions must be non-zero'; return
		if wellbore_radius+casing_width+pipe_width>xL: print 'Error: No room for reservoir rock.'; return
		#if (abs(feedzone_depth)+feedzone_width/2)>abs(zL): print 'Error: model not deep enough to include feedzone.'; return
		
		if gridfilename: meshname = gridfilename
		elif inputfilename: meshname = inputfilename.split('.')[0]+'grid'
		else: meshname = 'wellbore.grid'
		zL = abs(zL); xL = abs(xL)
		reservoir_width = xL - wellbore_radius - pipe_width - casing_width
		
		x = (list(np.linspace(0,wellbore_radius,wellbore_xdiv))+
			list(np.linspace(wellbore_radius,wellbore_radius+pipe_width,pipe_xdiv))[1:]+
			list(np.linspace(wellbore_radius+pipe_width,wellbore_radius+pipe_width+casing_width,casing_xdiv))[1:]+
			list(np.linspace(wellbore_radius+pipe_width+casing_width,xL,reservoir_xdiv))[1:])
		z = np.linspace(-zL,0,zdiv)
		
		self.grid.make(meshname,x=x,y=z,z=[0.],radial=True)
		
		# 3. create zones
		x0,x1 = self.grid.xmin,self.grid.xmax
		z0,z1 = self.grid.ymin,self.grid.ymax
		
		wb = fzone(index=1,name='wellbore')
		wb.rect([x0-0.01,z0-0.01],[x0-0.01+wellbore_radius,z1+0.01])
		self.add(wb)
		
		wb = fzone(index=2,name='wellhead')
		wb.rect([x0-0.01,z1-0.01],[x0-0.01+wellbore_radius,z1+0.01])
		self.add(wb)
		
		wb = fzone(index=3,name='wellbase')
		wb.rect([x0-0.01,z0-0.01],[x0-0.01+wellbore_radius,z0+0.01])
		self.add(wb)
		
		pp = fzone(index=10,name='pipe_upper')
		pp.rect([x0-0.01+wellbore_radius,z0-0.01],[x0-0.01+wellbore_radius+pipe_width,z1+0.01])
		self.add(pp)
		
		cs = fzone(index=20,name='casing_upper')
		cs.rect([x0-0.01+wellbore_radius+pipe_width,z0-0.01],
				[x0-0.01+wellbore_radius+pipe_width+casing_width,z1+0.01])
		self.add(cs)
		
		rk = fzone(index=100,name='reservoir')
		rk.rect([x0-0.01+wellbore_radius+pipe_width+casing_width,z0-0.01],[x1+0.01,z1+0.01])
		self.add(rk)
		
		sz = fzone(index=900, name = 'surface')
		sz.rect([x0-0.01+wellbore_radius,z1-0.01],[x1+0.01,z1+0.01])
		self.add(sz)
		
		ff = fzone(index=901, name = 'farfield')
		ff.rect([x1-0.01,z0-0.01],[x1+0.01,z1+0.01])
		self.add(ff)
		
		# 4. assign material properties
		# 	- conductivity
		self.add(fmacro('cond',zone=0,param=(('cond_x',1),('cond_y',1),('cond_z',1))))
		self.add(fmacro('cond',zone=10,param=(('cond_x',self.pipe_conductivity),('cond_y',self.pipe_conductivity),('cond_z',self.pipe_conductivity))))
		self.add(fmacro('cond',zone=20,param=(('cond_x',self.casing_conductivity),('cond_y',self.casing_conductivity),('cond_z',self.casing_conductivity))))
		self.add(fmacro('cond',zone=100,param=(('cond_x',self.reservoir_conductivity),('cond_y',self.reservoir_conductivity),('cond_z',self.reservoir_conductivity))))
		# 	- rock
		self.add(fmacro('rock',zone=10,param=(('density',self.pipe_density),('porosity',0),('specific_heat',self.pipe_specific_heat))))
		self.add(fmacro('rock',zone=20,param=(('density',self.casing_density),('porosity',0),('specific_heat',self.casing_specific_heat))))
		self.add(fmacro('rock',zone=100,param=(('density',self.reservoir_density),('porosity',self.reservoir_porosity),('specific_heat',self.reservoir_specific_heat))))
		self.add(fmacro('rock',zone=1,param=(('density',2500.),('porosity',1.),('specific_heat',1.e3))))
		# 	- permeability
		self.add(fmacro('perm',zone=0,param=(('kx',1.e-15),('ky',1.e-15),('kz',1.e-15))))
		self.add(fmacro('perm',zone=10,param=(('kx',1.e-20),('ky',1.e-20),('kz',1.e-20))))
		self.add(fmacro('perm',zone=20,param=(('kx',1.e-20),('ky',1.e-20),('kz',1.e-20))))
		
		reservoirPerm = fmacro('perm',zone=100)
		if (isinstance(reservoir_permeability,list) or isinstance(reservoir_permeability,tuple) or
			isinstance(reservoir_permeability,np.ndarray)):
			reservoirPerm.param['kx']=reservoir_permeability[0]
			reservoirPerm.param['ky']=reservoir_permeability[1]
			if len(reservoir_permeability) == 2:
				reservoirPerm.param['kz']=reservoir_permeability[1]
			elif len(reservoir_permeability) == 3:
				reservoirPerm.param['kz']=reservoir_permeability[2]
		else:
			reservoirPerm.param['kx']=reservoir_permeability
			reservoirPerm.param['ky']=reservoir_permeability
			reservoirPerm.param['kz']=reservoir_permeability
		self.add(reservoirPerm)
		
		wellborePerm = fmacro('perm',zone=1)
		wellborePerm.param['kx']=wellbore_permeability
		wellborePerm.param['ky']=wellbore_permeability
		wellborePerm.param['kz']=wellbore_permeability
		self.add(wellborePerm)
		
		
		# 5. assign injection - use boun and distributed flow source
		self.add(fboun(zone=[self.zone['wellhead']],type='ti',times=[0,self.simulation_time],
		variable=[['dsw',-self.injection_flow_rate,-self.injection_flow_rate],
		['t',self.injection_temperature,self.injection_temperature]]))
		
		# 6. assign initial conditions
		self.add(fmacro('grad',zone=0,param=(('reference_coord',z1),('direction',2),('variable',1),
			('reference_value',0.1),('gradient',-0.00981))))
		self._set_reservoir_temperature()
		
		# 7. assign boundary conditions
		self.add(fmacro('flow',zone=900,param=(('rate',self.surface_pressure),('energy',-self.surface_temperature),('impedance',100))))
		self.add(fmacro('flow',zone=901,param=(('rate',0),('energy',-self.surface_temperature),('impedance',-100))))
		self.add(fmacro('flow',zone=3,param=(('rate',0),('energy',-self.surface_temperature),('impedance',-100))))
		
		# 8. assign simulation parameters
		self.ctrl['geometry_ICNL'] = 4
		self.ctrl['gravity_direction_AGRAV'] = 2
		self.ctrl['stor_file_LDA'] = 0
		self.ctrl['min_timestep_DAYMIN'] = 1.e-8
		self.ctrl['max_timestep_DAYMAX'] = self.simulation_time/10
		
		self.time['initial_timestep_DAY'] = self.simulation_time/1.e2
		self.time['max_time_TIMS'] = self.simulation_time
		self.time['max_timestep_NSTEP'] = 2000
		
		# 9. assign history and contour output
		self.hist.variables.append(['pressure','temperature','flow'])
		self.hist.timestep_interval=1
		self.hist.format='tec'
		self.hist.time_interval=1.e20
		self.hist.nodelist = self.zone['wellbore'].nodelist

		self.cont.variables.append(['temperature','pressure','xyz','liquid'])
	def _set_reservoir_temperature(self):
		if isinstance(self._initial_temperature,str):
			if not os.path.isfile(self._initial_temperature): print 'ERROR: '+self._initial_temperature+' does not exist.'; return
			tempfile = open(self._initial_temperature,'r')
			lns = tempfile.readlines()
			commaFlag = False; spaceFlag = False
			if len(lns[0].split(',')) > 1: commaFlag = True
			elif len(lns[0].split()) > 1: spaceFlag = True
			if not commaFlag and not spaceFlag: print 'ERROR: incorrect formatting for '+self._initial_temperature+'. Expect first column depth (m) and second column temperature (degC), either comma or space separated.'; return
			zs = []; Ts = []
			for ln in lns:
				if commaFlag: ln = ln.split(',')
				elif spaceFlag: ln = ln.split()
				zs.append(float(ln[0])); Ts.append(float(ln[1]))
			allNeg = True 	#flag that info read in has positive z-coords (make negative)
			for z in zs:
				if z>0: allNeg = False; break
			if allNeg: zs = [-z for z in zs]
			
			zm = np.unique([nd.position[1] for nd in self.grid.nodelist])
			zn_ind = 0
			x0,x1 = self.grid.xmin,self.grid.xmax
			y0,y1 = self.grid.ymin,self.grid.ymax
			for z in zm:
				zn=fzone(index=200+zn_ind);zn.rect([x0-0.01,z-0.01],[x1+0.01,z+0.01]);self.add(zn)
				T = np.interp(-z,zs,Ts)
				self.add(fmacro('pres',zone=200+zn_ind,param=(('pressure',5),('temperature',T),('saturation',1))))
				zn_ind+=1			
		else:
			if self._initial_temperature > 0:
				temp = fmacro('pres',zone=0,param=(('pressure',1),('temperature',self._initial_temperature),('saturation',1)))
			elif self._initial_temperature < 0:
				temp = fmacro('grad',zone=0,param=(('reference_coord',self.grid.zmax),('direction',2),('variable',2),
					('reference_value',self.surface_temperature),('gradient',self._initial_temperature)))
			self.add(temp)
			self._remember_inittemp = temp
	def _clear_reservoir_temperature(self):
		if isinstance(self._remember_inittemp,fmacro):
			self.delete(self._remember_inittemp)
			self._remember_inittemp = None
	def plot(self,temperature_lims = [],pdf = '',combineString = 'gswin64',
	
		Tslice = True, Tslice_xlims=[],Tslice_ylims=[],	Tslice_divisions=[100,100], Tslice_method = 'nearest',
		
		Pslice = True, Pslice_xlims=[],Pslice_ylims=[],	Pslice_divisions=[100,100], Pslice_method = 'nearest',
		
		Ttime = True, Ttime_xlims=[], Ttime_ylims=[],
		
		Twell = True, Twell_times = [], Twell_xlims=[], Twell_ylims = [], Twell_initial=True, Twell_profiles=None,
		Twell_output = False,
		
		Pcorrection = True, Pcorrection_xlims=[], Pcorrection_ylims=[],
		
		imperial_units = False,
		
		write_out = False
		):
		'''Generate plots of wellbore simulation.
		
		:param temperature_lims: Limits on temperature axis of temperature vs. time plot.
		:type temperature_lims: lst[fl64,fl64]
		
		:param pdf: Name of pdf file to combine all output plots. If not specified, pdf will not be created.
		:type pdf: str
		
		:param combineString: Name of ghostscript executable.
		:type combineString: str
		'''
		# read in contour data
		if self.cont.format == 'surf':
			if self.work_dir:
				cont = fcontour(self.work_dir+os.sep+self.files.root+'.*_days_sca_node.csv',latest=True)
			else:
				cont = fcontour(self.files.root+'.*_days_sca_node.csv',latest=True)
		# read in history data
		if self.hist.format == 'tec':
			if self.work_dir:
				hist = fhistory(self.work_dir+os.sep+self.files.root+'_*_his.dat')
			else:
				hist = fhistory(self.files.root+'_*_his.dat')
		
		if pdf: 
			if self.work_dir: pdf = self.work_dir+os.sep+pdf
			ext = 'eps'
			mp = multi_pdf(save=pdf)
			mp.combineString = combineString
		else: ext = 'png'
		Tmax = np.max(cont[cont.times[0]]['T'])
		
		# 1. slice plot of temperature
		if Tslice:
			xlims = [self.grid.xmin,self.grid.xmax]
			ylims = [self.grid.ymin,self.grid.ymax]
			eqa = True
			if Tslice_xlims: xlims = Tslice_xlims; eqa = False
			if Tslice_ylims: ylims = Tslice_ylims; eqa = False
			scale = 1.
			title = 'temperature / $^o$C'
			if imperial_units: scale = [9/5.,+32.]; title = 'temperature / $^o$F'
			cont.slice_plot(variable = 'T',save=self.files.root+'_temperature.'+ext,cbar = True,
				ylabel='z / m',title = title,xlims = xlims,
				ylims = ylims,divisions=Tslice_divisions,method=Tslice_method,
				equal_axes = eqa,scale = scale)			
			if pdf: mp.add(self.files.root+'_temperature.eps')
		
		# 2. slice plot of pressure
		if Pslice:
			xlims = [self.grid.xmin,self.grid.xmax]
			ylims = [self.grid.ymin,self.grid.ymax]
			eqa = True
			if Pslice_xlims: xlims = Pslice_xlims; eqa = False
			if Pslice_ylims: ylims = Pslice_ylims; eqa = False
			scale = 1.
			title = 'pressure / MPa'
			if imperial_units: scale = 145.05; title = 'pressure / psi'
			cont.slice_plot(variable = 'P',save=self.files.root+'_pressure.'+ext,cbar = True,
				ylabel='z / m',title = title,xlims = xlims,
				ylims = ylims,divisions=Pslice_divisions,method=Pslice_method,
				equal_axes = eqa,scale=scale)			
			if pdf: mp.add(self.files.root+'_pressure.eps')
		
		# 3. time series plot of bottom hole temperature
		if Ttime:
			xlims = [hist.times[0],hist.times[-1]]
			if Ttime_xlims: xlims = Ttime_xlims
			nd = self.grid.node_nearest_point([self.grid.xmin,self.grid.ymin,0])
			Tnd = hist['T'][nd.index]
			y0,y1 = np.min(Tnd),np.max(Tnd); ymid = (y0+y1)/2; yrange = y1-y0
			ylims = [ymid-0.55*yrange,ymid+0.55*yrange]			
			if Ttime_ylims: ylims = Ttime_ylims
			scale = 1.
			ylabel = 'temperature / $^o$C'
			if imperial_units: 
				scale = [9/5.,+32.]
				ylabel = 'temperature / $^o$F'
				ylims = list(scale[0]*np.array(ylims)+scale[1])
			hist.time_plot(node = nd.index,save=self.files.root+'_temperatureWB.'+ext,
			variable='T',xlabel='t / days', ylabel=ylabel,title='temperature at well bore',
			var_lim=ylims,t_lim=xlims,scale = scale)
			mp.add(self.files.root+'_temperatureWB.'+ext)
			
			# write out text 
			if write_out:
				if self.work_dir:
					of = open(self.work_dir+os.sep+self.files.root+'_DH_temp.dat','w')
				else:
					of = open(self.files.root+'_DH_temp.dat','w')
				t = hist.times
				T = hist['T'][nd.index]
				for ti,Ti in zip(t,T):
					of.write(str(ti)+','+str(Ti)+'\n')
				of.close()
			
		# 4. down well plot of temperature, multiple times
		if Twell:
			xlabel = 'temperature / $^o$C'
			if imperial_units: 
				xlabel = 'temperature / $^o$F'
			cols = ['k-','r-','b-','k--','r--','b--','k:','r:','b:']
			if len(Twell_times)>9: Twell_times = Twell_times[:9]
			# get times
			if not Twell_times: ts = [hist.times[-1]]
			else:
				ts = []
				for t in Twell_times:
					dt = abs(np.array(hist.times) - t)
					ts.append(hist.times[np.where(dt == np.min(dt))[0][0]])
			# plot all profiles
			plt.figure(figsize=[8,8])
			ax = plt.axes([0.15,0.15,0.75,0.75])
			ax.set_title('downhole temperature profiles',size='medium')
			x0 = self.grid.xmin
			allNd = [nd for nd in self.hist.nodelist if nd.position[0] == x0]
			allNd.sort(key=lambda x: x.position[1])
			for t,col in zip(ts,cols[:len(ts)]):
				tind = np.where(hist.times==t)[0][0]
				if imperial_units:
					T = [hist['T'][nd.index][tind]*9/5.+32. for nd in allNd]
				else:
					T = [hist['T'][nd.index][tind] for nd in allNd]
				z = [nd.position[1] for nd in allNd]
				
				ax.plot(T,z,col)
				
				if Twell_output:
					filename = self.files.root+'_T_'+str(hist.times[tind])+'days'
					if imperial_units: filename+='_imperial'
					filename+='.txt'
					lns = []
					for zi,Ti in zip(z,T): lns.append(str(zi)+','+str(Ti)+'\n')
					file = open(filename,'w'); file.writelines(lns); file.close()
				
			if Twell_initial:
				if imperial_units:
					T = [hist['T'][nd.index][0]*9/5.+32. for nd in allNd]
				else:
					T = [hist['T'][nd.index][0] for nd in allNd]
				z = [nd.position[1] for nd in allNd]
				
				ax.plot(T,z,'g-')
				
				if Twell_output:
					filename = self.files.root+'_T_0days'
					if imperial_units: filename+='_imperial'
					filename+='.txt'
					lns = []
					for zi,Ti in zip(z,T): lns.append(str(zi)+','+str(Ti)+'\n')
					file = open(filename,'w'); file.writelines(lns); file.close()
				
				
				
			if Twell_xlims: ax.set_xlim(Twell_xlims)
			if Twell_ylims: ax.set_ylim(Twell_ylims)
			ax.set_xlabel(xlabel)
			ax.set_ylabel('z / m')
			xlims = ax.get_xlim()
			ylims = ax.get_ylim()

						
			if Twell_profiles:
				if isinstance(Twell_profiles,str): Twell_profiles = [Twell_profiles]
				
				ylims = ax.get_ylim()
				for profile,col in zip(Twell_profiles,cols[:len(Twell_profiles)]):
					if not os.path.isfile(profile): print 'WARNING: '+profile+' not found'; continue
					dat = np.loadtxt(profile)
					
					zs = dat[:,0]; Ts = dat[:,1]					
					allNeg = True 	#flag that info read in has positive z-coords (make negative)
					for z in zs:
						if z>0: allNeg = False; break
					if allNeg: zs = [-z for z in zs]
					
					N = int(abs(zs[-1]-zs[0])/abs(ylims[0]-ylims[1])*24)+3
					z = np.linspace(np.min(zs),np.max(zs),N)
					Ts = np.interp(z,zs,Ts)
					
					ax.plot(Ts,-1*np.array(z),col+'x')
					
			# set legend for flow rates	
			text_size = 'small'
			xc,yc = 0.7,0.94
			ln_len = 0.06
			txt_gap = 0.04
			dy = 0.04; 
			
			xr = xlims[1]-xlims[0]
			yr = ylims[1]-ylims[0]
			dyi = -dy*(ylims[1]-ylims[0])
			x1 = xlims[0] +xc*xr
			y1 = ylims[0] +yc*yr
			import math
			if Twell_initial:
				ax.plot([x1,x1+ln_len*xr],[y1,y1],'g-')
				ax.text(x1+(ln_len+txt_gap)*xr,y1,'t = 0 days',size=text_size,ha='left',va='center')
				y1 = y1 + dyi
			
			cnt = 0
			for t,col in zip(ts,cols[:len(ts)]):
				ax.plot([x1,x1+ln_len*xr],[y1,y1],col)
				if isinstance(t,list): t = t[0]
				if t>=1 and t<10: tstr = str(round(t*10)/10)
				elif t>=10: tstr = str(int(t))
				else: tstr = str(t)
				tstr
				ax.text(x1+(ln_len+txt_gap)*xr,y1,'t = '+tstr+' days',size=text_size,ha='left',va='center')
				y1 = y1 + dyi
				
				if Twell_profiles:
					if cnt == len(Twell_profiles): continue
					ax.plot([x1,x1+ln_len*xr],[y1,y1],col)
					ax.plot((2*x1+ln_len*xr)/2,y1,col+'x')
					if isinstance(t,list): t = t[0]
					if t>=1 and t<10: tstr = str(round(t*10)/10)
					elif t>=10: tstr = str(int(t))
					else: tstr = str(t)
					tstr
					ax.text(x1+(ln_len+txt_gap)*xr,y1,Twell_profiles[cnt].split('.')[0],size=text_size,ha='left',va='center')
					y1 = y1 + dyi
					
					cnt +=1
				
			ax.set_xlim(xlims)
			ax.set_ylim(ylims)
			
			plt.savefig(self.files.root+'_Twell.'+ext, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
			format=ext,transparent=True, bbox_inches=None, pad_inches=0.1)
			if pdf: mp.add(self.files.root+'_Twell.'+ext)
		
		# 5. plot pressure correction at base of hole
		if Pcorrection:
			x0 = self.grid.xmin
			allNd = [nd for nd in self.hist.nodelist if nd.position[0] == x0]
			allNd.sort(key=lambda x: x.position[1])			
			T = [hist['T'][nd.index][0] for nd in allNd] 	# temperature profile
			z = [nd.position[1] for nd in allNd]
			P = [abs(nd.position[1])*1e3*9.81/1e6 for nd in allNd]
			rho = dens(P,T)[0]
			dP0 = abs(np.trapz(rho,z))*9.81/1e6
			dP=[]
			for t in hist.times:
				tind = np.where(hist.times==t)[0][0]
				T = [hist['T'][nd.index][tind] for nd in allNd] 	# temperature profile
				rho = dens(P,T)[0]
				dP.append(abs(np.trapz(rho,z))*9.81/1e6 - dP0)
			plt.figure(figsize=[8,8])
			plt.clf()
			ax = plt.axes([0.15,0.15,0.75,0.75])
			ax.set_title('downhole density pressure correction',size='medium')
			scale = 1.
			if imperial_units: scale = 145.05
			ax.plot(hist.times,np.array(dP)*scale,'bx-')
			
			# write out text 
			if write_out:
				if self.work_dir:
					of = open(self.work_dir+os.sep+self.files.root+'_DH_dens_corr.dat','w')
					of2 = open(self.work_dir+os.sep+self.files.root+'_DH_pres.dat','w')
				else:
					of = open(self.files.root+'_DH_dens_corr.dat','w')
					of2 = open(self.files.root+'_DH_pres.dat','w')
				t = hist.times
				for ti,dp in zip(t,dP):
					of.write(str(ti)+','+str(dp)+'\n')
					of2.write(str(ti)+','+str(dp+dP0)+'\n')
				of.close()
				of2.close()
				
			ax.set_xlabel('time / days')
			if imperial_units:
				ax.set_ylabel('pressure correction / psi')
			else:
				ax.set_ylabel('pressure correction / MPa')
			
			ax.set_xlim([hist.times[0],hist.times[-1]])
			if Pcorrection_xlims: ax.set_xlim(Pcorrection_xlims)
			if Pcorrection_ylims: ax.set_ylim(Pcorrection_ylims)
			
			plt.savefig(self.files.root+'_Pcorrection.'+ext, dpi=100, facecolor='w', edgecolor='w',orientation='portrait', 
			format=ext,transparent=True, bbox_inches=None, pad_inches=0.1)
			if pdf: mp.add(self.files.root+'_Pcorrection.'+ext)
			
		if pdf: mp.make()
	# -------------------------------------- ATTRIBUTES ------------------------------------------------
	def _get_pipe_density(self): return self._pipe_density
	def _set_pipe_density(self,value): 
		self._pipe_density = value
		self.rock['pipe_upper'].param['density'] = value
	pipe_density = property(_get_pipe_density, _set_pipe_density) #: (*fl64*) Density of pipe.
	def _get_pipe_specific_heat(self): return self._pipe_specific_heat
	def _set_pipe_specific_heat(self,value): 
		self._pipe_specific_heat = value
		self.rock['pipe_upper'].param['specific_heat'] = value
	pipe_specific_heat = property(_get_pipe_specific_heat, _set_pipe_specific_heat) #: (*fl64*) Specific heat of pipe.
	def _get_casing_density(self): return self._casing_density
	def _set_casing_density(self,value): 
		self._casing_density = value
		self.rock['casing_upper'].param['density'] = value
	casing_density = property(_get_casing_density, _set_casing_density) #: (*fl64*) Density of casing.
	def _get_casing_specific_heat(self): return self._casing_specific_heat
	def _set_casing_specific_heat(self,value): 
		self._casing_specific_heat = value
		self.rock['casing_upper'].param['specific_heat'] = value
	casing_specific_heat = property(_get_casing_specific_heat, _set_casing_specific_heat) #: (*fl64*) Specific heat of casing
	def _get_reservoir_density(self): return self._reservoir_density
	def _set_reservoir_density(self,value): 
		self._reservoir_density = value
		self.rock['reservoir'].param['density'] = value
	reservoir_density = property(_get_reservoir_density, _set_reservoir_density) #: (*fl64*) Density of reservoir.
	def _get_reservoir_specific_heat(self): return self._reservoir_specific_heat
	def _set_reservoir_specific_heat(self,value): 
		self._reservoir_specific_heat = value
		self.rock['reservoir'].param['specific_heat'] = value
	reservoir_specific_heat = property(_get_reservoir_specific_heat, _set_reservoir_specific_heat) #: (*fl64*) Specific heat of reservoir
	def _get_reservoir_porosity(self): return self._reservoir_porosity
	def _set_reservoir_porosity(self,value): 
		self._reservoir_porosity = value
		self.rock['reservoir'].param['porosity'] = value
	reservoir_porosity = property(_get_reservoir_porosity, _set_reservoir_porosity) #: (*fl64*) Porosity of reservoir.
	def _get_reservoir_permeability(self): return self._reservoir_permeability
	def _set_reservoir_permeability(self,value): 
		self._reservoir_permeability = value
		if (isinstance(value,list) or isinstance(value,tuple) or
			isinstance(value,np.ndarray)):
			self.perm['reservoir'].param['kx']=value[0]
			self.perm['reservoir'].param['ky']=value[1]
			if len(value) == 2:
				self.perm['reservoir'].param['kz']=value[1]
			elif len(value) == 3:
				self.perm['reservoir'].param['kz']=value[2]
		else:
			self.perm['reservoir'].param['kx']=value
			self.perm['reservoir'].param['ky']=value
			self.perm['reservoir'].param['kz']=value
	reservoir_permeability = property(_get_reservoir_permeability, _set_reservoir_permeability) #: (*fl64*) Permeability of reservoir. If three element list, tuple or ndarray is passed, this is assumed to correspond to [kx,ky,kz].
	def _get_wellbore_permeability(self): return self._wellbore_permeability
	def _set_wellbore_permeability(self,value): 
		self._wellbore_permeability = value
		self.perm['wellbore'].param['kx'] = value
		self.perm['wellbore'].param['ky'] = value
		self.perm['wellbore'].param['kz'] = value
	wellbore_permeability = property(_get_wellbore_permeability, _set_wellbore_permeability) #: (*fl64*) Permeability in the wellbore. This should be set to a high number, representing free-flowing water.
	def _get_pipe_conductivity(self): return self._pipe_conductivity
	def _set_pipe_conductivity(self,value): 
		self._pipe_conductivity = value
		self.cond['pipe_upper'].param['cond_x'] = value
		self.cond['pipe_upper'].param['cond_y'] = value
		self.cond['pipe_upper'].param['cond_z'] = value
	pipe_conductivity = property(_get_pipe_conductivity, _set_pipe_conductivity) #: (*fl64*) Thermal conductivity of the steel pipe.
	def _get_casing_conductivity(self): return self._casing_conductivity
	def _set_casing_conductivity(self,value): 
		self._casing_conductivity = value
		self.cond['casing_upper'].param['cond_x'] = value
		self.cond['casing_upper'].param['cond_y'] = value
		self.cond['casing_upper'].param['cond_z'] = value
	casing_conductivity = property(_get_casing_conductivity, _set_casing_conductivity) #: (*fl64*) Thermal conductivity of the casing.
	def _get_reservoir_conductivity(self): return self._reservoir_conductivity
	def _set_reservoir_conductivity(self,value): 
		self._reservoir_conductivity = value
		self.cond['reservoir'].param['cond_x'] = value
		self.cond['reservoir'].param['cond_y'] = value
		self.cond['reservoir'].param['cond_z'] = value
	reservoir_conductivity = property(_get_reservoir_conductivity, _set_reservoir_conductivity) #: (*fl64*) Thermal conductivity of the reservoir.
	def _get_surface_pressure(self): return self._surface_pressure
	def _set_surface_pressure(self,value): 
		self._surface_pressure = value
		self.bounlist[0].variables = [['dsw',self.injection_flow_rate,self.injection_flow_rate],
		['t',self.injection_temperature,self.injection_temperature],
		['pw',value,value]]
	surface_pressure = property(_get_surface_pressure, _set_surface_pressure) #: (**) Pressure at top surface of model, default is atmospheric.
	def _get_surface_temperature(self): return self._surface_temperature
	def _set_surface_temperature(self,value): 
		self._surface_temperature = value
		self.flow['surface'].param['temperature'] = value
	surface_temperature = property(_get_surface_temperature, _set_surface_temperature) #: (**) Temperature at top surface of model, default is 25degC.
	def _get_injection_temperature(self): return self._injection_temperature
	def _set_injection_temperature(self,value): 
		self._injection_temperature = value
		self.bounlist[0].variables = [['dsw',self.injection_flow_rate,self.injection_flow_rate],['t',value,value],
		['pw',self.surface_pressure,self.surface_pressure]]
	injection_temperature = property(_get_injection_temperature, _set_injection_temperature) #: (*fl64*) Temperature of the fluid injected at the wellhead.
	def _get_injection_flow_rate(self): return self._injection_flow_rate
	def _set_injection_flow_rate(self,value): 
		self._injection_flow_rate = value
		self.bounlist[0].variables = [['dsw',value,value],['t',self.injection_temperature,self.injection_temperature],
		['pw',self.surface_pressure,self.surface_pressure]]
	injection_flow_rate = property(_get_injection_flow_rate, _set_injection_flow_rate) #: (*fl64*) Flow rate of fluid injected at the wellhead.
	def _get_simulation_time(self): return self._simulation_time
	def _set_simulation_time(self,value): 
		self._simulation_time = value
		self.bounlist[0].times = [0.,value]
		self.time['initial_timestep_DAY'] = value/1.e2
		self.time['max_time_TIMS'] = value
	simulation_time = property(_get_simulation_time, _set_simulation_time) #: (*fl64*) Length of simulation in days
	def _get_initial_temperature(self): return self._initial_temperature
	def _set_initial_temperature(self,value): 
		self._initial_temperature = value
		self._clear_reservoir_temperature()
		self._set_reservoir_temperature()
	initial_temperature = property(_get_initial_temperature, _set_initial_temperature) #: (*fl64*, *str*) Specifies initial temperature conditions in the reservoir. If positive, initial temperature is interpreted as isotropic. If negative, initial temperature is interpreted as a vertical gradient. If a string, initial temperature corresponds to a text file containing a temperature-depth profile and applies this as the initial reservoir temperatures.
