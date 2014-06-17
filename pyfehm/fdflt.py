"""Environment file for PyFEHM. Set default attribute values."""

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

import os,platform,pkgutil

from types import*

floatKeys = ['linear_converge_NRmult_G1','quadratic_converge_NRmult_G2','stop_criteria_NRmult_G3',
	'machine_tolerance_TMCH','overrelaxation_factor_OVERF','newton_cycle_tolerance_EPM',
	'upstream_weighting_UPWGT','timestep_multiplier_AIAA','min_timestep_DAYMIN','max_timestep_DAYMAX',
	'initial_timestep_DAY','max_time_TIMS','initial_year_YEAR','initial_month_MONTH','initial_day_INITTIME',
	'init_solute_conc_ANO','implicit_factor_AWC','tolerance_EPC','upstream_weight_UPWGTA','solute_start_DAYCS',
	'solute_end_DAYCF','flow_end_DAYHF','flow_start_DAYHS','max_iterations_IACCMX','timestep_multiplier_DAYCM',
	'initial_timestep_DAYCMM','max_timestep_DAYCMX','print_interval_NPRTTRC','alpha1_A1ADSF',
	'alpha2_A2ADSF','beta_BETADF']	
intKeys = ['reduced_dof_IRDOF','reordering_param_ISLORD','IRDOF_param_IBACK','number_SOR_iterations_ICOUPL',
	'max_machine_time_RNMAX','max_newton_iterations_MAXIT','number_orthogonalizations_NORTH',
	'max_solver_iterations_MAXSOLVE','JA','JB','JC','order_gauss_elim_NAR',	'max_multiply_iterations_IAMM',			
	'implicitness_factor_AAW','gravity_direction_AGRAV','geometry_ICNL','stor_file_LDA','max_timestep_NSTEP',
	'print_interval_IPRTOUT','coupling_NTT','element_integration_INTG','type_IADSF']
boolKeys = []
strKeys = ['acceleration_method_ACCM']
class fdflt(object):
	def __init__(self):			
		# material properties - these values will be assigned as defaults if not otherwise set
		self.permeability			=	1.e-15
		self.conductivity			=	2.2
		self.density				=	2500.
		self.porosity				=	0.1	
		self.specific_heat			=	1.e3	
		self.youngs_modulus			=	1.e4 	# MPa
		self.poissons_ratio			=	0.25
		self.pressure_coupling		=	1.
		self.thermal_expansion 		=	3.e-5	# / K
		
		# initial conditions
		self.Pi						=	1. 		# pressure
		self.Ti						=	30. 	# temperature 
		self.Si						=	1.		# saturation

		# output data formats
		self.hist_format 			= 	'tec'
		self.cont_format 			= 	'surf'
		self.parental_cont 			= 	True

		# set this to the fehm executable to be used if no default assigned
		self.fehm_path 				=	'c:\\path\\to\\fehm\\fehm.exe'
		if os.name is not 'posix':
			self.paraview_path 			=	'paraview.exe'
		else:
			self.paraview_path 			=	'paraview'
		self.lagrit_path			=	'c:\\path\\to\\lagrit\\lagrit.exe'
		self.files 					=	['outp','hist','check']
		self.co2_interp_path 		= 	'c:\\path\\to\\co2\\co2_interp_table.txt'
		self.co2_interp_path_2 		= 	'/alternate/path/to/co2/co2_interp_table.txt'
		if not os.path.isfile(self.co2_interp_path):
			self.co2_interp_path = self.co2_interp_path_2

		# fdata booleans
		self.associate 				= 	True		# associate macro, zone information with nodes
		self.sticky_zones 			= 	True		# print zone definitions immediately before use in input file
		self.full_connectivity 		=	True	
		self.sleep_time 			= 	1.
		self.keep_unknown 			= 	True 		# set true if PyFEHM should preserve unknown macros in future output files
		self.silent 				=	False		# turns off all PyFEHM verbiage
		
		# default values for mactro ITER (parameters controlling solver)
		self.iter = {
			'linear_converge_NRmult_G1':1.e-5, 			# convergence criteria
			'quadratic_converge_NRmult_G2':1.e-5,
			'stop_criteria_NRmult_G3':1.e-3,
			'machine_tolerance_TMCH':-1.e-5,
			'overrelaxation_factor_OVERF':1.1,
			
			'reduced_dof_IRDOF':0,
			'reordering_param_ISLORD':0,
			'IRDOF_param_IBACK':0,
			'number_SOR_iterations_ICOUPL':0,
			'max_machine_time_RNMAX':3600, 				# number of minutes at which FEHM will cut a simulation
			}
		# default values for macro CTRL (parameters controlling simulation)
		self.ctrl = {
			'max_newton_iterations_MAXIT':10,			# solver parameters
			'newton_cycle_tolerance_EPM':1.e-5,         # solver parameters
			'number_orthogonalizations_NORTH':8,        # solver parameters
			
			'max_solver_iterations_MAXSOLVE':24,
			'acceleration_method_ACCM':'gmre',
			
			'JA':1,'JB':0,'JC':0,
			'order_gauss_elim_NAR':2,
			
			'implicitness_factor_AAW':1,
			'gravity_direction_AGRAV':3, 				# direction of gravity
			'upstream_weighting_UPWGT':1.0,
			
			'max_multiply_iterations_IAMM':7,
			'timestep_multiplier_AIAA':1.5, 			# acceleration, time step multiplier
			'min_timestep_DAYMIN':1.e-5, 				# minimum allowable time step (days)
			'max_timestep_DAYMAX':30.,					# maximum allowable time step (days)
			
			'geometry_ICNL':0, 							# problem geometry (0 = 3-D)
			'stor_file_LDA':0 							# flag to use stor file
			}
		# default values for macro TIME
		self.time = {
			'initial_timestep_DAY':1., 					# initial time step size (days)
			'max_time_TIMS':365., 						# maximum simulation time (days)
			'max_timestep_NSTEP':200, 					# maximum number of time steps 
			'print_interval_IPRTOUT':1,					# for printing information to screen
			'initial_year_YEAR':None, 					# initial simulation time (years)	
			'initial_month_MONTH':None,					# (months)
			'initial_day_INITTIME':None					# (years)
			}
		# default values for macro SOL
		self.sol = {
			'coupling_NTT':1,
			'element_integration_INTG':-1
			}
		# default values for macro TRAC
		self.trac = {
			'init_solute_conc_ANO':0.,
			'implicit_factor_AWC':1.,
			'tolerance_EPC':1.e-7,
			'upstream_weight_UPWGTA':0.5,
			
			'solute_start_DAYCS':1.,
			'solute_end_DAYCF':2.,
			'flow_end_DAYHF':1.,
			'flow_start_DAYHS':2.,
			
			'max_iterations_IACCMX':50,
			'timestep_multiplier_DAYCM':1.2,
			'initial_timestep_DAYCMM':1.,
			'max_timestep_DAYCMX':1000.,
			'print_interval_NPRTTRC':1.
			}
		self.adsorption = {
			'type_IADSF':None,
			'alpha1_A1ADSF':None,
			'alpha2_A2ADSF':None,
			'beta_BETADF':None
			}
			
		# check to see if rc file exist, update defaults
		self._check_rc()
	def _check_rc(self):
		# check if pyfehmrc file exists		
		rc_lib = pkgutil.get_loader('fdflt').filename.split(os.sep)
		rc_lib1 = os.sep.join(rc_lib[:-1])+os.sep+'.pyfehmrc'
		rc_lib2 = os.sep.join(rc_lib[:-1])+os.sep+'pyfehmrc'
		rc_home1 = os.path.expanduser('~')+os.sep+'.pyfehmrc'
		rc_home2 = os.path.expanduser('~')+os.sep+'pyfehmrc'
		if os.path.isfile(rc_lib1): fp = open(rc_lib1)
		elif os.path.isfile(rc_lib2): fp = open(rc_lib2)
		elif os.path.isfile(rc_home1): fp = open(rc_home1)
		elif os.path.isfile(rc_home2): fp = open(rc_home2)
		else: return
		
		lns = fp.readlines()
		for ln in lns:
			ln = ln.split('#')[0]		# strip off the comment
			if ln.startswith('#'): continue
			elif ln.strip() == '': continue
			elif '&' in ln:
				if len(ln.split('&')) == 2:
					self._update_attribute(ln)
				elif len(ln.split('&')) == 3:
					self._update_dict(ln)
				else:
					print 'WARNING: unrecognized .pyfehmrc line \''+ln.strip()+'\''
			else:
				print 'WARNING: unrecognized .pyfehmrc line \''+ln.strip()+'\''
	def _update_attribute(self,ln):
		name,value = ln.split('&')
		name,value = name.strip(), value.strip()
		
		attributelist = self.__dict__.keys()
		
		if name not in attributelist:
			print 'ERROR: no attribute \''+name+'\''; return
			
		if isinstance(self.__dict__[name],dict):
			print 'ERROR: \''+name+'\' a dictionary. To set a dictionary value supply the dictionary key in format:'
			print 'dict_name : dict_key : value'
			return
		
		# translate None string
		if value in ['','None','none']: value = None
		
		if type(self.__dict__[name]) is IntType: 
			if value is not None: self.__setattr__(name,int(float(value)))
			else: self.__setattr__(name,None)
		elif type(self.__dict__[name]) is FloatType: 
			if value is not None: self.__setattr__(name,float(value))
			else: self.__setattr__(name,None)
		elif type(self.__dict__[name]) is StringType: 
			if value is not None: self.__setattr__(name,value)
			else:  self.__setattr__(name,None)
		elif type(self.__dict__[name]) is NoneType: 
			if value is not None: self.__setattr__(name,value)
			else: self.__setattr__(name,None)
		elif type(self.__dict__[name]) is BooleanType: 
			if value in ['True','1','1.']:
				self.__setattr__(name,True)
			elif value in ['False','0.','0'] or value == None:
				self.__setattr__(name,False)
			else:
				print 'ERROR: unrecognized boolean type \''+value+'\''; return
	def _update_dict(self,ln):
		name,key,value = ln.split('&')
		name,key,value = name.strip(), key.strip(), value.strip()
		
		dictlist = [k for k in self.__dict__.keys() if type(self.__dict__[k]) is dict]
		
		if name not in dictlist:
			print 'ERROR: no dictionary \''+name+'\''; return
			
		keys = self.__dict__[name].keys()
		if key not in keys:
			print 'ERROR: no such key \''+key+'\' in dictionary \''+name+'\''; return
		
		# translate None string
		if value in ['','None','none']: value = None
		
		if type(self.__dict__[name][key]) is IntType: 
			if value is not None: self.__dict__[name].__setitem__(key,int(float(value)))
			else: self.__setattr__(name,None)
		elif type(self.__dict__[name][key]) is FloatType: 
			if value is not None: self.__dict__[name].__setitem__(key,float(value))
			else: self.__setattr__(name,None)
		elif type(self.__dict__[name][key]) is StringType: 
			if value is not None: self.__dict__[name].__setitem__(key,value)
			else:  self.__setattr__(name,None)
		elif type(self.__dict__[name][key]) is NoneType: 
			if key in strKeys:
				if value is not None: self.__dict__[name].__setitem__(key,value)
				else:  self.__dict__[name].__setitem__(key,None)
			elif key in intKeys:
				if value is not None: self.__dict__[name].__setitem__(key,int(float(value)))
				else:  self.__dict__[name].__setitem__(key,None)
			elif key in floatKeys:
				if value is not None: self.__dict__[name].__setitem__(key,float(value))
				else:  self.__dict__[name].__setitem__(key,None)
			elif key in boolKeys:
				if value in ['True','1','1.']:
					self.__dict__[name].__setitem__(key,True)
				elif value in ['False','0.','0'] or value == None:
					self.__dict__[name].__setitem__(key,False)
			else: self.__setattr__(name,None)
		elif type(self.__dict__[name][key]) is BooleanType: 
			if value in ['True','1','1.']:
				self.__dict__[name].__setitem__(key,True)
			elif value in ['False','0.','0'] or value == None:
				self.__dict__[name].__setitem__(key,False)
			else:
				print 'ERROR: unrecognized boolean type \''+value+'\''; return
		
		
		
		
