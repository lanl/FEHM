"""For manipulating FEHM grids."""

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
import os,math,platform,string,subprocess,shutil
from subprocess import Popen

from time import time
from copy import deepcopy
from glob import glob

try:
	from matplotlib import pyplot as plt
	from mpl_toolkits.mplot3d import axes3d
except ImportError:
	'placeholder'
	
from fpost import*
from ftool import*

from itertools import izip

dflt = fdflt()

node_props = ('kx','ky','kz','cond_x','cond_y','cond_z','density','specific_heat','porosity','thermal_expansion','pressure_coupling',
'youngs_modulus','poissons_ratio')

nbr_update = dict((
	(0,[[0,2,4],[1,2,4]]),
	(1,[[1,2,4],[0,3,5]]),
	(2,[[0,3,4],[3,0,6]]),
	(3,[[1,3,4],[1,2,7]]),
	(4,[[0,2,5],[5,6,0]]),
	(5,[[1,2,5],[4,7,1]]),
	(6,[[0,3,5],[7,4,2]]),
	(7,[[1,3,5],[6,5,3]])
	))
def repair_grid(gridfilename):
	geo = fgrid()
	geo._read_inp(gridfilename)
	print 'The following duplicates were removed...'
	# with the repair flag, duplicate nodes will be deleted
	geo.add_nodetree(repair = geo)
	# now need to renumber nodes
	i0 = geo._nodelist[0].index
	for i,nd in enumerate(geo._nodelist[1:]):
		di = nd.index - 1 - i0
		if di != 0:
			for nd2 in geo._nodelist[i+1:]:
				nd2._index -= di
		i0 = nd.index
	
	newgridfilename = gridfilename.split('.')[0]+'.repaired' 
	print 'Writing repaired grid file '+newgridfilename
	geo.write(newgridfilename)
	geo.read(newgridfilename)
	return geo
class fnode(object):				#Node object.
	""" FEHM grid node object.
	"""
	__slots__ = ['_index','_position','_connections','_elements','_generator','_zone','_connected_nodes',
		'_permeability','_conductivity','_density','_specific_heat','_porosity','_youngs_modulus','_poissons_ratio','_thermal_expansion',
		'_pressure_coupling','_Pi','_Ti','_Si','_S_co2gi','_S_co2li','_co2aqi','_strsi','_dispi','_P','_T','_S',
		'_S_co2g','_S_co2l','_co2aq','_strs','_disp','_vol','_rlpmodel','_permmodel','_pormodel','_condmodel']
	def __init__(self,index,position):		
		self._index = index			
		self._position=position	
		self._connected_nodes = []
		self._connections = []	
		self._elements = []	
		self._generator = {}
		self._zone = {}
		self._permeability = None
		self._conductivity = None
		self._density = None
		self._specific_heat = None
		self._porosity = None
		self._youngs_modulus = None
		self._poissons_ratio = None
		self._thermal_expansion = None
		self._pressure_coupling = None
		self._Pi = None
		self._Ti = None
		self._Si = None
		self._S_co2gi = None
		self._S_co2li = None
		self._co2aqi = None
		self._strsi = None
		self._dispi = None	
		self._P = None
		self._T = None
		self._S = None
		self._S_co2g = None
		self._S_co2l = None
		self._co2aq = None
		self._strs = None
		self._disp = None		
		self._vol = None
		self._rlpmodel = None
		self._permmodel = None
		self._pormodel = None
		self._condmodel = None
	def __repr__(self): return 'nd'+str(self._index)
	def __getstate__(self):
		return dict((k, getattr(self, k)) for k in self.__slots__)
	def __setstate__(self, data_dict):
		for (name, value) in data_dict.iteritems():
			setattr(self, name, value)
	def _get_index(self): return self._index
	index = property(_get_index) #: (*int*) Integer number denoting the node.	
	def _get_position(self): return self._position
	position = property(_get_position) #: (*lst[fl64]*) List of the node's coordinates in 2- or 3-D space.
	def _get_zone(self): return self._zone
	def _set_zone(self,value): self._zone = value
	zone = property(_get_zone, _set_zone) #: (*dict*) Dictionary of zones to which the node belongs.
	def _get_generator(self): return self._generator
	def _set_generator(self,value): self._generator = value
	generator = property(_get_generator, _set_generator) #: (*dict*) Dictionary of generator properties associated with node.
	def _get_info(self):
		prntStr='\n Node number: '+str(self.index)+'\n'
		prntStr+='\nGeometric properties.......\n'
		prntStr+='          x = '+str(self.position[0])+'\n'
		prntStr+='          y = '+str(self.position[1])+'\n'
		prntStr+='          z = '+str(self.position[2])+'\n'
		if self.vol != None:
			prntStr+='     volume = '+str(self.vol)+'\n'
		prntStr+=' neighbours = ['
		for nd in self.connected_nodes:	prntStr+=str(nd.index)+','
		prntStr = prntStr[:-1]+']\n'
		if self.zone:
			prntStr+='\nZone membership............\n'
			for zn in self.zonelist:
				prntStr+='    '+str(zn.index)
				if zn.name:
					prntStr+=': '+zn.name
				prntStr+='\n'
		#if self.variable:
		prntStr+='\nState......................\n'
		if self.T or self.Ti: 
			if self.T and self.Ti:
				prntStr+='temperature = '+str(self.T)+' (init: '+str(self.Ti)+')\n'
			elif self.T:
				prntStr+='temperature = '+str(self.T)+'\n'
			elif self.Ti:
				prntStr+='temperature = '+str(self.Ti)+' (initial)\n'
		if self.P or self.Pi: 
			if self.P and self.Pi:
				prntStr+='   pressure = '+str(self.P)+' (init: '+str(self.Pi)+')\n'
			elif self.P:
				prntStr+='   pressure = '+str(self.P)+'\n'
			elif self.Pi:
				prntStr+='   pressure = '+str(self.Pi)+' (initial)\n'
		if self.S or self.Si: 
			if self.S and self.Si:
				prntStr+='  sat water = '+str(self.S)+' (init: '+str(self.Si)+')\n'
			elif self.S:
				prntStr+='  sat water = '+str(self.S)+'\n'
			elif self.Si:
				prntStr+='  sat water = '+str(self.Si)+' (initial)\n'
		prntStr+='\nMaterial properties........\n'
		if isinstance(self.permeability,(list,tuple,np.ndarray)):
			prntStr += '    kx = ' +str(self.permeability[0])+'\n'
			prntStr += '    ky = ' +str(self.permeability[1])+'\n'
			prntStr += '    kz = ' +str(self.permeability[2])+'\n'
		if isinstance(self.conductivity,(list,tuple,np.ndarray)):
			prntStr += '    cond_x = ' +str(self.conductivity[0])+'\n'
			prntStr += '    cond_y = ' +str(self.conductivity[1])+'\n'
			prntStr += '    cond_z = ' +str(self.conductivity[2])+'\n'
		if self.density:
			prntStr += '    density = ' +str(self.density)+'\n'
		if self.specific_heat:
			prntStr += '    spec. heat = ' +str(self.specific_heat)+'\n'
		if self.porosity:
			prntStr += '    porosity = ' +str(self.porosity)+'\n'
		if self.youngs_modulus:
			prntStr += '    Youngs modulus = ' +str(self.youngs_modulus)+'\n'
		if self.poissons_ratio:
			prntStr += '    Poissons ratio = ' +str(self.poissons_ratio)+'\n'
		if self.thermal_expansion:
			prntStr += '    Thermal expansion = ' +str(self.thermal_expansion)+'\n'
		if self.pressure_coupling:
			prntStr += '    Pressure coupling = ' +str(self.pressure_coupling)+'\n'
		if self.generator:
			prntStr+='\nGenerator properties.......\n'
			ks = self.generator.keys()
			for k in ks:
				prntStr += '    '+k + ' = ' +str(self.generator[k])+'\n'
		print prntStr
	what = property(_get_info)							#: Print to screen information about the node.
	def _get_connected_nodes(self):
		return self._connected_nodes
	connected_nodes = property(_get_connected_nodes)		#: (*lst[fnode]*) List of node objects connected to this node. This information only available if full_connectivity=True passed to fgrid.read()
	def _get_connections(self): return self._connections
	connections = property(_get_connections)#: (*lst[fconn]*) List of connection objects of which the node is a member.
	def _get_elements(self): return self._elements
	elements = property(_get_elements)#: (*lst[felem]*) List of element objects of which the node is a member.
	def _get_zonelist(self): return [self.zone[k] for k in self.zone.keys()]
	zonelist = property(_get_zonelist)	#: (*lst[fzone]*) List of zones of which the node is a member
	def _get_vol(self): return self._vol
	vol = property(_get_vol)	#: (*fl64*) Control volume associated with the node. This information only available if volumes() method called from grid attribute.
	def _get_permeability(self): return self._permeability
	permeability = property(_get_permeability) #: (*list*) permeability values at node.
	def _get_conductivity(self): return self._conductivity
	conductivity = property(_get_conductivity) #: (*list*) conductivity values at node.
	def _get_density(self): return self._density
	density = property(_get_density) #: (*fl64*) density at node.
	def _get_specific_heat(self): return self._specific_heat
	specific_heat = property(_get_specific_heat) #: (*fl64*) specific heat at node.
	def _get_porosity(self): return self._porosity
	porosity = property(_get_porosity) #: (*fl64*) porosity at node.
	def _get_youngs_modulus(self): return self._youngs_modulus
	youngs_modulus = property(_get_youngs_modulus) #: (*fl64*) Youngs modulus at node.
	def _get_poissons_ratio(self): return self._poissons_ratio
	poissons_ratio = property(_get_poissons_ratio) #: (*fl64*) Poissons ratio at node.
	def _get_thermal_expansion(self): return self._thermal_expansion
	thermal_expansion = property(_get_thermal_expansion) #: (*fl64*) Coefficient of thermal expansion at node.
	def _get_pressure_coupling(self): return self._pressure_coupling
	pressure_coupling = property(_get_pressure_coupling) #: (*fl64*) Biot pressure coupling coefficient at node.
	def _get_Pi(self): return self._Pi
	Pi = property(_get_Pi) #: (*fl64*) initial pressure at node.
	def _get_Ti(self): return self._Ti
	Ti = property(_get_Ti) #: (*fl64*) initial temperature at node.
	def _get_Si(self): return self._Si
	Si = property(_get_Si) #: (*fl64*) initial water saturation at node.
	def _get_S_co2gi(self): return self._S_co2gi
	S_co2gi = property(_get_S_co2gi) #: (*fl64*) initial gaseous CO2 saturation at node.
	def _get_S_co2li(self): return self._S_co2li
	S_co2li = property(_get_S_co2li) #: (*fl64*) initial liquid CO2 saturation at node.
	def _get_co2aqi(self): return self._co2aqi
	co2aqi = property(_get_co2aqi) #: (*fl64*) initial dissolved CO2 concentration at node.
	def _get_strsi(self): return self._strsi
	strsi = property(_get_strsi) #: (*fl64*) initial stresses at node.
	def _get_dispi(self): return self._dispi
	dispi = property(_get_dispi) #: (*fl64*) initial displacements at node.
	def _get_P(self): return self._P
	P = property(_get_P) #: (*fl64*) pressure at node during a simulation.
	def _get_T(self): return self._T
	T = property(_get_T) #: (*fl64*) temperature at node.
	def _get_S(self): return self._S
	S = property(_get_S) #: (*fl64*) water saturation at node.
	def _get_S_co2g(self): return self._S_co2g
	S_co2g = property(_get_S_co2g) #: (*fl64*) gaseous CO2 saturation at node.
	def _get_S_co2l(self): return self._S_co2l
	S_co2l = property(_get_S_co2l) #: (*fl64*) liquid CO2 saturation at node.
	def _get_co2aq(self): return self._co2aq
	co2aq = property(_get_co2aq) #: (*fl64*) dissolved CO2 concentration at node.
	def _get_strs(self): return self._strs
	strs = property(_get_strs) #: (*list*) stresses at node ([xx,yy,xy] for 2D, [xx,yy,zz,xy,yz,xz] for 3D).
	def _get_disp(self): return self._disp
	disp = property(_get_disp) #: (*list*) displacements at node ([x,y] for 2D, [x,y,z] for 3D).
	def _get_rlpmodel(self): return self._rlpmodel
	rlpmodel = property(_get_rlpmodel) #: (*int*) index of relative permeability model assigned to node
	def _get_permmodel(self): return self._permmodel
	permmodel = property(_get_permmodel) #: (*int*) index of stress-permeability model assigned to node.
	def _get_pormodel(self): return self._pormodel
	pormodel = property(_get_pormodel) #: (*int*) index of variable porosity model assigned to node.
	def _get_condmodel(self): return self._condmodel
	condmodel = property(_get_condmodel) #: (*int*) index of variable conductivity model assigned to node.
class fconn(object):				#Connection object.
	"""Connection object, comprising two connected nodes, separated by some distance.

	A connection is associated with a distance between the two nodes.
	"""
	__slots__=['_nodes','_distance','_geom_coef']
	def __init__(self,nodes):
		self._nodes = nodes		
		self._geom_coef = 0.
	def __repr__(self):	return 'n'+str(self._nodes[0]._index)+':n'+str(self._nodes[1]._index)	
	def __getstate__(self):
		return dict((k, getattr(self, k)) for k in self.__slots__)
	def __setstate__(self, data_dict):
		for (name, value) in data_dict.iteritems():
			setattr(self, name, value)
	def _get_distance(self): return self._distance
	distance = property(_get_distance)	#: (*fl64*) Distance between the two connected nodes.
	def _get_nodes(self): return self._nodes
	nodes = property(_get_nodes)#: (*lst[fnode]*) List of node objects (fnode()) that define the connection.
	def _get_geom_coef(self): return self._geom_coef
	geom_coef = property(_get_geom_coef)#: (*fl64*) Geometric coefficient associated with the connection (connected area divided by distance).
class felem(object):				#Element object.
	"""Finite element object, comprising a set of connected nodes.
	
	A finite element is associated with an element centre and an element volume.
	"""
	__slots__=['_index','_nodes']
	def __init__(self,index=None,nodes = []):
		self._index = index			
		self._nodes = nodes			
	def __repr__(self): 
		retStr = 'el'+str(self.index)+': '
		for nd in self.nodes: retStr+='nd'+str(nd.index)+', '
		return retStr
	def __getstate__(self):
		return dict((k, getattr(self, k)) for k in self.__slots__)
	def __setstate__(self, data_dict):
		for (name, value) in data_dict.iteritems():
			setattr(self, name, value)
	def _get_index(self): return self._index
	index = property(_get_index)#: (*int*) Integer number denoting the element.		
	def _get_centre(self):
		c = np.array([0,0,0])
		for nd in self.nodes: c += np.array(nd.position)
		c = c/len(self.nodes)
		return c
	centre = property(_get_centre)	#: (*ndarray*) Coordinates of the element centroid.
	def _get_info(self):
		print 'Element number: '+str(self.index)
		print '      Centroid: '
		print '            x = '+str(self.centre[0])
		print '            y = '+str(self.centre[1])
		print '            z = '+str(self.centre[2])
		prntStr =  'Contains nodes: ['
		for nd in self.nodes: prntStr += str(nd.index) +','
		prntStr = prntStr[:-1]+']'
		print prntStr
	what = property(_get_info)				#: Print to screen information about the element.
	def _get_nodes(self): return self._nodes
	nodes = property(_get_nodes)#: (*lst[fnode]*) List of node objects that define the element.		
	def _get_volume(self):
		return None
	vol = property(_get_volume)		#: (*fl64*) Volume of the element *** NOT DONE ***.		
class octree(object):				#Octree object.
	"""Octree for spatial searching in 3D grids.
	
	The octree is automatically constructed when a new grid file is parsed.
	
	Note: direct interaction with this method is not generally required (or advised!).
	
	The octree is generated recursively, by sub-dividing the model domain into smaller and smaller compartments.
	If a compartment contains more than one point then it continues to sub-divide until there is no more than one point
	per compartment.
	"""
	def __init__(self,bounds,elements,parent=None,repair=False):
		self.parent=parent		#: Parent compartment to which the compartment belongs.
		self.bounds=bounds		#: Bounding box defining the compartment.
		self.elements=elements	#: (*fnode*) Objects contained in the compartment.
		self.child=[] 			#: Children compartments contained within the compartment.
		self.neighbour = [None,None,None,None,None,None]	#: Six element list containing, in order, compartments +x, -x, +y, -y, +z, -z
		if self.parent: 		# if has a parent, is a child, iterate generation, inherit all elements
			self.generation=self.parent.generation+1	#: Generation of the compartment (a measure of how many sub-divisions have occurred.)
			self.all_elements=self.parent.all_elements
		else: 					# if not a parent, generation 0, elements only those existing
			self.generation=0
			self.all_elements=set(elements)
		if self.generation>50:
			if not repair:
				pos0 = self.elements[0].position
				multipleFlag = False
				for elt in self.elements[1:]:
					pos = elt.position
					dist = np.array(pos)-np.array(pos0)
					dist = np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
					if dist == 0: multipleFlag = True; break
				if multipleFlag:
					print 'ERROR: multiple nodes specified at same location. See below for details.'
					elt = self.elements[0]
					print ''
					print 'Node '+str(elt.index)+': '+str(elt.position)
					for elt in self.elements[1:]:
						pos = elt.position
						dist = np.array(pos)-np.array(pos0)
						dist = np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
						if dist == 0: 
							print 'Node '+str(elt.index)+': '+str(elt.position)
					print ''
					print 'Run repair_grid() to attempt fix.'
					return
			else:
				pos0 = self.elements[0].position
				multipleFlag = False
				for elt in self.elements[1:]:
					pos = elt.position
					dist = np.array(pos)-np.array(pos0)
					dist = np.sqrt(dist[0]**2+dist[1]**2+dist[2]**2)
					if dist == 0: multipleFlag = True; break
				if multipleFlag:
					nds = self.elements
					# find node with minimum index - this takes priority
					nd1 = nds[0]
					for nd2 in nds:
						if nd2.index < nd1.index: nd1 = nd2
					nds.remove(nd1)
					# for all other nodes, clean from parent grid
					nd1 = repair.node[nd1.index]					
					for nd2 in nds:
						nd2 = repair.node[nd2.index]
						for el in nd2._elements:
							el._nodes = [nd1 if nd.index == nd2.index else nd for nd in el._nodes]
							
						# remove node from nodelist
						ind = repair._nodelist.index(nd2)
						repair._nodelist.remove(nd2)
							
						print '  ',nd2,' duplicated ',nd1
					self.elements = [nd1]
		if self.num_elements>1: 	# if more than one element in a cube, sub-divide
			cubes=sub_cubes(self.bounds)
			# update neighbour locations
			cube_elements=[[],[],[],[],[],[],[],[]]
			for elt in self.elements:
				for icube,cube in enumerate(cubes):
					if in_cube(elt.position,cube):
						cube_elements[icube].append(elt)
						break
			nbrs=[]
			for cube,elts in zip(cubes,cube_elements):
				if len(elts)>0: 
					if repair:
						self.child.append(octree(cube,elts,self,repair))
					else:
						self.child.append(octree(cube,elts,self))
					nbrs.append(1)
				else: nbrs.append(0)
			self.assign_neighbours(nbrs)
	def __repr__(self): return self.bounds.__repr__()
	def assign_neighbours(self,neighbours):
		cube_inds = []
		for cube_ind,exist in enumerate(neighbours):
			if exist: cube_inds.append(cube_ind)
		for child,update in zip(self.child,cube_inds):
			child.neighbour = deepcopy(self.neighbour)
			posI,nbrs = nbr_update[update]
			for pos,nbr in zip(posI,nbrs):
				if nbr in cube_inds:
					child.neighbour[pos] = self.child[cube_inds.index(nbr)]
			
	def get_num_elements(self): return len(self.elements)
	num_elements=property(get_num_elements)	#: (*int*) Number of objects in the compartment.
	def get_num_children(self): return len(self.child)
	num_children=property(get_num_children) #: (*int*) Number of sub-compartments in the compartment.
	def leaf(self,pos): 							# find leaf corresponding to position
		if in_cube(pos,self.bounds):
			for child in self.child:
				childleaf=child.leaf(pos)
				if childleaf: return childleaf
			return self
		else: return None
	def min_dist(self,pos):
		min_dist = 1.e10
		if self.child:
			for child in self.child:
				dist,cube = child.min_dist(pos)
				if dist<min_dist: min_dist = dist; min_cube = cube
		else:
			min_dist = (np.sqrt((self.elements[0].position[0]-pos[0])**2+
				(self.elements[0].position[1]-pos[1])**2+
				(self.elements[0].position[2]-pos[2])**2))
			min_cube = self
		#print min_cube.elements[0].index, ' ', min_dist
		return min_dist,min_cube
class fgrid(object):				#Grid object.
	""" FEHM grid object.
	
	"""
	__slots__ = ['_nodelist','_node','_connlist','_conn','_elemlist','_elem','_octree','_dimensions','_parent',
		'_full_connectivity','_path','_stor','_pos_matrix']
	def __init__(self,full_connectivity=True):
		self._nodelist=[]			
		self._node={}				
		self._connlist=[]			
		self._conn={}				
		self._elemlist=[]			
		self._elem={}				
		self._octree=None			
		self._dimensions=3
		self._parent = None
		self._full_connectivity = full_connectivity
		self._path = fpath(parent=self)	
		self._stor = fpath(parent=self)		
		self._pos_matrix = None
	def __repr__(self): 
		if self.filename == None:
			return 'no grid'
		else:
			return self.filename			#Print out details
	def __getstate__(self):
		return dict((k, getattr(self, k)) for k in self.__slots__)
	def __setstate__(self, data_dict):
		for (name, value) in data_dict.iteritems():
			setattr(self, name, value)
	def read(self,gridfilename,full_connectivity=True,octree=False,storfilename=None): 
		"""Read data from an FEHM or AVS grid file. If an AVS grid is specified, PyFEHM will write out the corresponding FEHM grid file.

		:param gridfilename: name of grid file, including path specification.
		:type gridfilename: str
		:param full_connectivity: read element and connection data and construct corresponding objects. Defaults to False. Use if access to connectivity information will be useful.
		:type full_connectivity: bool
		:param octree: flag to use octree search algorithm for finding node locations (default = False).
		:type octree: bool
		:param storfilename: name of optional stor file, including path specification.
		:type storfilename: bool
		"""
		self._full_connectivity = full_connectivity
		self._path.filename = gridfilename 
		if not os.path.isfile(gridfilename):
			pyfehm_print('ERROR: file at '+self._path.full_path+' not found.')
			return
		if storfilename:
			self._stor.filename = storfilename 
			if not os.path.isfile(storfilename):
				pyfehm_print('ERROR: file at '+self._stor.full_path+' not found.')
				return
		# assess format of file	from first two lines
		infile = open(self._path.full_path)
		ln1 = infile.readline()
		ln2 = infile.readline()
		infile.close()
		isAvs = False
		try:
			a = int(ln2.split()[0]) 	# first number of second line is a 1, five numbers in header
			if a  == 1 and len(ln1.strip().split())==5: isAvs = True
		except ValueError: pass
		if ln1.strip().split()[0].startswith('coor'):
			self._read_fehm(storfilename=storfilename) 	# first line starts with 'coor'
		elif isAvs:
			self._read_avs()		
			newgridfilename = self._path.full_path.split('.')[:-1]
			newgridfilename = string.join(newgridfilename,'.')+'.inp'
			if os.path.isfile(newgridfilename):
				newgridfilename = newgridfilename[:-4] + '_new' + newgridfilename[-4:]
			self.write(newgridfilename, 'fehm')		# write out equivalent fehm grid
			if self._parent: self._parent.files.grid = self._path.filename
		else:
			pyfehm_print('ERROR: Unrecognized grid format.')
			
		if octree: self.add_nodetree()
		else: self._pos_matrix = np.array([nd._position for nd in self._nodelist])
		if self._parent: self._parent._add_boundary_zones()
	def _read_fehm(self,storfilename=None): 		#Read in fehm meshfile for node,element data .
		self._nodelist = []
		self._connlist = []
		self._elemlist = []
		infile = open(self._path.full_path)
		
		if self._parent: self._parent.files.grid = self._path.filename
			
		ln = infile.readline()
		N = int(infile.readline())
		self._nodelist = [None]*N
		for i in range(N):
			nd = infile.readline().strip().split()
			new_node = fnode(int(nd[0]),np.array([float(nd[1]),float(nd[2]),float(nd[3])]))
			self.add_node(new_node)	
		
		infile.readline()
		infile.readline()
		N = infile.readline()
		connectivity = int(N.strip().split()[0])
		N = int(N.strip().split()[1])
		self._elemlist = [None]*N
		
		lns = infile.readlines()
		
		for i,ln in enumerate(lns[:N]):
			el = [int(eli) for eli in ln.strip().split()]
			
			if self._full_connectivity:
				new_elem = felem(index = el[0], nodes = [self._nodelist[eli-1] for eli in el[1:]])
				self.add_elem(new_elem)
				if connectivity == 8:
					nds1 = [el[1],el[1],el[1],el[7],el[7],el[7],el[4],el[4],el[6],el[6],el[2],el[5]]
					nds2 = [el[5],el[2],el[4],el[6],el[3],el[8],el[3],el[8],el[5],el[2],el[3],el[8]]
				elif connectivity == 4:
					nds1 = [el[1],el[2],el[3],el[4]]
					nds2 = [el[2],el[3],el[4],el[1]]
				elif connectivity == 3:
					nds1 = [el[1],el[2],el[3]]
					nds2 = [el[2],el[3],el[1]]
				if storfilename is None:					
					for ndi1,ndi2 in izip(nds1,nds2):					
						# check if connection exists					
						if ndi1>ndi2: ndi = ndi2; ndi2 = ndi1; ndi1 = ndi
						nd1 = self._nodelist[ndi1-1]; nd2 = self._nodelist[ndi2-1]
						if nd2 in nd1._connected_nodes: continue
						
						nd1._connected_nodes.append(nd2)
						nd2._connected_nodes.append(nd1)
						
						self.add_conn(fconn([nd1,nd2]))
						nd1._connections.append(self.connlist[-1])
						nd2._connections.append(self.connlist[-1])
						
				# associate element nodes with element
				for nd in self._elemlist[i]._nodes: 
					self._nodelist[nd._index-1]._elements.append(self._elemlist[i])		
			else:
				self._elemlist[i] = el[1:]
				self._elem[el[0]] = self.elemlist[i]
		self._conn_dist()
		infile.close()		
		if not storfilename is None:
			self._read_stor(storfilename)
		if ((len(np.unique([nd._position[0] for nd in self._nodelist])) == 1) or
			 (len(np.unique([nd._position[1] for nd in self.nodelist])) == 1) or
			 (len(np.unique([nd._position[2] for nd in self.nodelist])) == 1)):
			self._dimensions = 2
	def _read_stor(self,filename):
		# Define function to parse values in file one at a time
		self._full_connectivity = True
		def valgen(fhandle):
			for line in fhandle:
				for v in line.split():
					yield v
		storfile = open(filename)
		storfile.readline()
		storfile.readline()
		# Ncons: number of written connections
		# Nnds: number of equations (nodes)
		# Nptrs: number of pointers
		# Nac: format of area coefficients: 1-scalar, 3-vector, 4-vector followed by scalar
		# NconMax: Max number of connections for a single node
		Ncons, Nnds, Nptrs, Nac, NconMax = map(int,storfile.readline().split())
		Ncoef = Nptrs - Nnds - 1 # Number of connections in mesh
		if Nac != 1:
			pyfehm_print('ERROR: Vector or vector/scalar coefficients are not currently supported - stor will not be read in')
			storfile.close()
			return
		# Read in volumes
		nextval = valgen(storfile)
		for i in range(Nnds):
			self.nodelist[i]._vol = float(nextval.next())
		sparse = np.zeros(Nnds+1)
		# Read in sparse - sparse[i+1] - sparse[i] is the number of connections for node in nodelist[i] 
		for i in range(Nnds+1):
			sparse[i] = int(nextval.next())
		# Calculate number of connections per node (nconns) using sparse
		s2 = np.delete(sparse,0).astype(int)	
		s1 = np.delete(sparse,-1).astype(int)	
		nodeconns = s2 - s1
		# Create connections betwween nodes
		self._connlist = []
		self._conn = {}
		cindexlist = [] # Create index of connections since not all connections in stor file are created 
		cindex = 0
		for i,nc in enumerate(nodeconns):
			for j in range(nc):
				ndid = int(nextval.next())
				if ndid > self.nodelist[i].index: 
					new_conn = fconn(nodes= [self.nodelist[i],self.nodelist[ndid-1]])
					self.add_conn(new_conn)
					self.nodelist[ndid-1].connections.append(new_conn)
					self.nodelist[i].connections.append(new_conn)
					cindexlist.append(cindex)
				cindex+=1
		# Collect pointers into coefficient array
		ptrs = np.zeros(Ncoef,dtype=np.int)
		for i in range(Ncoef):
			ptrs[i] = int(nextval.next())
		# Read past padded zeros
		for i in range(Nnds+1):
			dum = nextval.next()
		# Read diagonal pointers, not needed here
		dptrs = np.zeros(Nnds,dtype=np.int)
		for i in range(Nnds):
			dptrs[i] = nextval.next()
		# Read in coefficients
		coeffs = np.zeros(Ncons*Nac)
		for i in range(len(coeffs)):
			coeffs[i] = float(nextval.next())
		storfile.close()
		# Set coefficients in connections
		for ptr,c in zip(ptrs[cindexlist],self.connlist):
			c._geom_coef = coeffs[ptr-1]
	def _read_avs(self): 		#Read in avs meshfile for node, element data.
		self._nodelist = []
		self._connlist = []
		self._elemlist = []
		infile = open(self._path.full_path)
		
		if self._parent: self._parent.files.grid = self._path.filename
			
		ln = infile.readline()
		N = int(ln.strip().split()[0]) 			# number nodes
		N_el = int(ln.strip().split()[1])		# number elements
		for i in range(N): 						# FOR each node
			nd = infile.readline().strip().split()	# read line
			new_node = fnode(index=int(nd[0]),position=np.array([float(nd[1]),float(nd[2]),float(nd[3])]))
			self.add_node(new_node)	 				# add node object
		
		N = N_el
		connectivity = None
		for i in range(N): 						# FOR each element
			el = infile.readline().strip().split()
			el = [el[0],] + el[3:]
			el = [int(eli) for eli in el]
			if connectivity == None: connectivity = len(el) - 1
			
			if self._full_connectivity:
				new_elem = felem(index = el[0], nodes = [self.node[eli] for eli in el[1:]])
				self.add_elem(new_elem)
				if connectivity == 8:
					nds1 = [el[1],el[1],el[1],el[7],el[7],el[7],el[4],el[4],el[6],el[6],el[2],el[5]]
					nds2 = [el[5],el[2],el[4],el[6],el[3],el[8],el[3],el[8],el[5],el[2],el[3],el[8]]
				elif connectivity == 4:
					nds1 = [el[1],el[2],el[3],el[4]]
					nds2 = [el[2],el[3],el[4],el[1]]
				elif connectivity == 3:
					nds1 = [el[1],el[2],el[3]]
					nds2 = [el[2],el[3],el[1]]
				else:
					pyfehm_print('ERROR: unrecognized connectivity')
					return
				for nd1,nd2 in zip(nds1,nds2):
					if nd1>nd2: ndi = nd2; nd2 = nd1; nd1 = ndi
					nd1 = self.node[nd1]; nd2 = self.node[nd2]
					new_conn = fconn(nodes = [nd1,nd2])
					
					nd1inds = []
					for con in nd1.connections:
						for nd in con.nodes: nd1inds.append(nd.index)
					
					if nd2.index in nd1inds: continue
					
					self.add_conn(new_conn)
					nd1.connections.append(new_conn)
					nd2.connections.append(new_conn)
				# associate element nodes with element
				for nd in self.elemlist[-1].nodes: 
					self._node[nd.index].elements.append(self.elemlist[-1])		
			else:
				self._elemlist.append(el[1:])
				self._elem[el[0]] = self.elemlist[-1]
		infile.close()		
		if ((len(np.unique([nd.position[0] for nd in self.nodelist])) == 1) or
			 (len(np.unique([nd.position[1] for nd in self.nodelist])) == 1) or
			 (len(np.unique([nd.position[2] for nd in self.nodelist])) == 1)):
			self._dimensions = 2
		else: self._dimensions = 3
	def write(self,filename=None,format='fehm', remove_duplicates = True, recalculate_coefficients=True,compress_eps=1.e-5):
		"""Write grid object to a grid file (FEHM, AVS STOR file formats supported). Stor file support only for orthogonal hexahedral grids.

		:param filename: name of FEHM grid file to write to, including path specification, e.g. c:\\path\\file_out.inp
		:type filename: str
		:param format: FEHM grid file format ('fehm','avs','stor'). Defaults to 'fehm' unless filename is passed with extension '*.avs' or '*.stor'.
		:type format: str
		:param remove_duplicates: Remove duplicate entries in the stor file. Improves simulation run time, particularly for structured grids.
		:type remove_duplicates: bool
		:param compress_eps: Float used to determine geometric coefficients to remove (Coefficients less than max(geom_coef)*compress_eps will be removed)
		:type compress_eps: fl64
		:param keepcons: List of connections to keep in stor file no matter what
		:type keepcons: lst(connection keys)

		"""
		
		# autodetect format
		if filename:
			extension = filename.split('.')[-1]
			if extension == 'avs': format = 'avs'
			elif extension == 'stor': format = 'stor'
		
		if format == 'stor' and not self._full_connectivity:
			pyfehm_print('ERROR: STOR file cannot be written without full connectivity information. Read or create grid file with flag full_connectivity=True')
			return

		if format == 'stor': old_path = copy(self._path)		# save path
		if filename: self._path.filename=filename
		if filename == None: self._path.filename='default_GRID.inp'
		
		temp_path = fpath()
		temp_path.filename = filename
		
		path = temp_path.full_path
		try: os.makedirs(self._path.absolute_to_file)
		except: pass
		if self._parent:
			if self._parent.work_dir:
				try: os.makedirs(self._path.absolute_to_workdir)
				except:	pass
				path = self._path.absolute_to_workdir+os.sep+temp_path.filename
				
		outfile = open(path,'w')
		
		if format == 'fehm': self._write_fehm(outfile)
		elif format == 'avs': self._write_avs(outfile)
		elif format == 'stor': self._write_stor(outfile,remove_duplicates,recalculate_coefficients,compress_eps=compress_eps)
		else: pyfehm_print('ERROR: Unrecognized format '+format+'.');return
		
		if format == 'stor': 
			if not self._parent is None:
				self._path = old_path
				self._parent.ctrl['stor_file_LDA'] = 1
				self._parent.files.stor = path
	def _write_fehm(self,outfile):
			
		outfile.write('coor\n')
		outfile.write('   '+str(len(self.nodelist))+'\n')		
		for nd in self.nodelist: 
			outfile.write('%11d' % nd.index +'        ')
			outfile.write('%14.8f' % nd.position[0]+'        ')
			outfile.write('%14.8f' % nd.position[1]+'        ')
			outfile.write('%14.8f' % nd.position[2])
			outfile.write('\n')			
		outfile.write('\t0\n')
		outfile.write('elem\n')
		if self._full_connectivity:
			outfile.write(str(len(self.elemlist[0].nodes))+' '+str(len(self.elemlist))+'\n')	
			for el in self.elemlist:
				outfile.write(str(int(el.index))+'   ')
				for nd in el.nodes:
					outfile.write(str(nd.index)+'   ')
				outfile.write('\n')	
		else:
			outfile.write(str(len(self.elemlist[0]))+' '+str(len(self.elemlist))+'\n')	
			for i,el in enumerate(self.elemlist):
				outfile.write(str(i+1)+'   ')
				for nd in el:
					outfile.write(str(nd)+'   ')
				outfile.write('\n')		
			
		outfile.write('\nstop\n')
		outfile.close()
	def _write_avs(self,outfile):
		outfile.write(str(self.number_nodes)+' '+str(self.number_elems)+' 0 0 0\n')
		for nd in self.nodelist: 
			outfile.write('%11d' % nd.index +'        ')
			outfile.write('%14.8f' % nd.position[0]+'        ')
			outfile.write('%14.8f' % nd.position[1]+'        ')
			outfile.write('%14.8f' % nd.position[2])
			outfile.write('\n')
			
		if self._full_connectivity:			
			for el in self.elemlist:
				outfile.write(str(int(el.index))+'   1  ')
				if len(el.nodes) == 8: outfile.write('hex  ')
				for nd in el.nodes:
					outfile.write(str(nd.index)+'   ')
				outfile.write('\n')
		else:		
			for i,el in enumerate(self.elemlist):
				outfile.write(str(i+1)+'   1  ')
				if len(el) == 8: outfile.write('hex  ')
				for nd in el:
					outfile.write(str(nd)+'   ')
				outfile.write('\n')
		outfile.close()
	def _write_stor(self,outfile,remove_duplicates,recalculate_coefficients,compress_eps=1.e-5):
		# calculate volumes and geometric coefficients
		if recalculate_coefficients:
			for conn in self.connlist: conn._geom_coef = None
			for nd in self.nodelist: nd._vol = None
			for nd in self.nodelist:
				self._vorvol(nd) 		# calculate voronoi volume of node
		
		if remove_duplicates:
			tol = compress_eps
			indices1 = [] 			# for populating later
			geom_coefs = [con.geom_coef for con in self.connlist] 	# all coefficients
			gcU = np.ones((1,len(geom_coefs)))[0]*float('Inf')
			NgcU = 0
			for gc in geom_coefs:
				new_gcU = True
				if gc == 0:
					dgcu = abs(gc-gcU[:NgcU]) 		# distances between geom_coef and existing bins
				else:
					dgcu = abs((gc-gcU[:NgcU])/gc) 		# distances between geom_coef and existing bins
				if all(dgcu>tol): 
					gcU[NgcU] = gc
					NgcU += 1
			gcU[NgcU] = 0.
			gcU = gcU[:NgcU+1]
			gcU.sort()
			pyfehm_print("Number of compressed connections: "+ str(len(gcU)))
			
		# write file
		outfile.write('fehmstor ascir8i4 PyFEHM Sparse Matrix Voronoi Coefficients\n')
		import time
		lt = time.localtime()
		mon = str(lt.tm_mon)
		if len(mon) == 1: mon = '0'+mon
		day = str(lt.tm_mday)
		if len(day) == 1: day = '0'+day
		yr = str(lt.tm_year)
		hr = str(lt.tm_hour)
		if len(hr) == 1: hr = '0'+hr
		min = str(lt.tm_min)
		sec = str(lt.tm_sec)
		outfile.write(mon+'/'+day+'/'+yr+'    '+hr+':'+min+':'+sec+'\n')
		
		# write matrix parameters
		Nnds = len(self.nodelist)
		Ncons = 2*len(self.connlist)+Nnds
		if remove_duplicates:
			outfile.write('\t%5i'%len(gcU))
		else:
			outfile.write('\t%5i'%Ncons)
		outfile.write('\t%5i'%Nnds)
		outfile.write('\t%5i'%(Ncons+Nnds+1))
		outfile.write('\t%5i'%1)
		outfile.write('\t%5i'%7)
		outfile.write('\n')

		# write Voronoi Volumes
		cnt = 0
		for nd in self.nodelist:
			outfile.write('  %13.12E'%nd.vol)
			cnt +=1
			if cnt ==5: 
				outfile.write('\n')
				cnt = 0			
		if cnt != 0: outfile.write('\n')
				
		# write sparse matrix connectivity
		sparse = np.zeros((1,Ncons+Nnds+1))[0]
		gcStr = ''
		cnt = 0
		stride = Nnds+1
		diagonals = []
		for i, nd in enumerate(self.nodelist):
			# populate sparse
			# add the stride
			sparse[i] = stride
			diagonals.append(stride+1)
			# add the connections
			cons = []
			for con in nd.connections:
				if con.nodes[0].index == nd.index: cons.append([con.nodes[1].index,con.geom_coef])
				else: cons.append([con.nodes[0].index,con.geom_coef])
				if cons[-1][1] == None: cons[-1][1] = 0.
			cons.append([nd.index,0.])
			cons.sort(key=lambda x: x[0])
			for j,con in enumerate(cons):
				if con[0] == nd.index: diagonals[-1] += j
				# add connection to sparse
				sparse[stride] = con[0]
				stride +=1 	# increment the stride
				# write geometric coefficients
				gcStr +='  %13.12E'%con[1]
				if remove_duplicates: indices1.append(abs(gcU-con[1]).argmin()+1)
				cnt +=1
				if cnt ==5: 
					gcStr+='\n'
					cnt = 0			
		if cnt != 0: gcStr+='\n'
		sparse[i+1] = stride
		
		# row entries
		cnt = 0
		for i in sparse:
			outfile.write('  %6i'%i)
			cnt +=1
			if cnt ==5: 
				outfile.write('\n')
				cnt = 0			
		if cnt != 0: outfile.write('\n')
				
		# indices into coefficient list
		if not remove_duplicates: indices1 = range(1,Ncons+1)
		
		indices = indices1+list(np.zeros((1,Nnds+1))[0])
		cnt = 0
		for i in indices:
			outfile.write('  %6i'%i)
			cnt +=1
			if cnt ==5: 
				outfile.write('\n')
				cnt = 0			
		if cnt != 0: outfile.write('\n')
		cnt = 0
		for i in diagonals:
			outfile.write('  %6i'%i)
			cnt +=1
			if cnt ==5: 
				outfile.write('\n')
				cnt = 0			
		if cnt != 0: outfile.write('\n')
				
		# geometric area coefficient values
		if remove_duplicates:
			cnt = 0
			for i in gcU:
				outfile.write('  %13.12E'%i)
				cnt +=1
				if cnt ==5: 
					outfile.write('\n')
					cnt = 0			
			if cnt != 0: outfile.write('\n')
		else:
			outfile.write(gcStr)
		outfile.close()
	def _vorvol(self,nd):
		# array of connection mid-point positions
		pos = np.array([(con.nodes[0].position+con.nodes[1].position)/2. for con in nd.connections])
		# min and max of these mid-points, use to find bounding box of volume
		xm,ym,zm,xM,yM,zM = np.min(pos[:,0]),np.min(pos[:,1]),np.min(pos[:,2]),np.max(pos[:,0]),np.max(pos[:,1]),np.max(pos[:,2])
		# calculate volume
		nd._vol = np.prod([xM-xm,yM-ym,zM-zm])
		areas = [-(yM-ym)*(zM-zm),-(xM-xm)*(zM-zm),-(yM-ym)*(xM-xm)] 	# interface areas
		
		for con in nd.connections:
			if con.geom_coef != None: continue
			# establish direction of connection
			N = abs(con.nodes[0].position - con.nodes[1].position)
			if N[0]>N[1] and N[0]>N[2]:	
				#con._geom_coef = [areas[0]/(con.distance/2.),]
				con._geom_coef = areas[0]/con.distance
			elif N[1]>N[0] and N[1]>N[2]:	
				#con._geom_coef = [areas[1]/(con.distance/2.),]
				con._geom_coef = areas[1]/con.distance
			elif N[2]>N[1] and N[2]>N[0]:	
				#con._geom_coef = [areas[2]/(con.distance/2.),]
				con._geom_coef = areas[2]/con.distance
	def remove_zeros(self, tolerance = 1.e-8, keepcons = []):
		""" Removes node connections with geometric coefficients smaller than a supplied tolerance.
		
		:param tolerance: Relative tolerance below which connections will be removed (default = 1.e-8).
		:type tolerance: fl64
		:param keepcons: List of connections to keep in stor file ignoring tolerance.
		:type keepcons: lst(connection keys) or lst(connections)
		"""
		
		connlist = self.connlist # Save temporary connlist
		conndict = self.conn
		geom_coefs = [con.geom_coef for con in self.connlist] 	# all coefficients
		self._conn = {}
		self._connlist = []		 # Erase connlist
		mincoef = np.min(geom_coefs) * tolerance # Determine minimum allowed coef
		pyfehm_print("Minimum coefficient value to keep:"+str(mincoef))
		pyfehm_print("Number of connections (incl. zeros):"+str(len(geom_coefs)))
		# Collect connections to keep (less than because coef are given negative sign)
		for i,j in enumerate(keepcons):
			if isinstance(j,fconn):
				keepcons[i] = (j.nodes[0].index,j.nodes[1].index)
			elif isinstance(j,tuple): pass
			else:
				print "ERROR: keepconns contains unknown connection type"
				return
		keepcoef = [c for c in connlist if (c.geom_coef < mincoef or (c.nodes[0].index,c.nodes[1].index) in keepcons)]
		for con in keepcoef:
				self.add_conn(con)
		for nd in self.nodelist: nd._connections = [] 		# reassociate connections with nodes
		for con in self.connlist:
			con.nodes[0].connections.append(con)
			con.nodes[1].connections.append(con)
		pyfehm_print("Number of connections (zeros removed):"+str(len(self.connlist)))
	def make(self,gridfilename,x,y,z,full_connectivity=True,octree=False,radial=False):
		""" Generates an orthogonal mesh for input node positions. 
		
		The mesh is constructed using the ``fgrid.``\ **fmake** object and an FEHM grid file is written for the mesh.
		
		:param gridfilename: Name to which to save the grid file.
		:type gridfilename: str
		:param x: Unique set of x-coordinates.
		:type x: list[fl64]
		:param y: Unique set of y-coordinates.
		:type y: list[fl64]
		:param z: Unique set of z-coordinates.
		:type z: list[fl64]
		:param full_connectivity: read element and connection data and construct corresponding objects. Defaults to False. Use if access to connectivity information will be useful.
		:type full_connectivity: bool
		:param radial: Creates a radial grid (ignores z coordinate)
		:type radial: bool
		"""
		# ASSUMPTIONS FOR MAKE
		# if no parent - interpret string literally
		# if parent but NO work dir - interpret string literally
		# if parent and work dir but relative or absolute - interpret string relative to cwd
		# if parent and work dir and grid name ONLY - grid goes in work dir
		
		# PASS FULL PATH specification into fmake
		
		temp_path = fpath()
		temp_path.filename = gridfilename
		
		path = temp_path.full_path
		if self._parent:
			if self._parent.work_dir:
				path = self._path.absolute_to_workdir+os.sep+temp_path.filename
		
		fm = fmake(path,x,y,z)
		fm.write(radial=radial)		
		self.read(path,full_connectivity,octree)
	def lagrit_stor(self, grid = None, stor = None, exe = dflt.lagrit_path, overwrite = False):
		"""Uses LaGriT to create a stor file for the simulation, this will be used in subsequent runs.
		To create the stor file, LaGriT will convert a mesh comprised of hexahedral elements into one comprising
		only tetrahedrals. Therefore, a new FEHM grid file will be created and parsed, reflecting the modified 
		element structure.
		
		:param grid: Name of grid file to be created. Destination directory supported.
		:type grid: str
		:param stor: Name of stor file to be created. Destination directory supported.
		:type stor: str
		:param exe: Path to lagrit executable (default, fdflt.lagrit_path, in environment file 'fdflt.py').
		:type exe: str
		:param overwrite: Flag to request that the new FEHM grid file overwrites the old one.
		:type overwrite: bool
		
		"""
		if self.filename == None: 	
			pyfehm_print('ERROR: a grid must first be loaded/created before a stor file can be created.')	
			return
		if not os.path.isfile(exe): 
			pyfehm_print('ERROR: cannot find lagrit executable at '+exe)
			return
		
		# assign default stor file name if none given
		if stor == None:
			stor = self._path.full_path.split('.')[:-1]
			stor = string.join(stor,'.')+'.stor'
		# assign default grid file name, location if none given
		if grid == None:
			grid = self._path.full_path.split('.')[:-1]
			grid = string.join(grid,'.')+'_lg.'+self.filename.split('.')[-1]
		# write mesh to avs
		avs = self.filename.split('.')[:-1]
		avs = string.join(avs,'.')+'.avs'
		fname_save = self._path.full_path
		self.write(filename = avs, format = 'avs')
		self._path.filename = fname_save
		
		cwd = os.getcwd()
		if self._parent:
			if self._parent.work_dir:
				os.chdir(self._parent.work_dir)
				
		# create lagrit instruction file
		fp = open('lagrit_instructions.lgi','w')
		lns = []
		lns.append('# read the mesh in AVS format\n')
		lns.append('read avs '+avs+' cmohex\n')
		lns.append('\n')
		lns.append('# create a tet version of this point distribution\n')
		lns.append('cmo / create / cmotet / / / tet\n')
		lns.append('copypts / cmotet / cmohex\n')
		lns.append('  cmo / select cmotet\n')
		lns.append('  cmo / setatt / cmotet / imt / 1 0 0 / 1\n')
		lns.append('  cmo / setatt / cmotet / itp / 1 0 0 / 0\n')
		lns.append('\n')
		lns.append('# remove duplicates\n')
		lns.append('filter / 1 0 0\n')
		lns.append('rmpoint / compress\n')
		lns.append('\n')
		lns.append('# if sort is needed to reorder the nodes\n')
		lns.append('# otherwise keep same node ordering and comment out\n')
		lns.append('#sort / cmotet / index / ascending / ikey / zic yic xic\n')
		lns.append('#reorder / cmotet / ikey\n')
		lns.append('#cmo / DELATT / cmotet / ikey\n')
		lns.append('#cmo / status\n')
		lns.append('\n')
		lns.append('# connect into delaunay tet mesh \n')
		lns.append('connect noadd\n')
		lns.append('resetpts / itp\n')
		lns.append('quality\n')
		lns.append('\n')
		lns.append('# write tet and stor files\n')
		lns.append('dump fehm lagrit_out cmotet ascii \n')
		lns.append('\n')
		lns.append('\n')
		lns.append('finish\n')
		fp.writelines(lns)
		fp.close()
		
		os.system(exe+' < lagrit_instructions.lgi > nul')
		# isolate new grid and stor files, delete others
		files = glob('lagrit_out*.*')
		for file in files:
			if not file.endswith('.stor') and not file.endswith('.fehmn'):
				os.remove(file)
		os.remove('lagrit_instructions.lgi')
		os.remove(avs)
				
		# rename and move files, parse
		if overwrite: 	# overwriting old grid file
			shutil.move('lagrit_out.fehmn',self._path.full_path)
			self.read(self._path.full_path)
		elif grid: 			# grid name/location specified
			temp_path = fpath()
			temp_path.filename = grid
			try: os.makedirs(temp_path.absolute_to_file)
			except: pass
			shutil.move('lagrit_out.fehmn',grid)
			self.read(grid)
		
		temp_path = fpath()
		temp_path.filename = stor
		
		#path = temp_path.full_path
		try: os.makedirs(self._path.absolute_to_file)
		except: pass
		if self._parent:
			if self._parent.work_dir:
				try: os.makedirs(self._path.absolute_to_workdir)
				except:	pass
				stor = self._path.absolute_to_workdir+os.sep+temp_path.filename
		
		shutil.move('lagrit_out.stor',stor)
		if self._parent:
			self._parent.files.stor = stor
			self._parent.ctrl['stor_file_LDA'] = 1
			
		os.chdir(cwd)
	def volumes(self,volumefilename):
		""" Reads a lagrit generated file containing control volume information.
		
		:param volumefilename: Name of lagrit output file containing control volume information.
		:type volumefilename: str
		"""
		infile = open(volumefilename,'r')
		line = infile.readline() 
		line = infile.readline().strip().split()
		Ncol = int(line[0])
		colnames = []
		for i in range(Ncol):
			line = infile.readline().strip().split(',')
			colnames.append(line[0])
		keepReading = True
		for i,name in enumerate(colnames):
			if name == 'vorvol': break
		while keepReading:
			line = infile.readline().strip().split()
			if not line: keepReading = False; continue
			ndI = int(line[1])
			ndV = float(line[i+1])
			self.node[ndI]._vol = ndV
	def add_node(self,node):		#Add a node object.
		self._nodelist[node._index-1] = node
		self._node[node._index] = node
	def add_conn(self,conn):		#Add a connection object.
		self._connlist.append(conn)
		self._conn[(conn._nodes[0]._index,conn._nodes[1]._index)] = self._connlist[-1]
	def _conn_dist(self):
		if not self._full_connectivity: return
		# populate distance attribute for each conn object
		dist = np.array([c._nodes[0]._position-c._nodes[1]._position for c in self._connlist])
		dist = np.sqrt((dist[:,0])**2+(dist[:,1])**2+(dist[:,2])**2)			
		for i,c in enumerate(self.connlist): c._distance = dist[i]			
	def add_elem(self,elem):		#Add an element object.
		self._elemlist[elem._index-1] = elem
		self._elem[elem._index] = elem
	def node_nearest_point(self,pos = []):
		"""Return node object nearest to position in space. Method uses octree structure for speed up if available.

		:param pos: Coordinates, e.g. [2300., -134.8, 0.].
		:type pos: list
		:returns:  fnode() -- node object closest to position.

		"""
		if self.octree != None:
			min_dist = 1.e10
			nd = None
			if self.octree.leaf(pos) == None:
				pyfehm_print('point outside domain')
				return None
			lf = self.octree.leaf(pos)
			min_dist= 1.e10
			
			for leaf in ([lf,]+lf.neighbour):
				if leaf == None: continue
				dist,cube = leaf.min_dist(pos)
				if dist<min_dist: min_dist = dist; min_cube = cube
			return min_cube.elements[0]
		else:
			idx = np.abs(np.sum((self._pos_matrix - pos)**2,axis=1)).argmin()
			return self.nodelist[idx]
	def plot(self,save='',angle=[45,45],color='k',connections=False,equal_axes=True,
		xlabel='x / m',ylabel='y / m',zlabel='z / m',title='',font_size='small',cutaway=[],
		zones = []): 		#generates a 3-D plot of the zone.
		"""Generates and saves a 3-D plot of the grid.
		
		:param save: Name of saved zone image.
		:type save: str
		:param angle: 	View angle of zone. First number is tilt angle in degrees, second number is azimuth. Alternatively, if angle is 'x', 'y', 'z', view is aligned along the corresponding axis.
		:type angle: [fl64,fl64], str
		:param color: Colour of zone.
		:type color: str, [fl64,fl64,fl64]
		:param connections: Plot connections. If ``True`` all connections plotted. If between 0 and 1, random proportion plotted. If greater than 1, specified number plotted.
		:type connections: bool
		:param equal_axes: Force plotting with equal aspect ratios for all axes.
		:type equal_axes: bool
		
		:param xlabel: Label on x-axis.
		:type xlabel: str
		:param ylabel: Label on y-axis.
		:type ylabel: str
		:param zlabel: Label on z-axis.
		:type zlabel: str
		:param title: Title of plot.
		:type title: str
		
		:param font_size: Size of text on plot.
		:type font_size: str, int
		
		:param cutaway: Coordinate from which cutaway begins. Alternatively, specifying 'middle','centre' with choose the centre of the grid as the cutaway point.
		:type cutaway: [fl64,fl64,fl64], str
		
		:param zones: List of zone indices or names. If these are defined then nodes contained in those zones are highlighted in the output plot. Maximum of six zones.
		:type zones: int, str
		
		"""
		if not save: save = self.filename.split('.')[0]+'.png'
		if cutaway in ['middle','center','centre','mid']:			
			cutaway = [(self.xmin+self.xmax)/2,(self.ymin+self.ymax)/2,(self.zmin+self.zmax)/2]
		if isinstance(angle,str):
			if angle == 'x': angle = [0,0]
			elif angle == 'y': angle = [0,90]
			elif angle == 'z': angle = [90,90]
			else: return
			face1 = True; face2 = True; face3 = True; face4 = True; face5 = True; face6 = True
		else:
			while angle[0]<-90: angle[0]+=180
			while angle[0]>90: angle[0]-=180
			while angle[1]<0: angle[1]+=180
			while angle[1]>360: angle[1]-=180
			if angle[0]>0: face1 = True; face2 = False
			else: face1 = False; face2 = True
			if angle[1]>270 or angle[1]<=90: face3 = True; face4 = False
			else: face3 = False; face4 = True
			if angle[1]>0 and angle[1]<=180: face5 = True; face6 = False
			else: face5 = False; face6 = True
		# plot bounding box
		plt.clf()
		fig = plt.figure(figsize=[8.275,11.7])
		ax = plt.axes(projection='3d')
		ax.set_aspect('equal', 'datalim')
		
		ax.set_xlabel(xlabel,size=font_size)
		ax.set_ylabel(ylabel,size=font_size)
		ax.set_zlabel(zlabel,size=font_size)
		ax.set_title(title,size=font_size)
		
		for t in ax.get_xticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_yticklabels():
			t.set_fontsize(font_size)
		for t in ax.get_zticklabels():
			t.set_fontsize(font_size)

		xmin,xmax = self.xmin, self.xmax
		ymin,ymax = self.ymin, self.ymax
		zmin,zmax = self.zmin, self.zmax
		
		if equal_axes:
			MAX = np.max([xmax-xmin,ymax-ymin,zmax-zmin])/2
			C = np.array([xmin+xmax,ymin+ymax,zmin+zmax])/2
			for direction in (-1, 1):
				for point in np.diag(direction * MAX * np.array([1,1,1])):
					ax.plot([point[0]+C[0]], [point[1]+C[1]], [point[2]+C[2]], 'w')
		ax.view_init(angle[0],angle[1])
		
		if cutaway: xmid,ymid,zmid = cutaway
		else: 
			if face1:
				if face5:
					if face3: xmid,ymid,zmid = xmax,ymax,zmax
					else: xmid,ymid,zmid = xmin,ymax,zmax
				else:
					if face3:xmid,ymid,zmid = xmax,ymin,zmax
					else:xmid,ymid,zmid = xmin,ymin,zmax
			else:
				if face5:
					if face3: xmid,ymid,zmid = xmax,ymax,zmin
					else: xmid,ymid,zmid = xmin,ymax,zmin
				else:
					if face3:xmid,ymid,zmid = xmax,ymin,zmin
					else:xmid,ymid,zmid = xmin,ymin,zmin
					
		p13 = [xmid,ymid,zmid]
		if face1:
			if face5:
				if face3:					
					p1=[xmin,ymin,zmax]; p2=[xmin,ymax,zmax]; p3=[xmin,ymax,zmin]; p4=[xmax,ymax,zmin]
					p5=[xmax,ymin,zmin]; p6=[xmax,ymin,zmax]; p7=[xmax,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymax,zmax]; p10=[xmid,ymax,zmid]; p11=[xmax,ymax,zmid];p12=[xmax,ymid,zmid]
				else:						
					p1=[xmax,ymin,zmax]; p2=[xmax,ymax,zmax]; p3=[xmax,ymax,zmin]; p4=[xmin,ymax,zmin]
					p5=[xmin,ymin,zmin]; p6=[xmin,ymin,zmax]; p7=[xmin,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymax,zmax]; p10=[xmid,ymax,zmid]; p11=[xmin,ymax,zmid];p12=[xmin,ymid,zmid]
			else:
				if face3:						
					p1=[xmin,ymax,zmax]; p2=[xmin,ymin,zmax]; p3=[xmin,ymin,zmin]; p4=[xmax,ymin,zmin]
					p5=[xmax,ymax,zmin]; p6=[xmax,ymax,zmax]; p7=[xmax,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymin,zmax]; p10=[xmid,ymin,zmid]; p11=[xmax,ymin,zmid];p12=[xmax,ymid,zmid]
				else:		
					p1=[xmax,ymax,zmax]; p2=[xmax,ymin,zmax]; p3=[xmax,ymin,zmin]; p4=[xmin,ymin,zmin]
					p5=[xmin,ymax,zmin]; p6=[xmin,ymax,zmax]; p7=[xmin,ymid,zmax]; p8=[xmid,ymid,zmax]
					p9=[xmid,ymin,zmax]; p10=[xmid,ymin,zmid]; p11=[xmin,ymin,zmid];p12=[xmin,ymid,zmid]
		else:
			if face5:
				if face3:					
					p1=[xmin,ymin,zmin]; p2=[xmin,ymax,zmin]; p3=[xmin,ymax,zmax]; p4=[xmax,ymax,zmax]
					p5=[xmax,ymin,zmax]; p6=[xmax,ymin,zmin]; p7=[xmax,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymax,zmin]; p10=[xmid,ymax,zmid]; p11=[xmax,ymax,zmid];p12=[xmax,ymid,zmid]
				else:						
					p1=[xmax,ymin,zmin]; p2=[xmax,ymax,zmin]; p3=[xmax,ymax,zmax]; p4=[xmin,ymax,zmax]
					p5=[xmin,ymin,zmax]; p6=[xmin,ymin,zmin]; p7=[xmin,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymax,zmin]; p10=[xmid,ymax,zmid]; p11=[xmin,ymax,zmid];p12=[xmin,ymid,zmid]
			else:
				if face3:						
					p1=[xmin,ymax,zmin]; p2=[xmin,ymin,zmin]; p3=[xmin,ymin,zmax]; p4=[xmax,ymin,zmax]
					p5=[xmax,ymax,zmax]; p6=[xmax,ymax,zmin]; p7=[xmax,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymin,zmin]; p10=[xmid,ymin,zmid]; p11=[xmax,ymin,zmid];p12=[xmax,ymid,zmid]
				else:		
					p1=[xmax,ymax,zmin]; p2=[xmax,ymin,zmin]; p3=[xmax,ymin,zmax]; p4=[xmin,ymin,zmax]
					p5=[xmin,ymax,zmax]; p6=[xmin,ymax,zmin]; p7=[xmin,ymid,zmin]; p8=[xmid,ymid,zmin]
					p9=[xmid,ymin,zmin]; p10=[xmid,ymin,zmid]; p11=[xmin,ymin,zmid];p12=[xmin,ymid,zmid]
		pt1=p1;pt2=p2;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p9;pt2=p2;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p3;pt2=p2;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p3;pt2=p4;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p11;pt2=p4;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p5;pt2=p4;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p5;pt2=p6;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p1;pt2=p6;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p7;pt2=p6;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p7;pt2=p8;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p7;pt2=p12;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p11;pt2=p12;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p13;pt2=p12;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p13;pt2=p8;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p13;pt2=p10;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p9;pt2=p8;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p9;pt2=p10;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
		pt1=p11;pt2=p10;ax.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],[pt1[2],pt2[2]],'k-')
					
		# minor lines
		xs = np.unique([nd.position[0] for nd in self.nodelist])
		ys = np.unique([nd.position[1] for nd in self.nodelist])
		zs = np.unique([nd.position[2] for nd in self.nodelist])
		
		for x in xs:
			if x>=np.min([p2[0],p9[0]]) and x<=np.max([p2[0],p9[0]]):
				ax.plot([x,x],[p1[1],p2[1]],[p2[2],p2[2]],color=color,linewidth=0.5)
				ax.plot([x,x],[p2[1],p2[1]],[p2[2],p3[2]],color=color,linewidth=0.5)
			else:
				ax.plot([x,x],[p12[1],p2[1]],[p10[2],p10[2]],color=color,linewidth=0.5)
				ax.plot([x,x],[p12[1],p6[1]],[p2[2],p2[2]],color=color,linewidth=0.5)				
				ax.plot([x,x],[p7[1],p7[1]],[p7[2],p11[2]],color=color,linewidth=0.5)
				ax.plot([x,x],[p11[1],p11[1]],[p11[2],p3[2]],color=color,linewidth=0.5)
		for y in ys:
			if y>=np.min([p6[1],p7[1]]) and y<=np.max([p6[1],p7[1]]):
				ax.plot([p6[0],p1[0]],[y,y],[p2[2],p2[2]],color=color,linewidth=0.5)				
				ax.plot([p6[0],p6[0]],[y,y],[p2[2],p3[2]],color=color,linewidth=0.5)
			else:             
				ax.plot([p10[0],p11[0]],[y,y],[p10[2],p10[2]],color=color,linewidth=0.5)
				ax.plot([p9[0],p2[0]],[y,y],[p2[2],p2[2]],color=color,linewidth=0.5)							
				ax.plot([p4[0],p4[0]],[y,y],[p4[2],p11[2]],color=color,linewidth=0.5)
				ax.plot([p10[0],p10[0]],[y,y],[p10[2],p9[2]],color=color,linewidth=0.5)
		for z in zs:
			if z>=np.min([p4[2],p11[2]]) and z<=np.max([p4[2],p11[2]]):
				ax.plot([p4[0],p4[0]],[p5[1],p4[1]],[z,z],color=color,linewidth=0.5)				
				ax.plot([p4[0],p3[0]],[p4[1],p4[1]],[z,z],color=color,linewidth=0.5)
			else:
				ax.plot([p4[0],p4[0]],[p6[1],p7[1]],[z,z],color=color,linewidth=0.5)			
				ax.plot([p10[0],p10[0]],[p7[1],p11[1]],[z,z],color=color,linewidth=0.5)						
				ax.plot([p2[0],p8[0]],[p4[1],p4[1]],[z,z],color=color,linewidth=0.5)				
				ax.plot([p7[0],p8[0]],[p7[1],p7[1]],[z,z],color=color,linewidth=0.5)				
		
		if zones and self._parent:
			colors = ['r','g','b','k','m','c','y']
			colors.remove(color)
			if isinstance(zones,(int,str)): zones = [zones] 		#make an iterable
			for zn,col in zip(zones,colors): 	# iterate through zones
				if zn in self._parent.zone.keys():
					for nd in self._parent.zone[zn].nodelist:
						ax.plot([nd.position[0]],[nd.position[1]],[nd.position[2]],col+'.',ms=3,alpha=0.5)
		
		extension, save_fname, pdf = save_name(save,variable='grid',time=1)
		if self._parent:
			if self._parent.work_dir and not os.path.isdir(self._parent.work_dir): 
				os.makedirs(self._parent.work_dir)
			if self._parent.work_dir:
				plt.savefig(self._parent.work_dir+os.sep+save_fname, dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
					format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
			else:
				plt.savefig(save_fname, dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
					format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		else:
			plt.savefig(save_fname, dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
				format=extension,transparent=True, bbox_inches=None, pad_inches=0.1)
		if pdf: 
			os.system('epstopdf ' + save_fname)
			os.remove(save_fname)			
	def geom_coef_hist(self,log=True,plot=True):
		""" Plot histogram of geometric coefficients, sign is changed so they are positive

			:param log: Log10 scale coefficients if True, coefficients less than or equal to zero will not be included
			:type log: bool
			:param plot: Show plot if True
			:type plot: bool
			:returns: (lst(int),lst(fl64)) - tuple of count array and bin array
		"""
		conns = np.array([-c.geom_coef for c in self.connlist])
		if log: 
			conns = conns[np.where(conns!=0)[0]]
			conns = np.log10(conns)
			plt.xlabel('Log10(geom_coef)')
		else:
			plt.xlabel('geom_coef')
		h = plt.hist(conns)
		plt.ylabel('Count')
		plt.show()
		return h
	def add_nodetree(self,repair=False):
		""" Constuct octree for node positions. Call to update if changes made to grid.
		"""
		self._octree=octree(self.bounding_box,self.nodelist,repair=self)	
	def _summary(self):		
		L = 62
		s = ['']
		s.append(' ####---------------------------------------------------------####')
		line = ' #### FEHM grid file \''+self.filename+'\' summary.'
		for i in range(L-len(line)): line += ' '
		s.append(line+'####')
		s.append(' ####---------------------------------------------------------####')
		
		lines = []
		lines.append(' #### Domain extent:')
		lines.append('x = ['+str(self.xmin) + ', ' + str(self.xmax) + ']')		
		lines.append('y = ['+str(self.ymin) + ', ' + str(self.ymax) + ']')
		lines.append('z = ['+str(self.zmin) + ', ' + str(self.zmax) + ']')
		lines.append(' #### Statistics:')
		lines.append('nodes = ' +str(self.number_nodes))
		lines.append('elements = ' +str(self.number_elems))
		
		for line in lines:
			if line.startswith(' ##'):
				for i in range(L-len(line)): line += ' '
				s.append(line+'####')
			else:
				prntStr = ' #### -'
				keepGoing = True
				line = line.split()
				while keepGoing:
					if not line: 
						for i in range(L-len(prntStr)): prntStr += ' '
						s.append(prntStr+'####')
						prntStr = ' #### '
						break
					if len(prntStr)<(L-len(line[0])): 
						prntStr += ' '+line[0]
						line = line[1:]
					else:
						for i in range(L-len(prntStr)): prntStr += ' '
						s.append(prntStr+'####')
						prntStr = ' ####   '
		s.append(' ####---------------------------------------------------------####')
		s.append('')
		s = '\n'.join(s)
		pyfehm_print(s)
	def rotate(self,angle=0.,centre=[0.,0.]):
		'''Rotates the grid by some angle about a specified vertical axis.
		
		:param angle: Clockwise angle by which to rotate grid.
		:type angle: fl64
		:param centre: x and y coordinates of vertical axis about which to rotate. Alternatively, the centre of the computational domain can be specified by passing 'mid','middle','centre', or 'center'.
		:type centre: [fl64,fl64], str
		'''
		if centre in ['middle','mid','centre','center']:
			centre = [(self.xmin+self.xmax)/2.,(self.ymin+self.ymax)/2.]
		for nd in self.nodelist:
			old_pos = np.array(nd.position[0:2]) - np.array(centre) 	# position relative to centre of rotation
			theta_f = math.atan2(old_pos[1],old_pos[0]) + angle/180.*math.pi
			dist = np.sqrt(np.dot(old_pos,old_pos))
			new_pos = [dist*math.cos(theta_f),dist*math.sin(theta_f)]
			nd.position[0] = new_pos[0]+centre[0]
			nd.position[1] = new_pos[1]+centre[1]
	def get_bounding_box(self):
		minx = self.nodelist[0].position[0]
		maxx = self.nodelist[0].position[0]
		miny = self.nodelist[0].position[1]
		maxy = self.nodelist[0].position[1]
		minz = self.nodelist[0].position[2]
		maxz = self.nodelist[0].position[2]
		for node in self.nodelist:
			pos = node.position
			if pos[0]<minx: minx = pos[0]
			if pos[1]<miny: miny = pos[1]
			if pos[2]<minz: minz = pos[2]
			if pos[0]>maxx: maxx = pos[0]
			if pos[1]>maxy: maxy = pos[1]
			if pos[2]>maxz: maxz = pos[2]
		xr = maxx-minx+1.; yr = maxy-miny+1.; zr = maxz-minz+1.;
		
		c = np.array([minx+maxx,miny+maxy,minz+maxz])/2.
		
		r = np.max([xr,yr,zr])
		
		if minx == maxx:
			c[0] = minx; r = np.array([1,r,r])
			return [c-0.51*r,c+0.51*r]
		elif miny == maxy:
			c[1] = miny; r = np.array([r,1,r])
			return [c-0.51*r,c+0.51*r]
		elif minz == maxz:
			c[2] = minz; r = np.array([r,r,1])
			return [c-0.51*r,c+0.51*r]
		else:			
			return [c-0.51*r,c+0.51*r]
	bounding_box = property(get_bounding_box)	
	def _get_filename(self): return self._path.filename
	filename = property(_get_filename)#: (*str*) Name of FEHM grid file.
	def _get_node(self): return self._node
	node = property(_get_node)#: (*dict[fnode]*) Dictionary of grid nodes, indexed by node integer.
	def _get_nodelist(self): return self._nodelist
	nodelist = property(_get_nodelist)#: (*lst[fnode]*) List of all node objects in the grid.
	def _get_elem(self): return self._elem
	elem = property(_get_elem)#: (*dict[felem]*) Dictionary of elements, indexed by element integer.
	def _get_elemlist(self): return self._elemlist
	elemlist = property(_get_elemlist)#: (*lst[felem]*) List of all element objects in the grid.
	def _get_conn(self): return self._conn
	conn = property(_get_conn)#: (*dict[fconn]*) Dictionary of connections, indexed by a two element tuple of the member node integers.
	def _get_connlist(self): return self._connlist
	connlist = property(_get_connlist)#: (*lst[fconn]*) List of all connection objects in the grid.
	def _get_dimensions(self): return self._dimensions
	dimensions = property(_get_dimensions) #: (*int*) Dimensions of the grid.
	def get_xmin(self): return np.min([nd._position[0] for nd in self._nodelist])
	xmin = property(get_xmin) 				#: Minimum x-coordinate for all nodes.
	def get_xmax(self): return np.max([nd._position[0] for nd in self._nodelist])
	xmax = property(get_xmax)				#: Maximum x-coordinate for all nodes.
	def get_ymin(self): return np.min([nd._position[1] for nd in self._nodelist])
	ymin = property(get_ymin)				#: Minimum y-coordinate for all nodes.
	def get_ymax(self): return np.max([nd._position[1] for nd in self._nodelist])
	ymax = property(get_ymax)				#: Maximum y-coordinate for all nodes.
	def get_zmin(self): return np.min([nd._position[2] for nd in self._nodelist])
	zmin = property(get_zmin)				#: Minimum z-coordinate for all nodes.
	def get_zmax(self): return np.max([nd._position[2] for nd in self._nodelist])
	zmax = property(get_zmax)				#: Maximum z-coordinate for all nodes.
	def get_node_number(self): return len(self.nodelist)
	number_nodes = property(get_node_number)#: Number of nodes in grid.
	def get_element_number(self): return len(self.elemlist)
	number_elems = property(get_element_number)#: Number of elements in grid.
	def _get_stor(self): return self._stor.full_path
	stor = property(_get_stor) #: (*str*) Path to stor file.
	def _get_octree(self): return self._octree
	def _set_octree(self,value): self._octree = value
	octree = property(_get_octree, _set_octree) #: (*octree*) Octree object associated with the grid.
	def get_info(self):
		print 'FEHM grid file \''+self.filename+'\' summary.'
		print 'Model domain: x = ['+str(self.xmin) + ', ' + str(self.xmax) + ']'
		print '              y = ['+str(self.ymin) + ', ' + str(self.ymax) + ']'
		print '              z = ['+str(self.zmin) + ', ' + str(self.zmax) + ']'
		print '          nodes = ' +str(self.number_nodes)
		print '       elements = ' +str(self.number_elems)
		print ' '
	what = property(get_info) 				#: Print to screen information about the grid.
class fmake(object): 				#Rectilinear grid constructor.
	"""Generate an orthogonal mesh corresponding to vectors of nodal positions.
	"""
	def __init__(self,meshname,x=None,y=None,z=None):
		self._x = list(np.unique(x))
		self._y = list(np.unique(y))
		self._z = list(np.unique(z))
		self._dimension = None
		self._meshname = ''
		if meshname: self._meshname = meshname
	def write(self,meshname='',radial=False):
		"""Write out the grid file.
		
		:param meshname: Name of the grid file.
		:type meshname: str
		"""
		self.refresh()
		if meshname: self._meshname = meshname
		if self.meshname=='': self._meshname='default_GRID.inp'
		
		temp_path = fpath()
		temp_path.filename = self.meshname
		try:
			os.makedirs(temp_path.absolute_to_file)
		except:
			pass		
		outfile = open(temp_path.full_path,'w')
		
		outfile.write('coor\n')
		outfile.write('   '+str(len(self.nodelist))+'\n')
		for nd in self.nodelist: 
			outfile.write('%11d' % nd.index +'        ')
			outfile.write('%14.8f' % nd.position[0]+'        ')
			outfile.write('%14.8f' % nd.position[1]+'        ')
			outfile.write('%14.8f' % nd.position[2])
			outfile.write('\n')
		outfile.write('\t0\n')
		outfile.write('elem\n')
		outfile.write(str(len(self.elemlist[0]))+' '+str(len(self.elemlist))+'\n')		
		for i,el in enumerate(self.elemlist):
			outfile.write(str(i+1)+'   ')
			for nd in el:
				outfile.write(str(nd._index)+'   ')
			outfile.write('\n')
		outfile.write('\nstop\n')
		outfile.close()
	def refresh(self):
		"""Generate grid node and element objects corresponding to x y and z seeds.
		"""
		# first determine dimension of grid
		xF = len(self.x) > 1; yF = len(self.y) > 1; zF = len(self.z) > 1
		if xF and yF and zF: self.dimension = 3
		elif (xF and yF and not zF) or (xF and zF and not yF) or (zF and yF and not xF): self.dimension = 2
		else: 
			pyfehm_print('ERROR: not enough dimensions specified')
			return
		if self.dimension == 2: 
			# create nodes
			self._nodelist = []
			ind = 1
			self._z = 0.
			self._x = list(np.sort(self._x))
			self._y = list(np.sort(self._y))
			for yi in self._y:
				for xi in self._x:
					self._nodelist.append(fnode(index=ind,position=[xi,yi,0.]))
					ind +=1
			
			# create elements
			self._elemlist = []
			ind = 1
			xL = len(self._x); yL = len(self._y)
			for j in range(1,len(self._y)):
				for k in range(1,len(self._x)):
					nodes=[
						self._nodelist[(j-1)*xL + k-1],
						self._nodelist[(j-1)*xL + k],
						self._nodelist[j*xL + k],
						self._nodelist[j*xL + k-1],
						]
					self._elemlist.append(nodes)
		if self.dimension == 3:
			# create nodes
			self._nodelist = []
			ind = 1
			self._z = list(np.sort(self._z))
			self._x = list(np.sort(self._x))
			self._y = list(np.sort(self._y))
			for zi in self._z:
				for yi in self._y:
					for xi in self._x:
						self._nodelist.append(fnode(index=ind,position=[xi,yi,zi]))
						ind +=1
			
			# create elements
			self._elemlist = []
			ind = 1
			xL = len(self._x); yL = len(self._y); zL = len(self._z)
			for i in range(1,len(self._z)):
				for j in range(1,len(self._y)):
					for k in range(1,len(self._x)):
						nodes=[
							self._nodelist[i*xL*yL + (j-1)*xL + k-1],
							self._nodelist[i*xL*yL + (j-1)*xL + k],
							self._nodelist[i*xL*yL + j*xL + k],
							self._nodelist[i*xL*yL + j*xL + k-1],
							self._nodelist[(i-1)*xL*yL + (j-1)*xL + k-1],
							self._nodelist[(i-1)*xL*yL + (j-1)*xL + k],
							self._nodelist[(i-1)*xL*yL + j*xL + k],
							self._nodelist[(i-1)*xL*yL + j*xL + k-1],
							]
						self._elemlist.append(nodes)
	def _get_x(self): return self._x
	def _set_x(self,value): self._x = value
	x = property(_get_x, _set_x) #: (*lst[fl64]*) x coordinates of nodes.
	def _get_y(self): return self._y
	def _set_y(self,value): self._y = value
	y = property(_get_y, _set_y) #: (*lst[fl64]*) y coordinates of nodes.
	def _get_z(self): return self._z
	def _set_z(self,value): self._z = value
	z = property(_get_z, _set_z) #: (*lst[fl64]*) z coordinates of nodes.
	def _get_meshname(self): return self._meshname
	def _set_meshname(self,value): self._meshname = value
	meshname = property(_get_meshname, _set_meshname) #: (*str*) Name of grid file to write out.
	def _get_nodelist(self): return self._nodelist
	def _set_nodelist(self,value): self._nodelist = value
	nodelist = property(_get_nodelist, _set_nodelist) #: (*lst[fnode]*) List of node objects in the grid.
	def _get_elemlist(self): return self._elemlist
	def _set_elemlist(self,value): self._elemlist = value
	elemlist = property(_get_elemlist, _set_elemlist) #: (*lst[felem]*) List of element objects in the grid.
