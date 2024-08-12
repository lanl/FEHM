"""For reading FEHM output files."""

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
import os
from copy import copy,deepcopy
import platform
WINDOWS = platform.system()=='Windows'
if WINDOWS: slash = '\\'
else: slash = '/'


if True:    # output variable dictionaries defined in here, indented for code collapse
    cont_var_names_avs=dict([
    ('X coordinate (m)','x'),
    ('Y coordinate (m)','y'),
    ('Z coordinate (m)','z'),
    ('node','n'),
    ('Liquid Pressure (MPa)','P'),
    ('Vapor Pressure (MPa)','P_vap'),
    ('Capillary Pressure (MPa)','P_cap'),
    ('Saturation','saturation'),
    ('Temperature (deg C)','T'),
    ('Porosity','por'),
    ('X Permeability (log m**2)','perm_x'),
    ('Y Permeability (log m**2)','perm_y'),
    ('Z Permeability (log m**2)','perm_z'),
    ('X displacement (m)','disp_x'),
    ('Y displacement (m)','disp_y'),
    ('Z displacement (m)','disp_z'),
    ('X stress (MPa)','strs_xx'),
    ('Y stress (MPa)','strs_yy'),
    ('Z stress (MPa)','strs_zz'),
    ('XY stress (MPa)','strs_xy'),
    ('XZ stress (MPa)','strs_xz'),
    ('YZ stress (MPa)','strs_yz'),
    ('Youngs Mod (MPa)','E'),
    ('Excess Shear (MPa)','tau_ex'),
    ('Shear Angle (deg)','phi_dil'),
    ('Zone','zone'),
    ('Liquid Density (kg/m**3)','density'),
    ('Vapor Density (kg/m**3)','density_vap'),
    ('Source (kg/s)','flow'),
    ('Liquid Flux (kg/s)','flux'),
    ('Vapor Flux (kg/s)','flux_vap'),
    ('Volume Strain','strain'),
    ('Vapor X Volume Flux (m3/[m2 s])','flux_x_vap'),
    ('Vapor Y Volume Flux (m3/[m2 s])','flux_y_vap'),
    ('Vapor Z Volume Flux (m3/[m2 s])','flux_z_vap'),
    ('Liquid X Volume Flux (m3/[m2 s])','flux_x'),
    ('Liquid Y Volume Flux (m3/[m2 s])','flux_y'),
    ('Liquid Z Volume Flux (m3/[m2 s])','flux_z'),
    ]) 

    cont_var_names_tec=dict([
    ('X coordinate (m)','x'),
    ('Y coordinate (m)','y'),
    ('Z coordinate (m)','z'),
    ('X Coordinate (m)','x'),
    ('Y Coordinate (m)','y'),
    ('Z Coordinate (m)','z'),
    ('node','n'),
    ('Node','n'),
    ('Liquid Pressure (MPa)','P'),
    ('Vapor Pressure (MPa)','P_vap'),
    ('Capillary Pressure (MPa)','P_cap'),
    ('Saturation','saturation'),
    ('Water Saturation','water'),
    ('Super-Critical/Liquid CO2 Saturation','co2_sc_liquid'),
    ('Gaseous CO2 Saturation','co2_gas'),
    ('Dissolved CO2 Mass Fraction','co2_aq'),
    ('CO2 Phase State','co2_phase'),
    ('CO2 Gas Density (kg/m**3)','density_co2_gas'),
    ('CO2 Liquid Density (kg/m**3)','density_co2_sc_liquid'),
    ('Temperature (<sup>o</sup>C)','T'),
    ('Temperature (deg C)','T'),
    ('Porosity','por'),
    ('X Permeability (log m**2)','perm_x'),
    ('Y Permeability (log m**2)','perm_y'),
    ('Z Permeability (log m**2)','perm_z'),
    ('X displacement (m)','disp_x'),
    ('Y displacement (m)','disp_y'),
    ('Z displacement (m)','disp_z'),
    ('X stress (MPa)','strs_xx'),
    ('Y stress (MPa)','strs_yy'),
    ('Z stress (MPa)','strs_zz'),
    ('XY stress (MPa)','strs_xy'),
    ('XZ stress (MPa)','strs_xz'),
    ('YZ stress (MPa)','strs_yz'),
    ('Youngs Mod (MPa)','E'),
    ('Excess Shear (MPa)','tau_ex'),
    ('Shear Angle (deg)','phi_dil'),
    ('Zone','zone'),
    ('Liquid Density (kg/m**3)','density'),
    ('Vapor Density (kg/m**3)','density_vap'),
    ('Source (kg/s)','flow'),
    ('Liquid Flux (kg/s)','flux'),
    ('Vapor Flux (kg/s)','flux_vap'),
    ('Volume Strain','strain'),
    ('Vapor X Volume Flux (m3/[m2 s])','flux_x_vap'),
    ('Vapor Y Volume Flux (m3/[m2 s])','flux_y_vap'),
    ('Vapor Z Volume Flux (m3/[m2 s])','flux_z_vap'),
    ('Liquid X Volume Flux (m3/[m2 s])','flux_x'),
    ('Liquid Y Volume Flux (m3/[m2 s])','flux_y'),
    ('Liquid Z Volume Flux (m3/[m2 s])','flux_z'),
    ])

    cont_var_names_surf=dict([
    ('X coordinate (m)','x'), 
    ('X Coordinate (m)','x'), 
    ('X (m)','x'), 
    ('Y coordinate (m)','y'), 
    ('Y Coordinate (m)','y'), 
    ('Y (m)','y'), 
    ('Z coordinate (m)','z'), 
    ('Z Coordinate (m)','z'), 
    ('Z (m)','z'),
    ('node','n'),
    ('Node','n'),
    ('Liquid Pressure (MPa)','P'),
    ('Vapor Pressure (MPa)','P_vap'),
    ('Capillary Pressure (MPa)','P_cap'),
    ('Saturation','saturation'),
    ('Temperature (deg C)','T'),
    ('Porosity','por'),
    ('X Permeability (log m**2)','perm_x'),
    ('Y Permeability (log m**2)','perm_y'),
    ('Z Permeability (log m**2)','perm_z'),
    ('X displacement (m)','disp_x'),
    ('Y displacement (m)','disp_y'),
    ('Z displacement (m)','disp_z'),
    ('X stress (MPa)','strs_xx'),
    ('Y stress (MPa)','strs_yy'),
    ('Z stress (MPa)','strs_zz'),
    ('XY stress (MPa)','strs_xy'),
    ('XZ stress (MPa)','strs_xz'),
    ('YZ stress (MPa)','strs_yz'),
    ('Youngs Mod (MPa)','E'),
    ('Excess Shear (MPa)','tau_ex'),
    ('Shear Angle (deg)','phi_dil'),
    ('Zone','zone'),
    ('Liquid Density (kg/m**3)','density'),
    ('Vapor Density (kg/m**3)','density_vap'),
    ('Source (kg/s)','flow'),
    ('Liquid Flux (kg/s)','flux'),
    ('Vapor Flux (kg/s)','flux_vap'),
    ('Volume Strain','strain'),
    ('Vapor X Volume Flux (m3/[m2 s])','flux_x_vap'),
    ('Vapor Y Volume Flux (m3/[m2 s])','flux_y_vap'),
    ('Vapor Z Volume Flux (m3/[m2 s])','flux_z_vap'),
    ('Liquid X Volume Flux (m3/[m2 s])','flux_x'),
    ('Liquid Y Volume Flux (m3/[m2 s])','flux_y'),
    ('Liquid Z Volume Flux (m3/[m2 s])','flux_z'),
    ('Water Saturation','saturation'),
    ('Super-Critical/Liquid CO2 Saturation','co2_liquid'),
    ('Gaseous CO2 Saturation','co2_gas'),
    ('Dissolved CO2 Mass Fraction','co2_aq'),
    ('CO2 Phase State','co2_phase'),
    ('Aqueous_Species_001','species001_aq'),
    ('Aqueous_Species_002','species002_aq'),
    ('Aqueous_Species_003','species003_aq'),
    ('Aqueous_Species_004','species004_aq'),
    ('Aqueous_Species_005','species005_aq'),
    ('Aqueous_Species_006','species006_aq'),
    ('Aqueous_Species_007','species007_aq'),
    ('Aqueous_Species_008','species008_aq'),
    ('Aqueous_Species_009','species009_aq'),
    ('Aqueous_Species_010','species010_aq'),
    ('Aqueous_Species_011','species011_aq'),
    ('Aqueous_Species_012','species012_aq'),
    ('Aqueous_Species_013','species013_aq'),
    ('Aqueous_Species_014','species014_aq'),
    ('Aqueous_Species_015','species015_aq'),
    ('Aqueous_Species_016','species016_aq'),
    ('Aqueous_Species_017','species017_aq'),
    ('Aqueous_Species_018','species018_aq'),
    ('Aqueous_Species_019','species019_aq'),
    ('Aqueous_Species_020','species020_aq'),
    ])

    hist_var_names=dict([
    ('denAIR','density_air'),
    ('disx','disp_x'),
    ('disy','disp_y'),
    ('disz','disp_z'),
    ('enth','enthalpy'),
    ('glob','global'),
    ('humd','humidity'),
    ('satr','saturation'),
    ('strain','strain'),
    ('strx','strs_xx'),
    ('stry','strs_yy'),
    ('strz','strs_zz'),
    ('strxy','strs_xy'),
    ('strxz','strs_xz'),
    ('stryz','strs_yz'),
    ('wcon','water_content'),
    ('denWAT','density'),
    ('flow','flow'),
    ('visAIR','viscosity_air'),
    ('visWAT','viscosity'),
    ('wt','water_table'),
    ('presCAP','P_cap'),
    ('presVAP','P_vap'),
    ('presWAT','P'),
    ('presCO2','P_co2'),
    ('temp','T'),
    ('co2md','massfrac_co2_aq'),
    ('co2mf','massfrac_co2_free'),
    ('co2mt','mass_co2'),
    ('co2sg','saturation_co2g'),
    ('co2sl','saturation_co2l'),
    ])
    
    flxz_water_names = [
    'water_source',
    'water_sink',
    'water_net',
    'water_boundary',]
    
    flxz_vapor_names = [
    'vapor_source',
    'vapor_sink',
    'vapor_net',
    'vapor_boundary',]
    
    flxz_co2_names = [
    'co2_source',
    'co2_sink',
    'co2_in',
    'co2_out',
    'co2_boundary',
    'co2_sourceG',
    'co2_sinkG',
    'co2_inG',
    'co2_outG']
class fcontour(object):                     
    '''
    ****************************************************************
    Reading and plotting methods associated with contour output data.
    This includes file types avs, tec, & surf.
    ****************************************************************
    '''
    def __init__(self,filename=None,latest=False,first=False,nearest=None):
        if not isinstance(filename,list):
            self._filename=os_path(filename)
        self._silent = True#dflt.silent
        self._times=[]   
        self._format = ''
        self._data={}
        self._material = {}
        self._material_properties = []
        self._row=None
        self._variables=[]  
        self._user_variables = []
        self.key_name=[]
        self._keyrows={}
        self.column_name=[]
        self.num_columns=0
        self._x = []
        self._y = []
        self._z = []
        self._xmin = None
        self._ymin = None
        self._zmin = None
        self._xmax = None
        self._ymax = None
        self._zmax = None
        self._latest = latest
        self._first = first
        self._nearest = nearest
        if isinstance(self._nearest,(float,int)): self._nearest = [self._nearest]
        self._nkeys=1
        if filename is not None: self.read(filename,self._latest,self._first,self._nearest)
    def __getitem__(self,key):
        if key in self.times:
            return self._data[key]
        elif np.min(abs(self.times-key)/self.times)<.01:
            ind = np.argmin(abs(self.times-key))
            return self._data[self.times[ind]]
        else: return None
    def read(self,filename,latest=False,first=False,nearest=[]):                        # read contents of file
        '''
        Read in FEHM contour output information.
        
        :param filename: File name for output data, can include wildcards to define multiple output files.
        :type filename: str
        :param latest:  Boolean indicating PyFEHM should read the latest entry in a wildcard search.
        :type latest: bool
        :param first: Boolean indicating PyFEHM should read the first entry in a wildcard search.
        :type first: bool
        :param nearest: Read in the file with date closest to the day supplied. List input will parse multiple output files.
        :type nearest: fl64,list
        '''
        from glob import glob
        if isinstance(filename,list):
            files = filename
        else:
            filename = os_path(filename)
            files=glob(filename)
        if len(files)==0:   
            print('ERROR: '+filename+' not found')
            return
        # decision-making
        mat_file = None
        multi_type = None
        # are there multiple file types? e.g., _con_ and _sca_?
        # is there a material properties file? e.g., 'mat_nodes'?
        file_types = []
        for file in files:
            if '_sca_node' in file and 'sca' not in file_types: file_types.append('sca')
            if '_vec_node' in file and 'vec' not in file_types: file_types.append('vec')
            if '_con_node' in file and 'con' not in file_types: file_types.append('con')
            if '_hf_node' in file and 'hf' not in file_types: file_types.append('hf')
            if 'mat_node' in file: mat_file = file
        
        if self._nearest or latest or first:
            files = list(filter(os.path.isfile, glob(filename)))
            if mat_file: files.remove(mat_file)
            files.sort(key=lambda x: os.path.getmtime(x))
            files2 = []
            
            # retrieve first created and same time in group
            if first:
                files2.append(files[0])
                for file_type in file_types:
                    tag = '_'+file_type+'_node'
                    if tag in files2[-1]:
                        prefix = files2[-1].split(tag)[0]
                        break                       
                for file in files:
                    if file.startswith(prefix) and tag not in file: files2.append(file)
                    
            # retrieve files nearest in time to given (and same time in group)
            if self._nearest:
                ts = []
                for file in files:
                    file = file.split('_node')[0]
                    file = file.split('_sca')[0]
                    file = file.split('_con')[0]
                    file = file.split('_vec')[0]
                    file = file.split('_hf')[0]
                    file = file.split('_days')[0]                       
                    file = file.split('.')
                    file = file[-2:]
                    ts.append(float('.'.join(file)))
                ts = np.unique(ts)
                
                for near in self._nearest:
                    tsi = min(enumerate(ts), key=lambda x: abs(x[1]-near))[0]
                    files2.append(files[tsi])
                    for file_type in file_types:
                        tag = '_'+file_type+'_node'
                        if tag in files2[-1]:
                            prefix = files2[-1].split(tag)[0]
                            break
                    for file in files:
                        if file.startswith(prefix) and tag not in file: files2.append(file)                     
                        
            # retrieve last created and same time in group
            if latest:
                files2.append(files[-1])
                for file_type in file_types:
                    tag = '_'+file_type+'_node'
                    if tag in files2[-1]:
                        prefix = files2[-1].split(tag)[0]
                        break
                for file in files:
                    if file.startswith(prefix) and tag not in file: files2.append(file)
                    
            # removes duplicates
            files = []
            for file in files2:
                if file not in files: files.append(file)
        
        # group files into their types
        FILES = []
        for file_type in file_types:
            tag = '_'+file_type+'_node'
            FILES.append(sort_tec_files([file for file in files if tag in file]))
        FILES = np.array(FILES)
        
        # determine headers for 'tec' output
        for i in range(FILES.shape[1]):
            if not self._variables:
                files = FILES[:,i]
                headers = []
                for file in sort_tec_files(files):
                    with open(file, 'r') as fp:
                        headers.append(fp.readline())
                    #fp.close()
                firstFile = self._detect_format(headers)
                if self._format=='tec' and firstFile: 
                    headers = []
                    for file in sort_tec_files(files):
                        with open(file, 'r') as fp:
                            fp.readline()
                            headers.append(fp.readline())
                        #fp.close()
                    self._setup_headers_tec(headers)        
        
        # read in output data
        for i in range(FILES.shape[1]):
            files = FILES[:,i]
            # Skip -1 file if present
            if '-1' in files[0]: continue
            #for file in sort_tec_files(files): print(file)
            if not self._variables:
                headers = []
                for file in sort_tec_files(files):
                    with open(file, 'r') as fp:
                        headers.append(fp.readline())
                    #fp.close()
                self._detect_format(headers)
                #if self._format=='tec': self._setup_headers_tec(headers)
                if self._format=='avs': self._setup_headers_avs(headers,files)
                elif self._format=='avsx': self._setup_headers_avsx(headers)
                elif self._format=='surf': self._setup_headers_surf(headers)
                else: print('ERROR: Unrecognised format');return
                self.num_columns = len(self.variables)+1
            if self._format == 'tec': self._read_data_tec(files,mat_file)
            elif self._format == 'surf': self._read_data_surf(files,mat_file)
            elif self._format == 'avs': self._read_data_avs(files,mat_file)
            elif self._format == 'avsx': self._read_data_avsx(files,mat_file)
        
        # assemble grid information
        if 'x' in self.variables:
            self._x = np.unique(self[self.times[0]]['x'])
            self._xmin,self._xmax = np.min(self.x), np.max(self.x)
        if 'y' in self.variables:
            self._y = np.unique(self[self.times[0]]['y'])
            self._ymin,self._ymax = np.min(self.y), np.max(self.y)
        if 'z' in self.variables:
            self._z = np.unique(self[self.times[0]]['z'])
            self._zmin,self._zmax = np.min(self.z), np.max(self.z)
        #if dflt.parental_cont:
        #   print('')
        #   print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        #   print('WARNING:')
        #   print('')
        #   print('Contour data is indexed using the Pythonic convention in which the first index is 0. FEHM node numbering convention begins at 1.')
        #   print('')
        #   print('THEREFORE, to get the correct contour value for a particular node, you need to pass the node index MINUS 1. Using node index to access contour data will return incorrect values.')
        #   print('')
        #   print('For example:')
        #   print('>>> node10 = dat.grid.node[10]')
        #   print('>>> c = fcontour(\'*.csv\')')
        #   print('>>> T_node10 = c[c.times[-1]][\'T\'][node10.index - 1]')
        #   print('  or')
        #   print('>>> T_node10 = c[c.times[-1]][\'T\'][9]')
        #   print('will return the correct value for node 10.')
        #   print('')
        #   print('Do not turn off this message unless you understand how to correctly access nodal values from contour data.') 
        #   print('To turn off this message, open the environment file \'fdflt.py\' and set self.parental_cont = False')
        #   print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        #   print('')
    def _detect_format(self,headers):
        if headers[0].startswith('TITLE ='):        # check for TEC output
            self._format = 'tec'
        if headers[0].startswith('ZONE '):      # check for TEC output
            self._format = 'tec'
            return False
        elif headers[0].startswith('node, '):       # check for SURF output
            self._format = 'surf'
        elif headers[0].startswith('nodes at '):    # check for AVSX output
            self._format = 'avsx'
        elif headers[0].split()[0].isdigit():           # check for AVS output
            self._format = 'avs'
        return True
    def _setup_headers_avsx(self,headers):      # headers for the AVSX output format
        self._variables.append('n')
        for header in headers:
            header = header.strip().split(' : ')
            for key in header[1:]: 
                if key in list(cont_var_names_avs.keys()):
                    var = cont_var_names_avs[key]
                else: var = key
                self._variables.append(var)
    def _read_data_avsx(self,files,mat_file):               # read data in AVSX format
        datas = []
        for file in sorted(files):      # use alphabetical sorting
            with open(file, 'r') as fp:
                header = fp.readline()
                if file == sorted(files)[0]:
                    header = header.split('nodes at ')[1]
                    header = header.split('days')[0]
                    time = float(header)*24*2600
                    self._times.append(time)        
                lns = fp.readlines()
            #fp.close()
            datas.append(np.array([[float(d) for d in ln.strip().split(':')] for ln in lns]))
        data = np.concatenate(datas,1)
        self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
        
        if mat_file and not self._material_properties:
            with open(mat_file, 'r') as fp:
                header = fp.readline()
                self._material_properties = header.split(':')[1:]
            lns = fp.readlines()
            #fp.close()
            data = np.array([[float(d) for d in ln.strip().split(':')[1:]] for ln in lns])
            self._material= dict([(var,data[:,icol]) for icol,var in enumerate(self._material_properties)])
    def _setup_headers_avs(self,headers,files):         # headers for the AVS output format
        for header,file in zip(headers,files):
            lns_num = int(header.strip().split()[0])
            fp = open(file)
            lns = [fp.readline() for i in range(lns_num+1)][1:]
            #fp.close()
            self._variables.append('n')
            for ln in lns:
                varname = ln.strip().split(',')[0]
                if varname in list(cont_var_names_avs.keys()):
                    var = cont_var_names_avs[varname]
                else: var = varname
                if var not in self._variables: self._variables.append(var)
    def _read_data_avs(self,files,mat_file):        # read data in AVS format
        datas = []
        for file in sorted(files):
            first = (file == sorted(files)[0])
            with open(file, 'r') as fp:
                lns = fp.readlines()
                lns = lns[int(float(lns[0].split()[0]))+1:]
            
            if first: 
                file = file.split('_node')[0]
                file = file.split('_sca')[0]
                file = file.split('_con')[0]
                file = file.split('_vec')[0]
                file = file.split('_hf')[0]
                file = file.split('_days')[0]                       
                file = file.split('.')
                file = [fl for fl in file if fl.isdigit() or 'E-' in fl]
                time = float('.'.join(file))
                self._times.append(time)
            
            if first: 
                datas.append(np.array([[float(d) for d in ln.strip().split()] for ln in lns]))
            else:
                datas.append(np.array([[float(d) for d in ln.strip().split()[4:]] for ln in lns]))
            
        data = np.concatenate(datas,1)
        self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
    def _setup_headers_surf(self,headers):      # headers for the SURF output format
        for header in headers:
            header = header.strip().split(', ')
            for key in header: 
                varname = key.split('"')[0]
                varname = varname.strip()
                if varname in list(cont_var_names_surf.keys()):
                    var = cont_var_names_surf[varname]
                else: var = varname
                if var not in self._variables: self._variables.append(var)
    def _read_data_surf(self,files,mat_file):       # read data in SURF format
        datas = []
        for file in sorted(files):
            first = (file == sorted(files)[0])
            with open(file, 'r') as fp:
                lni = file.split('.',1)[1]
            
                if first: 
                    file = file.split('_node')[0]
                    file = file.split('_sca')[0]
                    file = file.split('_con')[0]
                    file = file.split('_vec')[0]
                    file = file.split('_hf')[0]
                    file = file.split('_days')[0]                       
                    file = file.split('.')
                    file = [fl for fl in file if fl.isdigit() or 'E-' in fl]
                    time = float('.'.join(file))
                    self._times.append(time)
            
                lni=fp.readline()           
                lns = fp.readlines()
                #fp.close()
            
                if first: 
                    datas.append(np.array([[float(d) for d in ln.strip().split(',')] for ln in lns]))
                else:
                    datas.append(np.array([[float(d) for d in ln.strip().split(',')[4:]] for ln in lns]))
            
            data = np.concatenate(datas,1)
            self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
        
            if mat_file and not self._material_properties:
                with open(mat_file, 'r') as fp:
                    header = fp.readline()
                for mat_prop in header.split(',')[1:]:
                    if 'specific heat' not in mat_prop:
                        self._material_properties.append(mat_prop.strip())
                lns = fp.readlines()
                #fp.close()
                data = np.array([[float(d) for d in ln.strip().split(',')[1:]] for ln in lns])
                self._material= dict([(var,data[:,icol]) for icol,var in enumerate(self._material_properties)])
    def _setup_headers_tec(self,headers):       # headers for the TEC output format
        for header in headers:
            header = header.split(' "')
            for key in header[1:]: 
                varname = key.split('"')[0].strip()
                if varname in list(cont_var_names_tec.keys()):
                    var = cont_var_names_tec[varname]
                else: var = varname
                
                if var not in self._variables: self._variables.append(var)
    def _read_data_tec(self,files,mat_file):                        # read data in TEC format
        datas = []
        for file in sorted(files):
            first = (file == sorted(files)[0])
            with open(file, 'r') as fp:
                ln = fp.readline()
                has_xyz = False
                while not ln.startswith('ZONE'):
                    ln = fp.readline()
                    has_xyz = True
                
                if first: 
                    lni = ln.split('"')[1]
                    time = lni.split('days')[0].strip()
                    time = float(time.split()[-1].strip())
                    try:
                        if time<self._times[0]: return
                    except: pass
                    self._times.append(time)
                    nds = None
                    if 'N =' in ln: 
                        nds = int(ln.split('N = ')[-1].strip().split(',')[0].strip())
            
                lns = fp.readlines()
                #fp.close()
                if nds: lns = lns[:nds]         # truncate to remove connectivity information
                
                if has_xyz:
                    if first: 
                        datas.append(np.array([[float(d) for d in ln.strip().split()] for ln in lns]))
                    else:
                        datas.append(np.array([[float(d) for d in ln.strip().split()[4:]] for ln in lns]))
                else:
                    if first: 
                        datas.append(np.array([[float(d) for d in ln.strip().split()] for ln in lns]))
                    else:
                        datas.append(np.array([[float(d) for d in ln.strip().split()[1:]] for ln in lns]))
                        
        data = np.concatenate(datas,1)
        #print('data ', data)
        #print('times: ', self._times)
        if data.shape[1]< len(self.variables):      # insert xyz data from previous read
            data2 = []
            j = 0
            for var in self.variables:
                if var == 'x': 
                    data2.append(self._data[self._times[0]]['x'])
                elif var == 'y': 
                    data2.append(self._data[self._times[0]]['y'])
                elif var == 'z': 
                    data2.append(self._data[self._times[0]]['z'])
                elif var == 'zone': 
                    data2.append(self._data[self._times[0]]['zone'])
                else: 
                    data2.append(data[:,j]); j +=1
            data = np.transpose(np.array(data2))
        self._data[time] = dict([(var,data[:,icol]) for icol,var in enumerate(self.variables)])
        if mat_file and not self._material_properties:
            with open(mat_file, 'r') as fp:
                fp.readline()
                header = fp.readline()
            for mat_prop in header.split(' "')[5:]:
                self._material_properties.append(mat_prop.split('"')[0].strip())
            lns = fp.readlines()
            if lns[0].startswith('ZONE'): lns = lns[1:]
            #fp.close()
            if nds: lns = lns[:nds]         # truncate to remove connectivity information
            data = np.array([[float(d) for d in ln.strip().split()[4:]] for ln in lns[:-1]])
            self._material= dict([(var,data[:,icol]) for icol,var in enumerate(self._material_properties)])
    def _check_inputs(self,variable, time, slice):  # assesses whether sufficient input information for slice plot
        if not variable: 
            s = ['ERROR: no plot variable specified.']
            s.append('Options are')
            for var in self.variables: s.append(var)
            s = '\n'.join(s)
            print(s)
            return True
        if time==None: 
            s = ['ERROR: no plot time specified.']
            s.append('Options are')
            for time in self.times: s.append(time)
            s = '\n'.join(s)
            print(s)
            return True
        if not slice: 
            s = ['Error: slice orientation undefined.']
            s.append('Options are')
            s.append('[\'x\',float] - slice parallel to y-axis at x=float')
            s.append('[\'y\',float] - slice parallel to x-axis at y=float')
            s.append('[\'theta\',float] - angle measured anti-clockwise from +x')
            s.append('[[float,float],[float,float]] - point to point')
            s = '\n'.join(s)
            print(s)
            return True
        return False
    def new_variable(self,name,time,data):  
        '''Creates a new variable, which is some combination of the available variables.
        
        :param name: Name for the variable.
        :type name: str
        :param time: Time key which the variable should be associated with. Must be one of the existing keys, i.e., an item in fcontour.times.
        :type time: fl64
        :param data: Variable data, most likely some combination of the available parameters, e.g., pressure*temperature, pressure[t=10] - pressure[t=5]
        :type data: lst[fl64]
        '''
        if time not in self.times: 
            print('ERROR: supplied time must correspond to an existing time in fcontour.times')
            return
        if name in self.variables:
            print('ERROR: there is already a variable called \''+name+'\', please choose a different name')
            return
        self._data[time][name] = data
        if name not in self._user_variables:
            self._user_variables.append(name)
    def slice(self, variable, slice, divisions, time=None, method='nearest'):
        '''Returns mesh data for a specified slice orientation from 3-D contour output data.
        
        :param variable: Output data variable, for example 'P' = pressure. Alternatively, variable can be a five element list, first element 'cfs', remaining elements fault azimuth (relative to x), dip, friction coefficient and cohesion. Will return coulomb failure stress.
        :type variable: str
        :param time: Time for which output data is requested. Can be supplied via ``fcontour.times`` list. Default is most recently available data.
        :type time: fl64
        :param slice: List specifying orientation of output slice, e.g., ['x',200.] is a vertical slice at ``x = 200``, ['z',-500.] is a horizontal slice at ``z = -500.``, [point1, point2] is a fixed limit vertical or horizontal domain corresponding to the bounding box defined by point1 and point2.
        :type slice: lst[str,fl64]
        :param divisions: Resolution to supply mesh data.
        :type divisions: [int,int]
        :param method: Method of interpolation, options are 'nearest', 'linear'.
        :type method: str:  
        :returns: X -- x-coordinates of mesh data.
        
        '''
        if time==None: 
            if np.min(self.times)<0: time = self.times[0]
            else: time = self.times[-1]
        from scipy.interpolate import griddata      
        delta = False
        if isinstance(time,list) or isinstance(time,np.ndarray):
            if len(time)>1: 
                time0 = np.min(time)
                time = np.max(time)
                delta=True      
        dat = self[time]
        
        # check to see if cfs plot requested
        cfs = False
        if isinstance(variable,list):
            if variable[0] in ['cfs','CFS']: cfs = True
        
        if not cfs:
            if delta: dat0 = self[time0]
            if isinstance(slice[0],str):
                if slice[0].startswith('y'):
                    xmin = np.min(dat['x']);xmax = np.max(dat['x'])
                    ymin = np.min(dat['z']);ymax = np.max(dat['z'])     
                    if slice[1] is None:
                        points = np.transpose(np.array([dat['x'],dat['z'],np.ones((1,len(dat['z'])))[0]]))
                        slice[1] = 1
                    else:
                        points = np.transpose(np.array([dat['x'],dat['z'],dat['y']]))
                elif slice[0].startswith('x'):
                    xmin = np.min(dat['y']);xmax = np.max(dat['y'])
                    ymin = np.min(dat['z']);ymax = np.max(dat['z'])     
                    if slice[1] is None:
                        points = np.transpose(np.array([dat['y'],dat['z'],np.ones((1,len(dat['z'])))[0]]))
                        slice[1] = 1
                    else:
                        points = np.transpose(np.array([dat['y'],dat['z'],dat['x']]))
                elif slice[0].startswith('z'):
                    xmin = np.min(dat['x']);xmax = np.max(dat['x'])
                    ymin = np.min(dat['y']);ymax = np.max(dat['y'])     
                    if slice[1] is None:
                        points = np.transpose(np.array([dat['x'],dat['y'],np.ones((1,len(dat['y'])))[0]]))
                        slice[1] = 1
                    else:
                        points = np.transpose(np.array([dat['x'],dat['y'],dat['z']]))
                elif slice[0].startswith('theta'): 
                    print('ERROR: theta slicing not supported yet') 
                    return
                xrange = np.linspace(xmin,xmax,divisions[0])
                yrange = np.linspace(ymin,ymax,divisions[1])
                X,Y = np.meshgrid(xrange,yrange)
                Z = (X+np.sqrt(1.757))/(X+np.sqrt(1.757))*slice[1]
                pointsI = np.transpose(np.reshape((X,Y,Z),(3,X.size)))
                vals = np.transpose(np.array(dat[variable]))
                valsI = griddata(points,vals,pointsI,method=method)
                valsI =  np.reshape(valsI,(X.shape[0],X.shape[1]))
                if delta:
                    vals = np.transpose(np.array(dat0[variable]))
                    valsI0 = griddata(points,vals,pointsI,method=method)
                    valsI0 =  np.reshape(valsI0,(X.shape[0],X.shape[1]))
                    valsI = valsI - valsI0
            elif isinstance(slice[0],list):
                # check if horizontal or vertical slice
                dx,dy,dz = abs(slice[0][0]-slice[1][0]),abs(slice[0][1]-slice[1][1]),abs(slice[0][2]-slice[1][2])
                if 100*dz<dx and 100*dz<dy:     #horizontal
                    xmin,xmax = np.min([slice[0][0],slice[1][0]]),np.max([slice[0][0],slice[1][0]])
                    ymin,ymax = np.min([slice[0][1],slice[1][1]]),np.max([slice[0][1],slice[1][1]])
                    xrange = np.linspace(xmin,xmax,divisions[0])
                    yrange = np.linspace(ymin,ymax,divisions[1])
                    X,Y = np.meshgrid(xrange,yrange)
                    Z = (X+np.sqrt(1.757))/(X+np.sqrt(1.757))*(slice[0][2]+slice[1][2])/2
                else:                           #vertical 
                    xmin,xmax = 0,np.sqrt((slice[0][0]-slice[1][0])**2+(slice[0][1]-slice[1][1])**2)
                    ymin,ymax = np.min([slice[0][2],slice[1][2]]),np.max([slice[0][2],slice[1][2]])
                    xrange = np.linspace(xmin,xmax,divisions[0])
                    yrange = np.linspace(ymin,ymax,divisions[1])
                    X,Z = np.meshgrid(xrange,yrange)
                    Y = X/xmax*abs(slice[0][1]-slice[1][1]) + slice[0][1]
                    X = X/xmax*abs(slice[0][0]-slice[1][0]) + slice[0][0]
                points = np.transpose(np.array([dat['x'],dat['y'],dat['z']]))
                pointsI = np.transpose(np.reshape((X,Y,Z),(3,X.size)))
                vals = np.transpose(np.array(dat[variable]))
                valsI = griddata(points,vals,pointsI,method=method)
                valsI =  np.reshape(valsI,(X.shape[0],X.shape[1]))
                if delta:
                    vals = np.transpose(np.array(dat0[variable]))
                    valsI0 = griddata(points,vals,pointsI,method=method)
                    valsI0 =  np.reshape(valsI0,(X.shape[0],X.shape[1]))
                    valsI = valsI - valsI0
            
        else:
            if delta: time0 = time[0]; time = time[-1]
            X,Y,Z,sxx = self.slice('strs_xx', slice, divisions, time, method)
            X,Y,Z,syy = self.slice('strs_yy', slice, divisions, time, method)
            X,Y,Z,szz = self.slice('strs_zz', slice, divisions, time, method)
            X,Y,Z,sxy = self.slice('strs_xy', slice, divisions, time, method)
            X,Y,Z,sxz = self.slice('strs_xz', slice, divisions, time, method)
            X,Y,Z,syz = self.slice('strs_yz', slice, divisions, time, method)
            X,Y,Z,sp  = self.slice('P',    slice, divisions, time, method)
            
            dip = variable[2]/180.*math.pi
            azi = variable[1]/180.*math.pi+3.14159/2.
            nhat = np.array([np.cos(azi)*np.sin(dip),np.sin(azi)*np.sin(dip),np.cos(dip)])
            mu = variable[3]
            cohesion = variable[4]
            
            px = sxx*nhat[0]+sxy*nhat[1]+sxz*nhat[2]
            py = sxy*nhat[0]+syy*nhat[1]+syz*nhat[2]
            pz = sxz*nhat[0]+syz*nhat[1]+szz*nhat[2]
            
            sig = px*nhat[0]+py*nhat[1]+pz*nhat[2]
            tau = np.sqrt(px**2+py**2+pz**2 - sig**2)
            valsI = tau - mu*(sig-sp) - cohesion
            if delta:
                X,Y,Z,sxx = self.slice('strs_xx', slice, divisions, time0, method)
                X,Y,Z,syy = self.slice('strs_yy', slice, divisions, time0, method)
                X,Y,Z,szz = self.slice('strs_zz', slice, divisions, time0, method)
                X,Y,Z,sxy = self.slice('strs_xy', slice, divisions, time0, method)
                X,Y,Z,sxz = self.slice('strs_xz', slice, divisions, time0, method)
                X,Y,Z,syz = self.slice('strs_yz', slice, divisions, time0, method)
                X,Y,Z,sp  = self.slice('P',    slice, divisions, time0, method)
                
                px = sxx*nhat[0]+sxy*nhat[1]+sxz*nhat[2]
                py = sxy*nhat[0]+syy*nhat[1]+syz*nhat[2]
                pz = sxz*nhat[0]+syz*nhat[1]+szz*nhat[2]
                
                sig = px*nhat[0]+py*nhat[1]+pz*nhat[2]
                tau = np.sqrt(px**2+py**2+pz**2 - sig**2)
                valsI = valsI - (tau - mu*(sig-sp) - cohesion)
                
        return X, Y, Z, valsI

    def node(self,node,time=None,variable=None):
        '''Returns all information for a specific node.
        
        If time and variable not specified, a dictionary of time series is returned with variables as the dictionary keys.
        
        If only time is specified, a dictionary of variable values at that time is returned, with variables as dictionary keys.
        
        If only variable is specified, a time series vector is returned for that variable.
        
        If both time and variable are specified, a single value is returned, corresponding to the variable value at that time, at that node.
        
        :param node: Node index for which variable information required.
        :type node: int 
        :param time: Time at which variable information required. If not specified, all output.
        :type time: fl64
        :param variable: Variable for which information requested. If not specified, all output.
        :type variable: str
        
        '''
        if 'n' not in self.variables: 
            print('Node information not available')
            return
        nd = np.where(self[self.times[0]]['n']==node)[0][0]
        if time is None and variable is None:
            ks = copy(self.variables); ks.remove('n')
            outdat = dict([(k,[]) for k in ks])
            for t in self.times:
                dat = self[t]
                for k in list(outdat.keys()):
                    outdat[k].append(dat[k][nd])
        elif time is None:
            if variable not in self.variables: 
                print('ERROR: no variable by that name')
                return
            outdat = []
            for t in self.times:
                dat = self[t]
                outdat.append(dat[variable][nd])
            outdat = np.array(outdat)
        elif variable is None:
            ks = copy(self.variables); ks.remove('n')
            outdat = dict([(k,self[time][k][nd]) for k in ks])          
        else:
            outdat = self[time][variable][nd]
        return outdat
    def _get_variables(self): return self._variables
    variables = property(_get_variables)#: (*lst[str]*) List of variables for which output data are available.
    def _get_user_variables(self): return self._user_variables
    user_variables = property(_get_user_variables) #: (*lst[str]*) List of user-defined variables for which output data are available.
    def _get_format(self): return self._format
    format = property(_get_format) #: (*str*) Format of output file, options are 'tec', 'surf', 'avs' and 'avsx'.
    def _get_filename(self): return self._filename
    filename = property(_get_filename)  #: (*str*) Name of FEHM contour output file. Wildcards can be used to define multiple input files.
    def _get_times(self): return np.sort(self._times)
    times = property(_get_times)    #: (*lst[fl64]*) List of times (in seconds) for which output data are available.
    def _get_material_properties(self): return self._material_properties
    def _set_material_properties(self,value): self._material_properties = value
    material_properties = property(_get_material_properties, _set_material_properties) #: (*lst[str]*) List of material properties, keys for the material attribute.
    def _get_material(self): return self._material
    def _set_material(self,value): self._material = value
    material = property(_get_material, _set_material) #: (*dict[str]*) Dictionary of material properties, keyed by property name, items indexed by node_number - 1. This attribute is empty if no material property file supplied.
    def _get_x(self): return self._x
    def _set_x(self,value): self._x = value
    x = property(_get_x, _set_x) #: (*lst[fl64]*) Unique list of nodal x-coordinates for grid.
    def _get_y(self): return self._y
    def _set_y(self,value): self._y = value
    y = property(_get_y, _set_y) #: (*lst[fl64]*) Unique list of nodal y-coordinates for grid.
    def _get_z(self): return self._z
    def _set_z(self,value): self._z = value
    z = property(_get_z, _set_z) #: (*lst[fl64]*) Unique list of nodal z-coordinates for grid.
    def _get_xmin(self): return self._xmin
    def _set_xmin(self,value): self._xmin = value
    xmin = property(_get_xmin, _set_xmin) #: (*fl64*) Minimum nodal x-coordinate for grid.
    def _get_xmax(self): return self._xmax
    def _set_xmax(self,value): self._xmax = value
    xmax = property(_get_xmax, _set_xmax) #: (*fl64*) Maximum nodal x-coordinate for grid.
    def _get_ymin(self): return self._ymin
    def _set_ymin(self,value): self._ymin = value
    ymin = property(_get_ymin, _set_ymin) #: (*fl64*) Minimum nodal y-coordinate for grid.
    def _get_ymax(self): return self._ymax
    def _set_ymax(self,value): self._ymax = value
    ymax = property(_get_ymax, _set_ymax) #: (*fl64*) Maximum nodal y-coordinate for grid.
    def _get_zmin(self): return self._zmin
    def _set_zmin(self,value): self._zmin = value
    zmin = property(_get_zmin, _set_zmin) #: (*fl64*) Minimum nodal z-coordinate for grid.
    def _get_zmax(self): return self._zmax
    def _set_zmax(self,value): self._zmax = value
    zmax = property(_get_zmax, _set_zmax) #: (*fl64*) Maximum nodal z-coordinate for grid.
    #def _get_information(self):
    #    print('FEHM contour output - format '+self._format)
    #    print(' call format: fcontour[time][variable][node_index-1]')
    #    prntStr =  '    times ('+str(len(self.times))+'): '
    #    for time in self.times: prntStr += str(time)+', '
    #    print(prntStr[:-2]+' days')
    #    prntStr = ' variables: '
    #    for var in self.variables: prntStr += str(var)+', '
    #    for var in self.user_variables: prntStr += str(var)+', '
    #    print(prntStr)
    #what = property(_get_information) #:(*str*) Print out information about the fcontour object.
class fhistory(object): 
    '''
    ****************************************************************
    Reading and plotting methods associated with history output data.
    This includes file types his & temp.his.
    ****************************************************************
    '''
    def __init__(self,filename=None,verbose=True):
        self._filename=None 
        #self._silent = dflt.silent
        self._silent = True
        self._format = ''
        self._times=[]  
        self._verbose = verbose
        self._data={}
        self._row=None
        self._nodes=[]  
        self._zones = []
        self._variables=[] 
        self._user_variables = []
        self._keyrows={}
        self.column_name=[]
        self.num_columns=0
        self._nkeys=1
        filename = os_path(filename)
        if filename: self._filename=filename; self.read(filename)
    def __getitem__(self,key):
        if key in self.variables or key in self.user_variables:
            return self._data[key]
        else: return None
    #def __repr__(self): 
    #    retStr =  'History output for variables '
    #    for var in self.variables:
    #        retStr += var+', '
    #    retStr = retStr[:-2] + ' at '
    #    if len(self.nodes)>10:
    #        retStr += str(len(self.nodes)) + ' nodes.'
    #    else:
    #        if len(self.nodes)==1:
    #            retStr += 'node '
    #        else:
    #            retStr += 'nodes '
    #        for nd in self.nodes:
    #            retStr += str(nd) + ', '
    #        retStr = retStr[:-2] + '.'
    #    return retStr
    def read(self,filename):                        # read contents of file
        '''
        Read in FEHM history output information.
        
        :param filename: File name for output data, can include wildcards to define multiple output files.
        :type filename: str
        '''
        from glob import glob
        import re
        glob_pattern = re.sub(r'\[','[[]',filename)
        glob_pattern = re.sub(r'(?<!\[)\]','[]]', glob_pattern)
        files=glob(glob_pattern)
        configured=False
        for i,fname in enumerate(files):
            #print('Names of files: \n', files)
            #if self._verbose:
            if '..' in fname:
                #print('fname is in compare: ', fname)
                if os.name == 'nt':  # For Windows
                    tmp=fname.split('\\')[-1]
                else:
                    tmp = fname.split('/')[-1]
                #print('tmp: ', tmp)
                if os.path.exists(tmp):
                    #print('valid comparison of: ', fname, ' and ', tmp)
                    pass
                else:
                    print('\n    **WARNING**   The file ', tmp , ' is in compare, but not in output. Skipping file...')
                    continue
            elif '..' not in fname:
                #print('fname is in output: ', fname)
                path=os.path.join('..', 'compare', '')+fname
                #print('PATH', path)
                if os.path.exists(path):
                    #print('valid comparison of: ', fname, ' and ', path)
                    pass
                else:
                    #print('\n    **WARNING**   The file ', fname , 'is in output, but doesn''t have a valid compare file. Skipping file...')
                    continue

            with open(fname, 'r') as self._file:
                header=self._file.readline()
                if header.strip()=='': continue             # empty file
                self._detect_format(header)
                if self._format=='tec': 
                    header=self._file.readline()
                    if header.strip()=='': continue         # empty file
                    i = 0; sum_file = False
                    while not header.startswith('variables'): 
                        header=self._file.readline()
                        i = i+1
                        if i==10: sum_file=True; break
                    if sum_file: continue
                    self._setup_headers_tec(header)
                elif self._format=='surf': 
                    self._setup_headers_surf(header)
                elif self._format=='default': 
                    header=self._file.readline()
                    header=self._file.readline()
                    if header.strip()=='': continue         # empty file
                    i = 0; sum_file = False
                    while 'Time' not in header: 
                        header=self._file.readline()
                        i = i+1
                        if i==10: sum_file=True; break
                    if sum_file: 
                        continue
                    #print('header is: ', header)
                    self._setup_headers_default(header)
                else: print('Unrecognised format');return
                if not configured:
                    #print('NOT CONFIGURED', len(self.nodes)+1)
                    self.num_columns = len(self.nodes)+1
                if self.num_columns>0: configured=True
                if self._format=='tec':
                    self._read_data_tec(fname.split('_')[-2])
                elif self._format=='surf':
                    self._read_data_surf(fname.split('_')[-2])
                elif self._format=='default':
                    if 'temp' in fname:
                        self._read_data_default(fname.split('_')[1].split('.')[0])
                    else: 
                        self._read_data_default(fname.split('_')[-1].split('.')[0])
                
                #self._get_information()    # Uncomment to print information about each test.
                #self._file.close()
    def _detect_format(self,header):
        if header.startswith('TITLE'):
            self._format = 'tec'
        elif header.startswith('Time '):
            self._format = 'surf'
        else:
            self._format = 'default'
    def _setup_headers_tec(self,header):
        header=header.split('" "Node')
        if self.nodes: return
        for key in header[1:-1]: self._nodes.append(int(key))
        self._nodes.append(int(header[-1].split('"')[0]))
    def _setup_headers_surf(self,header):
        header=header.split(', Node')
        if self.nodes: return
        for key in header[1:]: self._nodes.append(int(key))
    def _setup_headers_default(self,header):
        #print('\n_setup_headers_default', header)
        header=header.split(' Node')
        #print('\nheaders after split', header)
        if self.nodes: return
        for key in header[1:]: self._nodes.append(int(key))
    def _read_data_tec(self,var_key):
        self._variables.append(hist_var_names[var_key])
        lns = self._file.readlines()
        i = 0
        while lns[i].startswith('text'): i+=1
        data = []
        for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
        data = np.array(data)
        if data[-1,0]<data[-2,0]: data = data[:-1,:]
        self._times = np.array(data[:,0])
        self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
    def _read_data_surf(self,var_key):
        self._variables.append(hist_var_names[var_key])
        lns = self._file.readlines()
        data = []
        for ln in lns: data.append([float(d) for d in ln.strip().split(',')])
        data = np.array(data)
        if data[-1,0]<data[-2,0]: data = data[:-1,:]
        self._times = np.array(data[:,0])
        self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
    def _read_data_default(self,var_key):
        self._variables.append(hist_var_names[var_key])
        #print('These are the variables: ', self._variables)
        lns = self._file.readlines()
        if lns == []:
            #print('No time values in file')
            data = []
            #print('self._node', self.nodes)
            i = 0
            while i <= len(self._nodes):
                data.append(np.zeros([1] , dtype=float))
                i+=1
            #print('data', data)
            self._times = np.array(data[0],  dtype=float)
            #print('times', self._times)
            #print('dict([(node,data[icol+1]) for icol,node in enumerate(self.nodes)])', dict([(node,data[icol+1]) for icol,node in enumerate(self.nodes)]))
            self._data[hist_var_names[var_key]] = dict([(node,data[icol+1]) for icol,node in enumerate(self.nodes)])
            #print('self._data[hist_var_names[var_key]]', self._data[hist_var_names[var_key]])
        else:
            #print('lns: ', lns)
            data = []
            for ln in lns: data.append([float(d) for d in ln.strip().split()])
            data = np.array(data)
            #print('data [-1,0]', data[-1,0])
            #print('data [-2,0]', data[-2,0])
            if data[-1,0]<data[-2,0]: data = data[:-1,:]
            self._times = np.array(data[:,0])
            #print(' TIMES ', self._times)
            #print('dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])', dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)]))
            self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
    def _get_variables(self): return self._variables
    variables = property(_get_variables)#: (*lst[str]*) List of variables for which output data are available.
    def _get_user_variables(self): return self._user_variables
    def _set_user_variables(self,value): self._user_variables = value
    user_variables = property(_get_user_variables, _set_user_variables) #: (*lst[str]*) List of user-defined variables for which output data are available.
    def _get_format(self): return self._format
    format = property(_get_format) #: (*str*) Format of output file, options are 'tec', 'surf', 'avs' and 'avsx'.
    def _get_filename(self): return self._filename
    filename = property(_get_filename)  #: (*str*) Name of FEHM contour output file. Wildcards can be used to define multiple input files.
    def _get_times(self): return np.sort(self._times)
    times = property(_get_times)    #: (*lst[fl64]*) List of times (in seconds) for which output data are available.
    def _get_nodes(self): return self._nodes
    nodes = property(_get_nodes)    #: (*lst[fl64]*) List of node indices for which output data are available.
    def _get_information(self):
        print('FEHM history output - format '+self._format)
        print('call format: fhistory[variable][node][time_index]')
        prntStr = ' nodes: '
        for nd in self.nodes: prntStr += str(nd)+', '
        print(prntStr)
        prntStr =  'times \n('+str(len(self.times))+'): '
        for time in self.times: prntStr += str(time)+', '
        print(prntStr[:-2]+' days')
        prntStr = 'variables: '
        for var in self.variables: prntStr += str(var)+', '
        print(prntStr)
    what = property(_get_information) #:(*str*) Print out information about the fhistory object.
class fzoneflux(fhistory):
    '''
    Derived class of fhistory, for zoneflux output
    Zone flux history output information object.
    '''
#   __slots__ = ['_filename','_times','_verbose','_data','_row','_zones','_variables','_keyrows','column_name','num_columns','_nkeys']
    def __init__(self,filename=None,verbose=True):
        super(fzoneflux,self).__init__(filename, verbose)
        self._filename=None 
        self._times=[]  
        self._verbose = verbose
        self._data={}
        self._row=None
        self._zones=[]  
        self._variables=[] 
        self._keyrows={}
        self.column_name=[]
        self.num_columns=0
        self._nkeys=1
        if filename: self._filename=filename; self.read(filename)
    def _setup_headers_tec(self,header):
        'placeholder'
    def _read_data_tec(self,var_key):
        zn = int(var_key[-5:])
        if var_key.startswith('c'):
            if zn not in self._zones: self._zones.append(zn)
            if 'co2_source' not in self._variables:
                self._variables += flxz_co2_names
                for var in flxz_co2_names: self._data[var] = {}
                
            lns = self._file.readlines()
            i = 0
            while lns[i].startswith('text'): i+=1
            data = []
            for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
            data = np.array(data)
            if data[-1,0]<data[-2,0]: data = data[:-1,:]
            self._times = np.array(data[:,0])
            for j,var_key in enumerate(flxz_co2_names):
                self._data[var_key].update(dict([(zn,data[:,j+1])]))
        elif var_key.startswith('w'):
            if zn not in self._zones: self._zones.append(zn)
            if 'water_source' not in self._variables:
                self._variables += flxz_water_names
                for var in flxz_water_names: self._data[var] = {}
                
            lns = self._file.readlines()
            i = 0
            while lns[i].startswith('text'): i+=1
            data = []
            for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
            data = np.array(data)           
            if data[-1,0]<data[-2,0]: data = data[:-1,:]
            self._times = np.array(data[:,0])
            for j,var_key in enumerate(flxz_water_names):
                self._data[var_key].update(dict([(zn,data[:,j+1])]))
        elif var_key.startswith('v'):
            if zn not in self._zones: self._zones.append(zn)
            if 'vapor_source' not in self._variables:
                self._variables += flxz_vapor_names
                for var in flxz_vapor_names: self._data[var] = {}
                
            lns = self._file.readlines()
            i = 0
            while lns[i].startswith('text'): i+=1
            data = []
            for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
            data = np.array(data)           
            if data[-1,0]<data[-2,0]: data = data[:-1,:]
            self._times = np.array(data[:,0])
            for j,var_key in enumerate(flxz_vapor_names):
                self._data[var_key].update(dict([(zn,data[:,j+1])]))
    def _read_data_surf(self,var_key):
        self._variables.append(hist_var_names[var_key])
        lns = self._file.readlines()
        data = []
        for ln in lns[i:]: data.append([float(d) for d in ln.strip().split(',')])
        data = np.array(data)           
        if data[-1,0]<data[-2,0]: data = data[:-1,:]
        self._times = np.array(data[:,0])
        self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
    def _read_data_default(self,var_key):
        self._variables.append(hist_var_names[var_key])
        lns = self._file.readlines()
        data = []
        for ln in lns[i:]: data.append([float(d) for d in ln.strip().split()])
        data = np.array(data)
        if data[-1,0]<data[-2,0]: data = data[:-1,:]
        self._times = np.array(data[:,0])
        self._data[hist_var_names[var_key]] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])
    def _get_zones(self): return self._zones
    def _set_zones(self,value): self._zones = value
    zones = property(_get_zones, _set_zones) #: (*lst[int]*) List of zone indices for which output data are available.
class fnodeflux(object):                    
    '''
        Reading and plotting methods associated with internode flux files.
        Internode flux information.
        
        Can read either water or CO2 internode flux files.
        
        The fnodeflux object is indexed first by node pair - represented as a tuple of node indices - and then
        by either the string 'liquid' or 'vapor'. Data values are in time order, as given in the 'times' attribute.
    '''
    def __init__(self,filename=None):
        self._filename = filename
        self._silent = True #dflt.silent
        self._nodepairs = []
        self._times = []
        self._timesteps = []
        self._data = {}
        if self._filename: self.read(self._filename)
    def __getitem__(self,key):
        if key in self.nodepairs:
            return self._data[key]
        else: return None
    def read(self,filename):
        '''Read in FEHM contour output information.
        
        :param filename: File name for output data, can include wildcards to define multiple output files.
        :type filename: str
        '''
        if not os.path.isfile(filename):
            print('ERROR: cannot find file at '+filename)
            return
        fp = open(filename)
        lns = fp.readlines()
        N = int(lns[0].split()[1])
        
        data = np.zeros((N,len(lns)/(N+1),2))       # temporary data storage struc
        
        for ln in lns[1:N+1]:
            ln = ln.split()
            self._nodepairs.append((int(float(ln[0])),int(float(ln[1]))))   # append node pair
        
        for i in range(len(lns)/(N+1)):
            ln = lns[(N+1)*i:(N+1)*(i+1)]
            
            nums = ln[0].split()
            self._timesteps.append(float(nums[2]))
            self._times.append(float(nums[3]))
            for j,lni in enumerate(ln[1:]):
                lnis = lni.split()
                data[j,i,0] = float(lnis[2])
                data[j,i,1] = float(lnis[3])
                
        for i,nodepair in enumerate(self.nodepairs):
            self._data[nodepair] = dict([(var,data[i,:,icol]) for icol,var in enumerate(['vapor','liquid'])])
            
    def _get_filename(self): return self._filename
    def _set_filename(self,value): self._filename = value
    filename = property(_get_filename, _set_filename) #: (*str*) filename target for internode flux file.
    def _get_timesteps(self): return np.sort(self._timesteps)
    def _set_timesteps(self,value): self._timesteps = value
    timesteps = property(_get_timesteps, _set_timesteps) #: (*lst*) timestep for which node flux information is reported.
    def _get_times(self): return np.sort(self._times)
    def _set_times(self,value): self._times = value
    times = property(_get_times, _set_times) #: (*lst*) times for which node flux information is reported.
    def _get_nodepairs(self): return self._nodepairs
    def _set_nodepairs(self,value): self._nodepairs = value
    nodepairs = property(_get_nodepairs, _set_nodepairs) #: (*lst*) node pairs for which node flux information is available. Each node pair is represented as a two item tuple of node indices.
class ftracer(fhistory):                    
    '''
    Derived class of fhistory, for tracer output
    Tracer history output information object.
    '''
    def __init__(self,filename=None,verbose=True):
        super(ftracer,self).__init__(filename, verbose)
        self._filename=None 
        self._times=[]  
        self._verbose = verbose
        self._data={}
        self._row=None
        self._nodes=[]  
        self._variables=[] 
        self._keyrows={}
        self.column_name=[]
        self.num_columns=0
        self._nkeys=1
        if filename: self._filename=filename; self.read(filename)
    def _read_data_default(self,var_key):
        try: var_key = hist_var_names[var_key]
        except: pass
        self._variables.append(var_key)
        lns = self._file.readlines()
        data = []
        for ln in lns: data.append([float(d) for d in ln.strip().split()])
        data = np.array(data)
        if data[-1,0]<data[-2,0]: data = data[:-1,:]
        self._times = np.array(data[:,0])
        self._data[var_key] = dict([(node,data[:,icol+1]) for icol,node in enumerate(self.nodes)])

class fptrk(fhistory):                      
    '''
    Derived class of fhistory, for particle tracking output
    Tracer history output information object.
    '''
    def __init__(self,filename=None,verbose=True):
        super(fptrk,self).__init__(filename, verbose)
        self._filename=None 
        self._silent = True#dflt.silent
        self._times=[]  
        self._verbose = verbose
        self._data={}
        self._row=None
        self._nodes=[0] 
        self._variables=[] 
        self._keyrows={}
        self.column_name=[]
        self.num_columns=0
        self._nkeys=1
        if filename: self._filename=filename; self.read(filename)
    def read(self,filename):                        # read contents of file
        '''
        Read in FEHM particle tracking output information. Index by variable name.
        
        :param filename: File name for output data, can include wildcards to define multiple output files.
        :type filename: str
        '''
        from glob import glob
        import re
        glob_pattern = re.sub(r'\[','[[]',filename)
        glob_pattern = re.sub(r'(?<!\[)\]','[]]', glob_pattern)
        files=glob(glob_pattern)
        configured=False
        for i,fname in enumerate(files):
            #if self._verbose:
            #   print(fname)
            with open(fname, 'r') as self._file:
                header=self._file.readline()
                if header.strip()=='': continue             # empty file
                header=self._file.readline()
                self._setup_headers_default(header)
                self._read_data_default()
                #self._file.close()
    def _setup_headers_default(self,header):
        header=header.strip().split('"')[3:-1]
        header = [h for h in header if h != ' ']
        for var in header: self._variables.append(var)
    def _read_data_default(self):
        lns = self._file.readlines()
        data = []
        for ln in lns: data.append([float(d) for d in ln.strip().split()])
        data = np.array(data)
        if data[-1,0]<data[-2,0]: data = data[:-1,:]
        self._times = np.array(data[:,0])
        self._data = dict([(var,data[:,icol+1]) for icol,var in enumerate(self.variables)])

#### ADDING NEW CLASS FOR COMPARISONS ####
class fcomparison(object): 
    '''
    ****************************************************************
    Reading and plotting methods associated with history output data.
    This includes file types his & temp.his.
    ****************************************************************
    '''
    def __init__(self, filename=None, verbose=True):
        self._filename = None 
        self._silent = True
        self._format = ''
        self._times = [] 
        self._info = [] 
        self._verbose = verbose
        self._data = {}
        self._row = None
        self._nodes = []  
        self._zones = []
        self._variables = [] 
        self._user_variables = []
        self._keyrows = {}
        self.column_name = []
        self.num_columns = 0
        self._nkeys = 1
        self.files_info = [] 
        filename = os_path(filename)
        if filename: 
            self._filename = filename
            self.read(filename)

    def __getitem__(self, key):
        if key in self._variables or key in self._user_variables:
            return self._data.get(key, None)
        else:
            return None

    def read(self, filename):  # read contents of file
        '''
        Read in FEHM history output information.
        
        :param filename: File name for output data, can include wildcards to define multiple output files.
        :type filename: str
        '''
        from glob import glob
        import re
        glob_pattern = re.sub(r'\[','[[]',filename)
        glob_pattern = re.sub(r'(?<!\[)\]','[]]', glob_pattern)
        files = glob(glob_pattern)

        # Sort the files to ensure consistent order across platforms
        files.sort()

        for fname in files:
            if '..' in fname:
                if os.name == 'nt':  # For Windows
                    tmp = fname.split('\\')[-1]
                else:
                    tmp = fname.split('/')[-1]
                if os.path.exists(tmp):
                    pass
                else:
                    continue
            elif '..' not in fname:
                path = os.path.join('..', 'compare', '') + fname
                if os.path.exists(path):
                    pass
                else:
                    continue
            
            with open(fname, 'r') as file:
                self._read_data(file, fname)
                
    def _read_data(self, file, filename):
        lines = file.readlines()[2:]
        data = [line.strip().split() for line in lines]
        #print('DATA', data)

        max_length = max(len(sublist) for sublist in data)
        data = [sublist + [''] * (max_length - len(sublist)) for sublist in data]

        def to_float(value):
            try:
                return float(value)
            except ValueError:
                return 0

        data = np.array([[to_float(item) for item in sublist] for sublist in data], dtype=float)
        self.files_info.append((filename, data))  # Store the data as a tuple (filename, data)
        #print('Processed file: {filename}')
        #print(f'Data for {filename}: {data}')


    #########################################################################################                    


def fdiff( in1, in2, format='diff', times=[], variables=[], components=[], nodes=[], info=[]):
    '''
    ****************************************************************

    Fdiff Takes the difference of two fpost objects

    ****************************************************************
    
    :param in1: First fpost object
    :type filename: fpost object (fcontour)
    :param in2: First fpost object
    :type filename: fpost object (fcontour)
    :param format: Format of diff: diff->in1-in2 relative->(in1-in2)/abs(in2) percent->100*abs((in1-in2)/in2)
    :type format: str
    :param times: Times to diff
    :type times: lst(fl64)
    :param variables: Variables to diff
    :type variables: lst(str)
    :param components: Components to diff (foutput objects)
    :type components: lst(str)
    :returns: fpost object of same type as in1 and in2
    '''
    
    # Copy in1 and in2 in case they get modified below
    in1 = copy(in1)
    in2 = copy(in2)

    if type(in1) is not type(in2):
        print("ERROR: fpost objects are not of the same type: "+str(type(in1))+" and "+str(type(in2)))
        return
    if isinstance(in1, fcontour) or isinstance(in1, fhistory) or 'foutput' in str(in1.__class__):
        # Find common timesclear
        t = np.intersect1d(in1.times,in2.times)
        #print('lenth of t',len(t), 'length of times', len(times) )
        if len(t) == 0:
            print("ERROR: fpost object times do not have any matching values")
            return
        if len(times) > 0:
            times = np.intersect1d(times,t)
            if len(times) == 0:
                print("ERROR: provided times are not coincident with fpost object times")
                return
        else:
            #print('Time reached else.')
            times = t
            
    if isinstance(in1, fcomparison):
        atol = 1e-10  # Absolute tolerance
        rtol = 1e-05  # Relative tolerance

        out = copy(in1)
        out._info = []

        if len(in1.files_info) != len(in2.files_info):
            print('ERROR: The number of files in in1 and in2 do not match')
            return


        #print('in1.files_info: {in1.files_info}')
        #print('in2.files_info: {in2.files_info}')

        for (file1, data1), (file2, data2) in zip(in1.files_info, in2.files_info):
            print(f'\nComparing Files... {file1} with {file2}')
            #print(f'Type of data1: {type(data1)}, Type of data2: {type(data2)}')
            #print(f'Data1: {data1}')
            #print(f'Data2: {data2}')

            if data1.size == 0 or data2.size == 0:
                print(f'Missing data in comparison for files: {file1} and {file2}')
                continue

            data1_flat = data1.flatten()
            data2_flat = data2.flatten()

            #print('Data1 Flat:', data1_flat, len(data1_flat))
            #print('Data2 Flat:', data2_flat, len(data2_flat))

            if len(data1_flat) == 0 or len(data2_flat) == 0:
                print(f'No data to compare in files: {file1} and {file2}')
                continue

            # Check if arrays are close enough within the tolerance
            if np.allclose(data1_flat, data2_flat, rtol=rtol, atol=atol, equal_nan=True):
                print(f'Files are EQUAL')
                out._info.append(0)
            else:
                # Calculate the differences
                difference = np.abs(data1_flat - data2_flat)
                max_difference = np.max(difference)
                
                if max_difference <= atol + rtol * np.abs(data1_flat).max():
                    print(f'Files are considered EQUAL within a relative tolerance of:{rtol}, and an absolute tolerance of:{atol}')
                    #print(f'Max difference: {max_difference}\n')
                    out._info.append(0)
                else:
                    differences = np.setdiff1d(data1_flat, data2_flat)
                    differences2 = np.setdiff1d(data2_flat, data1_flat)

                    # print(f'Differences: {differences}')
                    # print(f'Differences2: {differences2}')

                    print(f'Files are UNEQUAL.')
                    print(f'Found {len(differences)} differences in {file1} and {file2}.\nDifferences are {differences} and {differences2}')
                    out._info.append(1)
                    
        return out

    elif isinstance(in1, fcontour):
        # Find common variables
        v = np.intersect1d(in1.variables,in2.variables)
        if len(v) == 0:
            print("ERROR: fcontour object variables do not have any matching values")
            return
        if len(variables) > 0:
            variables = np.intersect1d(variables,v)
            if len(variables) == 0:
                print("ERROR: provided variables are not coincident with fcontour object variables")
                return
        else:
            variables = v
        out = copy(in1)
        out._times = times
        out._variables = variables
        out._data = {}
        for t in times:
            if format == 'diff':
                out._data[t] = dict([(v,in1[t][v] - in2[t][v]) for v in variables])
            elif format == 'relative':
                out._data[t] = dict([(v,(in1[t][v] - in2[t][v])/np.abs(in2[t][v])) for v in variables])
            elif format == 'percent':
                out._data[t] = dict([(v,100*np.abs((in1[t][v] - in2[t][v])/in2[t][v])) for v in variables])
        return out
    
    #Takes the difference of two fhistory objects.  
    elif isinstance(in1, fhistory):
        # Find common variables
        v = np.intersect1d(in1.variables, in2.variables)
        if len(v) == 0:
            print("ERROR: fhistory object variables do not have any matching values")
            return
        if len(variables) > 0:
            variables = np.intersect1d(variables,v)
            #print('variables: ', variables)
            if len(variables) == 0:
                print("ERROR: provided variables are not coincident with fhistory object variables")
                return
        else:
            variables = v

        #Find common nodes.
        n = np.intersect1d(in1.nodes, in2.nodes)
        #print('int.nodes', in1.nodes, 'in2.nodes', in2.nodes)
        #print('n', len(n))
        if len(n) == 0:
            print("ERROR: fhistory object nodes do not have any matching values")
            return
        if len(nodes) > 0:
            nodes = np.intersect1d(nodes,n)
            if len(nodes) == 0:
                print("ERROR: provided nodes are not coincident with fhistory object nodes")
                return
        else:
            nodes = n
            #print('nodes from else: ', nodes)
        
        #Set up the out object.
        out = copy(in1)
        out._times = times
        out._variables = variables
        out._nodes = nodes
        out._data = {}
        
        #Find the difference at each time index for a variable and node.
        for v in variables:
            for n in nodes:
                i = 0
                diff = []
                for i in range(len(times)):
                    if format == 'diff':
                        #Quick fix to handle ptrk files.
                        if isinstance(in1, fptrk):
                            diff.append(in1[v][n]-in2[v][n])
                        else:
                            try:
                                diff.append(in1[v][n][i] - in2[v][n][i])
                            except IndexError:
                                #print(f"File is truncated at line {i}")
                                continue  
                    elif format == 'relative':
                        diff.append((in1[t][v] - in2[t][v])/np.abs(in2[t][v])) 
                    elif format == 'percent':
                        diff.append(100*np.abs((in1[t][v] - in2[t][v])/in2[t][v]))  
                    i = i + 1
                if isinstance(in1, fptrk):
                    out._data[v] = np.array(diff)
                    #print('out', out._data[v])
                else:
                    out._data[v] = dict([(n, diff)])
                    #print('out', out._data[v])
                
        #Return the difference.
        #print('out', out._data[v])     
        return out
            
    elif 'foutput' in str(in1.__class__):
        # Find common components
        c = np.intersect1d(in1.components,in2.components)
        if len(c) == 0:
            print("ERROR: foutput object components do not have any matching values")
            return
        if len(components) > 0:
            components = np.intersect1d(components,c)
            if len(components) == 0:
                print("ERROR: provided components are not coincident with foutput object components")
                return
        else:
            components = c      
        out = deepcopy(in1)
        out._times = times
        out._node = {}
        for tp in ['water','gas','tracer1','tracer2']:      
            out._node[tp] = None
        for cp in components:
            if format == 'diff':
                if len(variables):
                    out._node[cp] = dict([(n,dict([(v,np.array(in1._node[cp][n][v]) - np.array(in2._node[cp][n][v])) for v in list(in1._node[cp][n].keys()) if v in variables])) for n in in1.nodes])
                else:
                    out._node[cp] = dict([(n,dict([(v,np.array(in1._node[cp][n][v]) - np.array(in2._node[cp][n][v])) for v in list(in1._node[cp][n].keys())])) for n in in1.nodes])
            elif format == 'relative':
                if len(variables):
                    out._node[cp] = dict([(n,dict([(v,(np.array(in1._node[cp][n][v]) - np.array(in2._node[cp][n][v]))/np.abs(in2._node[cp][n][v])) for v in list(in1._node[cp][n].keys()) if v in variables])) for n in in1.nodes])
                else:
                    out._node[cp] = dict([(n,dict([(v,(np.array(in1._node[cp][n][v]) - np.array(in2._node[cp][n][v]))/np.abs(in2._node[cp][n][v])) for v in list(in1._node[cp][n].keys())])) for n in in1.nodes])
            elif format == 'percent':
                if len(variables):
                    out._node[cp] = dict([(n,dict([(v,100*np.abs((np.array(in1._node[cp][n][v]) - np.array(in2._node[cp][n][v]))/in2._node[cp][n][v])) for v in list(in1._node[cp][n].keys()) if v in variables])) for n in in1.nodes])
                else:
                    out._node[cp] = dict([(n,dict([(v,100*np.abs((np.array(in1._node[cp][n][v]) - np.array(in2._node[cp][n][v]))/in2._node[cp][n][v])) for v in list(in1._node[cp][n].keys())])) for n in in1.nodes])
                        
        return out
        
def sort_tec_files(files):
    # sort first by number, then by type
    #from string import join
    for file in files:
        if not file.endswith('.dat'): return files
    paths = [os.sep.join(file.split(os.sep)[:-1]) for file in files]
    files = [file.split(os.sep)[-1] for file in files]
    
    times = []
    for file in files: 
        for type in ['_days_sca_node','_days_vec_node','_days_hf_node','_days_con_node','_sca_node','_vec_node','_hf_node','_con_node']:
            if type in file:                
                times.append(float('.'.join(file.split(type)[0].split('.')[1:])))
                break
    times = sorted(enumerate(times), key=lambda x: x[1])
    
    paths = [paths[ind] for ind,time in times]
    files = [files[ind] for ind,time in times]

    return [path+os.sep+file if path else file for path,file in zip(paths,files)]
    
class foutput(object):
    def __init__(self,filename = None, input=None, grid = None, hide = True, silent=True, write = False):
        self._filename = filename
        self._silent = True#dflt.silent
        if self._filename:
            if input and grid:
                diag = process_output(filename, hide = hide,silent = silent,input=input,grid=grid,write=write)
            else:
                diag = process_output(filename, hide = hide,silent=silent,write=write)      
        self._node = deepcopy(diag.node)
        self._times = deepcopy(diag.time.data[1:])
    def _get_node(self): return self._node
    node = property(_get_node) #: (*dict*) Dictionary of node output, keyed first on component ('water','gas','tracer1'), then on node number, then on variable.
    def _get_nodes(self): 
        for type in ['water','gas','tracer1']:
            if self._node[type] == None: continue
            nds = list(self._node[type].keys())
            nds.sort()
            return nds
        return None
    nodes = property(_get_nodes)
    
    def _get_variables(self):
        """ 
        Get Variables
        
        Returns the variables in the foutput object or returns None if no 
        variables. 
        """

        for type in ['water','gas','tracer1']:
            for node in self.nodes:
                if self._node[type][node] == None: 
                    continue
                    
                vrbls = list(self._node[type][node].keys())
                vrbls.sort()
                
                return vrbls
            
        return None

    variables = property(_get_variables) #: (*lst*) List of variables available for the various components.
    
    def _get_components(self): 
        cpts = []
        for type in ['water','gas','tracer1','tracer2']:        
            if self._node[type] != None: cpts.append(type)
        return cpts
    components = property(_get_components) #: (*lst*) Component names for which nodal information available
    
    def _get_times(self): return self._times
    times = property(_get_times) #: (*ndarray*) Vector of output times.
    def _get_filename(self): return self._filename
    def _set_filename(self,value): self._filename = value
    filename = property(_get_filename, _set_filename) #: (*str*) Name of output file
    def _get_information(self):
        print('FEHM output file \''+self.filename+'\'')
        print('    call format: foutput.node[component][node][variable]')
        prntStr =  '    components: '
        for cpt in self.components: prntStr += str(cpt)+', '
        print(prntStr[:-2])
        prntStr = '    nodes: '
        for nd in self.nodes: prntStr += str(nd)+', '
        print(prntStr[:-2])
    what = property(_get_information) #:(*str*) Print out information about the fcontour object.
    
def os_path(path):
    if WINDOWS: path = path.replace('/','\\')
    else: path = path.replace('\\','/')
    return path
