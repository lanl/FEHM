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

"""Classes and methods for interactive assistance with PyFEHM."""

def textline(symbol,N):
	out = ''
	for i in range(N): out += symbol
	return out
class fhelp(object):
	def __init__(self, parent = None):
		self._parent = parent
	def __repr__(self): 
		outStr = 'PyFEHM PaperClip assistant: "It looks like you\'re building a model..."\n'
		outStr+='\n' 
		outStr+='		        _.-;:q=._                                  \n'
		outStr+='                      .\' j=""^k;:\\.                                 \n'
		outStr+='                     ; .F       ";`Y                                \n'
		outStr+='                    ,;.J_        ;\'j                                \n'
		outStr+='                  ,-;"^7F       : .F           _________________    \n'
		outStr+='                 ,-\'-_<.        ;gj. _.,---""\'\'               .\'    \n'
		outStr+='                ;  _,._`\.     : `T"5,                       ;      \n'
		outStr+='                : `?8w7 `J  ,-\'" -^q. `                     ;       \n'
		outStr+='                 \;._ _,=\' ;   n58L Y.                     .\'       \n'
		outStr+='                   F;";  .\' k_ `^\'  j\'                     ;        \n'
		outStr+='                   J;:: ;     "y:-=\'                      ;         \n'
		outStr+='                    L;;==      |:;   jT\                  ;         \n'
		outStr+='                    L;:;J      J:L  7:;\'       _         ;          \n'
		outStr+='                    I;|:.L     |:k J:.\' ,  \'       .     ;          \n'
		outStr+='                    |;J:.|     ;.I F.:      .           :           \n'
		outStr+='                   ;J;:L::     |.| |.J  , \'   `    ;    ;           \n'
		outStr+='                 .\' J:`J.`.    :.J |. L .    ;         ;            \n'
		outStr+='                ;    L :k:`._ ,\',j J; |  ` ,        ; ;             \n'
		outStr+='              .\'     I :`=.:."_".\'  L J             `.\'             \n'
		outStr+='            .\'       |.:  `"-=-\'    |.J              ;              \n'
		outStr+='        _.-\'         `: :           ;:;           _ ;               \n'
		outStr+='    _.-\'"             J: :         /.;\'       ;    ;                \n'
		outStr+='  =\'_                  k;.\.    _.;:Y\'     ,     .\'                 \n'
		outStr+='     `"---..__          `Y;."-=\';:=\'     ,      .\'                  \n'
		outStr+='              `""--..__   `"==="\'    -        .\'                    \n'
		outStr+='                       ``""---...__    itz .-\'                      \n'
		outStr+='                                   ``""---\'                         \n'
		outStr+='Available help topics: \n'
		outStr += '   .help.permmodel()'
		return outStr
	def permmodel(self,index=None):
		print ''
		if not index: # summary of available models
			ws = 'FEHM Stress-permeability models'
			print ws
			print textline('-',len(ws))
			print 'Permeability depends dynamically on the evolving stress state, generally changing to reflect some mechanism of mechanical failure.\n'
			ws = 'PyFEHM input'
			print ws
			print textline('-',len(ws))
			print 'Permmodels are defined using the fmodel() class and added via the fdata.add() method.\n'
			print 'Inputs:'
			print 'index - denotes the permmodel to be used, see below for a summary'
			print 'zone - the zone object (or index or name) to which the permmodel is applied'
			print 'param - a dictionary of parameters for the specific permmodel\n'

			ws = 'Summary of available permmodels'
			print ws
			print textline('-',len(ws))
			
			# model 1
			print '1 - no permeability dependence on stress\n'
			
			# model 21
			ws = '21 - permeability and Youngs modulus are modified following Mohr-Coulomb failure on a user specified fracture orientation.\n'			
			print ws
			
			# model 22
			ws = '22 - permeability and Youngs modulus are modified following Mohr-Coulomb failure on the plane of maximum shear.\n'			
			print ws
			
			# model 24
			ws = '24 - permeability is modified due to Mohr-Coulomb failure on an ensemble of fracture planes.\n'			
			print ws
			
			# model 25
			ws = '25 - modification of permmodel 25 to allow user specification of the fracture orientation distribution. \n'
			print ws
			
			ws = 'Example'
			print ws
			print textline('-',len(ws))
			
			ws = 'Assign permmodel 24 to the zone called \'damage\'.\n\n'
			ws += '>>> pm = fmodel(\'permmodel\',index=24)  # create the permmodel object\n'
			ws += '>>> pm.zone = dat.zone[\'damage\']       # assign a zone\n'
			ws += '>>> pm.param[\'shear_frac_tough\'] = 1.e6# assign parameters\n'
			ws += '>>> pm.param[\'static_frict_coef\'] = 0.65\n'
			ws += '>>> ...\n'
			ws += '>>> pm.param[\'frac_dens\'] = 1.\n'
			ws += '>>> dat.add(pm)\n'
			print ws
			
		elif index in [1,21,22,24,25]: 	# summary of specific models
			ws = 'Permmodel '+str(index)
			print ws
			print textline('-',len(ws))
			if index == 21:
				ws = 'Permeability and Youngs modulus are modified following Mohr-Coulomb failure on a user-specified fracture orientation.'			
				ws += ' The user supplies the components of the fracture normal in the model coordinate system. This fracture is assessed for' 
				ws += ' Mohr-Coulomb failure. If failure occurs, permeability multipliers are applied in the fracture frame and then rotated back into the model frame.\n'
				print ws
				ws = 'Parameters'
				print ws
				print textline('-',len(ws))
				ws = 'nx - x-component of fracture normal.\n\n'
				ws += 'ny - y-component of fracture normal.\n\n'
				ws += 'nx - z-component of fracture normal.\n\n'
				ws += 'frict_coef - friction coefficient used in calculation of Mohr-Coulomb failure.\n\n'
				ws += 'cohesion - cohesion of the fracture (MPa).\n\n'
				ws += 'pore_press_factor - multiplying factor for pressure in effective stress calculation, e.g., if 0, pressure in fault is not accounted for in failure.\n\n'
				ws += 'tau_ex_ramp_range - range of excess shear stress over which to ramp permeability.\n\n'
				ws += 'yngs_mod_mult_x - multiplier for Youngs modulus in x-direction.\n\n'
				ws += 'yngs_mod_mult_y - as above in y-direction.\n\n'
				ws += 'yngs_mod_mult_z - as above in z-direction.\n\n'
				ws += 'perm_mult_x - multiplier for permeability in x-direction.\n\n'
				ws += 'perm_mult_y - as above for y-direction.\n\n'				
				ws += 'perm_mult_z - as above for z-direction.\n\n'
				print ws
				
				ws = 'Example'
				print ws
				print textline('-',len(ws))
				
				ws = 'Assign permmodel 21 to the zone called \'damage\'. Youngs modulus is unchanged and permeability'
				ws += ' in the fracture plane (x and y directions) can increase by a factor of 100.\n\n'
				ws += '>>> pm = fmodel(\'permmodel\',index=21)  # create the permmodel object\n'
				ws += '>>> pm.zone = dat.zone[\'damage\']       # assign a zone\n'
				ws += '>>> pm.param[\'nx\'] = 0.86\n'
				ws += '>>> pm.param[\'ny\'] = 0\n'
				ws += '>>> pm.param[\'nz\'] = 0.5\n'
				ws += '>>> pm.param[\'frict_coef\'] = 0.65      # assign parameters\n'
				ws += '>>> pm.param[\'cohesion\'] = 1.\n'
				ws += '>>> pm.param[\'pore_press_factor\'] = 1.\n'
				ws += '>>> pm.param[\'tau_ex_ramp_range\'] = 2.\n'
				ws += '>>> pm.param[\'yngs_mod_mult_x\'] = 1.\n'
				ws += '>>> pm.param[\'yngs_mod_mult_y\'] = 1.\n'
				ws += '>>> pm.param[\'yngs_mod_mult_z\'] = 1.\n'
				ws += '>>> pm.param[\'perm_mult_x\'] = 100.\n'
				ws += '>>> pm.param[\'perm_mult_y\'] = 100.\n'
				ws += '>>> pm.param[\'perm_mult_z\'] = 1.\n'
				ws += '>>> dat.add(pm)\n'
				print ws
			if index == 22:
				ws = 'Permeability and Youngs modulus are modified following Mohr-Coulomb failure on the plane of maximum shear.'			
				ws += ' Permeability enhancement is specified for three axes in a fracture coordinate frame, with the modified permeability' 
				ws += ' tensor subsequently rotated back into the global frame.\n'
				print ws
				ws = 'Parameters'
				print ws
				print textline('-',len(ws))
				ws = 'frict_coef - friction coefficient used in calculation of Mohr-Coulomb failure.\n\n'
				ws += 'cohesion - cohesion of the fracture (MPa).\n\n'
				ws += 'pore_press_factor - multiplying factor for pressure in effective stress calculation, e.g., if 0, pressure in fault is not accounted for in failure.\n\n'
				ws += 'tau_ex_ramp_range - range of excess shear stress over which to ramp permeability.\n\n'
				ws += 'yngs_mod_mult_x - multiplier for Youngs modulus in x-direction.\n\n'
				ws += 'yngs_mod_mult_y - as above in y-direction.\n\n'
				ws += 'yngs_mod_mult_z - as above in z-direction.\n\n'
				ws += 'perm_mult_x - multiplier for permeability in x-direction.\n\n'
				ws += 'perm_mult_y - as above for y-direction.\n\n'				
				ws += 'perm_mult_z - as above for z-direction.\n\n'
				ws += 'por_mult - multiplier for porosity.\n\n'
				ws += 'perm_incremental - boolean, permeability a functino of incremental shear stress.\n\n'
				ws += 'tau_ex_init - excess shear stress that must be exceeded before damage occurs.\n\n'
				print ws
				
				ws = 'Example'
				print ws
				print textline('-',len(ws))
				
				ws = 'Assign permmodel 22 to the zone called \'damage\'. Youngs modulus is unchanged and permeability'
				ws += ' in the fracture plane (x and y directions) can increase by a factor of 100.\n\n'
				ws += '>>> pm = fmodel(\'permmodel\',index=22)  # create the permmodel object\n'
				ws += '>>> pm.zone = dat.zone[\'damage\']       # assign a zone\n'
				ws += '>>> pm.param[\'frict_coef\'] = 0.65      # assign parameters\n'
				ws += '>>> pm.param[\'cohesion\'] = 1.\n'
				ws += '>>> pm.param[\'pore_press_factor\'] = 1.\n'
				ws += '>>> pm.param[\'tau_ex_ramp_range\'] = 2.\n'
				ws += '>>> pm.param[\'yngs_mod_mult_x\'] = 1.\n'
				ws += '>>> pm.param[\'yngs_mod_mult_y\'] = 1.\n'
				ws += '>>> pm.param[\'yngs_mod_mult_z\'] = 1.\n'
				ws += '>>> pm.param[\'perm_mult_x\'] = 100.\n'
				ws += '>>> pm.param[\'perm_mult_y\'] = 100.\n'
				ws += '>>> pm.param[\'perm_mult_z\'] = 1.\n'
				ws += '>>> pm.param[\'por_mult\'] = 1.\n'
				ws += '>>> pm.param[\'perm_incremental\'] = 0\n'
				ws += '>>> pm.param[\'tau_ex_init\'] = 0.5\n'
				ws += '>>> dat.add(pm)\n'
				print ws
			elif index == 24:
				ws = 'Permeability is modified due to Mohr-Coulomb failure on an ensemble of fracture planes.'			
				ws += ' Fracture orientations are uniformly distributed on a steronet. Permeability of each fracture is enhanced in its own plane' 
				ws += ' with the total enhancement given by ensemble average of all fractures rotated back into the global frame. This '
				ws += 'permmodel exhibits anisotropic enhancement due to preferential activation of stress-sensitive fractures.\n'
				print ws
				ws = 'Parameters'
				print ws
				print textline('-',len(ws))
				ws = 'shear_frac_tough - parameter converting excess shear stress to displacement.\n\n'
				ws += 'static_frict_coef - static coefficient of friction on the fractures.\n\n'
				ws += 'dynamic_frict_coef - dynamic coefficient that friction drops to during failure.\n\n'
				ws += 'frac_num - number of fractures per control volume.\n\n'
				ws += 'onset_disp - shear displacement of the fracture at which permeability modification begins.\n\n'
				ws += 'disp_interval - displacement interval over which full permeability enhancement is realised.\n\n'
				ws += 'max_perm_change - maximum change in permeability in log-space (log(multiplier)).\n\n'
				ws += 'frac_cohesion - facture cohesion, used in Mohr-Coulomb failure calculation.\n\n'
				ws += 'frac_dens - density of fractures in the control volume. This number multiplies permeability.\n\n'				
				print ws
				
				ws = 'Example'
				print ws
				print textline('-',len(ws))
				
				ws = 'Assign permmodel 24 to the zone called \'damage\'.\n\n'
				ws += '>>> pm = fmodel(\'permmodel\',index=24)  # create the permmodel object\n'
				ws += '>>> pm.zone = dat.zone[\'damage\']       # assign a zone\n'
				ws += '>>> pm.param[\'shear_frac_tough\'] = 1.e6    # assign parameters\n'
				ws += '>>> pm.param[\'static_frict_coef\'] = 0.65\n'
				ws += '>>> pm.param[\'dynamic_frict_coef\'] = 0.55\n'
				ws += '>>> pm.param[\'frac_num\'] = 100\n'
				ws += '>>> pm.param[\'onset_disp\'] = 1.5e-3\n'
				ws += '>>> pm.param[\'disp_interval\'] = 5.e-3\n'
				ws += '>>> pm.param[\'max_perm_change\'] = 1.8\n'
				ws += '>>> pm.param[\'frac_cohesion\'] = 1.\n'
				ws += '>>> pm.param[\'frac_dens\'] = 1.\n'
				ws += '>>> dat.add(pm)\n'
				print ws
			elif index == 25:
				ws = 'Modification of permmodel 25 to allow user specification of the fracture orientation distribution. When using this permmodel, '			
				ws += 'FEHM looks for an auxilliary file called fracture_orientations.dat, which contains a two column list of fracture dips ' 
				ws += 'and azimuths. These data define a population from which synthetic fracture distributions are constructed for each node. '
				ws += 'See Dempsey et al. (2013, Int. J. Rock Mech. Min. Sci.) for more details.\n'
				print ws
				ws = 'Parameters'
				print ws
				print textline('-',len(ws))
				ws = 'shear_frac_tough - parameter converting excess shear stress to displacement.\n\n'
				ws += 'static_frict_coef - static coefficient of friction on the fractures.\n\n'
				ws += 'dynamic_frict_coef - dynamic coefficient that friction drops to during failure.\n\n'
				ws += 'frac_num - number of fractures per control volume.\n\n'
				ws += 'onset_disp - shear displacement of the fracture at which permeability modification begins.\n\n'
				ws += 'disp_interval - displacement interval over which full permeability enhancement is realised.\n\n'
				ws += 'max_perm_change - maximum change in permeability in log-space (log(multiplier)).\n\n'
				ws += 'frac_cohesion - facture cohesion, used in Mohr-Coulomb failure calculation.\n\n'
				print ws
				
				ws = 'Example'
				print ws
				print textline('-',len(ws))
				
				ws = 'Assign permmodel 25 to the zone called \'damage\'.\n\n'
				ws += '>>> pm = fmodel(\'permmodel\',index=25)  # create the permmodel object\n'
				ws += '>>> pm.zone = dat.zone[\'damage\']       # assign a zone\n'
				ws += '>>> pm.param[\'shear_frac_tough\'] = 1.e6    # assign parameters\n'
				ws += '>>> pm.param[\'static_frict_coef\'] = 0.65\n'
				ws += '>>> pm.param[\'dynamic_frict_coef\'] = 0.55\n'
				ws += '>>> pm.param[\'frac_num\'] = 100\n'
				ws += '>>> pm.param[\'onset_disp\'] = 1.5e-3\n'
				ws += '>>> pm.param[\'disp_interval\'] = 5.e-3\n'
				ws += '>>> pm.param[\'max_perm_change\'] = 1.8\n'
				ws += '>>> pm.param[\'frac_cohesion\'] = 1.\n'
				ws += '>>> dat.add(pm)\n'
				print ws
		else:
			print 'No information available for that permmodel.'