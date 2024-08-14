#***********************************************************************
# Copyright 2014 Los Alamos National Security, LLC All rights reserved
# Unless otherwise indicated, this information has been authored by an
# employee or employees of the Los Alamos National Security, LLC (LANS),
# operator of the Los Alamos National Laboratory under Contract No.
# DE-AC52-06NA25396 with the U.S. Department of Energy. The U.S.
# Government has rights to use, reproduce, and distribute this
# information. The public may copy and use this information without
# charge, provided that this  Notice and any statement of authorship are
# reproduced on all copies. Neither the Government nor LANS makes any
# warranty, express or  implied, or assumes any liability or
# responsibility for the use of this information.      
#***********************************************************************
import unittest
import os
import sys
import argparse
import re
import numpy as np
import glob
import itertools
from subprocess import call
from subprocess import PIPE
from contextlib import contextmanager
import shutil
from tpl_write import tpl_write

try:
    #sys.path.insert(0,os.path.join('pyfehm'))
    import fpost
except ImportError as err:
    print('ERROR: Unable to import pyfehm fpost module')
    print(err)
    os._exit(0)

#Suppresses tracebacks
__unittest = True 

class fehmTest(unittest.TestCase):
    """
    Represents a FEHM test-case. 
    
    Initialize with the name of the test-case folder as a string argument.
    Example: test_case = Tests('saltvcon')
    
    To add the test-case to a suite, initialize a suite object, and call
    addTest() with the test-case as an argument.
    Example: suite = unittest.TestSuite()
             suite.addTest(test_case)
              
    Authors: Dylan Harp, Mark Lange
    Updated: June 2014 
    """    
        
    def __init__(self, testname, log):
        """
        Call the unittest constructor then initialize the main directory and
        log switch values and create a fail log if the log switch is turned on.
        
        :param testname: The name of the test method.
        :type name: str
        
        :param log: Switch that determines if a fail log will be generated.
        :type name: bool
        
        .. Authors: Mark Lange
        .. Updated: July 2014
        """
    
        super(fehmTest, self).__init__(testname)
        self.log = log
        
        #If log switch is on, create the fail log file.
        if self.log:
            self.fail_log = open('fail_log.txt', 'a')
        
        self.maindir = os.getcwd()
    
    # TESTS ######################################################### 
        
    def rad_decay(self):
        """
        **Test radioactive decay option in rxn macro**

        The simulation is a batch reactor without flow and comparison is
        made to the Bateman equation. Decay of 135I->135Xe->135Cs is modeled.
        The test ensures that FEHM matches the bateman equation within 10% 
        relative error for all concentrations greater than 1e-6 moles/kg vapor.

        H. Bateman. "Solution of a System of Differential Equations Occurring in the 
            Theory of Radio-active Transformations," Proc. Cambridge Phil. Soc. IS, 
            423 (1910) https://archive.org/details/cbarchive_122715_solutionofasystemofdifferentia1843

        .. Authors: Dylan Harp, Michelle Bourret
        .. Updated May 2016 by Dylan Harp 
        """

        # Change to test directory
        os.chdir('rad_decay')
        # Import python Bateman equation module
        sys.path.append('.')
        from bateman import bateman

        # Create parameter dictionary with half lives and initial concentations for I, Xe and Cs
        pars = {'thalf_I': 0.00075, 
                'thalf_Xe': 0.001043379,
                'C0_I': 3.48e-4,
                'C0_Xe': 1e-30,
                'C0_Cs': 1e-30} 

        # Create simulation run directory
        output_dir = '_output'
        if os.path.exists( output_dir ): shutil.rmtree(output_dir)
        os.mkdir( output_dir )
        os.chdir(output_dir)

        # Write simulation input file using parameter dictionary
        tpl_write(pars,'../input/run.tpl','run.dat')
        # Run fehm
        self._run_fehm('')

        # Collect results
        CI_fehm = np.genfromtxt('run_135iodine.dat',skip_header=4)
        CXe_fehm = np.genfromtxt('run_135xenon.dat',skip_header=4)
        CCs_fehm = np.genfromtxt('run_135cesium.dat',skip_header=4)
        times = CI_fehm[:,0]/365.

        # Run bateman equation
        thalf = np.array([pars['thalf_I'], pars['thalf_Xe'], 1.])
        lmbda = np.log(2)/thalf
        lmbda[2] = 0
        C0 = np.array([pars['C0_I'], pars['C0_Xe'], pars['C0_Cs']])

        CI_b = bateman(times,[C0[0]],[lmbda[0]])
        CXe_b = bateman(times,C0[0:2],lmbda[0:2])
        CCs_b = bateman(times,C0,lmbda)
        
        for t,f,b in np.column_stack([CI_fehm, CI_b]):
            if b > 1e-6:
                self.assertTrue(abs(f-b)/f<0.1, "Concentration mismatch for Iodine at time %g, FEHM: %g, Bateman: %f"%(t,f,b))

        for t,f,b in np.column_stack([CXe_fehm, CXe_b]):
            if b > 1e-6:
                self.assertTrue(abs(f-b)/f<0.1, "Concentration mismatch for Xenon at time %g, FEHM: %g, Bateman: %f"%(t,f,b))

        for t,f,b in np.column_stack([CCs_fehm, CCs_b]):
            if b > 1e-6:
                self.assertTrue(abs(f-b)/f<0.1, "Concentration mismatch for Cesium at time %g, FEHM: %g, Bateman: %f"%(t,f,b))

        os.chdir(self.maindir)

    def saltvcon(self):
        """
        **Test the Salt Variable Conductivity Macro**
         
        Tests the calculations of thermal conductivity of crushed and intact 
        salt.

        Intact salt:
            ``kxi = k_{t-300}(300/T)^1.14``
            *Munson et al. (1990) Overtest for Simulate Defense High-Level*
            *Waste (Room B): In Situ Data Report. WIPP. Sandia National*
            *Laboratories, SAND89-2671*

        Thermal conductivity of crushed salt from Asse mine:
            ``kx_asse = -270*phi^4+370*phi^3-136*phi^2+1.5*phi+5`` 
            *Bechtold et al. (2004) Backfilling and sealing of underground*
            *respositories for radioactive waste in salt 
            *(BAMBUS II project), EUR 20621, ISBN 92-894-7767-9*

            
        ``kx = (k_{t-300}/kx_asse)*(300/T)^1.14`` if kx is less then 1.e-6, set 
        to 1.e-6. If porosity is greater than 0.4, it is truncated as 0.4 since 
        the Kx relationship is only valid within this range.

        The excel spreadsheet /information/saltvcon.xlsx contains the
        associated calculations. 
        
        .. Authors: Dylan Harp, Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['components'] = ['water']
        #arguments['variables']  = ['Kx','dens','P','sat','T']
        arguments['variables']  = ['Kx']
        arguments['test_measure'] = 'perc_difference' 
        
        #test_case() does not actually display this custom error message yet.
        arguments['err_msg'] = \
            '\nIncorrect intact salt thermal conductivity for node %s \
             \nRelative error: %s, Threshold: %s \
             \nExpected = %s Simulated = %s'
        
        self.test_case('saltvcon', arguments)
        
        
    def dissolution(self):
        """ 
        **Test the Dissoultion Macro**
        
        A one-dimensional transport simulation of calcite (CaC03(s)) 
        dissolution is tested. Profiles of concentration versus reactor 
        length, at selected times, will be compared against the analytical 
        solution.

        Details of this test are described in the FEHM V2.21 Validation Test 
        Plan on pages 93-95 *(STN: 10086-2.21-00, Rev.No. 00, Document ID: 
        10086-VTP-2.21-00, August 2003)* 
        
        .. Authors: Dylan Harp, Mark Lange    
        .. Updated: June 2014 by Mark Lange          
        """
    
        arguments = {}
        arguments['variables'] = ['Np[aq] (Moles/kg H20)']
        
        #test_case() does not actually display this custom error message yet.
        arguments['err_msg'] = '\nIncorrect concentration at time %s'
        
        self.test_case('dissolution', arguments)
    
    def evaporation(self):
        """ 
        **Test the evaporation Macro**

        Evaporation test                                                                
        Comparison of Changes in Mass with Time

        .. Authors: Mark Lange
        .. Updated: 2024 by Erica Hinrichs
                 
        """
    
        #arguments = {}
        #arguments['variables'] = ['Np[aq] (Moles/kg H20)']
        
        #test_case() does not actually display this custom error message yet.
        #arguments['err_msg'] = '\nIncorrect concentration at time %s'
        
        self.test_case('evaporation')

    def salt_perm_poro(self):
        """ 
        **Test the Salt Permeability and Porosity Macro**

        The porosity-permeability function for compacted salt from *Cinar et 
        at. (2006)* is tested using a six node problem with porosities from 
        0.01 to 0.2. The excel spreadsheet in information/salt-perm-poro.xlsx 
        contains calculations of the perm-poro function.

        *Cinar, Y, G Pusch and V Reitenbach (2006) Petrophysical and 
        capillary properties of compacted salt. Transport in Porous Media. 
        64, p. 199-228, doi: 10.1007/s11242-005-2848-1* 
        
        .. Authors: Dylan Harp, Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
    
        arguments = {}
        arguments['variables'] = ['n', 'perm_x']
        
        #test_case() does not actually display this custom error message yet.
        arguments['err_msg'] = \
            '\nIncorrect permeability at node %s. Expected %s, Simulated '    
        
        self.test_case('salt_perm_poro', arguments)
  
    def avdonin(self):
        """
        **Test the Radial Heat and Mass Transfer Problem**
        
        Compares the generated contour and history files to old contour and 
        history files that are known to be correct. For contour files, only the
        temperature values at time 2 are tested. For history files, all 
        temperature values are tested. 
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange    
        """
        
        arguments = {}
        arguments['times'] = [2.0]
        arguments['variables'] = ['T']
        
        self.test_case('avdonin', arguments)
        
    def boun(self): 
        """
        **Test the Boundry Functionality**
        
        Compares the generated contour files to old contour files that are known
        to be correct. Only the pressure and hydraulic head values at time 2 are 
        tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['times'] = [2.0]
        arguments['variables'] = ['P', 'Hydraulic Head (m)']

        self.test_case('boun', arguments) 
        
    def cden(self):
        """
        **Test the Concentration Dependent Brine Density Functionality**
        
        Compares generated history files to old history files that are known to 
        be correct. Only the density values are tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['variables'] = ['density']
         
        self.test_case('cden') 

    def cflxz(self):
        """
        **cflxz_test**
        
        Concentration Zone Flux Test                                                    
        Comparison of FEHm with 1-D Difusion Model  
        
        .. Authors: Mark Lange
        .. Updated: 2024 by Erica hinrichs
        """
        
        #arguments = {}
        #arguments['variables'] = ['density']
         
        self.test_case('cflxz_test')

    def heat3d(self):
        """
        **heat3d**
        
        3-D Heat Conduction Problem                                                     
        Comparison of Model and Analytical Solution for Temperature vs Time  
        
        .. Authors: Mark Lange
        .. Updated: 2024 by Erica hinrichs
        """
        
        #arguments = {}
        #arguments['variables'] = ['density']
         
        self.test_case('heat3d') 
        
    def doe(self):
        """
        **Test the DOE Code Comparison Project, Problem 5, Case A**
        
        Compares the generated contour and history files to old contour and
        history files that are known to be correct. For contour files, only the
        pressure, temperature, and saturation values at time 3 are tested. For 
        history files, all pressure, temperature, and saturation values are 
        tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['times'] = [3.0]
        arguments['variables'] = ['P', 'T', 'saturation']
        arguments['maxerr'] = 1.0
        
        self.test_case('doe', arguments)
        
    def head(self):
        """
        **Test Head Pressure Problem**
        
        Compares the generated contour files to old contour files that are known
        to be correct. Only the pressure values at day 2 are tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange 
        """
        
        arguments = {}
        arguments['times'] = [2.0]
        arguments['variables'] = ['P']
        
        self.test_case('head', arguments)

    def heat2d(self):
        """
        **Test heat2d**
         
        """
        
        #arguments = {}
        #arguments['times'] = [2.0]
        #arguments['variables'] = ['P']
        
        self.test_case('heat2d')
                    
    def heat2d_quad(self):
        """
        **Test the Heat 2D Quad Problem**

        Input files copied from heat2d_quad which writes output to heat_flux
        Compares heat2d_quad 00003 121 Node same as VV compare script

        .. Authors:  Terry Miller
        .. Updatd: July 2024
        """

        self.test_case('heat2d_quad')

    def henrys_law(self):
        """
        **Test the Henrys law Problem**

        From henry1.comparein
        It looks like trc values for time and concentration are compared to input/*.analyt
        after converting values by some number.

        We should be able to compare trc to trc as the numbers are close to equal
        Precision is outside of normal math ie 0.999477720105420 - 0.99947772010542046 = 0

        .. Authors:  Terry Miller tam
        .. Updatd: July 2024
        """

        self.test_case('henrys_law')

    def ramey(self):
        """
        **Test Temperature in a Wellbore Problem**
        
        Compares the generated contour and history files to old contour and
        history file that are known to be correct. For the contour files, only 
        the temperature values at time 2 are tested. For the history files, all 
        temperature values are tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['times'] = [2.0]
        arguments['variables'] = ['T']
        
        self.test_case('ramey', arguments)

    def sptr_btc(self):
        """
        **Test sptr_btc**
        
        .. Authors: Mark Lange
        .. Updated: July 2024 by Erica Hinrichs
        """
        
        self.test_case('sptr_btc')

    def perm_test(self):
        """
        **Test perm_test**

        Comparison of Ratio of Permeabilty Change to Ratio of Change in Flow Rate
        
        .. Authors: Mark Lange
        .. Updated: July 2024 by Erica Hinrichs
        """
        
        self.test_case('perm_test')

    def transport3d(self):
        """
        **Test transport3d**

        Three-Dimensional Radionuclide Transport Problem - trac_rlp                     
        Comparison of FEHM and TRACRN for Concentration vs Time
        
        .. Authors: Mark Lange
        .. Updated: July 2024 by Erica Hinrichs
        """
        
        self.test_case('transport3d')

    def uz_test(self):
        """
        **Test uz_test**

        
        .. Authors: Mark Lange
        .. Updated: July 2024 by Erica Hinrichs
        """
        
        self.test_case('uz_test')
   
    def theis(self):
        """
        **Test Pressure Transient Analysis Problem**
        
        Compares the generated contour files to old contour files known to be 
        correct. Only the pressure values at time 2 are tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['times'] = [2.0]
        arguments['variables'] = ['P']
        
        self.test_case('theis', arguments)
        
    def dryout(self):
        """
        **Test Dry-Out of a Partially Saturated Medium**
        
        Compares the generated contour files to old contour files known to be 
        correct. The saturation is tested for all times.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange 
        """
        
        arguments = {}
        arguments['variables'] = ['saturation']
        arguments['maxerr'] = 0.01
        
        self.test_case('dryout', arguments)
        
    def multi_solute(self):
        """
        **Test Multi-Solute Transport with Chemical Reaction**
        
        Compares the generated tracer files to old tracer files known to be 
        correct. All concentraction values are tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        self.test_case('multi_solute')
        
    def sorption(self):
        """
        **Test One Dimensional Reactive Solute Transport**
        
        Compares the generated tracer files to old tracer files known to be
        correct. All concentraction values are tested.
        
        .. Authors: Mark Lange
        .. Updated: June  2014 by Mark Lange
        """
        
        self.test_case('sorption')
        
    def baro_vel(self):
        """
        **Test Pore-Scale Velocity in a Homogeneous Media**
        
        Compares the generated contour files with the old files known to be
        correct. Tests times 3-6 for a root mean square difference of less than 
        0.01.
        
        .. Authors: Mark Lange
        .. Updated: July 2014 by Mark Lange
        """
        
        args = {}
        args['test_measure'] = 'rms_difference'
        args['maxerr'] = 0.01
        args['times'] = [3.0, 4.0, 5.0, 6.0]
        
        self.test_case('baro_vel', args)
        
    def cellbased(self):
        """
        **Test the Cell-Based Particle Tracking Model**
        
        Compares the generated output files with the old files known be correct.
        All values are tested for a root mean square difference of less than
        0.05.
        
        *TODO: This test does not regenerated chain.ini! This is something that 
         needs to be fixed in the future.*
        
        .. Authors: Mark Lange
        .. Updated: July 2014 by Mark Lange
        """
        
        args = {}
        args['test_measure'] = 'rms_difference'
        args['maxerr'] = 0.05
        
        self.test_case('cellbased',args)
        
    def heat_pipe(self):      
        """
        **Test the Heat Pipe Problem**
        
        Compares the generated output files with the old files known to be
        correct. All values are tested.
        
        RLPM subcases were excluded because they caused FEHM to produce too many
        negative volumes.
        
        .. Authors: Mark Lange
        .. Updatd: July 2014 by Mark Lange
        """
        
        self.test_case('heat_pipe')

    def toronyi(self):
        """
        **Test the Toronyi Two-Phase Problem**
        
        Compares the generated contour files with the old contour files known to
        be correct. Tests all values at time 2.0.
        
        .. Authors: Mark Lange
        .. Updated: July 2014 by Mark Lange
        """
    
        args = {}
        args['times'] = [2.0]
        
        self.test_case('toronyi', args)
        
        
    def colloid_filtration(self):
        '''
        Test Colloid Filtration
        
        Compares the generated ptrk files with the old ptrk files known to be 
        correct.
        
        .. Authors: Mark Lange
        .. Updated: July 2024 by Erica Hinrichs
        '''
        
        self.test_case('colloid_filtration')
        
    def mptr(self):
        '''
        Test Multi-Species Particle Tracking
        
        Compares the generated ptrk files with the old ptrk files knonw to be 
        correct.
        
        .. Authors: Mark Lange
        .. Updated: July 2014 by Mark Lange
        '''
        
        self.test_case('mptr')

    def bodyforce(self):
        """
        **Test greater flexibility for specifying body forces**
        
        Compares results generated by three available syntax inputs to an old output files known to be 
        correct.
        
        .. Authors: David Dempsey 
        .. Updated: May 2014 by Shaoping Chu 
        """
        
        arguments = {}
        arguments['variables'] = ['strs_xx','strs_yy','strs_zz']
        arguments['maxerr'] = 0.5

        self.test_case('bodyforce',arguments)
        
    def richards(self):
        """
        **Richards equation test**
        
        Comparison of Richards equation with 2-phase solution 
        Node by node comparison.
        
        .. Authors: Shaoping Chu, converted from Zora Dash's Perl-based test cases 
        .. Updated: May 2014 by Shaoping Chu 
        """
        
        arguments = {}
        arguments['variables'] = ['P','S']
        arguments['test_measure'] = 'perc_difference'
        arguments['maxerr'] = 0.5

        self.test_case('richards',arguments)
        
    def ppor_read(self):
        """
        **Test test ppor and porosity read/write**

        Compares results generated by model output to an old output files known to be
        correct.

        .. Authors: Shaoping, converted from Zhiming Lu's test case 
        .. Updated: Oct 2016 by Shaoping Chu
        """

        arguments = {}

        self.test_case('ppor_read',arguments)

    def wvtest(self):
        """
        **Test wvtest**

        Unsaturated Diffusion test                                                      
        Comparison of FEHM Saturation vs Time with original run

        .. Updated: July 2024 by Erica Hinrichs

        """

        #arguments = {}

        self.test_case('wvtest')

    def fracture_aperture(self):
        """
        **Test fracture_aperture**

        Comparison of 3D Wellbore Thermal Stress test results

        .. Updated: July 2024 by Erica Hinrichs

        """

        #arguments = {}

        self.test_case('fracture_aperture')

    def potential_energy(self):
        """
        **Test potential_energy**

        Comparison of Ratio of Permeabilty Change to Ratio of Change in Flow Rate

        .. Updated: July 2024 by Erica Hinrichs

        """

        #arguments = {}

        self.test_case('potential_energy')

    def vapor_extraction(self):
        """
        **Test vapor_extraction**

        Comparison of Model and Analytical Solution for Vapor Pressure vs Position

        .. Updated: July 2024 by Erica Hinrichs

        """

        #arguments = {}

        self.test_case('vapor_extraction')

    def heatflux_1DConvection(self):
        '''
        Test heat flux output for a 1D convection only case

        Compares the generated *_hf_node.csv files with a pre-established output

        ..Authors: Sharad Kelkar, Zora Dash
        ..Updated: 6 August 2014 by Sharad Kelkar
        '''

        args = {}
        args['times'] = [10.0]
        args['test_measure'] = 'perc_difference'
        args['maxerr'] = 2.0

        self.test_case('heatflux_1DConvection')
                    
    # Test Developer Functionality ############################################
        
    def test_case(self, name, parameters={}):
        """ 
        Performs a test on a FEHM simulation and raises an AssertError if it 
        fails the test. 
         
        :param name: The name of the test-case folder.
        :type name: str
        
        :param parameters: Attribute values that override default values.
                               Key Value Choices
                                   'variables': list[str, str, ...]
                                   'times': list[float, float, ...]
                                   'nodes': list[int, int, ...]
                                   'components': list[str, str, ...]
                                   'maxerr': float
                                   'test_measure': str - 'max_difference' or
                                                         'rms_difference' or
                                                         'perc_difference'
                                   
        :type parameters: dict
            
        The folder 'name' in fehmpytests must exist with correct structure.
        If parameters are not passed into this method, all simulated attributes 
        will be checked for a relative difference of less than 1.e-4. 
        
        Authors: Mark Lange
        Updated: June 2014 by Mark Lange                                
        """ 
         
        os.chdir(name)
        
        #Search for fehmn control files and extract subcases.
        filenames = glob.glob(os.path.join('input','control','*.files'))
        #print('filenames: ', filenames)
        subcases  = []
        for filename in filenames:
            #print('\nfilename: ', filename)

            path = os.path.join('input', 'control', ' ')
            if os.name == 'nt':  # For Windows
                subcase = re.sub(r'input[\\/]', '', filename)
                subcase = re.sub(r'control[\\/]', '', subcase)
                subcase = re.sub(r'\.files$', '', subcase)
            else:  # For Unix-based systems
                subcase = re.sub(os.path.join('input','control',''), '', filename)
                subcase = re.sub('.files', '', subcase)
            #print('subcase: ', subcase)
            
            #File named 'fehmn.files' to be used for tests with single case.
            if subcase != 'fehmn':
                subcases.append(subcase)
                #print('append subases', subcases)
            else:
                subcases = ['']
                #print('broken subcases')
                break
                
        try:
            #Test the new files generated with each subcase.
            for subcase in subcases:
                parameters['subcase'] = subcase
                # CD into run directory
                output_dir = subcase+'_output'
                #print('Output directory: ', output_dir)
                if os.path.exists( output_dir ): shutil.rmtree(output_dir)
                os.mkdir( output_dir )
                os.chdir( output_dir )
                filetypes = ['*.avs','*.csv','*.his','*.out','*.trc','*.ptrk','*.dat','*.sptr3', '*.cflx']
                test_flag = False
                
                for filetype in filetypes:
                    parameters['filetype'] = filetype
                    compare_pattern = (os.path.join('..', 'compare', '*' )+ subcase + filetype)
                    found_files = glob.glob(compare_pattern)

                    #print('Checking for files with pattern: ', compare_pattern)
                    #print('Found files: ', found_files)
                    #print('Number of found files: ', len(found_files))

                    if len(found_files) > 0:
                        #print(f'\n\nfound {len(found_files)} files. Continuing to test template')
                        test_method = self._test_template(filetype, subcase, parameters)
                        test_method()
                        #print(f'Test method executed for filetype: {filetype} on files: {found_files}')
                        test_flag = True
                    else:
                        #print('Test method NOT executed for filetype: ', filetype)
                        pass

                os.chdir('..')
                if not test_flag:
                    if self.log:
                        line = f'\nFailed at subcase: {subcase} filetype: {filetype}'
                        self.fail_log.write(line)
                    self.fail("Missing any valid comparison files, no test performed")

        finally:
            # Allows other tests to be performed after exception.
            os.chdir(self.maindir)
            
    def _test_template(self, filetype, subcase, parameters={}):
        """
        **Test Template**
        
        Calling this function with the filename and subcase will return the 
        correct test method. 
        
        :param parameters: Stores optional prespecified variable, time, node, 
                           component, and format values to override defaults.
        :type filesfile:   dict 
        
        The intent behind creating this method was to make it easier to add new 
        tests for other outputs in the future. If you need to modify this
        function and need help understanding it, look into 'closure' and
        functions as 'first class objects'.
        """
            
        #Get pre-specified parameters from call.
        #print('\n---- at test template ----')
        keys = ['variables', 'times', 'nodes', 'components', 'info']
        values = dict.fromkeys(keys, [])
        values['maxerr'] = 1.e-4
        values['test_measure'] = 'max_difference'
        for key in values:
            if key in parameters:
                values[key] = parameters[key]                      
        mxerr = values['maxerr']
        components = values['components']
        test_measure = values['test_measure']
        #print('\nsubcase: ', subcase, ' filetype: ', filetype, ' parameters: ', parameters)

        self._run_fehm(subcase)
       
        def contour_case():      
            #Find the difference between the old and new
            #print('Contour Case')
            f_old = fpost.fcontour(os.path.join('..','compare','*')+subcase+'.'+filetype)
            f_new = fpost.fcontour('*'+subcase+'.'+filetype)
            f_dif = fpost.fdiff(f_new, f_old)
                
            msg = 'Incorrect %s at time %s.'

            #If no pre-specified times, grab them from f_dif.
            if len(values['times']) == 0:
                times = f_dif.times
            else:
                times = values['times']
            
            #If no pre-specified variables, grab them from f_dif.         
            if len(values['variables']) == 0:
                variables = f_dif.variables
            else:
                variables = values['variables']
            
            #Check the variables at each time for any significant differences.
            test_flag = False
            for t in times: 
                #print(subcase)
                #Its possible some times do not have all variables in f_dif.
                for v in np.intersect1d(variables, list(f_dif[t].keys())):
                    #Measure the difference into a single quantity.
                    f_dif[t][v] = list(map(abs, f_dif[t][v]))
                    f_old[t][v] = list(map(abs, f_old[t][v]))
                    difference = { 
                        'max_difference': max(f_dif[t][v]),
                        'rms_difference': np.sqrt(np.mean(f_dif[t][v])), 
                        'perc_difference':
                            sum([x/float(y) 
                                 for x,y in zip(f_dif[t][v], f_old[t][v])
                                 if y != 0]) /
                            float(len(f_dif[t][v]))
                    }[test_measure]
                    try:
                        #print('true? ', difference, '<', mxerr, '=', difference<mxerr, 'for ', (os.path.join('..','compare','*')+subcase+'.'+filetype), ' and ', ('*'+subcase+'.'+filetype))
                        self.assertTrue(difference<mxerr, msg%(v, t))
                        test_flag = True
                    except AssertionError as e:
                        #Write to fail log if switch is on.
                        if self.log:
                            kvpairs = {'variable':str(v), 'time':str(t)}
                            line = '\nFailed at subcase:'+subcase
                            line = line+' filetype:'+filetype
                            for key in kvpairs:        
                                line = line+' '+key+':'+kvpairs[key]
                            self.fail_log.write(line)   
                        raise e
            if not test_flag:
                self.fail("Missing common nodes in compare and output contour files, no test performed")
        def history_case():
            #Find the difference between the old and new
            f_old = fpost.fhistory(os.path.join('..','compare','*')+subcase+filetype)
            f_new = fpost.fhistory('*'+subcase+filetype)
            f_dif = fpost.fdiff(f_new, f_old)

            #If no pre-specified variables, grab them from f_dif.         
            if len(values['variables']) == 0:
                variables = f_dif.variables
            else:
                variables = values['variables']
            
            #If no pre-specifed nodes, grab them from f_dif.
            if len(values['nodes']) == 0:
                nodes = f_dif.nodes
            else:
                nodes = values['nodes']   
            msg = 'Incorrect %s at node %s.'  
            
            #Check the nodes at each variable for any significant differences.   
            test_flag = False
            for v in variables:
                #Its possible some variables do not have all nodes in f_dif.
                for n in np.intersect1d(nodes, list(f_dif[v].keys())):  
                    f_dif[v][n] = list(map(abs, f_dif[v][n]))
                    difference = { 
                        'max_difference': max(f_dif[v][n]),
                        'rms_difference': np.sqrt(np.mean(f_dif[v][n])), 
                        'perc_difference':
                            sum([x/float(y) 
                                 for x,y in zip(f_dif[v][n], f_old[v][n])
                                 if y != 0]) /
                            float(len(f_dif[v][n]))
                    }[test_measure]
                    try:   
                        self.assertTrue(difference<mxerr, msg%(v, n))
                        test_flag = True
                    except AssertionError as e:
                        #Write to fail log if switch is on.
                        if self.log:
                            kvpairs = {'variable':str(v), 'node':str(n)}
                            line = '\nFailed at subcase:'+subcase
                            line = line+' filetype:'+filetype
                            for key in kvpairs:        
                                line = line+' '+key+':'+kvpairs[key]
                            print(line)
                            self.fail_log.write(line)   
                        raise e
            if not test_flag:
                self.fail("Missing common nodes in compare and output history files, no test performed")

        #### ADDING FOR COMPARISONS ####
        def comparison_case():
            #Find the difference between the old and new
            #print('---- at comparison in fhmpytests ----')
            f_old = fpost.fcomparison(os.path.join('..','compare','*')+subcase+filetype)
            f_new = fpost.fcomparison('*'+subcase+filetype)
            f_dif = fpost.fdiff(f_new, f_old)

            total=0

            for ele in range(0, len(f_dif._info)):
                total = total + f_dif._info[ele]

            test_flag = False
            try:   
                self.assertTrue(total == 0)
                test_flag = True
            except AssertionError as e:
                print("Significant differences found, no test performed.")
                print('fpost._info: ', f_dif._info)
                #Write to fail log if switch is on.
                if self.log:
                    line = '\nThere are significant differences between files, no test performed'
                    self.fail_log.write(line)   
                raise e
            
            if not test_flag:
                self.fail("There are significant differences between files, no test performed")

        ####################################################
                    
        def tracer_case():
            #Find the difference between the old and new
            f_old = fpost.ftracer(os.path.join('..','compare','*')+subcase+filetype) 
            f_new = fpost.ftracer('*'+subcase+filetype)
            f_dif = fpost.fdiff(f_new, f_old)
            
            #If no pre-specified variables, grab them from f_dif.         
            if len(values['variables']) == 0:
                variables = f_dif.variables
            else:
                variables = values['variables']
            
            #If no pre-specifed nodes, grab them from f_dif.
            if len(values['nodes']) == 0:
                nodes = f_dif.nodes
            else:
                nodes = values['nodes']    

            msg = 'Incorrect %s at node %s.'  
            
            #Check the nodes at each variable for any significant differences.
            test_flag = False
            for v in variables:
                #Its possible some variables do not have all nodes in f_dif.
                for n in np.intersect1d(nodes, list(f_dif[v].keys())):
                    #Measure the difference into a single quantity.
                    f_dif[v][n] = list(map(abs, f_dif[v][n]))
                    difference = { 
                        'max_difference': max(f_dif[v][n]),
                        'rms_difference': np.sqrt(np.mean(f_dif[v][n])), 
                        'perc_difference':
                            sum([x/float(y) 
                                 for x,y in zip(f_dif[v][n], f_old[v][n])
                                 if y != 0]) /
                            float(len(f_dif[v][n]))
                    }[test_measure]     
                    try:
                        self.assertTrue(difference<mxerr, msg%(v, n))
                        test_flag = True
                    except AssertionError as e:
                        #Write to fail log if switch is on.
                        if self.log:
                            kvpairs = {'variable':str(v), 'node':str(n)}
                            line = '\nFailed at subcase:'+subcase
                            line = line+' filetype:'+filetype
                            for key in kvpairs:        
                                line = line+' '+key+':'+kvpairs[key]
                            self.fail_log.write(line)   
                        raise e   
            if not test_flag:
                self.fail("Missing common nodes in compare and output tracer files, no test performed")
            
        def output_case():
            #Find difference between old and new file assume 1 file per subcase.
            old_filename = glob.glob(os.path.join('..','compare','*')+subcase+filetype)[0]
            new_filename = glob.glob('*'+subcase+filetype)[0]
            f_old = fpost.foutput(old_filename)
            f_new = fpost.foutput(new_filename)
            f_dif = fpost.fdiff(f_new, f_old)
            
            #If no pre-specified variables, grab them from f_dif.         
            if len(values['variables']) == 0:
                variables = f_dif.variables
            else:
                variables = values['variables']

            #If no pre-specifed nodes, grab them from f_dif.
            if len(values['nodes']) == 0:
                nodes = f_dif.nodes
            else:
                nodes = values['nodes']
                
            #If no pre-specifed components, grab them from f_dif.
            if len(values['components']) == 0:
                components = f_dif.components
            else:
                components = values['components']
                
            msg = 'Incorrect %s at %s node %s.'  
            
            #Check the node at each component for significant differences.   
            test_flag = False
            for c in components:
                for n in np.intersect1d(nodes, list(f_dif.node[c].keys())):
                    for v in np.intersect1d(variables, list(f_dif.node[c][n].keys())):
                        #Measure the difference into a single quantity.
                        fdiff_array = list(map(abs, f_dif.node[c][n][v]))
                        difference = { 
                            'max_difference': max(fdiff_array),
                            'rms_difference': np.sqrt(np.mean(fdiff_array)),
                            'perc_difference':
                                sum([x/float(y) 
                                     for x,y in zip(fdiff_array, fdiff_array)
                                     if y != 0]) /
                                float(len(fdiff_array)) 
                        }[test_measure]
                        try:
                            self.assertTrue(difference < mxerr, msg%(v,c,n))
                            test_flag = True
                        except AssertionError as e:
                            #Write to fail log if switch is on.
                            if self.log:
                                kvpairs = { 'component': str(c), 
                                            'node': str(n),
                                            'variable': str(v), }
                                line = '\nFailed at subcase:'+subcase
                                line = line+' filetype:'+filetype
                                for key in kvpairs:        
                                    line = line+' '+key+':'+kvpairs[key]
                                self.fail_log.write(line)   
                            raise e
            if not test_flag:
                self.fail("Missing common nodes and/or variables in compare and output out files, no test performed")
        
        def ptrack_case():
            #Find the difference between the old and new
            f_old = fpost.fptrk(os.path.join('..','compare','*')+subcase+filetype)
            f_new = fpost.fptrk('*'+subcase+filetype)  
            f_dif = fpost.fdiff(f_new, f_old)
            
            msg = 'Incorrect %s.'
            
            #If no pre-specified variables, grab them from f_dif.         
            if len(values['variables']) == 0:
                variables = f_dif.variables
            else:
                variables = values['variables']
                
            test_flag = False
            for v in variables:
                #Measure the difference into a single quantity.
                fdiff_array = list(map(abs, f_dif[v]))
                difference = { 
                    'max_difference': max(fdiff_array),
                    'rms_difference': np.sqrt(np.mean(fdiff_array)),
                    'perc_difference':
                            sum([x/float(y) 
                                 for x,y in zip(fdiff_array, fdiff_array)
                                 if y != 0]) /
                            float(len(fdiff_array)) 
                }[test_measure]
                #Perform test, if fail log switch is on, write a fail log.
                try:
                    self.assertTrue(difference < mxerr, msg%v)
                    test_flag = True
                except AssertionError as e:
                    #Write to fail log if switch is on.
                    if self.log:
                        kvpairs = {'variable': str(v)}
                        line = '\nFailed at subcase:'+subcase
                        line = line+' filetype:'+filetype
                        for key in kvpairs:        
                            line = line+' '+key+':'+kvpairs[key]
                        self.fail_log.write(line)   
                    raise e
            if not test_flag:
                self.fail("Missing common nodes in compare and output ptrk files, no test performed")

        #Returns the test method for filetype.
        return { '*.avs':  contour_case,
                 '*.csv':  contour_case,
                 '*.dat':  comparison_case,
                 '*.sptr3':  comparison_case,
                 '*.his':  history_case,
                 '*.trc':  comparison_case,
                 '*.out':  comparison_case, 
                 '*.cflx': comparison_case,
                 '*.ptrk': ptrack_case }[filetype]
                                    
    def _run_fehm(self, subcase):
        """ 
        **Utility function to run fehm**
        
        Asserts that fehm terminates successfully.

        :param filesfile: name of fehm files file
        :type filesfile: str 
        """
        
        #Find the control file for the test-case or for the subcase. 
        if subcase == '':
            #print('subcase empty assume fehm at: \n',  os.path.join('..','input','control','fehmn.files'))
            filesfile = os.path.join('..','input','control','fehmn.files')
        else:
            filesfile = os.path.join('..','input','control',subcase+'.files')
            
        evalstr = exe+' '+filesfile
        #print('evalstr: ', evalstr)
        
        with open(os.devnull, "w") as f:
            call(evalstr, shell=True, stdout=f)
        
        outfile = None
        errfile = 'fehmn.err'

        with open( filesfile, 'r' ) as f:
            #print('\nchecking for outp and error...')
            lines = f.readlines()
            # Check for new filesfile format
            for line in lines:
                if 'outp' in line:
                    outfile = line.split(':')[1].strip()
                elif 'error' in line:
                    errfile = line.split(':')[1].strip()
                           
            # Assume old format
            if outfile is None and ':' not in lines[0]:
                outfile=lines[3].strip()
                #print('\nNo outfile found...', outfile)
 
        complete = False
        if outfile:
            #print('\noutfile found...')
            with open(outfile, 'r' ) as f:
                for line in reversed(f.readlines()):
                    if 'End Date' in line:
                        complete = True
                        break
                        
        if os.path.exists(errfile): 
            errstr = open( errfile, 'r' ).read()
        else: 
            errstr = ''
        # Change to maindir in case assertTrue fails 
        #print('CWD: ',  os.getcwd() )   
        curdir = os.getcwd() 
        os.chdir(self.maindir)
        
        msg = 'Unsuccessful fehm simulation\nContents of '
        #print('!!!!!!',self.assertTrue(complete, msg+errfile+':\n\n'+errstr) )
        self.assertTrue(complete, msg+errfile+':\n\n'+errstr)
        os.chdir(curdir)
                     
def cleanup():
    """ 
    Utility function to remove files after test

    """
    dirs = glob.glob(os.path.join('*','*_output'))
    for d in dirs: shutil.rmtree( d )
                  
def suite(mode, test_case, log):
    suite = unittest.TestSuite()
    
    #Default mode is admin for now. Should it be different?
    if mode == 'admin' or mode == 'default':
        suite.addTest(fehmTest('avdonin', log))
        suite.addTest(fehmTest('baro_vel', log))
        suite.addTest(fehmTest('bodyforce', log))
        #suite.addTest(fehmTest('boun', log))
        suite.addTest(fehmTest('cden', log))
        #suite.addTest(fehmTest('cellbased', log))
        suite.addTest(fehmTest('cflxz', log))
        suite.addTest(fehmTest('colloid_filtration', log))
        #suite.addTest(fehmTest('dissolution', log))
        #suite.addTest(fehmTest('doe', log))
        suite.addTest(fehmTest('dryout', log))
        suite.addTest(fehmTest('evaporation', log))
        suite.addTest(fehmTest('fracture_aperture', log))
        suite.addTest(fehmTest('head', log))
        suite.addTest(fehmTest('heat2d', log))
        suite.addTest(fehmTest('heat2d_quad', log))
        suite.addTest(fehmTest('heat3d', log))
        suite.addTest(fehmTest('heat_pipe', log))
        suite.addTest(fehmTest('henrys_law', log))
        suite.addTest(fehmTest('mptr', log))
        suite.addTest(fehmTest('multi_solute', log))
        suite.addTest(fehmTest('perm_test', log))
        suite.addTest(fehmTest('potential_energy', log))
        suite.addTest(fehmTest('ramey', log))
        suite.addTest(fehmTest('richards', log))
        suite.addTest(fehmTest('salt_perm_poro', log))
        #suite.addTest(fehmTest('saltvcon', log))
        suite.addTest(fehmTest('sorption', log))
        suite.addTest(fehmTest('sptr_btc', log))
        suite.addTest(fehmTest('theis', log))
        suite.addTest(fehmTest('toronyi', log))
        #suite.addTest(fehmTest('transport3d', log))
        #suite.addTest(fehmTest('uz_test', log))
        suite.addTest(fehmTest('vapor_extraction', log))
        suite.addTest(fehmTest('wvtest', log))
        suite.addTest(fehmTest('rad_decay', log))   
    
    elif mode == 'developer':
        #This mode will be a reduced set that runs faster.
        pass
             
    elif mode == 'solo':
        suite.addTest(fehmTest(test_case, log))
            
    elif mode == 'admin':
        pass
        
    return suite
    
if __name__ == '__main__':
    
    #Unless the user specifies a single test-case, this isn't important.
    test_case = ''

    #Set up command-line interface.
    parser = argparse.ArgumentParser(description='FEHM Test-Suite')
    
    #Comand-line Options
    group = parser.add_mutually_exclusive_group()   
    a = 'store_true'
    h = 'Run the entire test-suite.'
    group.add_argument('-a', '--admin', help=h, action=a)
    h = 'Run a portion of the test-suite.'
    group.add_argument('-d', '--dev', help=h, action=a)
    h = 'Run a single test-case.'
    group.add_argument('-s', '--solo', help=h, action=a)
    h = "Create a fail statistics file 'fail_log.txt'"
    parser.add_argument('-l', '--log', help=h, action=a)
    h = "Clean up fehm output files"
    parser.add_argument( '--clean', help=h, action=a)
    #Positional Arguments
    h = 'Path to the FEHM executable.'
    parser.add_argument('exe', help=h)
    h = 'Single test-case to run.'
    parser.add_argument('testcase', nargs='?', help=h, default=None)
     
    args = vars(parser.parse_args())
    
    exe = os.path.abspath(args['exe'])
    
    if args['clean']:
        cleanup()
        sys.exit(0)

    #Determine the mode.    
    if args['solo']:
        #Make sure that the test-case was specified, otherwise show help.
        if args['testcase'] != None:
            mode = 'solo'
            test_case = args['testcase']
        else:
            parser.print_help()   
    else:
        #Make sure user did not attempt to specify a test-case, show help if so.
        if args['testcase'] == None:
            if args['admin']:
                mode = 'admin'
            elif args['dev']:
                mode = 'developer'
            else:
                mode = 'default'
        else:
            #If user didn't specify admin or dev, assume solo mode.
            if not args['admin'] and not args['dev']:
                mode = 'solo'
                test_case = args['testcase']
            else:       
                parser.print_help()
                
    #If the user wants a log, give them a log.
    log = False            
    if args['log']:
        log = True
    
    #Run the test suite.    
    runner = unittest.TextTestRunner(verbosity=2)
    test_suite = suite(mode, test_case, log)
    runner.run(test_suite)
    




