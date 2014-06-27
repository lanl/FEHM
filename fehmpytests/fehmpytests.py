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
from subprocess import call, PIPE
try:
    sys.path.insert(0,'../pyfehm')
    from fdata import *
except ImportError as err:
    print 'ERROR: Unable to import pyfehm fpost module'
    print err
    os._exit(0)
from glob import glob

#Suppresses tracebacks
__unittest = True 


class Tests(unittest.TestCase):
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
    
    # TESTS ######################################################### 
        
    def test_saltvcon(self):
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
        to 1.e-6.

        The excel spreadsheet /information/saltvcon.xlsx contains the
        associated calculations. 
        
        .. Authors: Dylan Harp, Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['components'] = ['water']
        arguments['variables']  = ['Kx']
        arguments['format'] = 'relative' 
        
        #test_case() does not actually display this custom error message yet.
        arguments['err_msg'] = \
            '\nIncorrect intact salt thermal conductivity for node %s \
             \nRelative error: %s, Threshold: %s \
             \nExpected = %s Simulated = %s'
        
        self.test_case('saltvcon', arguments)
        
        
    def test_dissolution(self):
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
        
    def test_salt_perm_poro(self):
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
  
    def test_avdonin(self):
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
        
    def test_boun(self): 
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

        self.test_case('boun_test', arguments) 
        
    def test_cden(self):
        """
        **Test the Concentration Dependent Brine Density Functionality**
        
        Compares generated history files to old history files that are known to 
        be correct. Only the density values are tested.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange
        """
        
        arguments = {}
        arguments['variables'] = ['density']
         
        self.test_case('cden_test') 
        
    def test_doe(self):
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
        
        self.test_case('doe', arguments)
        
    def test_head(self):
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
                    
    def test_ramey(self):
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
   
    def test_theis(self):
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
        
    def test_dryout(self):
        """
        Test Dry-Out of a Partially Saturated Medium
        
        Compares the generated contour files to old contour files known to be 
        correct. The saturation is tested for all times.
        
        .. Authors: Mark Lange
        .. Updated: June 2014 by Mark Lange 
        """
        
        arguments = {}
        arguments['variables'] = ['saturation']
        
        self.test_case('dryout')
        
    def test_multi_solute(self):
        self.test_case('multi_solute')
        
    # Test Developer Functionality ############################################
        
    def test_case(self, name, parameters={}):
        """ 
        Performs a test on a FEHM simulation and raises an AssertError if it 
        fails the test. 
         
        :param name: The name of the test-case folder.
        :type name: str
        
        :param parameters: Attribute values that override default values.
        :type parameters: dict
            
        The folder 'name' in fehmpytests must exist with correct structure.
        If parameters are not passed into this method, all simulated attributes 
        will be checked for a relative difference of less than 1.e-4. 
        
        Authors: Mark Lange
        Updated: June 2014 by Mark Lange                                
        """ 
         
        os.chdir(name)
        
        #Search for fehmn control files and extract subcases.
        filenames = glob('input/control/*.files')
        subcases  = []
        for filename in filenames:
            subcase = re.sub('input/control/', '', filename)
            subcase = re.sub('.files', '', subcase)
            
            #File named 'fehmn.files' to be used for tests with single case.
            if subcase != 'fehmn':
                subcases.append(subcase)
            else:
                subcases = ['']
                break
                
        try:
            #Test the new files generated with each subcase.
            for subcase in subcases:
                filetypes = ['*.avs', '*.csv', '*.his', '*.out', '*.trc']
                for filetype in filetypes:
                    #Check to make sure there are files of this type.
                    if len(glob('compare/'+'*'+subcase+filetype)) > 0:
                        #Check for any significant differences.
                        self._checkDifferences(filetype, subcase, parameters)
                 
        finally:
            #Allows other tests to be performed after exception.
            #self._cleanup(['*.*'])
            os.chdir(self.maindir)
            
                   
    # UTILITIES ######################################################
         
    def _checkDifferences(self, filetype, subcase, parameters={}):
        """ 
        **Check Difference**
        
        Checks the difference between old and new output. Fails if the
        difference is greater than maxerr, defaulted to 1.e-4. 
        
        .. Authors: Mark Lange
        .. Updated: June 2014
        """

        #Get pre-specified parameters from call.
        keys = ['variables', 'times', 'nodes', 'components', 'format']
        values = dict.fromkeys(keys, [])
        values['maxerr'] = 1.e-4
        format = 'relative'
        for key in values:
            if key in parameters:
                values[key] = parameters[key]
                
        variables = values['variables']
        times = values['times']
        nodes = values['nodes']
        mxerr = values['maxerr']
        components = values['components']
        format = values['format']
        
        #Read in the old comparison files. 
        f_old = self._fgeneral('compare/*'+subcase+filetype) 

        #Read in new files.
        self._run_fehm(subcase)
        f_new = self._fgeneral('*'+subcase+filetype)
        
        #Find the difference between the two files.
        f_dif = fdiff(f_new, f_old)
        
        #If no pre-specified variables, grab them from f_dif.
        if len(variables) == 0:
            variables = f_dif.variables

        #Choose if testing contour files.        
        condition1 = '.avs' in filetype
        condition2 = '.csv' in filetype  
        if condition1 or condition2:
            #If no pre-specified times, grab them from f_dif.
            if len(times) == 0:    
                times = f_dif.times
        
            msg = 'Incorrect %s at time %s.'
            
            #Check the variables at each time for any significant differences.
            for t in times: 
                #Its possible some times do not have all variables in f_dif.
                for v in np.intersect1d(variables, f_dif[t]):
                    f_dif[t][v] = map(abs, f_dif[t][v])   
            	    self.assertTrue(max(f_dif[t][v]) < mxerr, msg%(v, t))
        
        #Choose if testing history files.       
        elif '.his' in filetype:
            #If no pre-specifed nodes, grab them from f_dif.
            if len(nodes) == 0:
                nodes = f_dif.nodes
        
            msg = 'Incorrect %s at node %s.'  
            
            #Check the nodes at each variable for any significant differences.   
            for v in variables:
                #Its possible some variables do not have all nodes in f_dif.
                for n in np.intersect1d(nodes, f_dif[v]):     
            	    self.assertTrue(max(f_dif[v][n])<mxerr, msg%(v, n))

        #Choose if testing output files.
        elif '.out' in filetype:
            #If no pre-specified componets, grab them from f_dif.
            if len(components) == 0:
                components = f_dif.components
                
            #If no pre-specifed nodes, grab them from f_dif.
            if len(nodes) == 0:
                nodes = f_dif.nodes
                
            msg = 'Incorrect %s at %s node %s.'  
            
            #Check the node at each component for significant differences.   
            for c in components:
                for n in nodes:
                    for v in variables:
                        difference = max(f_dif.node[c][n][v])
                        self.assertTrue(difference < mxerr, msg%(v,c,n))
                        
    def _fgeneral(self, file_pattern):
        """
        Process General fdata
        
        Chooses the correct fdata method to process FEHM calculations.
        
        :param file_pattern: File pattern to processed.
        :type file_pattern:  str
        
        Author: Mark Lange
        Updated: June 2014 
        """
    
        #Is a history file even if it ends in .avs or .csv if _his in pattern.
        if '_his.' in file_pattern:
            return fhistory(file_pattern)  
        elif '.his' in file_pattern:
            return fhistory(file_pattern)
        elif '.avs' in file_pattern:
            return fcontour(file_pattern)
        elif '.csv' in file_pattern:
            return fcontour(file_pattern)
        elif '.out' in file_pattern:
            #Assuming each subcase has only 1 file, use glob to find it.
            filename = glob(file_pattern)[0]
            return foutput(filename)
            	    
    def setUp(self):
        self.maindir = os.getcwd()

    def _cleanup(self,files):
        """ 
        Utility function to remove files after test

        :param files: list of file names to remove
        :type files: lst(str) 
        """
            
        for g in files:
            for f in glob(g):
                if os.path.exists(f): os.remove(f)

    def _run_fehm(self, subcase):
        """ 
        **Utility function to run fehm**
        
        Asserts that fehm terminates successfully.

        :param filesfile: name of fehm files file
        :type filesfile: str 
        """
        
        #Find the control file for the test-case or for the subcase. 
        if subcase == '':
            filesfile = 'input/control/fehmn.files'
        else:
            filesfile = 'input/control/'+subcase+'.files'
            
        call(exe+' '+filesfile, shell=True, stdout=PIPE)
        outfile = None
        errfile = 'fehmn.err'

        with open( filesfile, 'r' ) as f:
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
 
        complete = False
        if outfile:
            with open( outfile, 'r' ) as f:
                for line in reversed(f.readlines()):
                    if 'End Date' in line:
                        complete = True
                        break
                        
        if os.path.exists(errfile): 
            errstr = open( errfile, 'r' ).read()
        else: 
            errstr = ''
            
        # Change to maindir in case assertTrue fails    
        curdir = os.getcwd() 
        os.chdir(self.maindir)
        
        msg = 'Unsuccessful fehm simulation\nContents of '
        self.assertTrue(complete, msg+errfile+':\n\n'+errstr)
        os.chdir(curdir)
      
def suite(mode, test_case):
    suite = unittest.TestSuite()
    
    #Default mode is admin for now. Should it be different?
    if mode == 'admin' or mode == 'default':
        suite.addTest(Tests('test_saltvcon'))
        suite.addTest(Tests('test_dissolution'))
        suite.addTest(Tests('test_salt_perm_poro'))
        suite.addTest(Tests('test_avdonin'))
        suite.addTest(Tests('test_boun'))
        suite.addTest(Tests('test_cden'))
        suite.addTest(Tests('test_doe'))       
        suite.addTest(Tests('test_head'))
        suite.addTest(Tests('test_ramey'))
        suite.addTest(Tests('test_theis'))
        suite.addTest(Tests('test_dryout'))
        
        #TODO - Look into why this test take so long.
        #suite.addTest(Tests('test_evaporation'))
        
        #TODO - Figure out how to read some other formats.
        #suite.addTest(Tests('test_sptr_btc'))
        #suite.addTest(Tests('test_sorption'))
        #suite.addTest(Tests('test_particle_capture'))
        #suite.addTest(Tests('test_multi_solute'))
        #suite.addTest(Tests('test_mptr'))
        #suite.addTest(Tests('test_lost_part'))
        #suite.addTest(Tests('test_chain'))
        #suite.addTest(Tests('test_co2test'))
        #suite.addTest(Tests('test_convection'))
        #suite.addTest(Tests('test_dpdp_rich'))
        #suite.addTest(Tests('test_erosion'))
        #suite.addTest(Tests('test_gdpm'))
        #suite.addTest(Tests('test_forward'))
    
    elif mode == 'developer':
        #This mode will be a reduced set that runs faster.
        pass
             
    elif mode == 'solo':
        suite.addTest(Tests(test_case))
            
    elif mode == 'admin':
        pass
        
    return suite

if __name__ == '__main__':
    
    #Unless the user specifies a single test-case, this isn't important.
    test_case = ''

    #Set up command-line interface.
    parser = argparse.ArgumentParser(description='FEHM Test-Suite')
    
    group = parser.add_mutually_exclusive_group()
    
    a = 'store_true'
    h = 'Run the entire test-suite.'
    group.add_argument('-a', '--admin', help=h, action=a)
    
    h = 'Run a portion of the test-suite.'
    group.add_argument('-d', '--dev', help=h, action=a)
    
    h = 'Run a single test-case.'
    group.add_argument('-s', '--solo', help=h, action=a)

    h = 'Path to the FEHM executable.'
    parser.add_argument('exe', help=h)
    
    h = 'Single test-case to run.'
    parser.add_argument('testcase', nargs='?', help=h, default=None)
     
    args = vars(parser.parse_args())
    
    exe = args['exe']
        
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
            else:       
                parser.print_help()
        
    runner = unittest.TextTestRunner(verbosity=2)
    test_suite = suite(mode, test_case)
    runner.run(test_suite)


