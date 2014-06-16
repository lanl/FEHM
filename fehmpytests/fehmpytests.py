import unittest
import os
import sys
import re
import numpy as np
from subprocess import call, PIPE
try:
    from fdata import *
except ImportError:
    try:
        sys.path.insert(0,'../pyfehm')
        from fdata import *
    except ImportError as err:
        print 'ERROR: Unable to import pyfehm fpost module'
        print err
        os._exit(0)
from glob import glob

__unittest = True # Suppresses tracebacks

class Tests(unittest.TestCase):

    # TESTS ######################################################### 

    def saltvcon(self):
        ''' Salt variable conductivity

            Modification: new Dec 19 2013 dharp@lanl.gov
            Modification: updated to use f_dif function Feb 6 2014

            This function tests the saltvcon macro which calculates thermal
            conductivity of crushed and intact salt.

            Intact salt:
            kxi = k_{t-300}(300/T)^1.14
            Munson et al. (1990) Overtest for Simulate Defense High-Level
                Waste (Room B): In Situ Data Report. WIPP. Sandia National
                Laboratories, SAND89-2671

            Crushed salt:
            Thermal conductivity of crushed salt from Asse mine
            kx_asse = -270*phi^4+370*phi^3-136*phi^2+1.5*phi+5  
            Bechtold et al. (2004) Backfilling and sealing of underground
                respositories for radioactive waste in salt (BAMBUS II
                project), EUR 20621, ISBN 92-894-7767-9
            kx = (k_{t-300}/kx_asse)*(300/T)^1.14
            if kx is less then 1.e-6, set to 1.e-6

            The excel spreadsheet in ./saltvcon/saltvcon.xlsx contains the
            associated calculations. '''
            
        # Relative error theshold
        maxerr = 1.e-4
           
        # Change directory to test directory
        os.chdir('saltvcon')
        
        # Read in comparison output file
        f_old = foutput('intact-compare.out')
        
        # Run intact salt model
        self.run_fehm('intact.files')
        
        # Read in simulation output file
        f_new = foutput('intact.out')
        
        #Calculate the difference between the new and the old data.
        f_dif = fdiff(f_new,f_old,components=['water'],variables=['Kx'],format='relative')
        
        for n in f_dif.nodes:
            self.assertTrue(f_dif.node['water'][n]['Kx'][0]<maxerr, '\nIncorrect intact salt thermal conductivity for node '+str(n)\
                                                                   +'\nRelative error: '+str(f_dif.node['water'][n]['Kx'][0])+\
                                                                    ', Threshold: '+str(maxerr)+'\nExpected='+\
                                                                    str(f_old.node['water'][n]['Kx'][0])+' Simulated='+\
                                                                    str(f_new.node['water'][n]['Kx'][0]))
                                                                                             
        del f_old, f_new, f_dif

        # Read in comparison output file
        f_old = foutput('crushed-compare.out')
        
        # Run intact salt model
        self.run_fehm('crushed.files')
        
        # Read in simulation output file
        f_new = foutput('crushed.out')
        
        #Calculate the difference between the new and the old data.
        f_dif = fdiff(f_new,f_old,components=['water'],variables=['Kx'],format='relative')
        
        for n in f_dif.nodes:
            self.assertTrue(f_dif.node['water'][n]['Kx'][0]<maxerr, '\nIncorrect intact salt thermal conductivity for node '+str(n)\
                                                                                             +'\nRelative error: '+str(f_dif.node['water'][n]['Kx'][0])+\
                                                                                             ', Threshold: '+str(maxerr)+'\nExpected='+\
                                                                                             str(f_old.node['water'][n]['Kx'][0])+' Simulated='\
                                                                                             +str(f_new.node['water'][n]['Kx'][0]))
                                                                                             
        # Remove files
        self.cleanup(['nop.temp','fehmn.err','*.avs*','*_head','intact.out','crushed.out'])
        
        # Return to main directory
        os.chdir(self.maindir)

    def dissolution(self):
        ''' Dissolution

            Modification: new Jan 7 2014 dharp@lanl.gov
            Modification: 

            A one-dimensional transport simulation of calcite (CaC03(s)) dissolution is
            tested. Profiles of concentration versus reactor length, at selected times,
            will be compared against the analytical solution.

            Details of this test are described in the FEHM V2.21 Validation Test Plan
            on pages 93-95 
            (STN: 10086-2.21-00, Rev.No. 00, Document ID: 10086-VTP-2.21-00, August 2003) 
            
        '''
             
        # Relative error theshold
        maxerr = 1.e-4
        # Test directory name
        dir = 'dissolution'
        # Change directory to test directory
        os.chdir(dir)

        # Read in comparison files
        f_old = fcontour('compare.*_con_node.csv')

        # Read in new output files
        self.run_fehm()
        f_new = fcontour('dissolution.*_con_node.csv')

        # Diff new and old files
        f_dif = fdiff(f_new,f_old,variables=['Np[aq] (Moles/kg H20)'])

        # Test for correct concentrations
        for t in f_dif.times:
            self.assertTrue(f_dif[t]['Np[aq] (Moles/kg H20)'].all()<maxerr, '\nIncorrect concentration at time '+str(t))
            
        # Remove created files
        self.cleanup(['nop.temp','fehmn.err','dissolution*.csv','*.avs_log','*geo','*.out','*.trc','*.his','*_head'])
        os.chdir(self.maindir)

    def salt_perm_poro(self):
        ''' Salt perm-poro function

            Modification: new Jan 29 2014 dharp@lanl.gov
            Modification: 

            The porosity-permeability function for compacted salt from Cinar et 
            at. (2006) is tested using a six node problem with porosities from 
            0.01 to 0.2. The excel spreadsheet in 
            ./salt_perm_poro/salt-perm-poro.xlsx contains calculations of the 
            perm-poro function.

            Cinar, Y, G Pusch and V Reitenbach (2006) Petrophysical and 
            capillary properties of compacted salt. Transport in Porous Media. 
            64, p. 199-228, doi: 10.1007/s11242-005-2848-1 '''
        
        cwd = os.getcwd()
        
        # Relative error theshold
        maxerr = 1.e-4
        
        # Test directory name
        dir = 'salt_perm_poro'
        
        # Change directory to test directory
        os.chdir(dir)
        # Read in comparison files
        f_old = fcontour('compare.00001_sca_node.csv')
        
        # Read in new output files
        self.run_fehm()
        f_new = fcontour('run.00001_sca_node.csv')
        
        # Diff new and old files
        f_dif = fdiff(f_new,f_old,variables=['n','perm_x'])
        
        # Test for correct permeabilities
        for node,dif,k_old,k_new in zip(f_dif[1]['n'],f_dif[1]['perm_x'],f_old[1]['perm_x'],f_new[1]['perm_x']):
            self.assertTrue(dif<maxerr, '\nIncorrect permeability at node '+str(node)+'. Expected '+str(k_old)+', Simulated '\
                                                      +str(k_new)) 
            
        # Remove created files
        self.cleanup(['nop.temp','fehmn.err','run*.csv','*.out','run.avs_log'])
        os.chdir(self.maindir)
        
    def test_saltvcon(self):
        """
        Tests that function calculates relatively correct Kx values in water. 
        Information about the function can be found in the folder 'information'.
        """
        
        arguments = {}
        arguments['components'] = ['water']
        arguments['variables']  = ['Kx']
        arguments['format'] = 'relative' 
        
        self._test_case('saltvcon', arguments)
        
        
    def test_dissolution(self):
        """
        Tests that the function calculates the correct concentrations. Because
        the folder structure is not set up yet to handle no subcases, the empty 
        list is used to find the folder 'subcase-' which has the control file. 
        Information about the function can be found in the folder 'information'.
        """
    
        arguments = {}
        arguments['variables'] = ['Np[aq] (Moles/kg H20)']
        arguments['subcases']  = ['']
        
        self._test_case('dissolution', arguments)
        
    def test_salt_perm_poro(self):
        """
        Tests that the function calculates the correct porosities. Because
        the folder structure is not set up yet to handle no subcases, the empty 
        list is used to find the folder 'subcase-' which has the control file. 
        Information about the function can be found in the folder 'information'.
        """
    
        arguments = {}
        arguments['variables'] = ['n', 'perm_x']    
        arguments['subcases']  = ['']
        
        self._test_case('salt_perm_poro', arguments)
  
    def test_avdonin(self):
        self._test_case('avdonin')
        
    def test_boun(self): 
        self._test_case('boun_test') 
        
    def test_cflxz(self):
        self._test_case('cflxz')
        
    def test_cden(self):
        self._test_case('cden_test')
        
    def test_chain(self):
        self._test_case('chain')
        
    def test_co2test(self):
        self._test_case('co2test')
        
    def test_convection(self):
        self._test_case('convection')
        
    def test_doe(self):
        self._test_case('doe')
    
    def test_dpdp_rich(self):
        self._test_case('dpdp_rich')
        
    def test_erosion(self):
        self._test_case('erosion_test')
    
    def test_evaporation(self):
        self._test_case('evaporation')
        
    def test_forward(self):
        self._test_case('forward')
        
    def test_gdpm(self):
        self._test_case('gdpm')
    
    def test_head(self):
        self._test_case('head')
        
    def test_lost_part(self):
        self._test_case('lost_part')
    
    def test_mptr(self):
        self._test_case('mptr_test')
        
    def test_multi_solute(self):
        self._test_case('multi_solute')
        
    def test_particle_capture(self):
        self._test_case('particle_capture')
    
    def test_ramey(self):
        self._test_case('ramey')
        
    def test_sorption(self):
        self._test_case('sorption')
        
    def test_sptr_btc(self):
        self._test_case('sptr_btc')
    
    def test_theis(self):
        self._test_case('theis')
                   
    # UTILITIES ######################################################
    
    def _test_case(self, name, parameters={}):
        """ General Test Case 
        Should be able to test any test case with the following folder set-up:
            test-case-name: [compare, input, subcase-N1, subcase-N2, etc.],
        where compare contains data known to be correct, input contains the
        needed input files, and subcase-X contains a control file for a subcase.
        
        Normally tests all values using a relative difference that must be less 
        than 1.e-4, but these values can be set by passing a dictionary with 
        keys variables, times, nodes, components and/or format assigned to the 
        new values you would like to use instead. """ 
         
        os.chdir(name)
        
        #Search for folders with 'subcase-' and extract subcase name.
        filenames = glob('subcase-*')
        subcases = []
        for filename in filenames:
            subcases.append(re.sub('subcase-', '', filename))

        #Test the new files generated with each subcase.
        for subcase in subcases:
            filetypes = ['*.avs', '*.csv', '*.his', '*.out']
            for filetype in filetypes:
                #Check to make sure there are files of this type.
                if len(glob('compare/'+'*'+subcase+filetype)) > 0:
                    #Check for any significant differences.
                    try:
                        self._checkDifferences(filetype, subcase, parameters)
                    finally:
                        #Allows other tests to be performed after exception.
                        self.cleanup(['*.*'])
                        os.chdir(self.maindir)
         
    def _checkDifferences(self, filetype, subcase, parameters={}):
        """ Check Difference
        Checks the difference between old and new output. Fails if the
        difference is greater than maxerr, defaulted to 1.e-4. """
        
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
        self.run_fehm(subcase)   
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
            #If no prespecified componets, grab them from f_dif.
            if len(components) == 0:
                components = f_dif.components
                
            #If no pre-specifed nodes, grab them from f_dif.
            if len(nodes) == 0:
                nodes = f_dif.nodes
                
            msg = 'Incorrect %s at %s node %s.'  
            
            #Check the node at each component for significant differences.   
            for c in components:
                #Its possible some componets do not have all nodes in f_dif.
                for n in nodes:
                    #In case variables are different per node, intersect.
                    for v in variables:
                        self.assertTrue(f_dif.node[c][n][v]<mxerr, msg%(v,c,n))
   
    def _fgeneral(self, file_pattern):
        #Chooses the correct fpost object to represent output files.
        if '.avs' in file_pattern:
            return fcontour(file_pattern)
        elif '.csv' in file_pattern:
            return fcontour(file_pattern)
        elif '.his' in file_pattern:
            return fhistory(file_pattern)
        elif '.out' in file_pattern:
            #Assuming each subcase has only 1 file, use glob to find it.
            filename = glob(file_pattern)[0]
            return foutput(filename)
            	    
    def setUp(self):
        self.maindir = os.getcwd()

    def cleanup(self,files):
        """ Utility function to remove files after test

            :param files: list of file names to remove
            :type files: lst(str) """
            
        for g in files:
            for f in glob(g):
                if os.path.exists(f): os.remove(f)

    def run_fehm(self, subcase):
        """ Utility function to run fehm
            Asserts that fehm terminates successfully

            :param filesfile: name of fehm files file
            :type filesfile: str """
        
        filesfile = 'subcase-'+subcase+'/fehmn.files'
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
            if outfile is None and ':' not in lines[0]: outfile=lines[3].strip()
 
        complete = False
        if outfile:
            with open( outfile, 'r' ) as f:
                for line in reversed(f.readlines()):
                    if 'End Date' in line:
                        complete = True
                        break
                        
        if os.path.exists(errfile): errstr = open( errfile, 'r' ).read()
        else: errstr = ''
        curdir = os.getcwd()
        # Change to maindir in case assertTrue fails
        os.chdir(self.maindir)
        self.assertTrue(complete, 'Unsuccessful fehm simulation\nContents of '+errfile+':\n\n'+errstr)
        os.chdir(curdir)
      
def suite(case, test_case):
    suite = unittest.TestSuite()
    
    if case == 'all':
        suite.addTest(Tests('test_saltvcon'))
        suite.addTest(Tests('test_dissolution'))
        suite.addTest(Tests('test_salt_perm_poro'))
        suite.addTest(Tests('test_avdonin'))
        suite.addTest(Tests('test_boun'))
        suite.addTest(Tests('test_cden'))
        
        #TODO - Find or generate some compare files.
        #suite.addTest(Tests('test_chain'))
        #suite.addTest(Tests('test_co2test'))
        #suite.addTest(Tests('test_convection'))
        #suite.addTest(Tests('test_doe'))
        #suite.addTest(Tests('test_dpdp_rich'))
        #suite.addTest(Tests('test_erosion'))
        #suite.addTest(Tests('test_evaporation'))
        #suite.addTest(Tests('test_forward'))
        #suite.addTest(Tests('test_gdpm'))
        #suite.addTest(Tests('test_head'))
        #suite.addTest(Tests('test_lost_part'))
        #suite.addTest(Tests('test_mptr'))
        #suite.addTest(Tests('test_multi_solute'))
        #suite.addTest(Tests('test_particle_capture'))
        #suite.addTest(Tests('test_ramey'))
        #suite.addTest(Tests('test_sorption'))
        #suite.addTest(Tests('test_sptr_btc'))
        #suite.addTest(Tests('test_theis'))
             
    elif case == 'single':
        suite.addTest(Tests(test_case))
        
    elif case == 'developer':
        pass
        
    elif case == 'admin':
        pass
        
    return suite

if __name__ == '__main__':
    #By default, run all test cases.
    if len(sys.argv) > 1:
        exe = os.path.abspath(sys.argv[1])
        mode = 'all'
        test_case = ''
    
        #Single mode can (currently) specify 1 test-case to run.    
        if len(sys.argv) > 2:
            mode = 'single'
            test_case = sys.argv[2]
                       
    else:
        print "Usage: python fehmpytests.py fehm-executable test(optional)"  
        os._exit(0)
    
    runner = unittest.TextTestRunner(verbosity=2)
    test_suite = suite(mode, test_case)
    runner.run(test_suite)


