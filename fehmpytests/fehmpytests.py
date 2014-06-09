import unittest
import os,sys
import re
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
        '''  Dissolution

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

            The porosity-permeability function for compacted salt from Cinar et at. (2006) is tested
            using a six node problem with porosities from 0.01 to 0.2. The excel spreadsheet in 
            ./salt_perm_poro/salt-perm-poro.xlsx contains calculations of the perm-poro function.

            Cinar, Y, G Pusch and V Reitenbach (2006) Petrophysical and capillary properties of compacted
                salt. Transport in Porous Media. 64, p. 199-228, doi: 10.1007/s11242-005-2848-1 '''
        
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
      
    def avdonin(self):
        ''' Avdonin Radial Heat and Mass Transfer Test 
        Tests that the temperatures at each time are correct.
        Modified by mlange806@gmail.com on June 3, 2014 '''
        
        #Error threshold - May need to change.
        maxerr = 1.e-4 
        
        #Change to test directory.
        os.chdir('avdonin')
        
        #Test three different cases.
        cases = ['84', '400', '800']
        for case in cases: 
            #Read in comparison files.
            f_old_con = fcontour('compare'+case+'.*_sca_node.csv')
            f_old_his = fhistory('compare'+case+'_temp.his')
            
            #Create fehmn.files for current case.
            generic_files = open('generic_fehmn.files')
            data = generic_files.read()
            data = re.sub('case', case, data)
            files = open('fehmn.files', 'w')
            files.write(data)
            generic_files.close()
            files.close()
            
            #Read in new output files.
            self.run_fehm()
            f_new_con = fcontour('avdonin'+case+'.*.csv')
            f_new_his = fhistory('avdonin'+case+'_temp.his')
            
            #Find the difference between the old and new. 
            f_dif_con = fdiff(f_new_con, f_old_con)
            f_dif_his = fdiff(f_new_his, f_old_his)
            
            #Get the contour information.
            times = f_dif_con.times
            variables = f_dif_con.variables
            
            #Error Message for Incorrect Contour    
            msg = 'Incorrect %s at time %s.'
            
            #Check that the new contour files are still the same.
            for t, v in [(t, v) for t in times for v in variables]:
            	self.assertTrue(max(f_dif_con[t][v])<maxerr, msg%(v, t))
            
            #Get the history information.
            variables = f_dif_his.variables
            nodes = f_dif_his.nodes
            
            #Error Message for Incorrect History    
            msg = 'Incorrect %s at node %s.'
            	
            #Check that the new history files are still the same.
            for v, n in [(v, n) for v in variables for n in nodes]:
            	self.assertTrue(max(f_dif_his[v][n])<maxerr, msg%(v, n))   
            
            #Remove created files.
            trash = ['avdonin'+case+'.*.csv', '*.avs_log', '*.con', '*.geo', 
                     'avdonin*.his', '*.out', '*.err', 'all', 'fehmn.files']
            self.cleanup(trash)
        
        #Return to the main directory.       
        os.chdir(self.maindir)
        
    def barometric(self):
    	""" Comming Soon """
    	pass
        
    def binmode(self):
        """ Comming Soon """
        pass
        
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
                   
    # UTILITIES ######################################################
    
    def _test_case(self, name):
        """ General Test Case 
        Should be able to test any test case.
        Modified by mlange806@gmail.com on June 4, 2014. """
        
        #Error Threshold
        maxerr = 1.e-4
        
        #Change to test directory.
        os.chdir(name)
        
        #Find the subcases.
        subcases = self.getSubcases()
        
        #Test each subcase.
        for subcase in subcases:
            #Test each output file type.
            file_types = ['*.csv', '*.his']
            for file_type in file_types:
                #Check to make sure there are files of this type.
                if len(glob('compare/'+file_type)) > 0:
                    #Read in the old comparison files.
                    f_old = self.fgeneral('compare/'+file_type)
                    
                    #Create fehmn.files for current case.
                    generic_files = open('input/generic_fehmn.files')
                    data = generic_files.read()
                    data = re.sub('N', subcase, data)
                    data = re.sub('UM', '', data)
                    files = open('fehmn.files', 'w')
                    files.write(data)
                    generic_files.close()
                    files.close()
                    
                    #Read in new files.
                    self.run_fehm()
                    f_new = self.fgeneral(file_type)
                    
                    #Find the difference between the two files.
                    f_dif = fdiff(f_new, f_old)
                    
                    #Get the history information.
                    variables = f_dif.variables
                    nodes = f_dif.nodes
                    
                    #Error Message for Incorrect History    
                    msg = 'Incorrect %s at node %s.'
                    	
                    #Check that the new history files are still the same.
                    for v, n in [(v, n) for v in variables for n in nodes]:
                        #If the key combination exists, test.
                        try:
                    	    self.assertTrue(max(f_dif[v][n])<maxerr, msg%(v, n))
                    	#Otherwise, ignore.
                    	except:  
                            pass
                            
                #There are no files of this type so ignore this test.    
                else:
                    pass
                                 
        #Remove all files created outside compare and input.
        self.cleanup(['*.*'])
        
        #Return to the main directory.       
        os.chdir(self.maindir)

    def setUp(self):
        # Set location of main directory
        self.maindir = os.getcwd()

    def cleanup(self,files):
        ''' Utility function to remove files after test

            :param files: list of file names to remove
            :type files: lst(str)
        '''
        for g in files:
            for f in glob(g):
                if os.path.exists(f): os.remove(f)

    def run_fehm(self, filesfile='fehmn.files'):
        """
            Utility function to run fehm
            Asserts that fehm terminates successfully

            :param filesfile: name of fehm files file
            :type filesfile: str
        """
        
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
        
    def getSubcases(self):
        """ Get Cases
        Assumming that subcases are numbers, returns a set of subcases using a 
        test-case's comparison files.
        
        *Must be inside the test-case folder.
        
        Developed by mlange806@gmail.com on June 4, 2014 """
    
        #Find the names of every comparison file.
        types = ['*.csv', '*.his']
        file_names = []
        for t in types:
            file_names = file_names+glob('compare/'+t)
        
        #Using the comparison files, extract the subcase numbers.    
        subcases = []
        for file_name in file_names:
            pattern = re.compile(r'\d+')
            subcase = pattern.findall(file_name)[0]
            if subcase not in subcases:
                subcases.append(subcase)
                
        return subcases
        
    def fgeneral(self, file_pattern):
        #Chooses the correct fpost object to represent output files.
        if '.csv' in file_pattern:
            return fcontour(file_pattern)
        elif '.his' in file_pattern:
            return fhistory(file_pattern)
             
def suite(case):
    suite = unittest.TestSuite()
    if case == 'all':
        suite.addTest(Tests('saltvcon'))
        suite.addTest(Tests('dissolution'))
        suite.addTest(Tests('salt_perm_poro'))
        suite.addTest(Tests('avdonin'))
        #suite.addTest(Tests('barometric'))
        #suite.addTest(Tests('binmode'))
        #suite.addTest(Tests('test_boun'))
        suite.addTest(Tests('test_cden'))
        #suite.addTest(Tests('test_cden'))
        suite.addTest(Tests('test_chain'))
        #suite.addTest(Tests('test_co2test'))
        #suite.addTest(Tests('test_convection'))
        #suite.addTest(Tests('test_doe'))
        suite.addTest(Tests('test_dpdp_rich'))
        #suite.addTest(Tests('test_erosion'))
        #suite.addTest(Tests('test_evaporation'))
        suite.addTest(Tests('test_forward'))
        suite.addTest(Tests('test_gdpm'))
        #suite.addTest(Tests('test_head'))
        #suite.addTest(Tests('test_lost_part'))
        suite.addTest(Tests('test_mptr'))
        
    elif case == 'developer':
        pass
    elif case == 'admin':
        pass
    return suite

if __name__ == '__main__':
    if len(sys.argv) > 1:
        exe = os.path.abspath(sys.argv[1])
    else:
        print "Usage: python fehmpytests.py 'fehm executable'"
        os._exit(0)
    # Case set to all for now
    case = 'all'
    runner = unittest.TextTestRunner(verbosity=2)
    test_suite = suite(case)
    runner.run (test_suite)


