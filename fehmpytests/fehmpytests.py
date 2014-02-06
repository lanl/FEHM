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

    # TESTS ######################################################

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
            associated calculations.
        '''
        # Relative error theshold
        maxerr = 1.e-4
        # Test directory name
        dir = 'saltvcon'
        # Change directory to test directory
        os.chdir(dir)
        # Read in comparison output file
        f_old = foutput('intact-compare.out')
        # Run intact salt model
        self.run_fehm('intact.files')
        # Read in simulation output file
        f_new = foutput('intact.out')
        f_dif = fdiff(f_new,f_old,components=['water'],variables=['Kx'],format='relative')
        for n in f_dif.nodes:
            self.assertTrue(f_dif.node['water'][n]['Kx'][0]<maxerr, '\nIncorrect intact salt thermal conductivity for node '+str(n)+'\nRelative error: '+str(f_dif.node['water'][n]['Kx'][0])+', Threshold: '+str(maxerr)+'\nExpected='+str(f_old.node['water'][n]['Kx'][0])+' Simulated='+str(f_new.node['water'][n]['Kx'][0]))
        del f_old, f_new, f_dif

        # Read in comparison output file
        f_old = foutput('crushed-compare.out')
        # Run intact salt model
        self.run_fehm('crushed.files')
        # Read in simulation output file
        f_new = foutput('crushed.out')
        f_dif = fdiff(f_new,f_old,components=['water'],variables=['Kx'],format='relative')
        for n in f_dif.nodes:
            self.assertTrue(f_dif.node['water'][n]['Kx'][0]<maxerr, '\nIncorrect intact salt thermal conductivity for node '+str(n)+'\nRelative error: '+str(f_dif.node['water'][n]['Kx'][0])+', Threshold: '+str(maxerr)+'\nExpected='+str(f_old.node['water'][n]['Kx'][0])+' Simulated='+str(f_new.node['water'][n]['Kx'][0]))
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
        #############################################################
        # Relative error theshold
        maxerr = 1.e-4
        # Test directory name
        dir = 'dissolution'
        # Change directory to test directory
        os.chdir(dir)
        #############################################################
        # Read in comparison files
        f_old = fcontour('compare.*_con_node.csv')
        # Run intact salt model
        self.run_fehm()
        # Read in new output files
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
                salt. Transport in Porous Media. 64, p. 199-228, doi: 10.1007/s11242-005-2848-1
        '''
        cwd = os.getcwd()
        # Relative error theshold
        maxerr = 1.e-4
        # Test directory name
        dir = 'salt_perm_poro'
        # Change directory to test directory
        os.chdir(dir)
        # Read in comparison files
        f_old = fcontour('compare.00001_sca_node.csv')
        # Run intact salt model
        self.run_fehm()
        # Read in new output files
        f_new = fcontour('run.00001_sca_node.csv')
        # Diff new and old files
        f_dif = fdiff(f_new,f_old,variables=['n','perm_x'])
        # Test for correct permeabilities
        for node,dif,k_old,k_new in zip(f_dif[1]['n'],f_dif[1]['perm_x'],f_old[1]['perm_x'],f_new[1]['perm_x']):
            self.assertTrue(dif<maxerr, '\nIncorrect permeability at node '+str(node)+'. Expected '+str(k_old)+', Simulated '+str(k_new)) 
        # Remove created files
        self.cleanup(['nop.temp','fehmn.err','run*.csv','*.out','run.avs_log'])
        os.chdir(self.maindir)

    # UTILITIES ######################################################

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
        ''' Utility function to run fehm
            Asserts that fehm terminates successfully

            :param filesfile: name of fehm files file
            :type filesfile: str
        '''
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
            if outfile is None and ':' not in lines[0]: outfile = lines[3].strip()
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



def suite(case):
    suite = unittest.TestSuite()
    if case == 'all':
        suite.addTest(Tests('saltvcon'))
        suite.addTest(Tests('dissolution'))
        suite.addTest(Tests('salt_perm_poro'))
    elif case == 'developer':
        pass
    elif case == 'admin':
        pass
    return suite

if __name__ == '__main__':
    if len(sys.argv) > 1:
        exe = os.path.abspath(sys.argv[1])
    else:
        print "Usage: python fehm_tests.py 'fehm executable'"
        os._exit(0)
    # Case set to all for now
    case = 'all'
    runner = unittest.TextTestRunner(verbosity=2)
    test_suite = suite(case)
    runner.run (test_suite) 


