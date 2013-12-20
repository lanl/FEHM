import unittest
import os,sys
import re
from subprocess import call, PIPE
__unittest = True # Suppresses tracebacks

class Tests(unittest.TestCase):

    def cleanup(self,files):
        ''' Utility function to remove files after test

            :param files: list of file names to remove
            :type files: lst(str)
        '''
        for f in files:
            if os.path.exists(f): os.remove(f)

    def saltvcon(self):
        ''' saltvcon macro test

            Modification: new Dec 19 2013 dharp@lanl.gov
            Modification: 

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
        '''
        cwd = os.getcwd()
        # Relative error theshold
        maxerr = 1.e-4
        # Test directory name
        dir = 'saltvcon'
        # Change directory to test directory
        os.chdir(dir)
        #############################################################
        # Calculate thermal conductivities
        T = [20,40,60,80,100,120] # Temperature [C] at nodes
        phi = [0.01,0.1,0.3,0.5,0.7,0.9] # Porosity
        Ai = [26.85,5.4,1.14] # Intact salt coefficients
        Ac = [26.85,1.08,-270,370,-136,1.5,5,1.14] # Crushed salt coefficients
        #Calculate intact salt thermal conductivities
        kxi = []
        for t in T:
            kxi.append(Ai[1]*((Ai[0]+273.15)/(t+273.15))**Ai[2])
        #Calculate crushed salt thermal conductivities
        kxc = []
        for t,p in zip(T,phi):
            kx_asse = Ac[2]*p**4+Ac[3]*p**3+Ac[4]*p**2+Ac[5]*p+Ac[6]
            kx_temp = Ac[1]*kx_asse*((Ac[0]+273.15)/(t+273.15))**Ac[7]
            if kx_temp < 1.e-6:
                kx_temp = 1.e-6
            kxc.append(kx_temp)
        #############################################################
        # Run intact salt model
        call(exe+' intact.files', shell=True, stdout=PIPE)
        # Open FEHM output file and read
        f = open('intact.out', 'r')
        data = f.readlines()
        f.close()
        os.remove('intact.out')
        # Search for thermal conductivity column heading 'Kx'
        for i in range(len(data)):
            if re.search('Kx',data[i]): break
        # Collect thermal conductivities
        kx_sim = []
        for j in range(6):
            i += 1
            kx_sim.append(float(data[i].split()[3])) # Grab values in 4th column for each node 
        # Compare true and simulated values
        nodeno = 1
        for true,sim in zip(kxi,kx_sim):   
            relerr = abs(true-sim)/true
            self.assertTrue(relerr<maxerr, '\nIncorrect intact salt thermal conductivity for node '+str(nodeno)+'\nRelative error: '+str(relerr)+', Threshold: '+str(maxerr)+'\nExpected='+str(true)+' Simulated='+str(sim))
            nodeno += 1
        #############################################################
        # Run crushed salt model
        call(exe+' crushed.files', shell=True, stdout=PIPE)
        # Open FEHM output file and read
        f = open('crushed.out', 'r')
        data = f.readlines()
        f.close()
        os.remove('crushed.out')
        # Search for thermal conductivity column heading 'Kx'
        for i in range(len(data)):
            if re.search('Kx',data[i]): break
        # Collect thermal conductivities
        kx_sim = []
        for j in range(6):
            i += 1
            kx_sim.append(float(data[i].split()[3])) # Grab values in 4th column for each node 
        # Compare true and simulated values
        nodeno = 1
        for true,sim in zip(kxc,kx_sim):   
            relerr = abs(true-sim)/true
            self.assertTrue(relerr<maxerr, '\nIncorrect crushed salt thermal conductivity for node '+str(nodeno)+'\nRelative error: '+str(relerr)+', Threshold: '+str(maxerr)+'\nExpected='+str(true)+' Simulated='+str(sim))
            nodeno += 1
        #############################################################
        # Return to main directory
        self.cleanup(['nop.temp','fehmn.err'])
        os.chdir(cwd)


def suite(case):
    suite = unittest.TestSuite()
    if case == 'all':
        suite.addTest(Tests('saltvcon'))
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


