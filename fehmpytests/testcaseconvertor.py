import os
import sys
import re
import glob
import distutils.core
import shutil

class TestCaseConvertor():
    """ Test Case Convertor 
    Converts old test-cases into new test-cases.
    
    Changed by mlange806@gmail.com on June 11, 2014. """
    
    def __init__(self, new_path):        
        self._old_path = \
          '/n/swdev/FEHM_dev/VERIFICATION_linux/VERIFICATION_V3.2linux/'        
        self._new_path = new_path  
    
    def convertCase(self, name):
        """ Convert Case
        Takes a test-case 'name' from the old path and copies files to the new 
        test directory. """
        
        #Find subcases.
        subcases = self._getSubcases(name)
     
        #Navigate to fehmpytest.
        os.chdir(self._new_path)
        
        #Create folder for test case, unless it already exists.
        try:
            self._createTestFolder(name, subcases)
        except:
            print 'ERROR: There is already a folder for this test case.'
            os._exit(0)
            
        #Copy the input files from old to new test suite.
        source = self._old_path+name+'/input'
        dest = self._new_path+name+'/input'
        distutils.dir_util.copy_tree(source, dest)  
        
        #Create a control file for each subcase.
        for subcase in subcases:
            self._primeControlFile(name, subcase)
  
        #Copy the comparison files over from the old test-case.
        self._copyCompareFiles(name)
       
        print 'Done.'
        
    def _getSubcases(self, name):
        """ Get Subcases
        Assumming that subcases are numbers, returns a set of subcases using a 
        test-case's output files. """
        
        #Get cwd and go to the old test-case's output folder.
        return_path = os.getcwd()
        os.chdir(self._old_path+name+'/output/')  

        #Find the names of every comparison file.
        types = ['*.avs', '*.csv', '*.his']
        file_names = []
        for t in types:
            file_names = file_names+glob.glob(t)
         
        #Using the comparison files, extract the subcase numbers.    
        subcases = []
        for file_name in file_names:
            pattern = re.compile(r'\d+')
            subcase = pattern.findall(file_name)[0]
            if subcase not in subcases:
                subcases.append(subcase)
                
        os.chdir(return_path)        
        return subcases
                                  
    def _primeControlFile(self, name, subcase):
        """ Prime Control File
        Makes the control file ready for FEHM execution by renaming it to 
        fehmn.files and removing 'output/' occurrences. """
        
        #Save the cwd and go to the new test-case.
        return_path = os.getcwd()
        os.chdir(self._new_path+name)
        
        #Copy the control file to the subcase folder.
        source = glob.glob('input/*.files')[0]
        dest = 'subcase-'+subcase+'/fehmn.files'
        shutil.copyfile(source, dest)
           
        #Open the file for editing.
        opened_file = open('subcase-'+subcase+'/fehmn.files', 'r+')
        
        #Remove all 'output/' from the control file.
        data = opened_file.read()
        data = re.sub('output/', '', data)
        
        #Subsitute N or NUM with the subcase number.
        data = re.sub('N', subcase, data)
        data = re.sub('UM', '', data)
        
        #Write the changes to the file.
        opened_file.seek(0)
        opened_file.write(data)
        opened_file.truncate()
        opened_file.close()
        
        os.chdir(return_path)
        
    def _copyCompareFiles(self, name):
        """ Copy Compare Files
        For test-case, 'name', copies all csv, avs, and history files from the 
        old test-suite to the new test-suite. """
        
        #Change directory to the old test-case's output folder.
        os.chdir(self._old_path+name+'/output')
        
        #Read all files in this list of types.
        types = ['*.avs', '*.csv' '*.his']
        compare_files = {}
        for t in types:
            compare_files.update(self._readToDictionary(t))
        
        #Write the files to the compare folder in the new test-case.
        os.chdir(self._new_path+name+'/compare')
        self._writeFromDictionary(compare_files)
        os.chdir('..')
                 
    def _readToDictionary(self, pattern):
        """ Read to Dictionary
        Given a glob file pattern, reads files into a dictionary. """
    
        #Determine the input file names.
        file_names = glob.glob(pattern)
        
        #Store the contents of each input file as a dictionary.
        files = {}
        for f in file_names:
            opened_file = open(f)
            files[f] = opened_file.read()
            opened_file.close()
            
        return files
        
    def _createTestFolder(self, name, subcases):
        """ Create Test Folder
        Creates a folder in the new test-suite for 'name' with a folder for 
        input, old comparison files, and a folder for each new subcase run. """
    
        os.makedirs(name)
        os.chdir(name)
        os.makedirs('input')
        os.makedirs('compare')
        for subcase in subcases:
            os.makedirs('subcase-'+subcase)
        
    def _writeFromDictionary(self, files):
        for key in files:  
            opened_file = open(key, 'w')
            opened_file.write(files[key])
            opened_file.close()
                                                 
if __name__ == '__main__':
    #Check that a directory name was specified.
    if len(sys.argv) == 3:
        name = sys.argv[1]
        new_path = sys.argv[2]
    else:
        print "Usage: python testcaseconvertor.py <test-name> <old-test-location>"
        os._exit(0)
     
    #Converts 'name' from old_test_suite to test case in fehmpytests.py. 
    convertor = TestCaseConvertor(new_path)   
    convertor.convertCase(name)
    
