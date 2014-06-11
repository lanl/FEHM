import os
import sys
import re
import glob
import distutils.core

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
     
        #Navigate to fehmpytest.
        os.chdir(self._new_path)
        
        #Create folder for test case, unless it already exists.
        try:
            self._createTestFolder(name)
        except:
            print 'ERROR: There is already a folder for this test case.'
            os._exit(0)
            
        #Copy the input files from old to new test suite.
        source = self._old_path+name+'/input'
        dest = self._new_path+name+'/input'
        distutils.dir_util.copy_tree(source, dest)  
        
        self._primeControlFile(name) 
             
        self._copyCompareFiles(name)
         
        print 'Done.'
                                  
    def _primeControlFile(self, name):
        """ Prime Control File
        Makes the control file ready for FEHM execution by renaming it to 
        fehmn.files and removing 'output/' occurrences. """
        
        #Change the control file's name.
        os.chdir(self._new_path+name+'/input')
        control_file = glob.glob('*.files')[0]
        os.rename(control_file, 'generic_fehmn.files')
        
        #Remove all 'output/' from the control file.
        opened_file = open('generic_fehmn.files', 'r+')
        data = opened_file.read()
        data = re.sub('output/', '', data)
        opened_file.seek(0)
        opened_file.write(data)
        opened_file.truncate()
        opened_file.close()
        
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
        
    def _createTestFolder(self, name):
        os.makedirs(name)
        os.chdir(name)
        os.makedirs('input')
        os.makedirs('compare')
        
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
        print "Usage: python old2new.py '<test-name>' '<old-test-location>'"
        os._exit(0)
     
    #Converts 'name' from old_test_suite to test case in fehmpytests.py. 
    convertor = TestCaseConvertor(new_path)   
    convertor.convertCase(name)
    
