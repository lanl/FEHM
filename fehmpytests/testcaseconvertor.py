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

import os
import sys
import re
import glob
import distutils.core
import shutil

class TestCaseConvertor():
    """ 
    Test Case Convertor
     
    Converts old test-cases into new test-cases.
    
    Authors: Mark Lange
    Updated: June 2014 by Mark Lange
    """
    
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
        if len(subcases) > 1:    
            for subcase in subcases:
                self._primeControlFile(name, subcase)
        else:
            self._primeControlFile(name, 'fehmn')
  
        #Copy the comparison files over from the old test-case.
        self._copyCompareFiles(name)
       
        print 'Done.'
        
    def _getSubcases(self, name):
        """ 
        Get Subcases
        
        Assumming that subcases are numbers, returns a set of subcases using a 
        test-case's output files. 
        """
        
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
            #Everything up to the first dot should be the unique subcase.
            subcase = file_name.split('.')[0]
            if subcase not in subcases:
                subcases.append(subcase)
                              
        os.chdir(return_path)        
        return subcases
                                  
    def _primeControlFile(self, name, subcase):
        """ 
        Prime Control File
        
        Makes the control file ready for FEHM execution by renaming it to 
        fehmn.files and removing 'output/' occurrences. 
        """
        
        #Save the cwd and go to the new test-case.
        return_path = os.getcwd()
        os.chdir(self._new_path+name+'/input')
        
        #Copy the control file to the control folder.
        source = glob.glob('*.files')[0]
        dest = 'control/'+subcase+'.files'
        shutil.copyfile(source, dest)
        
        #For a generic control file, it must be changed for the subcase.
        if subcase != 'fehmn':
            f = open(dest, 'r+')
            data = f.read()
            data = data.split('\n')
            
            edited_file = []
            
            for line in data:
                var_positions = ['N', 'base']
                for vp in var_positions:
                    #Search for the placeholder range.
                    start = line.find('/')
                    end = line.find(vp)
                    
                    #The end will be negative if it was not found.
                    if end > 0:                
                        #Replace the place holder with the subcase.
                        placeholder = line[start+1 : end+len(vp)]
                        line = re.sub(placeholder, subcase, line)
                        
                #Some files use NUM instead of N.
                line = re.sub('UM', '', line)
                
                #No output folder. Remove '/output/' occurences.
                line = re.sub('output/', '', line)
                
                #Append the line to edit list.
                edited_file.append(line)
            
            #Write the changes to the file.
            edited_file = '\n'.join(edited_file)    
            f.seek(0)
            f.write(edited_file)
            f.truncate()
            f.close()
  
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
        """ 
        Create Test Folder
        
        Creates a test folder for name with folder input and compare inside.
        Inside the input folder, creates a control folder. 
        """
    
        os.makedirs(name)
        os.chdir(name)
        os.makedirs('input')
        os.makedirs('input/control')
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
        print 'Usage: python testcaseconvertor.py <test-name> '+ \
              '<old-test-location>'
              
        os._exit(0)
     
    #Converts 'name' from old_test_suite to test case in fehmpytests.py. 
    convertor = TestCaseConvertor(new_path)   
    convertor.convertCase(name)
    
