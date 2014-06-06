import os,sys
import re
import glob
import distutils.core

class Convertor():
    """ Converter 
    Converts old test-cases into new tet-cases.
    Modified by mlange806@gmail.com on June 6, 2014."""
    
    def __init__(self, new_path):
        """ Default Constructor 
        Sets the default path for the old test suite. """
        
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
        
        #Prime the control file in the new test-case.
        self._primeControlFile(name) 
             
        #Get the compare files.
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
        old test-suite to the new test-suite. If the format is avs, they are 
        converted into csv and then copied. """
        
        #Change directory to the old test-case's output folder.
        os.chdir(self._old_path+name+'/output')
        
        #Read all files in this list of types.
        types = ['*.csv', '*.avs', '*.his']
        compare_files = {}
        for t in types:
            #If the file type is avs, convert to csv.
            if t == '*.avs':
                old_format = self._readFiles(t)
                new_format = {}
                for key in old_format:
                    new_contents = self._convertFormat(key)
                    new_key = re.sub('.avs', '.csv', key)
                    new_format[new_key] = new_contents     
                compare_files.update(new_format)          
            #Else, just read the files as they are.
            else:
                compare_files.update(self._readFiles(t))
        
        #Write the files to the compare folder in the new test-case.
        os.chdir(self._new_path+name+'/compare')
        self._writeFiles(compare_files)
        os.chdir('..')
                 
    def _readFiles(self, pattern):
        """ Read Files
        Given a glob pattern, reads files into a dictionary. """
    
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
        """ Create Test Folder
        Creates a folder called 'name' with subfolders for the files. """
        
        os.makedirs(name)
        os.chdir(name)
        os.makedirs('input')
        os.makedirs('compare')
        
    def _writeFiles(self, files):
        """ Write Files
        Given a dictionary of files, writes each file. """
        
        for key in files:  
            opened_file = open(key, 'w')
            opened_file.write(files[key])
            opened_file.close()
    
    def _convertFormat(self, filename):
        """ AVS to CSV Converter
        Takes in an AVS filename and converts it to CSV.
        Developed by mlange806@gmail.com May 29, 2014 """
        
        #Store the attributes and convert each line.
        data = []
        header = ''     
        with open(filename) as fp:
            #Read column numbers.
            first_line = fp.readline()
            first_line = first_line.split()
            column_num = int(first_line[0])
            
            #Read the column names.
            names = []
            i = 0
            while i < column_num:
                names.append(fp.readline())
                
                #Remove the extra dimension.
                names[i] = names[i].split(',')
                names[i] = names[i][0]
                
                i = i + 1
            
            #Create the header.        
            header = ', '.join(names)
                
            #For each row, change the spaces to commas.
            for line in fp:
                data.append(', '.join(line.split()))
    
        #Write the changes to the original file.
        new_format = header+'\n'  
        for line in data:
            new_format = new_format + '\n'+line
        
        return new_format
              
class ControlFile():
    """ Control File
    Stores file contents and the name of files by category. File contents and 
    the name of files in each category can get retrieved by get functions. """
    
    def __init__(self, data, files):
        #Initialize with file contents and dictionary of file information.
        self._data = data
        self._files = files
        
    def getData(self):
        #Returns the entire file contents.
        return self._data
        
    def getFiles(self, file_type):
        """ Get Files
        Returns glob pattern for category, i.e. 'input' or 'grid'. """
        
        #Work-around for different grid names in control file
        if file_type == 'grid':
            known_grids = ['grid', 'gridf']
            for grid in known_grids:
                if grid in self._files:
                    return self._files[grid]
                    
        #For now, there are not any other exceptions known.
        else:
            return self._files[file_type]               
                      
if __name__ == '__main__':
    #Check that a directory name was specified.
    if len(sys.argv) == 3:
        name = sys.argv[1]
        new_path = sys.argv[2]
    else:
        print "Usage: python old2new.py '<test-name>' '<old-test-location>'"
        os._exit(0)
     
    #Converts 'name' from old_test_suite to test case in fehmpytests.py. 
    convertor = Convertor(new_path)   
    convertor.convertCase(name)
    
