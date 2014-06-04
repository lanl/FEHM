import os,sys
import re
from glob import glob

def convert_test(name, new_path):
    """ Convert Test
    
    Automates part of old to new test-case conversion for a given test-case,
    'name', by creating a folder in the new test location, specified by 
    'new_path', and copying the control, input, grid, and comparison
    files. Code generation will be added soon.
    
    *Test-case must be in a single folder in the old suite.
    
    Modified by mlange806@gmail.com on June 4, 2014. """

    #Change directory to specified test case.
    old_path = '/n/swdev/FEHM_dev/VERIFICATION_linux/VERIFICATION_V3.2linux/'
    os.chdir(old_path+name+'/input')
    
    #Find the name of the file list.
    file_name = glob('*.files')[0]
    
    #Read in the contents of the file list.
    opened_file = open(file_name)
    file_list_whole = opened_file.read()
    opened_file.close()
    
    #Store the file list as a dictionary.
    file_list = {}
    with open(file_name) as opened_file:
        for line in opened_file:
            try:
                #Add line to file_list if it can be split into kvpairs.
                (key, value) = line.split(': ')
                
                #Tidy up the string.
                value = re.sub('N.', '*.', value)
                value = re.sub('\n', '', value)
                
                file_list[key] = value
            except:
                #If the line can not be split into kvpairs, ignore.
                pass
    
    os.chdir('..')
                
    #Determine the input file names.
    file_names = glob(file_list['input'])
    
    #Store the contents of each input file as a dictionary.
    input_files = {}
    for f in file_names:
        opened_file = open(f)
        input_files[f] = opened_file.read()
        opened_file.close()
    
    #Determine the grid file name.
    file_name = file_list['gridf']
    
    #Read in the contents of the grid file.
    opened_file = open(file_name)
    grid_file = opened_file.read()
    opened_file.close()
    
    #Get the names of the compare files.
    os.chdir('output')
    types = ['*.avs', '*.his']
    file_names = []
    for t in types:    
        file_names = file_names+glob(t)
         
    #Store the contents of each comparision file as a dictionary.
    comparison_files = {}
    for f in file_names:
        opened_file = open(f)
        comparison_files[f] = opened_file.read()
        opened_file.close()
    
    #Navigate to fehmpytest.
    os.chdir(new_path)
    
    #Create folder for test case, unless it already exists.
    try:
        os.makedirs(name)
        os.chdir(name)
        os.makedirs('input')
        os.makedirs('output')
    except:
        print 'ERROR: There is already a folder for this test case.'
        os._exit(0)
    
    #Write the file list from the old test suite.
    opened_file = open('fehmn.files', 'w')
    opened_file.write(file_list_whole)
    opened_file.close()
    
    #Write the input files from the old test suite.
    for key in input_files:
        opened_file = open(key, 'w')
        opened_file.write(input_files[key])
        opened_file.close()
    
    #Write the grid file from the old test suite.
    opened_file = open(file_list['gridf'], 'w')
    opened_file.write(grid_file)
    opened_file.close()
    
    #Write the comparison files from the old test suite.
    for key in comparison_files:
        opened_file = open('output/'+key, 'w')
        opened_file.write(comparison_files[key])
        opened_file.close()
        
    print 'Done.'
    
if __name__ == '__main__':
    #Check that a directory name was specified.
    if len(sys.argv) == 3:
        name = sys.argv[1]
        new_path = sys.argv[2]
    else:
        print "Usage: python old2new.py '<test-name>' '<test-location>'"
        os._exit(0)
     
    #Converts 'name' from old_test_suite to test case in fehmpytests.py.   
    convert_test(name, new_path)
    
