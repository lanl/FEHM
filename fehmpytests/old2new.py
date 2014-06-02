#IMPORTANT - This script is hardcoded to work with /home/mlange.

import os,sys
import re
from glob import glob

def convert_test(name):
    """ Convert Test
    Configured to be run in /home/mlange - Converts old test case to new by
    creating a folder in fehmpytests and copying files over.
    Developed by mlange806@gmail.com on June 2, 2014."""

    #Change directory to specified test case.
    os.chdir('old_test_suite/'+name+'/input')
    
    #Find the name of the file list.
    file_name = glob('*.files')[0]
    
    #Read in the contents of the file list.
    opened_file = open(file_name)
    file_list_whole = opened_file.read()
    opened_file.close()
    
    #Store the file list as a dictionary and extras as a list.
    file_list = {}
    with open(file_name) as opened_file:
        for line in opened_file:
            try:
                #Add line to file_list if it can be split into (key, value).
                (key, value) = line.split(' ')
                
                #Tidy up the string.
                value = re.sub('N.', '*.', value)
                value = re.sub('\n', '', value)
                
                #Set the entry.
                file_list[key] = value
            except:
                #Otherwise ignore.
                pass
    
    os.chdir('..')
                
    #Determine the input file names.
    file_names = glob(file_list['input:'])
    
    #Store the contents of each input file as a dictionary.
    input_files = {}
    for f in file_names:
        opened_file = open(f)
        input_files[f] = opened_file.read()
        opened_file.close()
    
    #Determine the grid file name.
    file_name = file_list['gridf:']
    
    #Read in the contents of the grid file.
    opened_file = open(file_name)
    grid_file = opened_file.read()
    opened_file.close()
    
    #Navigate to fehmpytest.
    path = '/home/mlange/sft-files/fehm-open/fehmpytests'
    os.chdir(path)
    
    #Create folder for test case.
    os.makedirs(name)
    os.chdir(name)
    os.makedirs('input')
    os.makedirs('output')
    
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
    opened_file = open(file_list['gridf:'], 'w')
    opened_file.write(grid_file)
    opened_file.close()
    
if __name__ == '__main__':
    #Check that a directory name was specified.
    if len(sys.argv) == 2:
        name = sys.argv[1]
    else:
        print "Usage: python old2new.py '<folder-name>'"
        os._exit(0)
     
    #Converts 'name' from old_test_suite to test case in fehmpytests.py.   
    convert_test(name)
    
