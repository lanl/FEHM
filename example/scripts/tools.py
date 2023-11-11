'''
Collection of useful tools for general use.
'''
import os,sys
import numpy as np
import tempfile
from math import floor,log10
import re
import pandas as pd

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific notation of the
    given number formatted for use with LaTeX or Mathtext, with
    specified number of significant decimal digits and precision (number
    of decimal digits to show). The exponent to be used can also be
    specified explicitly.
    (https://stackoverflow.com/questions/18311909/how-do-i-annotate-
    with-power-of-ten-formatting)

    Examples
    --------
    - sci_notation(3.45642e-12,1)
        Scientific notation using 10^() as separator.
    - sci_notation(3.45642e-12,1,2,exponent=-14)
        1 significant decimal digit and 2 decimal digits shown as well as a
        given exponent.
    """
    if exponent is None:
        exponent = int(floor(log10(abs(num))))
    coeff = round(num / float(10**exponent), decimal_digits)
    if precision is None:
        precision = decimal_digits
    #  return r"${0:.{2}f}\cdot10^{{{1:d}}}$".format(coeff, exponent, precision) # using the \cdot 
    return r"${0:.{2}f}\times10^{{{1:d}}}$".format(coeff, exponent, precision)

def sine_signal(period, ts, y_0, ampl):
    '''
    Create a sinusoidal wave by specifying the period, time series (in
    desired time units), y starting place, and amplitude of the disturbance.
    Useful for generating synthetic sine wave for FEHM *.boun files.
    '''
    w = 2.*np.pi/period             #angular freq
    ys = np.sin(w*ts) * ampl + y_0  #y values
    return ys

def cosine_signal(period, ts, y_0, ampl):
    '''
    Create a sinusoidal wave (COSINE) by specifying the period, time series (in
    desired time units), y starting place, and amplitude of the disturbance.
    Useful for generating synthetic cosine wave for FEHM *.boun files.
    '''
    w = 2.*np.pi/period             #angular freq
    ys = np.cos(w*ts) * ampl + y_0  #y values
    return ys

def symlink(target, link_name, overwrite=False):
    '''
    Create a symbolic link named link_name pointing to target.
    If link_name exists then FileExistsError is raised, unless overwrite=True.
    When trying to overwrite a directory, IsADirectoryError is raised.
    '''
    if not overwrite:
        os.symlink(target, link_name)
        return
    # os.replace() may fail if files are on different filesystems
    link_dir = os.path.dirname(link_name)
    # Create link to target with temporary filename
    while True:
        temp_link_name = tempfile.mktemp(dir=link_dir)
        # os.* functions mimic as closely as possible system functions
        # The POSIX symlink() returns EEXIST if link_name already exists
        # https://pubs.opengroup.org/onlinepubs/9699919799/functions/symlink.html
        try:
            os.symlink(target, temp_link_name)
            break
        except FileExistsError:
            pass
    # Replace link_name with temp_link_name
    try:
        # Pre-empt os.replace on a directory with a nicer message
        if not os.path.islink(link_name) and os.path.isdir(link_name):
            raise IsADirectoryError(f"Cannot symlink over existing directory: '{link_name}'")
        os.replace(temp_link_name, link_name)
    except:
        if os.path.islink(temp_link_name):
            os.remove(temp_link_name)
        raise


def find_first_following(filename,string1,string2,initial_final='initial'):
    '''
    Helper function for parsing data from FEHM run.out files.
    Gets the value following string1 and string2.

    Parameters
    ----------

    string1 : str
        String usually specifying what the species name is. e.g.,
        (' Solute output information, species number methane_l').
    string2 : str
        String specifying the quantity to grab, e.g.:
        '   current mass ='
    '''
    i=-1;ln=[]
    with open(filename,'r') as f:
        lines = f.readlines()
        for l in lines:
            i+=1
            #Find first line matching string1
            if re.search(string1,l):
                ln.append(i)
                if initial_final=='initial': break
    #Find the next line that matches string2
    #  for l in lines[i:]:
        if initial_final == 'initial': i = ln[0]
        else: i = ln[-1]
    for l in lines[i:]:
        imass_search = string2
        if re.search(imass_search,l):
            #  #Remove spaces, \n, and units ('mol') from string
            #  val.append(float(l.strip(imass_search).strip('\n').strip(' mol')))
            #Remove spaces, \n, =, and units ('mol') from string
            val = float(l.strip(imass_search).strip('\n').strip('=').strip(' mol'))
            break
    return val

def find_all_first_following(filename,string1,string2):
    i=0;ln=[];val=[];t=[]
    with open(filename,'r') as f:
        lines = f.readlines()
        for l in lines:
            i+=1
            #Find first line matching string1
            if re.search(string1,l):
                ln.append(i)
                #  if initial_final=='initial': break
    #Find the next line that matches string2
    #  for l in lines[i:]:
        #  if initial_final == 'initial': i = ln[0]
        #  else: i = ln[-1]
    #  for l in lines[i:]:
    #  print(ln)
    for i in ln:
        for l in lines[i:]:
            imass_search = string2
            if re.search(imass_search,l):
                #  #Remove spaces, \n, and units ('mol') from string
                #  val.append(float(l.strip(imass_search).strip('\n').strip(' mol')))
                #Remove spaces, \n, =, and units ('mol') from string
                val.append(float(l.strip(imass_search).strip('\n').strip('=').strip(' mol')))
                break
    #Get Timestamps
    for l in lines:
        solute_time_search = 'Solute information at time ='
        if re.search(solute_time_search,l):
            t.append(float(l.strip(solute_time_search).strip('\n').strip('days')))
    # Convert to two-column ndarray 
    data_array = np.column_stack( (t, val) )
    #  return [t,val]
    return data_array


def parse_outfile(filename, species_name, quantity='current mass'):
    '''
    Parse FEHM run.out file. Can be used to get time series data of each named
    solute tracer species in the model at each time step.
    Output is a 2-column Pandas DataFrame where the first column is time
    [days], and second column is whatever the specified quantity is.

    Parameters
    ----------

    species_name : str
        Tracer species name as displayed in *.out file (e.g., 'xenon_s', 'xenon_l', 'xenon_v', ...).
    quantity : str
        String to search for (will grab the value after the search term).
        For now, only works for solute tracer masses.
        Options: ['current mass']
    '''
    col1 = 'time'
    # Search term 1
    s1 = ' Solute output information, species number {}'.format(species_name)
    # Search term 2 (quantity following search term 1)
    if quantity=='current mass':
        col2 = 'moles'
        # Make sure correct number of spaces
        if '   ' not in quantity: s2 = '   {}'.format(quantity)
        else: s2=quantity
        #---- Get initial mass too (t=0 not listed under current mass)
        s2_initial = 'initial mass'
        initial = find_first_following(filename,s1,s2_initial,'initial')
    # [[ Can add more quantity terms here later as elif statements ]]
    else:
        try: s2
        except UnboundLocalError: raise NameError('quantity argument not recognized. Check spelling.')
    data_arr = find_all_first_following(filename, s1, s2)
    # Add in initial t=0 values
    data_arr = np.insert( data_arr, 0, np.column_stack( (0., initial) ), 0 )
    # Convert ndarray to Pandas DataFrame
    data_frame = pd.DataFrame(data_arr, columns=[col1,col2])
    return data_frame
