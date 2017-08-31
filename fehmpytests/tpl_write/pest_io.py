''' Utilities to handle reading and writing PEST files '''
from asteval import Interpreter
from glob import glob
import re
from numpy import recarray, array
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

def tpl_write( pardict, f, outflnm ):
    ''' Write model input file using PEST template file

        :param pardict: Dictionary of parameter values
        :type pardict: dict
        :param f: File handle or file name of PEST template file
        :type f: str or file handle
        :param outflnm: Name of model input file to be written
        :type outflnm: str
    '''
    # Check if f is a string or file and read in lines
    if isinstance( f, file ): 
        t = f.read()
        fnm = f.name
        f.close()
    elif isinstance( f, str ): 
        fnm = f
        with open( f, 'r') as fh:
            t = fh.read()
            fh.close()

    # Make sure file is PEST TPL file
    lh = t.split('\n')[0]
    k = lh.split()
    if k[0] != 'ptf':
        print fnm+" does not appear to be a PEST template file"
        return
    else: 
        tok = k[1] # Collect parameter identifier character
        t = re.sub( lh+'\n', '', t)

    # Complile regex pattern
    p = re.compile(tok+'[^'+tok+']*'+tok)
    # Find all parameter identifiers in file
    ms = p.findall(t)
    pd = {} # Pattern dictionary
    # Create evaluation interpreter
    aeval = Interpreter()
    for k,v in pardict.items():
        aeval.symtable[k] = v
    # Evaluate all unique expressions
    for m in ms:
        pstr = m.split(tok)[1].strip()
        if pstr not in pd:
            pd[m] = aeval(pstr)
    # Perform substitutions
    for k,v in pd.items():
        t = re.sub( re.escape(k), str(v), t)

    # Write output file
    fout = open( outflnm, 'w' )
    fout.write(t)
    fout.close()

def read_par_files( *files ):
    ''' Read in one or more PEST parameter files

        :param files: Strings of file names, glob characters are supported
        :type files: str
        :param output_format: Format to return in
        :type output_format: str ('numpy_array' or 'recarray')
        :returns: parameter names list and numpy array of parameters or recarray of parameters
    '''

    output_format = 'numpy_array' 

    ks = OrderedDict()
    ks_save = OrderedDict()
    pars = [] # List to collect parameter sets
    first = True # Flag to catch first file
    names = [] # List to collect parameter names
    # Loop through function arguments
    for fstr in files:
        fnms = glob( fstr )
        # Loop through files from glob of each function argument
        for fnm in fnms:
            with open(fnm, 'r') as f:
                l = f.readline() # Skip header
                # Store parameters in ks dictionary
                for l in f:
                    vs = l.split()[0:2]
                    ks[vs[0]] = float(vs[1])
                # If first file, just collect names and values
                if first:
                    ks_save = ks
                    names = ks.keys()
                    pars.append(ks.values())
                    first = False
                # Else, check that keys match and put values in correct order according to first file
                else:
                    if not ks_save.keys() == ks.keys():
                        print "Parameters in "+fnm+" differs from other files"
                        return 0
                    else:
                        ptemp = []
                        for v in ks_save.keys():
                            ptemp.append(ks[v])
                        pars.append(ptemp)

    # Create record array
    pars = array(pars)
    if output_format == 'numpy_array':
        return names, pars
    elif output_format == 'recarray':
        pars_rc = pars.view(dtype=zip(names,['float64']*len(names))).copy()
        return pars_rc



                       

                

                    

