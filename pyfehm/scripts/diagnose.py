
import os
from optparse import OptionParser

usage = 'Usage: diagnose.py FEHM_executable'
parser = OptionParser(usage=usage)

(options, args) = parser.parse_args()

from fdata import*

dat = fdata('fehmn.files')
dat.hist.variables=list(flatten(dat.hist.variables))
dat._diagnostic.refresh_nodes()
p = Popen(args[0],stdout=PIPE)
dat._diagnostic.stdout = p.stdout
dat._diagnostic.poll = p.poll
dat._diagnostic.construct_viewer() 
for line in iter(p.stdout.readline, b''): print line.rstrip() 