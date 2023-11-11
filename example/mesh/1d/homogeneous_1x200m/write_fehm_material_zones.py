'''
Get node numbers of each cell from mesh connectivity, then assign to the
correct zone number based on ``material_ids.txt``.
'''
import os,sys
from os.path import join
sys.path.append('/project/gas_seepage/jportiz/scripts')
import utils
import numpy as np
import numpy
from collections import Counter

fehm_mesh_file = 'grid.inp'


# Write the Top and Bottom zones to .zone files 
mesh = utils.grid_recarray(fehm_mesh_file)
topcoord = max(mesh['y'])
botcoord = min(mesh['y'])
# Find which nodes are on boundaries
topnodes=[]
botnodes=[]
for i in mesh:
    if i['y'] == topcoord: topnodes.append(i['node'])
    elif i['y'] == botcoord: botnodes.append(i['node'])
utils.write_zone('top.zone', topnodes,[100])
utils.write_zone('bot.zone', botnodes,[200])


