'''
Mars methane barometric pumping flow and transport simulations with
temperature-dependent (or constant) adsorption in the shallow regolith.

Uses subsurface temperatures from 1-D thermal simulations to assign temperature
boundary conditions in the subsurface where adsorption occurs (because the
thermal simulations require much more refined mesh than the flow & transport
simulations).
'''
import os,sys
from os.path import join
#  sys.path.append('/project/gas_seepage/jportiz/scripts/')
sys.path.append(join(os.getcwd(),'scripts/'))
from gasmigra import utils
import mars_tools
from mars_tools import sols2days, days2sols
from tools import symlink, find_nearest
from mesh_tools import mesh_files
from glob import glob
import baro_decomp.convert_boun as convert_boun
import pandas as pd
#import ipdb
import pickle
import numpy as np
import numpy
from matk import pest_io
import shutil
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import math
#  import logger
#-----------------------------------------------------
#MATPLOTLIBRC PLOTTING PARAMETERS
# Load up sansmath so that math --> helvetica font
# Also need to tell tex to turn on sansmath package
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{sansmath}',
    r'\sansmath']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
plt.rcParams['axes.labelweight']=u'normal'
#-----------------------------------------------------

try:
    RUNDIR = os.environ["RUNDIR"]
except KeyError:
    print("Set the RUNDIR evironment variable first, e.g.:")
    print("    export RUNDIR='<FEHM-MARS_REPO_DIRECTORY>/example/' ")
try:
    FEHM_MARS_SRC = os.environ["FEHM_MARS_SRC"]
except KeyError:
    print("Set the FEHM_MARS_SRC evironment variable first, e.g.:")
    print("    export FEHM_MARS_SRC='<FEHM-MARS_REPO_DIRECTORY>/src/' ")

#  main_mesh_dir = '/project/gas_seepage/jportiz/mars/mesh/2d_fracture_network/'
#  tpl_dir = '/project/gas_seepage/jportiz/mars/runs/tpls/thermalAdsorption_tpls'
#  exe = '/home/jportiz/software/FEHM/src/xfehm_v3.3.1_mars-debug'
main_mesh_dir = RUNDIR+'/mesh/2d_fracture_network/'
tpl_dir = RUNDIR+'/tpls/thermalAdsorption_tpls'
exe = FEHM_MARS_SRC+'/xfehm_v3.3.1_mars'

#DIRECTORY OF POTENTIAL MESH FILES 
mesh_dict = {
            # 50x200m, fracDen=0.000% (@b=1mm)
            'depth200_fracDen000': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen000/top.zone',
                'boun_file':'mars_pressTemp.boun',
                #  'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            # 50x200m, fracDen=0.035% (@b=1mm)
            'depth200_fracDen035': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            # 50x200m, fracDen=0.020% (@b=1mm)
            'depth200_fracDen020': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen020/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            # 50x200m, fracDen=0.010% (@b=1mm)
            'depth200_fracDen010': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen010/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            # 50x200m, fracDen=0.005% (@b=1mm)
            'depth200_fracDen005': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen005/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            # 50x200m, fracDen=0.0035% (@b=1mm)
            'depth200_fracDen003': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen003/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            # 50x200m, fracDen=0.001% (@b=1mm)
            'depth200_fracDen001': {
                'mesh_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x200m_fracDen001/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': RUNDIR+'/1d_heat_runs/heat_200m'},
            #  # 50x500m, fracDen=0.000% (@b=1mm)
            #  'depth500_fracDen000': {
                #  'mesh_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/grid.inp',
                #  'stor_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/grid.stor',
                #  'avs_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/fracture_2D.inp', #not used
                #  'matzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/fracture_material.zone',
                #  'matzonn_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/fracture_material.zonn',
                #  'outsidezone_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/fracture_2D_outside.zone', #not used
                #  'topzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen000/top.zone',
                #  'boun_file':'mars_pressTemp.boun',
                #  'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            # 50x500m, fracDen=0.035% (@b=1mm)
            'depth500_fracDen035': {
                'mesh_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen035/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            # 50x500m, fracDen=0.010% (@b=1mm)
            'depth500_fracDen010': {
                'mesh_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen010/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            # 50x500m, fracDen=0.001% (@b=1mm)
            'depth500_fracDen001': {
                'mesh_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x500m_fracDen001/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            # 50x1000, fracDen=0.035% (@b=1mm)
            'depth1000_fracDen035': {
                'mesh_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x1000_fracDen035/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            # 50x1000, fracDen=0.010% (@b=1mm)
            'depth1000_fracDen010': {
                'mesh_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x1000_fracDen010/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            # 50x1000, fracDen=0.001% (@b=1mm)
            'depth1000_fracDen001': {
                'mesh_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/grid.inp',
                'stor_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/grid.stor',
                'avs_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/fracture_2D.inp', #not used
                'matzone_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/fracture_material.zone',
                'matzonn_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/fracture_material.zonn',
                'outsidezone_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/fracture_2D_outside.zone', #not used
                'topzone_file':main_mesh_dir+'levyfractures_50x1000_fracDen001/top.zone',
                'boun_file':'mars_pressTemp.boun',
                'heat_dir': '/lclscratch/jportiz/projects/gas_seepage/mars/runs/1d_heat_runs/heat_200m'},
            }

def parse_hist(results_dir, param='temp', mesh_file=None):
    '''
    Parse FEHM history files, with node names, depths, etc. Assumes that the
    output nodes are in a vertical line on the left side of domain.
    '''
    cwd = os.getcwd()
    #  results_dir = join(cwd,'geo_gradient')
    results_dir = join(results_dir)
    #-------------------------------------------------- 
    #           GET MESH DATA
    #-------------------------------------------------- 
    # If mesh_file not provided, find automatically from results_dir by looking
    # at the fehmn.files file
    #   - Otherwise, you can specify a static fehmn.files to use w/ mesh_file 
    if mesh_file is None: mesh_path = join(results_dir,'fehmn.files')
    else:                 mesh_path = mesh_file
    #  with open(join(results_dir,'fehmn.files'),'r') as ff:
    with open(mesh_path,'r') as ff:
        i=0
        for line in ff.readlines():
            if line.startswith('grida'):
                i+=1
                os.chdir(results_dir); wd = os.getcwd()
                gf = line.split(": ",1)[1].rstrip()
                os.chdir(os.path.split(gf)[0]); meshdir = os.getcwd()
                meshfile = join(meshdir,os.path.split(gf)[-1])
                os.chdir(cwd)
                if i==1: break
    grid = utils.grid_recarray(meshfile)
    fracnodes = np.where(grid['x']==0.) #for 2D-planar, single fracture (left side)
    fracdepths = grid['y'][fracnodes]   #  ``     ``         ``
    mesh_dir1 = os.path.split(os.path.split(os.path.split(meshfile)[0])[0])[1]
    mesh_dir2 = os.path.split(os.path.split(meshfile)[0])[1]
    mesh_name = mesh_dir1+'/'+mesh_dir2
    mesh_name = mesh_name.replace('_','\_')
    #-------------------------------------------------- 
    #           FEHM  DATA 
    #-------------------------------------------------- 
    if   param == 'temp': hist_name='run_temp.his'
    elif param == 'pres': hist_name='run_presWAT.his' #may have to alter this
    elif param == 'visc': hist_name='run_visWAT.his'
    else                : sys.exit('MUST SPECIFY A VALID PARAMETER')
    tf = os.path.join(results_dir, hist_name) #e.g., 'run_temp.his'
    #---- Load the history file 
    tdata = np.loadtxt(tf,skiprows=4)
    t_cols = np.loadtxt(tf,skiprows=3,max_rows=1,dtype='str')[1::2]
    # Get depths of all nodes
    nodes =[]; ys = []
    for i in t_cols[1:]:
        nodes.append(int(i)-1)
        ys.append(grid[int(i)-1]['y'])
    nodes = np.asarray(nodes)
    ys    = np.asarray(ys)
    return tdata, t_cols, nodes, ys


def node_macro(grid,macro_dir='.',dirxn='vert',spacing='none'):
    '''
    Write ``nodeline.dat`` file for monitoring nodes in FEHM.

    Parameters
    ----------
    grid : str
        Mesh file.
    macro_dir : str
        Path to ``macros`` directory. Creates if does not exist.
    dirxn : str
        -``vert`` specifies a vertical line of nodes on left side of domain
        (used to be every 10th node).
        -``horz`` gives a horizontal line of nodes on the top of domain.
        -``both`` gives both (vertical then horizontal)
    spacing : float or 'none' (str)
        Specify the dy or dx spacing you want for the nodes [m].
        Only works for ``vert`` currently.
        - SPECIAL CASE: spacing:'top' will only have a single node at the top
    '''
    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir)
    fh = open(os.path.join(macro_dir,'nodeline.dat'),'w')
    #---- Vertical line of nodes on x=0
    #  get = grid[numpy.where(grid['x']==0)[0][0:-1:10]]
    get = grid[numpy.where(grid['x']==0)[0]]
    #  get = get[0::10] #<-- if you want only every 10th node
    get.sort(order='y')
    # if top node not included, add it
    if get[-1]['y'] != max(grid['y']):
        get = numpy.append(get, numpy.array(grid[numpy.where(grid['x']==0)[0]][-1]))
    #  get = get[0:-1]['node']

    #-- Spacing
    # Grab subset of nodes with the specified spacing (if used)
    if spacing is not 'none':
        if spacing == 'top':
            #  print(get['node'])
            topnode= max(get['node']) #get only the top node in this special case
            get = [topnode]
        else:
            first = get[0] #get first 
            sinds = [];ind=0
            for ys in enumerate(get):
                #  print(ys)
                if ys[0] == 0:
                    sinds.append(ys[1])
                else:
                    # Calculate spacing between current grid point and last recorded one
                    calcd_spacing = round(abs(sinds[-1]['y'] - ys[1]['y']),2)
                    if calcd_spacing == spacing:
                        sinds.append(ys[1])
            get = [] #get just the node numbers
            for i in sinds: get.append(i['node'])
    #If spacing not specified
    else: get = get['node']

    #---- Horizontal line of nodes on y max (aka top)
    gett = grid[numpy.where(grid['y']==numpy.max(grid['y']))[0]]
    gett.sort(order='x')
    gett = gett['node']
    #fh.write('node\n'+str(len(print_all))+'\n')
    #numpy.savetxt(fh,print_all,fmt=' %d')
    if dirxn=='vert':
        fh.write('node\n'+str(len(get))+'\n')
        numpy.savetxt(fh,get,fmt=' %d')
        fh.write('\nstop\n')
        fh.close()
    elif dirxn=='horz':
        fh.write('node\n'+str(len(gett))+'\n')
        numpy.savetxt(fh,gett,fmt=' %d')
        fh.write('\nstop\n')
        fh.close()
    elif dirxn =='both':
        print_all = numpy.concatenate([gett,get])
        fh.write('node\n'+str(len(print_all))+'\n')
        numpy.savetxt(fh,print_all,fmt=' %d')
        fh.write('\nstop\n')
        fh.close()

#  def node_macro(grid,macro_dir='.',dirxn='vert'):
    #  '''
    #  Write ``nodeline.dat`` file for monitoring nodes in FEHM.
#  
    #  Parameters
    #  ----------
    #  grid : str
        #  Mesh file.
    #  macro_dir : str
        #  Path to ``macros`` directory. Creates if does not exist.
    #  dirxn : str
        #  -``vert`` specifies a vertical line of nodes on left side of domain
        #  (every 10th node).
        #  -``horz`` gives a horizontal line of nodes on the top of domain.
        #  -``both`` gives both (vertical then horizontal)
    #  '''
    #  # Make sure macro_dir exists
    #  if not os.path.isdir(macro_dir): os.mkdir(macro_dir)
    #  fh = open(os.path.join(macro_dir,'nodeline.dat'),'w')
    #  # Vertical line of nodes (x=0)
    #  #  get = grid[numpy.where(grid['x']==0)[0][0:-1:10]]
    #  get = grid[numpy.where(grid['x']==0)[0]]
    #  get = get[0::10]
    #  get.sort(order='y')
    #  # if top node not included, add it
    #  if get[-1]['y'] != max(grid['y']):
        #  get = numpy.append(get, numpy.array(grid[numpy.where(grid['x']==0)[0]][-1]))
    #  #  get = get[0:-1]['node']
    #  get = get['node']
    #  # Horizontal line of nodes (y max) (aka top)
    #  gett = grid[numpy.where(grid['y']==numpy.max(grid['y']))[0]]
    #  gett.sort(order='x')
    #  gett = gett['node']
    #  print_all = numpy.concatenate([gett,get])
    #  #fh.write('node\n'+str(len(print_all))+'\n')
    #  #numpy.savetxt(fh,print_all,fmt=' %d')
    #  if dirxn=='vert':
        #  fh.write('node\n'+str(len(get))+'\n')
        #  numpy.savetxt(fh,get,fmt=' %d')
        #  fh.write('\nstop\n')
        #  fh.close()
    #  elif dirxn=='horz':
        #  fh.write('node\n'+str(len(gett))+'\n')
        #  numpy.savetxt(fh,gett,fmt=' %d')
        #  fh.write('\nstop\n')
        #  fh.close()
    #  elif dirxn =='both':
        #  fh.write('node\n'+str(len(print_all))+'\n')
        #  numpy.savetxt(fh,print_all,fmt=' %d')
        #  fh.write('\nstop\n')
        #  fh.close()

def perm_macro(fname, grid, k_bg, macro_dir='.',zones=[],zone_perms=[]):
    '''
    Write perm macro file for FEHM (no gdkm).

    Parameters
    ----------
    fname : str
        File name (e.g. 'perm.dat')
    grid : str
        Mesh file.
    k_bg : float
        Background (matrix) permeability [m^2].
    macro_dir : str
        Path to ``macros`` directory. Creates if does not exist.
    zones : list or ndarray
        Zone/zonn numbers to which non-background perms are assigned.
    zone_perms: list or ndarray
        Permeability values to assign to non-background [m^2].
    '''
    # Make sure macro_dir exists
    if not os.path.isdir(macro_dir): os.mkdir(macro_dir)
    #  fh = open(os.path.join(macro_dir,'perm.dat'),'w')
    fh = open(os.path.join(macro_dir,fname),'w')
    fh.write('perm\n  1 0 0 '+str(k_bg)+' '+str(k_bg)+' '+str(k_bg)+'\n')
    #  ps = numpy.column_stack([gdkm['node2'],gdkm['node2'],numpy.ones_like(gdkm['node2']),gdkm['k_f'],gdkm['k_f'],gdkm['k_f']])
    #  numpy.savetxt(fh,ps,fmt=' %d %d %d %16.15g %16.15g %16.15g')
    for zn,k in zip(zones,zone_perms):
        if zn < 0:
            fh.write(' {0} 0 0 {1} {1} {1}\n'.format( int(zn), k ))
        else:
            fh.write(' {0} {0} 1 {1} {1} {1}\n'.format( int(zn), k  ))
    fh.write('\nstop\n')
    fh.close()

def SSA_bet():
    '''
    Returns the measured BET specific surface area (SSA) for JSC-Mars-1 soil
    analog.
    Units : [m2/kg]
    Source: Gough2010
    '''
    return 1.00e5

def ML_ch4():
    '''
    Returns the monolayer coverage of CH4 as determined from size of an
    adsorbed CH4 molecule (19.18 Å^2 = 19.18e-20 m^2).
    Units : [molecules/m^2]
    Source: Gough2010
    '''
    return 5.21e18

def molar_mass_ch4():
    '''[kg/mol]'''
    return 0.01604

def v_ms(T,molar_mass):
    '''
    Returns the mean molecular speed of a gas at temperature T (K), given by
    Maxwell-Boltzmann distribution of gas velocities.
    Units : [m/s]
    Source: Gough2010
    '''
    R = 8.314      #[J/mol/K]
    #  # Convert molar mass to molecular weight (M_w) (might be wrong here...)
    #  M_w = molar_mass_ch4()/6.022e23         #[kg/molecule] molecular weight of CH4
    return ( (8.*R*T)/(np.pi*molar_mass) )**0.5
    #  return ( (8.*R*T)/(np.pi*M_w) )**0.5

def k_B():
    '''
    Returns Boltzmann's constant
    Units : [m^2/kg/s^2/K]
    Source: Gough2010
    '''
    return 1.380649e-23

def h():
    '''
    Returns Planck's constant.
    Units : [m^2/kg/s]
    '''
    return 6.62607015e-34

def maxconc():
    '''
    FOR JSC-Mars-1
    --------------
    Returns the MAXCONC parameter (maximum sorbed concentration) for kinetic
    Langmuir adsorption used by rxn macro .
    Units : [moles/kg_rock].
    Source: Gough2010
    '''
    mc = ML_ch4() * SSA_bet() / 6.022e23
    return mc

def K_eq(T, molar_mass, deltaE_obs):
    '''
    Returns the equilibrium constant for a simple, one-step adsorption process.
    Units : [1/concentration] or [kg/mol]
    Source: Gough2010 (Eqn 13)

    Parameters:
    ----------
    T : float
        Temperature [K].
    molar_mass : float
        Molar mass of the gas [kg/mol].
    deltaE_obs : float
        Adsorption energy (experimentally determined) [J/mol]
    '''
    R = 8.314     #[J/mol/K]
    Keq = ( v_ms(T,molar_mass)*h() )/(4.*k_B()*T) * np.exp( deltaE_obs /R/T )
    return Keq

#  # !!! DEBUG !!!
#  # 1st, load p dict from run.py, then:
#  macro_dir = 'macros'
#  init0_dir=None;init_dir=None;spinup_dir=None;m_air_dir=None;conc_init_dir=None; verbose=False
#  
def model(p,init0_dir=None,init_dir=None,spinup_dir=None,m_air_dir=None,conc_init_dir=None,
          exe='xfehm',tpl_dir='/project/gas_seepage/jportiz/mars/runs/tpls',verbose=False,
          macro_dir='macros', gdkm_top_nodes=False):
    '''
    Gas migration model for Mars baro pumping

    :param p: Dictionary of parameters values keyed by parameter names, so options below
    :type p: dict
    :param init0_dir: Directory of pre-initialization run with uniform properties, if None it will be generated
    :type init0_dir: string
    :param init_dir: Directory of initialization run, if None it will be generated
    :type init_dir: string
    :param spinup_dir: Directory of spinup run, if None it will be generated
    :type spinup_dir: string
    :param m_air_dir: Directory of run to determine mass of air in cavity/chimney in kg to determine initial concentration in moles/kg-air
    :type m_air_dir: string
    :param exe: FEHM executable, if None, environment variable FEHM_EXE or simply xfehm will be used
    :type exe: string
    :param tpl_dir: Directory containing FEHM template files
    :type tpl_dir: string
    :param macro_dir: Directory to put automatically generated FEHM macros
    :type macro_dir: string
    :param gdkm_top_nodes: remove top nodes from gdkm macro
    :type gdkm_top_nodes: bool

    Parameter dictionary 'p' options, defaults in ():
        mesh (101): mesh index, see gas_migration.mesh_dict for options
        dfmax (0.0009): maximum fracture aperture [m]
        tracer_start (0.0): Time to start tracer test, determines what point along baro pressure curve to start [d]
    '''

    # Uncomment to debug in ipython
    #ipdb.set_trace()
    # If parameters are missing, set defaults from Jordan et al. (2016) Scientific Reports
    if 'mesh' not in p: p['mesh']=50
    # Fracture aperture
    if 'df' not in p: p['df'] = 0.0001 #[m]
    if 'methane_prod_rate' not in p: p['methane_prod_rate'] = 0. #initial tracer source (no source strength)
    # Tracer start time, used to shift starting time along baro pressure signal
    #if 'tracer_start' not in p: p['tracer_start'] = 0.0
    # Perm cutoff for gdkm nodes
    if 'k_gdkm_cutoff' not in p: p['k_gdkm_cutoff'] = 1.e-12
    if 'sim_time' not in p: p['sim_time'] = mars_tools.sols2days(52000.) #mars_tools.sols2days(2712.41) #sols-->days (extent of the MSL baro record)
    if 'c_initial' not in p: p['c_initial'] = 0.0 #[mol/kg vapor] initial conc. of contaminated zone
    if 'phi_m' not in p: p['phi_m'] = 0.350 #default matrix porosity
    #  if 'phi_f' not in p: p['phi_f'] = 0.001 #default fracture porosity (upscaled if needed) 
    if 'k_m' not in p: p['k_m'] = 1.e-14 #default matrix perm
    if 'T_ref' not in p: p['T_ref'] = -46.93 #[°C] mean surf. temp @ Gale [Klusman2022]
    if 'no_adsorption' not in p: p['no_adsorption'] = False # "True" turns off adsorption rxn
    if 'noSorp_dir' not in p: p['noSorp_dir'] = None         #default name of noSorp dir to restart from
    if 'synthetic' not in p: p['synthetic'] = True           #default to use synthetic P and T

    #-------------------------------------------------- 
    # PRINT A WARNING IF ADSORPTION RUN 
    #-------------------------------------------------- 
    #  if p['no_adsorption'] is True:
    if (p['no_adsorption'] is False) and (p['noSorp_dir'] is not None):
        print('************************ WARNING *******************************')
        print('You are running an adsorption simulation. This will restart from')
        print('the following no-sorption simulation:')
        print('{}/noSorp_tracer'.format(p['noSorp_dir']))
        print()
        print('Please ensure that run is current and has the right properties.')
        print('Do you want to continue? (y/n)')
        proceed = input()
        if proceed == 'y' or proceed =='Y':
            print('Ok, proceeding with adsorption simulation')
        else: exit()
        print('****************************************************************')
        print()

    #-------------------------------------------------- 
    # KEEP TRACK OF TYPE(S) OF REGOLITH ADSORPTION 
    #-------------------------------------------------- 
    # Set defaults
    regolith_type = {}
    regolith_type['active']   = False
    regolith_type['constant'] = False
    # Now change based on regolith depths provided in run.py
    if 'active_regolith_base'   in p: regolith_type['active']   = True
    if 'constant_regolith_base' in p: regolith_type['constant'] = True
    #  # Now change depths to base of 'False' entries to bogus values 
    #  #  - Still want to create zone files
    #  if regolith_type['active']   is False: p['active_regolith_base'] = 9999
    #  if regolith_type['constant'] is False: p['constant_regolith_base'] = 9999
    # Ensure depths are negative  (unless bogus value)
    # ---- Default constant_regolith_base is bottom of domain
    if regolith_type['active'] is False: p['active_regolith_base'] = 9999
    else: p['active_regolith_base']   = -abs(p['active_regolith_base'])
    if regolith_type['constant'] is False: p['constant_regolith_base'] = 9999
    elif isinstance(p['constant_regolith_base'], str):  #if a string
        p['constant_regolith_base']     = -1e20
    else: p['constant_regolith_base']   = -abs(p['constant_regolith_base'])
    # Error checking
    if all(rg is True for rg in regolith_type.values()):
        if abs(p['constant_regolith_base']) <= abs(p['active_regolith_base']):
            print('B')
            print('************************ WARNING *******************************')
            print('Depth of constant_regolith_base must be greater than that of active_regolith_base')
            print('  Terminating simulation.')
            print('****************************************************************')
            exit()
        #
    elif (all(rg is False for rg in regolith_type.values()) and (p['no_adsorption'] is False)):
            print('C')
            print('************************ WARNING *******************************')
            print('You are running an adsorption simulation but have not specified')
            print('active_regolith_base and/or constant_regolith_base depths.')
            print('  Terminating simulation.')
            print('****************************************************************')
            exit()
    print('Type(s) of regolith adsorption in this simulation:')
    print(regolith_type)


    p['mesh_file'] = mesh_dict[p['mesh']]['mesh_file']
    p['stor_file'] = mesh_dict[p['mesh']]['stor_file']
    p['avs_file'] = mesh_dict[p['mesh']]['avs_file']
    p['matzone_file'] = mesh_dict[p['mesh']]['matzone_file']
    p['matzonn_file'] = mesh_dict[p['mesh']]['matzonn_file']
    p['outsidezone_file'] = mesh_dict[p['mesh']]['outsidezone_file']
    p['topzone_file'] = mesh_dict[p['mesh']]['topzone_file']
    p['boun_file'] = mesh_dict[p['mesh']]['boun_file']
    #  if p['no_adsorption'] is True:
        #  tpl_dir = '/project/gas_seepage/jportiz/mars/runs/tpls/'

    #-------------------------------------------------- 
    # THROW A WARNING IF EXE DOES NOT HAVE 'mars' IN IT
    #-------------------------------------------------- 
    if 'mars' not in exe:
        print('************************ WARNING *******************************')
        print('You did not provide an FEHM_mars executable for parameter "exe".')
        proceed = input('Do you want to continue? (y/n)')
        if proceed == 'y' or proceed =='Y':
            print('Ok, proceeding without FEHM_mars executable')
        else: exit()
        print('****************************************************************')


    g = utils.grid_recarray(p['mesh_file'])

    zones = [-9] #9 = fracture(s) zone
    k_f = p['df']**2./12 #[m2]
    #--------------------------------------------- 
    # WRITE perm.dat AND nodeline.dat 
    #--------------------------------------------- 
    #---- Write different perm.dat and nodeline.dat macros depending on mesh (2D-single fracture vs 2D-network) ((or not...)
    if '2d_planar' in p['mesh_file']:
        zone_perms = [k_f]
        node_macro(g,macro_dir=macro_dir,dirxn='vert')
    #  elif 'levy' in p['mesh_file']:
    elif '2d_fracture_network' in p['mesh_file']:
        #calculate upscaled fracture porosity with phi_f = b/dx 
        dx = g[1]['x']-g[0]['x']
        phi_f = p['df']/dx
        p['phi_f'] = phi_f
        #calculate effective permeability of the fractures using k1A1 = k2A2
        zone_perms = [(k_f*p['df']*1)/(1*1)]
        #  node_macro(g,macro_dir=macro_dir,dirxn='horz')
        #  node_macro(g,macro_dir=macro_dir,dirxn='vert')
        # If deep/shallow run, only have the top node in nodeline.dat
        spc = 'top'
        #  if 'deep' in p['noSorp_dir'] or p['no_adsorption']==True: spc = 'top'
        #  elif 'shallow' in p['noSorp_dir']                       : spc = 'top'
        #  else: spc = 1.0
        node_macro(g,macro_dir=macro_dir,dirxn='vert',spacing=spc)
    #  perm_macro(g, k_bg=1.e-14, macro_dir=macro_dir,zones=zones,zone_perms=zone_perms)
    #  perm_macro(g, k_bg=p['k_m'], macro_dir=macro_dir,zones=zones,zone_perms=zone_perms)
    perm_macro('perm.dat',g, k_bg=p['k_m'], macro_dir=macro_dir,zones=zones,zone_perms=zone_perms)

    #--------------------------------------------- 
    # WRITE ZONE FILE FOR CONTAMINATED ZONE [-8]
    #--------------------------------------------- 
    contamzonen = -8
    contam_thickness = 5. #[m]   (maybe make this 10m later)
    with open(join(macro_dir,'contaminated_layer.zone'),'w+') as fz:
        fz.write('zone\n{}   contaminated/source layer\n'.format(str(abs(contamzonen))))
        fz.write('-0.1 1.e4 1.e4 -0.1\n')
        # Get depth range of source/contam zone
        if 'source_depth_top' not in p:     # default to bottom of domain
            bb=min(g['y'])-0.1; tc=min(g['y'])+contam_thickness+0.1
        else:                               # if depth is specified in run.py
            # Force source_depth_top to be negative
            p['source_depth_top'] = -abs(p['source_depth_top'])
            #  bb=p['source_depth_top']-contam_thickness-0.1; tc=p['source_depth_top']+0.1
            bb=p['source_depth_top']-contam_thickness-0.1; tc=p['source_depth_top']+0.1
            # Check if there is enough space in domain to 
            if p['source_depth_top']-contam_thickness <=  min(g['y']):
                print()
                print('**** WARNING *****')
                print('Source zone top {}m is too close to bottom of domain ({}m).'.format(p['source_depth_top'], min(g['y'])) )
                proceed = input('Do you want to continue? (y/n)')
                if proceed.lower() == 'y': pass
                else: sys.exit('Cancelling simulation.')
        # Write to zone/zonn file
        fz.write('{:0.2f} {:0.2f} {:0.2f} {:0.2f}\n\n'.format(bb,bb,tc,tc))
        fz.write('stop')
    mesh_files.zone_to_zonn(join(macro_dir,'contaminated_layer.zone'))
    #--------------------------------------------- 
    # WRITE ZONE FILE FOR ACTIVE REGOLITH [-10] 
    #--------------------------------------------- 
    active_regolith_zonen = -10
    with open(join(macro_dir,'active_regolith.zone'),'w+') as fz:
        fz.write('zone\n{}     regolith (thermally active)\n'.format(str(abs(active_regolith_zonen))))
        fz.write('-0.1 1.e4 1.e4 -0.1\n')
        # Get depth range of regolith (specified in run.py file)
        top = max(g['y'])
        bb=(p['active_regolith_base'])-0.1; tc=top+0.1
        #  if 'regolith_base' not in p:     # default to 30 m depth 
            #  bb=-30.1; tc=top+0.1
        #  else:                               # if depth is specified in run.py
            #  bb=-abs(p['regolith_base'])-0.1; tc=top+0.1
        # Write to zone/zonn file
        if regolith_type['active'] is True:
            fz.write('{:0.2f} {:0.2f} {:0.2f} {:0.2f}\n\n'.format(bb,bb,tc,tc))
        # If no active regolith, write bogus values
        else:
            fz.write('{:0.2f} {:0.2f} {:0.2f} {:0.2f}\n\n'.format(bb,bb,bb,bb))
        fz.write('stop')
    mesh_files.zone_to_zonn(join(macro_dir,'active_regolith.zone'))
    #--------------------------------------------- 
    # WRITE ZONE FILE FOR CONSTANT REGOLITH [-11] 
    #--------------------------------------------- 
    constant_regolith_zonen = -11
    p['constant_regolith_zonen'] = constant_regolith_zonen
    with open(join(macro_dir,'constant_regolith.zone'),'w+') as fz:
        fz.write('zone\n{}     regolith (constant adsorption coeff)\n'.format(str(abs(constant_regolith_zonen))))
        fz.write('-0.1 1.e4 1.e4 -0.1\n')
        # Get depth range of constant regolith (specified in run.py file)
        if 'active_regolith_base' not in p:
            top = max(g['y'])
        else:
            #  top = -abs(p['active_regolith_base']);
            top = p['active_regolith_base'];
        #  bb=-abs(p['constant_regolith_base'])-0.1; tc=top-0.2
        bb=p['constant_regolith_base']-0.1; tc=top-0.2
        # Write to zone/zonn file
        if regolith_type['constant'] is True:
            fz.write('{:0.2f} {:0.2f} {:0.2f} {:0.2f}\n\n'.format(bb,bb,tc,tc))
        # If no constant regolith, write bogus values 
        else:
            fz.write('{:0.2f} {:0.2f} {:0.2f} {:0.2f}\n\n'.format(bb,bb,bb,bb))
        fz.write('stop')
    mesh_files.zone_to_zonn(join(macro_dir,'constant_regolith.zone'))
    #--------------------------------------------- 
    # WRITE ZONE FILE FOR BOTTOM ZONE [-200] 
    #--------------------------------------------- 
    bottomzonen = -200
    botnodes = g['node'][g['y']==min(g['y'])]
    utils.write_zone(join(macro_dir,'bottom.zone'),botnodes,[bottomzonen])
    mesh_files.zone_to_zonn(join(macro_dir,'bottom.zone'))
    #--------------------------------------------- 
    # ALSO WRITE TOP ZONE AS A ZONN [-100] 
    #--------------------------------------------- 
    #  symlink(join(tpl_dir,p['topzone_file']),join(macro_dir,os.path.split(p['topzone_file'])[1]),overwrite=True)
    shutil.copyfile(join(tpl_dir,p['topzone_file']),join(macro_dir,os.path.split(p['topzone_file'])[1]))
    mesh_files.zone_to_zonn(join(macro_dir,'top.zone'))
    #--------------------------------------------- 
    # WRITE ZONE FILE FOR TOP AND BOTTOM OBS NODES  [-12]
    #--------------------------------------------- 
    topbot_obszonen = -12
    topnode = min(g['node'][g['y']==max(g['y'])]) #get leftmost (min)
    botnode = min(g['node'][g['y']==min(g['y'])]) #get leftmost (min)
    utils.write_zone(join(macro_dir,'obs_topbot.zone'),[botnode,topnode],[topbot_obszonen])
    mesh_files.zone_to_zonn(join(macro_dir,'obs_topbot.zone'))

    #  #--------------------------------------------- 
    #  #--------------------------------------------- 
    #  # WRITE THE MAIN BOUN FILE TO BE USED 
    #  #--------------------------------------------- 
    #  #--------------------------------------------- 
    #  print('writing the main boun file to be used...')
    #  #---- (1) Get the pressures for the surface -------------------------------
    #  datadir='/project/gas_seepage/jportiz/mars/data/pressure/rems/processed_data'
    #  if mars_tools.days2sols(p['sim_time'])<=2712.41: #[sols] (extent of actual pressure record)
        #  pres_df = pd.read_csv(join(datadir,'pchip_press_elev_correction.csv'))
    #  else:
        #  #  pres_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_press_50.csv'))
        #  pres_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_press_50_Ls.csv'))
    #  actual_press_soltime_start = pres_df['SOL_TIME'].iat[0]       #[sols]
    #  actual_press_Ls_start = pres_df['L_s'].iat[0]                 #[°]
    #  ts = np.array(pres_df['SOL_TIME']-pres_df['SOL_TIME'].iat[0]) #[sols] set first time to 'zero'
    #  ts = mars_tools.sols2days(ts) #[d] convert sols to days for FEHM boun file
    #  ps = np.array(pres_df['PRESS_ADJUSTED']/1e6) #[MPa]
    #  #  temps = [0.05]*len(ps) #[deg C]
    #  #  #--Shift the pressures up to 100x the initial pressure value for FEHM
    #  #  fehm_upshift = ps[0]*100. - ps[0] #MPa
    #  #--Shift the pressures up to 1000x the initial pressure value for FEHM 
    #  #  pupshift_factor = 1000 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
    #  pupshift_factor = 100 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
    #  fehm_upshift = ps[0]*pupshift_factor - ps[0] #[MPa]
    #  ps = ps + fehm_upshift #this shifts up to FEHM range, but doesn't distort amplitudes
    #  P_initial = ps[0] #[MPa]
    #  p['P_initial'] = P_initial
    #  topzonen = -100 #zone to apply BC changes to
    #  outfile=p['boun_file']#'mars_press.boun'
#  
    #  #---- (2) Now get the subsurface temperatures from the heat/thermal sim -----------
    #  # Assumes that the T_upshift was the same as it is in current simulation
    #  tdata, t_cols, nodes, ys = parse_hist(join(mesh_dict[p['mesh']]['heat_dir'],'thermal_spinup'), param='temp')
    #  # Only need temperatures at depths < regolith_base  (plus bottom temp)
    #  T_inds   = ys >= -abs(p['regolith_base'])  #indices
    #  T_depths = ys[T_inds]
    #  # Subset the temperature data to those indices
    #  td = tdata[:,1:]
    #  td = td[:,T_inds]
    #  tdata_ss = np.column_stack( (tdata[:,0], td) )
    #  #----------------------------
    #  # MATCH UP STARTING TIMES WITH PRESSURE FORCING
    #  #----------------------------
    #  # Look at original stitched temperature timeseries to interpolate thermal
    #  # data to solar longitude (L_s)
    #  temp_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_amb_temp_80_Ls.csv'))
    #  # Get only the times from late in the simulation (@steady-state) that align
    #  # with the actual start time of the pressure forcing (start=1.551587 sols,
    #  # gets set to 0 for FEHM)
    #  # Grab the indices of the last 2.0 years of thermal sim data
    #  tfinal = tdata_ss[-1,0] #[days]
    #  latetime_inds = np.where(tdata_ss[:,0]>=(tfinal-mars_tools.sols2days(2.0*668.6)))[0]
    #  latetime_days = tdata_ss[latetime_inds,0]  #[days] start FEHM @ 0.0 days
    #  latetime_daysActual = latetime_days+mars_tools.sols2days(actual_press_soltime_start) #[days] actual start time is a couple of days (1.55 sols)
    #  latetime_solsActual = mars_tools.days2sols(latetime_daysActual)
    #  latetime_data = tdata_ss[latetime_inds,:]
    #  #---- Get the L_s of the latetime data from thermal sims
    #  f0 = interp1d(temp_df['SOL_TIME'],temp_df['L_s'], fill_value='extrapolate')
    #  latetime_LsActual = f0(latetime_solsActual)
    #  #  nearest_Ls = find_nearest(latetime_LsActual, actual_press_Ls_start)
    #  #  nearest_Ls_ind = np.where(latetime_LsActual==nearest_Ls)[0][0]
    #  # Subset only latetime data with at least 1 MY remaining    
    #  i = 0
    #  while tfinal-latetime_days[i]>=mars_tools.sols2days(668.6): i+=1
    #  nearest_Ls = find_nearest(latetime_LsActual[:i], actual_press_Ls_start)
    #  nearest_Ls_ind = np.where(latetime_LsActual[:i]==nearest_Ls)[0][0]
    #  # Use nearest_Ls_ind as the approximate start time for the thermal data
    #  # Get ind of late start time + 1 MY
    #  nearest_my = find_nearest(latetime_solsActual,latetime_solsActual[nearest_Ls_ind]+668.59)
    #  my_ind = np.where(latetime_solsActual==nearest_my)[0][0]
    #  #---- Finally, interpolate the times of temps to match pressure dt
    #  # (NOT SURE THIS WORKS)
    #  tinterp = ts[ts<=mars_tools.sols2days(668.59)]
    #  new_latetime_data = np.zeros( (len(tinterp), np.shape(latetime_data)[1]) )
    #  for col in np.arange(1,np.shape(latetime_data)[1]):
        #  f1 = interp1d(latetime_days-latetime_days[0], latetime_data[:,col])
        #  #  new_latetime_data[col] = f1(tinterpo)
        #  new_latetime_data[:,col] = f1(tinterp)
    #  new_latetime_data[:,0] = tinterp  #[days]
    #  # Make sure the temperatures are continuous – shave off a bit
    #  # MIGHT NOT NEED TO DO THIS....
    #  #-- Get "shave" inds by looking at surface temperatures
    #  new_t = new_latetime_data[:,0]   #times
    #  new_p = new_latetime_data[:,-1]  #surface temps
    #  w = 100  #window size
    #  # LOOK FOR MAXIMUM
    #  # Get last temp maximum before the added record (blue in plot)
    #  e_max = max(new_p[-w:])                   #temp
    #  e_max_ind = np.where(new_p==e_max)[0][0]  #index
    #  # Get first temp maximum in the added record    (red in plot)
    #  b_max = max(new_p[:w])                    #temp
    #  b_max_ind = np.where(new_p==b_max)[0][0]  #index
    #  # Get number of data points before maximum of beginning-record (=b_max_ind)
    #  #  #Plot 1 year repeated 1x for visual inspection 
    #  #  # Remove data after the end-record maximum (peak)
    #  #  # Also remove b_max_ind # of points before end-record maximum
    #  #  new_t = new_t[:(e_max_ind-b_max_ind)]
    #  #  new_p = new_p[:(e_max_ind-b_max_ind)]
    #  #  w = 100  #window size
    #  #  fig, (ax1,ax2) = plt.subplots(2)
    #  #  #Full 2 years
    #  #  ax1.plot(new_t,new_p,color='b')
    #  #  ax1.plot(new_t+(new_t[-1]-new_t[0]),new_p,color='r')
    #  #  #Zoom in to where they meet
    #  #  ax2.plot(new_t,new_p,color='b')
    #  #  ax2.plot(new_t+(new_t[-1]-new_t[0]),new_p,color='r')
    #  #  ax2.set_xlim(new_t[-1]-1.5, new_t[-1]+1.5)
    #  #  plt.savefig('test_post-shave.pdf')
    #  #-- Now apply the "shaving" to entire latetime temp data subset 
    #  new_latetime_data = new_latetime_data[:(e_max_ind-b_max_ind),:]
#  
    #  # Get nodes in F&T mesh (NOT heatflow mesh -- 1D) coinciding with those depths
    #  inds = np.where(g['y'] >= -abs(p['regolith_base']))[0]
    #  Tmeshdepths = g[inds]
    #  Tdepth_dict = {}
    #  for depth in np.unique(Tmeshdepths['y']):
        #  #  nn = g['node'][np.where(Tmeshdepths['y']==depth)[0]]
        #  nn = g['node'][np.where(g['y']==depth)[0]]
        #  #  minnode = min(nn); maxnode = max(nn)
        #  minnode = min(nn); maxnode = max(nn)
        #  Tdepth_dict[depth] = [minnode,maxnode]
#  
    #  #-------------------------------------------------- 
    #  #   NOW WRITE PRESSURES AND TEMPERATURES TO .boun FILE 
    #  #-------------------------------------------------- 
    #  print('writing T and P to boun file...')
    #  #---- GATHER BC CHANGES BY WHERE THEY TAKE PLACE -----
    #  bc_dict = {}  # fill with sub-directories for each set of depths
    #  #---- SURFACE MODEL
    #  y = 0
    #  topzonen = -100 #zone to apply BC changes to
    #  sm_dict = {}
    #  times_surf = ts[::4]
    #  press_surf = ps[::4]
    #  #  times_surf = ts[::8]  # JPO - try less temporal resolution
    #  #  press_surf = ps[::8]
    #  # Interpolate temperatures to the pressure forcing times
    #  temp_col = np.where(T_depths == y)[0][0]+1
    #  f0 = interp1d(tdata_ss[:,0], tdata_ss[:,temp_col], fill_value='extrapolate')
    #  #  f0 = interp1d(new_latetime_data[:,0], new_latetime_data[:,temp_col], fill_value='extrapolate')
    #  temp_surf = f0(times_surf)
    #  # Fill the sub-dictionary
    #  sm_dict['times'] = times_surf
    #  #  sm_dict['times'] = times_surf[:100]  #DEBUG
    #  sm_dict['type']  = 'ti_linear'
    #  sm_dict['vars']  = ['pw', 't']
    #  sm_dict['vals']  = [press_surf, temp_surf]
    #  #  sm_dict['vals']  = [press_surf[:100], temp_surf[:100]] #DEBUG
    #  sm_dict['loc']   = [topzonen]
    #  # Add this to main BC directory
    #  bc_dict['top'] = sm_dict
#  
    #  #---- DEPTHS -1 - regolith_base
    #  counter = 0
    #  for y in T_depths[:-1]:
        #  counter+=1
        #  sm_dict={}
        #  temp_col = np.where(T_depths == y)[0][0]+1
        #  #  # Interpolate temperatures to the pressure forcing times
        #  #  #  f0 = interp1d(tdata_ss[:,0], tdata_ss[:,temp_col], fill_value='extrapolate')
        #  #  f0 = interp1d(new_latetime_data[:,0], new_latetime_data[:,temp_col], fill_value='extrapolate')
        #  #  temp_sub = f0(times_surf)
        #  #  times_sub = new_latetime_data[:,0]
        #  #  temp_sub = new_latetime_data[:,temp_col]
        #  times_sub = new_latetime_data[::4,0]        #hourly
        #  temp_sub  = new_latetime_data[::4,temp_col] #hourly
        #  #  times_sub = new_latetime_data[::8,0]        #  -JPO try more temporal resolution
        #  #  temp_sub  = new_latetime_data[::8,temp_col] #
        #  #  temp_sub = f0(times_surf)[:100]  #DEBUG
        #  # Fill the sub-dictionary
        #  sm_dict['times'] = times_sub
        #  #  sm_dict['times'] = times_surf[:100]    #DEBUG
        #  sm_dict['type']  = 'cy_linear'
        #  sm_dict['vars']  = 't'
        #  sm_dict['vals']  = temp_sub
        #  #  sm_dict['vals']  = temp_sub[:100]     #DEBUG
        #  sm_dict['loc']   = Tdepth_dict[y]  #JA JB JC node specification (1 node apart)
        #  # Add this to main BC directory
        #  kname = 'depth'+str(counter)
        #  #  bc_dict[y] = sm_dict
        #  bc_dict[kname] = sm_dict
#  
    #  #  #---- BOTTOM OF DOMAIN (y=-200m) constant temperature
    #  #  sm_dict={}
    #  #  y = min(g['y'])
    #  #  # Fill the sub-dictionary
    #  #  sm_dict['times'] = [0., 1e20]
    #  #  sm_dict['type']  = 'ti'
    #  #  sm_dict['vars']  = 't'
    #  #  #  sm_dict['vals']  = [T_bot]*2       #calculated earlier
    #  #  sm_dict['vals']  = np.array([T_bot]*2)       #calculated earlier
    #  #  sm_dict['loc']   = [int(-200)]     #zone -200 is bottom 
    #  #  # Add this to main BC directory    #(maybe don't add)
    #  #  #  bc_dict[y] = sm_dict
    #  #  bc_dict['bottom'] = sm_dict
#  
    #  #---- "AIR" INJECTION IN CONTAMINATED ZONE (METHANE PRODUCTION ZONE) 
    #  sm_dict={}
    #  # Mars microbe cells distribution
    #  #  cells=1e6 #[cells/m^3]
    #  #--FOR REFERENCE--
    #  #  lowprod = 1.7e-17 #[mol/cell/day] 
    #  #  highprod = 4.88e-14 #[mol/cell/day] 
    #  #  serpentinizaton = 1e-17 #[kg/s/m3 solid rock] (max.)
    #  #-----------------
    #  #---- Air injection rate (distributed 'water' source in contaminated zone)
    #  inj_rate = -1e-15 #[kg/s] source is negative ((may need to make higher or lower so no pressure increase)) (-1e-3 seems too high)
    #  p['inj_rate'] = inj_rate
    #  #---- Total (bulk) volume of contaminated zone
    #  domain_width = max(g['x']) - min(g['x']) #[m]
    #  vol_bulk = domain_width * contam_thickness * 1. #[m^3]
    #  # Fill the sub-dictionary
    #  sm_dict['times'] = [0., 1e20]
    #  sm_dict['type']  = 'ti'
    #  sm_dict['vars']  = 'dsw'          #distributed "water" source
    #  #  sm_dict['vals']  = [inj_rate]*2
    #  sm_dict['vals']  = np.array([inj_rate]*2)
    #  sm_dict['loc']   = [contamzonen]   #zone -8 is contam zone 
    #  # Add this to main BC directory    #(maybe don't add)
    #  bc_dict['contam_zone'] = sm_dict
#  
    #  # Write the dictionary to lists to be read into write_boun() function
    #  timelist = []; varslist = []; valslist = []; typelist=[]; zonelist=[]
    #  for key in bc_dict:
        #  #  print(key)
        #  timelist.append(bc_dict[key]['times'])
        #  varslist.append(bc_dict[key]['vars'])
        #  valslist.append(bc_dict[key]['vals'])
        #  typelist.append(bc_dict[key]['type'])
        #  zonelist.append(bc_dict[key]['loc'])
    #  #DEBUG
    #  #  out_file = join(macro_dir,'test')
    #  #  t=timelist
    #  #  ys=valslist
    #  #  step=typelist
    #  #  varname=varslist
    #  #  zones=zonelist
#  
    #  # Write the .boun
    #  utils.write_boun(join(macro_dir,outfile),timelist,valslist,typelist,varslist,zonelist)
    #  print('done.')
    #  #-------------------------------------------------- 


    #--------------------------------------------- 
    #--------------------------------------------- 
    # WRITE THE MAIN BOUN FILE TO BE USED 
    #--------------------------------------------- 
    #--------------------------------------------- 
    # NOTE: This section of code reads in the surface pressures for the boun macro. 
    #       It also reads in the subsurface temperatures from the thermal simulations
    #       even though those are no longer used in the boun file (these sims are
    #       constant temperature. The subsurface temperatures will be used to calculate
    #       the distribution coefficients that will be written to the userr_data.dat 
    #       file for the new rxn subroutine called 'userr'.
    print('writing the main boun file to be used...')
    datadir='/project/gas_seepage/jportiz/mars/data/pressure/rems/processed_data'
    #--------------------------------------------- 
    # SYNTHETIC PRESSURE AND TEMPERATURE
    #--------------------------------------------- 
    if p['synthetic'] == True:
        #-------------------------------------------------- 
        #---- (1) Generate synthetic pressures 
        #-------------------------------------------------- 
        df = pd.read_csv(join(datadir,'pchip_press_elev_correction.csv'))
        mean_pas = df['PRESS_ADJUSTED'].mean()         #[Pa]
        # Load synthetic temperature components generated previously
        synth_components_dir = '/project/gas_seepage/jportiz/mars/analytical/fourier_analysis/generate_synthetic_data/output'
        with open(join(synth_components_dir, 'vars_pressFullSignal.pkl'), 'rb') as f:
            periods,amps,shifts = pickle.load(f)
        # Create times for BC changes
        # Convert sols to days for FEHM
        dt = 0.04                                      #[days]
        ts = np.arange(0., p['sim_time']+dt, dt)       #[days]     #full record

        # Generate synthetic pressures (and convert to MPa)
        ps = mars_tools.synth_fourier_series(periods+amps+shifts, mean_pas, days2sols(ts), num_modulations=1)/1e6 #[MPa]
        #--Shift the pressures up to 1000x the initial pressure value for FEHM 
        #  pupshift_factor = 1000 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
        pupshift_factor = 100 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
        fehm_upshift = ps[0]*pupshift_factor - ps[0] #[MPa]
        ps = ps + fehm_upshift #this shifts up to FEHM range, but doesn't distort amplitudes
        P_initial = ps[0] #[MPa]
        p['P_initial'] = P_initial
        topzonen = -100 #zone to apply BC changes to
        outfile=p['boun_file']#'mars_press.boun'


        #-------------------------------------------------- 
        #---- (2) Now get the synthetic subsurface temperatures from the heat/thermal sim -----------
        #-------------------------------------------------- 
        # Assumes that the T_upshift was the same as it is in current simulation
        mean_temp = p['T_ref']                         #[°C]

        tdata, t_cols, nodes, ys = parse_hist(join(mesh_dict[p['mesh']]['heat_dir']+'_synth','thermal_spinup'), param='temp')
        if regolith_type['active'] is True:
            # Only need temperatures at depths < regolith_base  (plus bottom temp)
            #  T_inds   = ys >= -abs(p['active_regolith_base'])  #indices
            T_inds   = ys >= p['active_regolith_base']  #indices
        # Even if not active, still want something for T_depths for other dependencies
        else:
            T_inds   = ys == max(g['y'])
        T_depths = ys[T_inds]
        # Subset the temperature data to those indices
        td = tdata[:,1:]
        td = td[:,T_inds]
        tdata_ss = np.column_stack( (tdata[:,0], td) )
        if regolith_type['constant'] is True:
            # Also get the temperatures within the constant_regolith_base 
            # - constant_regolith_base: no temperature changes, but still adsorption
            #  if regolith_type['constant'] is True:
            #  constT_inds = (ys >= -abs(p['constant_regolith_base'])) & (ys < -abs(p['active_regolith_base']))
            constT_inds = (ys >= p['constant_regolith_base']) & (ys < p['active_regolith_base'])
        # Even if no constant, still want something for constT_depths for other dependencies
        else:
            constT_inds = T_inds
        constT_depths = ys[constT_inds]
        # Now also subset the constant regolith temperature data
        const_td = tdata[:,1:]
        const_td = const_td[:,constT_inds]
        const_tdata_ss =  np.column_stack( (tdata[:,0], const_td) )

        # Get nodes in F&T mesh (NOT heatflow mesh -- 1D) coinciding with the active depths
        inds = np.where(g['y'] >= p['active_regolith_base'])[0]
        Tmeshdepths = g[inds]
        Tdepth_dict = {}
        for depth in np.unique(Tmeshdepths['y']):
            nn = g['node'][np.where(g['y']==depth)[0]]
            minnode = min(nn); maxnode = max(nn)
            Tdepth_dict[depth] = [minnode,maxnode]

        #-------------------------------------------------- 
        #   NOW WRITE PRESSURES AND (CONSTANT) TEMPERATURES TO .boun FILE 
        #-------------------------------------------------- 
        print('writing synthetic T and P to boun file...')
        #---- GATHER BC CHANGES BY WHERE THEY TAKE PLACE -----
        bc_dict = {}  # fill with sub-directories for each set of depths
        #---- SURFACE MODEL
        y = 0
        topzonen = -100 #zone to apply BC changes to
        sm_dict = {}

        #  times_surf = ts[::4]
        #  press_surf = ps[::4]
        times_surf = ts
        press_surf = ps
        # Interpolate temperatures to the pressure forcing times
        temp_col = np.where(T_depths == y)[0][0]+1
        f0 = interp1d(tdata_ss[:,0], tdata_ss[:,temp_col], fill_value='extrapolate')
        #  temp_surf = f0(times_surf)
        temp_surf = f0(ts)
        # Fill the sub-dictionary
        sm_dict['times'] = times_surf
        sm_dict['type']  = 'ti_linear'
        sm_dict['vars']  = ['pw', 't']
        sm_dict['vals']  = [press_surf, temp_surf]
        sm_dict['loc']   = [topzonen]
        # Add this to main BC directory
        bc_dict['top'] = sm_dict

        #Write the surface temperatures to a csv file for use elsewhere
        twocol = np.column_stack( (mars_tools.days2sols(sm_dict['times']), temp_surf) )
        np.savetxt('surface_temperatures.csv', twocol, delimiter=',')

        #---- DEPTHS -1 - active_regolith_base
        counter = 0
        for y in T_depths[:-1]:
            counter+=1
            sm_dict={}
            temp_col = np.where(T_depths == y)[0][0]+1
            # Interpolate temperatures to the pressure forcing times
            f0 = interp1d(tdata_ss[:,0], tdata_ss[:,temp_col], fill_value='extrapolate')
            temp_sub = f0(times_surf)
            # Fill the sub-dictionary
            #  sm_dict['times'] = times_sub
            sm_dict['times'] = times_surf
            #  sm_dict['times'] = times_surf[:100]    #DEBUG
            sm_dict['type']  = 'cy_linear'
            sm_dict['vars']  = ['t']
            sm_dict['vals']  = temp_sub
            sm_dict['loc']   = Tdepth_dict[y]  #JA JB JC node specification (1 node apart)
            # Add this to main BC directory
            kname = 'depth'+str(counter)
            print(kname)
            bc_dict[kname] = sm_dict

        #---- DEPTHS active_regolith_base - constant_regolith_base
        # NOTE: we are getting average temperatures for the "constant" adsorption part of regolith 
        counter = 0
        sm_dict={}
        constT_avgs = []
        kname = 'constant_regolith'
        for y in constT_depths:
            counter+=1
            temp_col = np.where(constT_depths == y)[0][0]+1
            #  #  print('y = {}'.format(y))
            #  times_sub = latetime_const_tdata[::4,0]         #hourly
            #  temp_sub  = latetime_const_tdata[::4,temp_col]  #hourly
            times_sub = times_surf.copy()
            temp_sub = const_tdata_ss[:,temp_col]
            constT_avgs.append(np.mean(temp_sub))
        constT_average = np.mean(constT_avgs) #average up all the averages
        sm_dict['times'] = [0., 1e20]
        sm_dict['type'] = 'ti'
        sm_dict['vars']  = 't'
        sm_dict['vals']  = np.array([constT_average]*2)
        sm_dict['loc']   = [constant_regolith_zonen]   #zone -11 is the constant regolith zone 
        # Add this to main BC directory    #(maybe don't add)
        bc_dict['constant_regolith'] = sm_dict

        #  #---- BOTTOM OF DOMAIN (y=-200m) constant temperature
        #  sm_dict={}
        #  y = min(g['y'])
        #  # Fill the sub-dictionary
        #  sm_dict['times'] = [0., 1e20]
        #  sm_dict['type']  = 'ti'
        #  sm_dict['vars']  = 't'
        #  #  sm_dict['vals']  = [T_bot]*2       #calculated earlier
        #  sm_dict['vals']  = np.array([T_bot]*2)       #calculated earlier
        #  sm_dict['loc']   = [int(-200)]     #zone -200 is bottom 
        #  # Add this to main BC directory    #(maybe don't add)
        #  #  bc_dict[y] = sm_dict
        #  bc_dict['bottom'] = sm_dict

        #---- "AIR" INJECTION IN CONTAMINATED ZONE (METHANE PRODUCTION ZONE) 
        sm_dict={}
        # Mars microbe cells distribution
        #  cells=1e6 #[cells/m^3]
        #--FOR REFERENCE--
        #  lowprod = 1.7e-17 #[mol/cell/day] 
        #  highprod = 4.88e-14 #[mol/cell/day] 
        #  serpentinizaton = 1e-17 #[kg/s/m3 solid rock] (max.)
        #-----------------
        #---- Air injection rate (distributed 'water' source in contaminated zone)
        inj_rate = -1e-15 #[kg/s] source is negative ((may need to make higher or lower so no pressure increase)) (-1e-3 seems too high)
        p['inj_rate'] = inj_rate
        #---- Total (bulk) volume of contaminated zone
        domain_width = max(g['x']) - min(g['x']) #[m]
        vol_bulk = domain_width * contam_thickness * 1. #[m^3]
        # Fill the sub-dictionary
        sm_dict['times'] = [0., 1e20]
        sm_dict['type']  = 'ti'
        sm_dict['vars']  = 'dsw'          #distributed "water" source
        #  sm_dict['vals']  = [inj_rate]*2
        sm_dict['vals']  = np.array([inj_rate]*2)
        sm_dict['loc']   = [contamzonen]   #zone -8 is contam zone 
        # Add this to main BC directory    #(maybe don't add)
        bc_dict['contam_zone'] = sm_dict

        #  # Write the dictionary to lists to be read into write_bounp['sim_time']/times_sub[-1]() function
        #  timelist = []; varslist = []; valslist = []; typelist=[]; zonelist=[]
        #  for key in bc_dict:
            # Calculate distcoeff using the temperatures
            #  #  print(key)
            #  timelist.append(bc_dict[key]['times'])
            #  varslist.append(bc_dict[key]['vars'])
            #  valslist.append(bc_dict[key]['vals'])
            #  typelist.append(bc_dict[key]['type'])
            #  zonelist.append(bc_dict[key]['loc'])
#  
        #  #DEBUG
        #  #  out_file = join(macro_dir,'test')
        #  #  t=timelist
        #  #  ys=valslist
        #  #  step=typelist
        #  #  varname=varslist
        #  #  zones=zonelist
        #---- Fill these out manually since boun only has a couple models now
        t = [bc_dict['top']['times'], bc_dict['contam_zone']['times']]
        valslist = [ [bc_dict['top']['vals'][0], [temp_surf[0]]*len(press_surf)],
                    bc_dict['contam_zone']['vals'] ]
        varslist = [ ['pw','t'], 'dsw' ]
        typelist = ['ti_linear', 'ti']
        zonelist = [ bc_dict['top']['loc'], bc_dict['contam_zone']['loc'] ]

        #---- Write the .boun
        #  utils.write_boun(join(macro_dir,outfile),timelist,valslist,typelist,varslist,zonelist)
        utils.write_boun(join(macro_dir,outfile),t,valslist,typelist,varslist,zonelist)
        print('done.')
        #-------------------------------------------------- 


    #--------------------------------------------- 
    # REAL PRESSURE AND TEMPERATURE
    #--------------------------------------------- 
    elif p['synthetic'] == False:
        #---- (1) Get the pressures for the surface -------------------------------
        if mars_tools.days2sols(p['sim_time'])<=2712.41: #[sols] (extent of actual pressure record)
            pres_df = pd.read_csv(join(datadir,'pchip_press_elev_correction.csv'))
        else:
            #  pres_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_press_50.csv'))
            pres_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_press_50_Ls.csv'))
        actual_press_soltime_start = pres_df['SOL_TIME'].iat[0]       #[sols]
        actual_press_Ls_start = pres_df['L_s'].iat[0]                 #[°]
        ts = np.array(pres_df['SOL_TIME']-pres_df['SOL_TIME'].iat[0]) #[sols] set first time to 'zero'
        ts = mars_tools.sols2days(ts) #[d] convert sols to days for FEHM boun file
        st_ind = np.where(ts<=p['sim_time'])[0][-1]+1 #last index of times <= sim_time
        ts = ts[:st_ind]  #new
        ps = np.array(pres_df['PRESS_ADJUSTED']/1e6) #[MPa]
        ps = ps[:st_ind]  #new
        #  temps = [0.05]*len(ps) #[deg C]
        #  #--Shift the pressures up to 100x the initial pressure value for FEHM
        #  fehm_upshift = ps[0]*100. - ps[0] #MPa
        #--Shift the pressures up to 1000x the initial pressure value for FEHM 
        #  pupshift_factor = 1000 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
        pupshift_factor = 100 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
        fehm_upshift = ps[0]*pupshift_factor - ps[0] #[MPa]
        ps = ps + fehm_upshift #this shifts up to FEHM range, but doesn't distort amplitudes
        P_initial = ps[0] #[MPa]
        p['P_initial'] = P_initial
        topzonen = -100 #zone to apply BC changes to
        outfile=p['boun_file']#'mars_press.boun'

        #---- (2) Now get the subsurface temperatures from the heat/thermal sim -----------
        # Assumes that the T_upshift was the same as it is in current simulation
        tdata, t_cols, nodes, ys = parse_hist(join(mesh_dict[p['mesh']]['heat_dir'],'thermal_spinup'), param='temp')
        if regolith_type['active'] is True:
            # Only need temperatures at depths < regolith_base  (plus bottom temp)
            #  T_inds   = ys >= -abs(p['active_regolith_base'])  #indices
            T_inds   = ys >= p['active_regolith_base']  #indices
        # Even if not active, still want something for T_depths for other dependencies
        else:
            T_inds   = ys == max(g['y'])
        T_depths = ys[T_inds]
        # Subset the temperature data to those indices
        td = tdata[:,1:]
        td = td[:,T_inds]
        tdata_ss = np.column_stack( (tdata[:,0], td) )
        if regolith_type['constant'] is True:
            # Also get the temperatures within the constant_regolith_base 
            # - constant_regolith_base: no temperature changes, but still adsorption
            #  if regolith_type['constant'] is True:
            #  constT_inds = (ys >= -abs(p['constant_regolith_base'])) & (ys < -abs(p['active_regolith_base']))
            constT_inds = (ys >= p['constant_regolith_base']) & (ys < p['active_regolith_base'])
        # Even if no constant, still want something for constT_depths for other dependencies
        else:
            constT_inds = T_inds
        constT_depths = ys[constT_inds]
        # Now also subset the constant regolith temperature data
        const_td = tdata[:,1:]
        const_td = const_td[:,constT_inds]
        const_tdata_ss =  np.column_stack( (tdata[:,0], const_td) )
        #----------------------------
        # MATCH UP STARTING TIMES WITH PRESSURE FORCING
        #----------------------------
        # Look at original stitched temperature timeseries to interpolate thermal
        # data to solar longitude (L_s)
        temp_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_amb_temp_80_Ls.csv'))
        # Get only the times from late in the simulation (@steady-state) that align
        # with the actual start time of the pressure forcing (start=1.551587 sols,
        # gets set to 0 for FEHM)
        # Grab the indices of the last 2.0 years of thermal sim data
        tfinal = tdata_ss[-1,0] #[days]
        latetime_inds = np.where(tdata_ss[:,0]>=(tfinal-mars_tools.sols2days(2.0*668.6)))[0]
        latetime_days = tdata_ss[latetime_inds,0]  #[days] start FEHM @ 0.0 days
        latetime_daysActual = latetime_days+mars_tools.sols2days(actual_press_soltime_start) #[days] actual start time is a couple of days (1.55 sols)
        latetime_solsActual = mars_tools.days2sols(latetime_daysActual)
        latetime_data = tdata_ss[latetime_inds,:]
        # Also get the latetime data of the constant-temperature regolith
        latetime_const_tdata = const_tdata_ss[latetime_inds,:]
        #---- Get the L_s of the latetime data from thermal sims
        f0 = interp1d(temp_df['SOL_TIME'],temp_df['L_s'], fill_value='extrapolate')
        latetime_LsActual = f0(latetime_solsActual)
        #  nearest_Ls = find_nearest(latetime_LsActual, actual_press_Ls_start)
        #  nearest_Ls_ind = np.where(latetime_LsActual==nearest_Ls)[0][0]
        # Subset only latetime data with at least 1 MY remaining    
        i = 0
        while tfinal-latetime_days[i]>=mars_tools.sols2days(668.6): i+=1
        nearest_Ls = find_nearest(latetime_LsActual[:i], actual_press_Ls_start)
        nearest_Ls_ind = np.where(latetime_LsActual[:i]==nearest_Ls)[0][0]
        # Use nearest_Ls_ind as the approximate start time for the thermal data
        # Get ind of late start time + 1 MY
        nearest_my = find_nearest(latetime_solsActual,latetime_solsActual[nearest_Ls_ind]+668.59)
        my_ind = np.where(latetime_solsActual==nearest_my)[0][0]
        #---- Finally, interpolate the times of temps to match pressure dt
        # (NOT SURE THIS WORKS)
        tinterp = ts[ts<=mars_tools.sols2days(668.59)]
        new_latetime_data = np.zeros( (len(tinterp), np.shape(latetime_data)[1]) )
        for col in np.arange(1,np.shape(latetime_data)[1]):
            f1 = interp1d(latetime_days-latetime_days[0], latetime_data[:,col])
            #  new_latetime_data[col] = f1(tinterpo)
            new_latetime_data[:,col] = f1(tinterp)
        new_latetime_data[:,0] = tinterp  #[days]
        # Make sure the temperatures are continuous – shave off a bit
        # MIGHT NOT NEED TO DO THIS....
        #-- Get "shave" inds by looking at surface temperatures
        new_t = new_latetime_data[:,0]   #times
        new_p = new_latetime_data[:,-1]  #surface temps
        w = 100  #window size
        # LOOK FOR MAXIMUM
        # Get last temp maximum before the added record (blue in plot)
        e_max = max(new_p[-w:])                   #temp
        e_max_ind = np.where(new_p==e_max)[0][0]  #index
        # Get first temp maximum in the added record    (red in plot)
        b_max = max(new_p[:w])                    #temp
        b_max_ind = np.where(new_p==b_max)[0][0]  #index
        # Get number of data points before maximum of beginning-record (=b_max_ind)
        #  #Plot 1 year repeated 1x for visual inspection 
        #  # Remove data after the end-record maximum (peak)
        #  # Also remove b_max_ind # of points before end-record maximum
        #  new_t = new_t[:(e_max_ind-b_max_ind)]
        #  new_p = new_p[:(e_max_ind-b_max_ind)]
        #  w = 100  #window size
        #  fig, (ax1,ax2) = plt.subplots(2)
        #  #Full 2 years
        #  ax1.plot(new_t,new_p,color='b')
        #  ax1.plot(new_t+(new_t[-1]-new_t[0]),new_p,color='r')
        #  #Zoom in to where they meet
        #  ax2.plot(new_t,new_p,color='b')
        #  ax2.plot(new_t+(new_t[-1]-new_t[0]),new_p,color='r')
        #  ax2.set_xlim(new_t[-1]-1.5, new_t[-1]+1.5)
        #  plt.savefig('test_post-shave.pdf')
        #-- Now apply the "shaving" to entire latetime temp data subset 
        new_latetime_data = new_latetime_data[:(e_max_ind-b_max_ind),:]

        # Get nodes in F&T mesh (NOT heatflow mesh -- 1D) coinciding with the active depths
        #  inds = np.where(g['y'] >= -abs(p['regolith_base']))[0]
        inds = np.where(g['y'] >= p['active_regolith_base'])[0]
        Tmeshdepths = g[inds]
        Tdepth_dict = {}
        for depth in np.unique(Tmeshdepths['y']):
            #  nn = g['node'][np.where(Tmeshdepths['y']==depth)[0]]
            nn = g['node'][np.where(g['y']==depth)[0]]
            #  minnode = min(nn); maxnode = max(nn)
            minnode = min(nn); maxnode = max(nn)
            Tdepth_dict[depth] = [minnode,maxnode]

        #-------------------------------------------------- 
        #   NOW WRITE PRESSURES AND (CONSTANT) TEMPERATURES TO .boun FILE 
        #-------------------------------------------------- 
        print('writing real T and P to boun file...')
        #---- GATHER BC CHANGES BY WHERE THEY TAKE PLACE -----
        bc_dict = {}  # fill with sub-directories for each set of depths
        #---- SURFACE MODEL
        y = 0
        topzonen = -100 #zone to apply BC changes to
        sm_dict = {}

        times_surf = ts[::4]
        press_surf = ps[::4]
        # Interpolate temperatures to the pressure forcing times
        temp_col = np.where(T_depths == y)[0][0]+1
        f0 = interp1d(tdata_ss[:,0], tdata_ss[:,temp_col], fill_value='extrapolate')
        #  f0 = interp1d(new_latetime_data[:,0], new_latetime_data[:,temp_col], fill_value='extrapolate')
        temp_surf = f0(times_surf)
        # Fill the sub-dictionary
        sm_dict['times'] = times_surf
        #  sm_dict['times'] = times_surf[:100]  #DEBUG
        sm_dict['type']  = 'ti_linear'
        #  sm_dict['vars']  = ['pw', 't','en']
        #  sm_dict['vars']  = ['pw', 't','s']
        sm_dict['vars']  = ['pw', 't']
        #  sm_dict['vals']  = [press_surf, temp_surf, np.zeros(len(temp_surf))]
        #  sm_dict['vals']  = [press_surf, temp_surf, np.ones(len(temp_surf))]
        sm_dict['vals']  = [press_surf, temp_surf]
        sm_dict['loc']   = [topzonen]
        # Add this to main BC directory
        bc_dict['top'] = sm_dict

        #Write the surface temperatures to a csv file for use elsewhere
        twocol = np.column_stack( (mars_tools.days2sols(sm_dict['times']), temp_surf) )
        np.savetxt('surface_temperatures.csv', twocol, delimiter=',')

        #---- DEPTHS -1 - active_regolith_base
        counter = 0
        for y in T_depths[:-1]:
            counter+=1
            sm_dict={}
            temp_col = np.where(T_depths == y)[0][0]+1
            #  # Interpolate temperatures to the pressure forcing times
            #  #  f0 = interp1d(tdata_ss[:,0], tdata_ss[:,temp_col], fill_value='extrapolate')
            #  f0 = interp1d(new_latetime_data[:,0], new_latetime_data[:,temp_col], fill_value='extrapolate')
            #  temp_sub = f0(times_surf)
            #  times_sub = new_latetime_data[:,0]
            #  temp_sub = new_latetime_data[:,temp_col]
            times_sub = new_latetime_data[::4,0]        #hourly
            temp_sub  = new_latetime_data[::4,temp_col] #hourly
            #  temp_sub = f0(times_surf)[:100]  #DEBUG
            # Fill the sub-dictionary
            sm_dict['times'] = times_sub
            #  sm_dict['times'] = times_surf[:100]    #DEBUG
            sm_dict['type']  = 'cy_linear'
            #  sm_dict['vars']  = ['t', 'en']
            #  sm_dict['vars']  = ['t', 's']
            sm_dict['vars']  = ['t']
            #  sm_dict['vals']  = [temp_sub, np.ones(len(temp_sub))]
            sm_dict['vals']  = temp_sub
            #  sm_dict['vals']  = temp_sub[:100]     #DEBUG
            sm_dict['loc']   = Tdepth_dict[y]  #JA JB JC node specification (1 node apart)
            # Add this to main BC directory
            kname = 'depth'+str(counter)
            print(kname)
            #  bc_dict[y] = sm_dict
            bc_dict[kname] = sm_dict

        #---- DEPTHS active_regolith_base - constant_regolith_base
        # NOTE: we are getting average temperatures for the "constant" adsorption part of regolith 
        counter = 0
        sm_dict={}
        constT_avgs = []
        kname = 'constant_regolith'
        for y in constT_depths:
            counter+=1
            temp_col = np.where(constT_depths == y)[0][0]+1
            #  print('y = {}'.format(y))
            times_sub = latetime_const_tdata[::4,0]         #hourly
            temp_sub  = latetime_const_tdata[::4,temp_col]  #hourly
            constT_avgs.append(np.mean(temp_sub))
        constT_average = np.mean(constT_avgs) #average up all the averages
        sm_dict['times'] = [0., 1e20]
        sm_dict['type'] = 'ti'
        sm_dict['vars']  = 't'
        sm_dict['vals']  = np.array([constT_average]*2)
        sm_dict['loc']   = [constant_regolith_zonen]   #zone -11 is the constant regolith zone 
        # Add this to main BC directory    #(maybe don't add)
        bc_dict['constant_regolith'] = sm_dict

        #  #---- BOTTOM OF DOMAIN (y=-200m) constant temperature
        #  sm_dict={}
        #  y = min(g['y'])
        #  # Fill the sub-dictionary
        #  sm_dict['times'] = [0., 1e20]
        #  sm_dict['type']  = 'ti'
        #  sm_dict['vars']  = 't'
        #  #  sm_dict['vals']  = [T_bot]*2       #calculated earlier
        #  sm_dict['vals']  = np.array([T_bot]*2)       #calculated earlier
        #  sm_dict['loc']   = [int(-200)]     #zone -200 is bottom 
        #  # Add this to main BC directory    #(maybe don't add)
        #  #  bc_dict[y] = sm_dict
        #  bc_dict['bottom'] = sm_dict

        #---- "AIR" INJECTION IN CONTAMINATED ZONE (METHANE PRODUCTION ZONE) 
        sm_dict={}
        # Mars microbe cells distribution
        #  cells=1e6 #[cells/m^3]
        #--FOR REFERENCE--
        #  lowprod = 1.7e-17 #[mol/cell/day] 
        #  highprod = 4.88e-14 #[mol/cell/day] 
        #  serpentinizaton = 1e-17 #[kg/s/m3 solid rock] (max.)
        #-----------------
        #---- Air injection rate (distributed 'water' source in contaminated zone)
        inj_rate = -1e-15 #[kg/s] source is negative ((may need to make higher or lower so no pressure increase)) (-1e-3 seems too high)
        p['inj_rate'] = inj_rate
        #---- Total (bulk) volume of contaminated zone
        domain_width = max(g['x']) - min(g['x']) #[m]
        vol_bulk = domain_width * contam_thickness * 1. #[m^3]
        # Fill the sub-dictionary
        sm_dict['times'] = [0., 1e20]
        sm_dict['type']  = 'ti'
        sm_dict['vars']  = 'dsw'          #distributed "water" source
        #  sm_dict['vals']  = [inj_rate]*2
        sm_dict['vals']  = np.array([inj_rate]*2)
        sm_dict['loc']   = [contamzonen]   #zone -8 is contam zone 
        # Add this to main BC directory    #(maybe don't add)
        bc_dict['contam_zone'] = sm_dict

        #  # Write the dictionary to lists to be read into write_bounp['sim_time']/times_sub[-1]() function
        #  timelist = []; varslist = []; valslist = []; typelist=[]; zonelist=[]
        #  for key in bc_dict:
            # Calculate distcoeff using the temperatures
            #  #  print(key)
            #  timelist.append(bc_dict[key]['times'])
            #  varslist.append(bc_dict[key]['vars'])
            #  valslist.append(bc_dict[key]['vals'])
            #  typelist.append(bc_dict[key]['type'])
            #  zonelist.append(bc_dict[key]['loc'])
#  
        #  #DEBUG
        #  #  out_file = join(macro_dir,'test')
        #  #  t=timelist
        #  #  ys=valslist
        #  #  step=typelist
        #  #  varname=varslist
        #  #  zones=zonelist
        #---- Fill these out manually since boun only has a couple models now
        t = [bc_dict['top']['times'], bc_dict['contam_zone']['times']]
        valslist = [ [bc_dict['top']['vals'][0], [temp_surf[0]]*len(press_surf)],
                    bc_dict['contam_zone']['vals'] ]
        varslist = [ ['pw','t'], 'dsw' ]
        typelist = ['ti_linear', 'ti']
        zonelist = [ bc_dict['top']['loc'], bc_dict['contam_zone']['loc'] ]

        #---- Write the .boun
        #  utils.write_boun(join(macro_dir,outfile),timelist,valslist,typelist,varslist,zonelist)
        utils.write_boun(join(macro_dir,outfile),t,valslist,typelist,varslist,zonelist)
        print('done.')
        #-------------------------------------------------- 

    #--------------------------------------------- 
    # WRITE THE EOS MACRO FOR CO2 "AIR" PROPERTIES
    #--------------------------------------------- 
    print('writing the eos macro...')
    #  P_ref = 0.0006      # [MPa] equivalent to 600 Pa  (for some reason, setting P_ref to 0.0008 makes negative densities...)
    P_ref = 0.0007     # [MPa] equivalent to 700 Pa  (for some reason, setting P_ref to 0.0008 makes negative densities...)
    T_ref_eos = -50   # [°C]  the EOS T_ref can be -50, but T_ref in init should be mean surf. temp. @ Gale 
    #  P_upshift = P_ref * 100 - P_ref #shift pressures up to FEHM range without distorting amplitudes
    #  P_upshift = P_ref * 100 #shift pressures up to FEHM range without distorting amplitudes
    P_upshift = P_ref * pupshift_factor
    T_upshift = 95.    # [°C] trick FEHM to think it's above 0° (and keep T between 0-100°C)
    co2_props_file = '/project/gas_seepage/jportiz/mars/data/co2_properties/co2_isothermal_allTemps.csv'
    d_co2 = pd.read_csv(co2_props_file)
    # Calculate and set the eos parameters
    # EW1-3 -----------------------------
    # (EW1=P_ref; EW2=T_ref; EW3=rho_ref)
    #  EW1 = P_upshift
    # Make the approximation that P_ref is same as starting Pressure
    #  EW1 = round(press_surf[0],5)
    EW1 = round(P_initial,6)
    #  EW2 = T_ref_eos+T_upshift
    #  EW2 = 48.07  # needs to be same as T_initial (surface)
    EW2 = temp_surf[0]  # needs to be same as T_initial (surface)
    EW3 = float(d_co2.loc[(d_co2['Temperature (K)']==T_ref_eos+273)&(round(d_co2['Pressure (MPa)'],5)==P_ref)]['Density (kg/m3)'])
    # EW4 --------------------------------
    # (EW4=d(rho)/dP)
    P0 = P_ref
    # next higher pressure at same temp
    index = min(d_co2.loc[(d_co2['Pressure (MPa)']>P_ref)&(d_co2['Temperature (K)']==T_ref_eos+273)].index)
    P1 = d_co2.loc[index]['Pressure (MPa)']
    rho0 = EW3
    rho1 = d_co2.loc[index]['Density (kg/m3)']
    EW4 = float((rho1-rho0)/(P1-P0))
    # EW5 --------------------------------
    # (EW5=d(rho)/dT)
    T0 = T_ref_eos+273
    # next higher temp at same pressure 
    index = min(d_co2.loc[(d_co2['Temperature (K)']>T0)&(round(d_co2['Pressure (MPa)'],5)==P_ref)].index)
    T1 = d_co2.loc[index]['Pressure (MPa)']
    rho0 = EW3
    rho1 = d_co2.loc[index]['Density (kg/m3)']
    EW5 = float((rho1-rho0)/(T1-T0))
    # EW6 --------------------------------
    # (EW6=E_ref)
    EW6 = float(d_co2.loc[(d_co2['Temperature (K)']==T_ref_eos+273)&(round(d_co2['Pressure (MPa)'],5)==P_ref)]['Enthalpy (kJ/kg)'] / 1000) # [MJ/kg]
    # EW7 --------------------------------
    # (EW7=dE/dP)
    P0 = P_ref
    index = min(d_co2.loc[(d_co2['Pressure (MPa)']>P_ref)&(d_co2['Temperature (K)']==T_ref_eos+273)].index)
    P1 = d_co2.loc[index]['Pressure (MPa)']
    H0 = EW6
    H1 = d_co2.loc[index]['Enthalpy (kJ/kg)']/1000 # [MJ/kg]
    EW7 = (H1-H0)/(P1-P0)    # [(MJ/kg)/MPa]
    # EW8 --------------------------------
    # (EW8=dE/dT)
    T0 = T_ref_eos+273
    index = min(d_co2.loc[(d_co2['Temperature (K)']>T0)&(round(d_co2['Pressure (MPa)'],5)==P_ref)].index)
    T1 = d_co2.loc[index]['Temperature (K)']
    H0 = EW6
    H1 = d_co2.loc[index]['Enthalpy (kJ/kg)']/1000 # [MJ/kg]
    EW8 = (H1-H0)/(T1-T0)    # [(MJ/kg)/°C]
    # EW9 --------------------------------
    # (EW9=mu_ref)
    EW9 = float(d_co2.loc[(d_co2['Temperature (K)']==T_ref_eos+273)&(round(d_co2['Pressure (MPa)'],5)==P_ref)]['Viscosity (Pa*s)']) # [Pa s]
    # EW10 -------------------------------
    # (EW10=d(mu)/dP)
    P0 = P_ref
    index = min(d_co2.loc[(d_co2['Pressure (MPa)']>P_ref)&(d_co2['Temperature (K)']==T_ref_eos+273)].index)
    P1 = d_co2.loc[index]['Pressure (MPa)']
    mu0 = EW9
    mu1 = d_co2.loc[index]['Viscosity (Pa*s)']
    EW10 = (mu1-mu0)/(P1-P0)    # [(Pa*s)/MPa]
    # EW11 -------------------------------
    # (EW11=d(mu)/dT)
    T0 = T_ref_eos+273
    index = min(d_co2.loc[(d_co2['Temperature (K)']>T0)&(round(d_co2['Pressure (MPa)'],5)==P_ref)].index)
    T1 = d_co2.loc[index]['Temperature (K)']
    mu0 = EW9
    mu1 = d_co2.loc[index]['Viscosity (Pa*s)']
    EW11 = (mu1-mu0)/(T1-T0)    # [(Pa*s)/°C]
    # EV1-11 -----------------------------
    EV1=EW1; EV2=EW2; EV3=1.0; EV4=0; EV5=0;EV6=-0.020773
    EV7=0.827; EV8=11.398; EV9=1.372e-4; EV10=0.; EV11=0.
    #---------------------------------------------------- 
    # Write the EOS file
    eos_filename ='eos_mars_co2.dat'
    mars_tools.write_eos(join(macro_dir,eos_filename),EW1,EW2,EW3,EW4,EW5,EW6,EW7,EW8,EW9,EW10,EW11,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11)
    print('done.')

    #  #--------------------------------------------- 
    #  #    (NOT USING THIS FOR TRACER RUNS)
    #  # WRITE THE geo-grad_air-static.boun FILE 
    #  #     - USED FOR INITIALIZATION OF GEOTHERMAL AND AIR-STATIC GRADIENT
    #  #--------------------------------------------- 
    #  print('writing geo-grad_air_static.boun file...')
    #  T_ref = p['T_ref']    #[°C] mean surf. temp. @ Gale (different than T_ref_eos)
    #  T_0 = T_ref+T_upshift #50.#0.#-40.                    #[°C] surf temp upshifted for FEHM
    #  p['T_0'] = T_0
    #  y_bot = min(g['y'])
    #  T_bot = T_0 + p['geothermal_gradient']*abs(y_bot)#calculate T at bottom of domain
    #  p['T_bot'] = T_bot
    #  #  pw = [0.070618]*2
    #  pw = [0.70618]*2
    #  surf_temps = [T_0]*2
    #  bot_temps = [T_bot]*2
    #  #  outfile = 'air_static.boun'
    #  outfile = 'geo-grad_air-static.boun'
    #  utils.write_boun(join(macro_dir,outfile),([0.,1e20],[0.,1e20]),[[pw,surf_temps],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]])
    #  print('done.')

    #--------------------------------------------- 
    # WRITE THE isothermal.boun FILE 
    #     - USED FOR INITIALIZATION AIR-STATIC PRESSURE GRADIENT
    #--------------------------------------------- 
    print('writing isothermal.boun file...')
    T_ref = p['T_ref']    #[°C] mean surf. temp. @ Gale (different than T_ref_eos)
    T_0 = T_ref+T_upshift #50.#0.#-40.                    #[°C] surf temp upshifted for FEHM
    p['T_0'] = T_0
    #  y_bot = min(g['y'])
    #  T_bot = T_0 + p['geothermal_gradient']*abs(y_bot)#calculate T at bottom of domain
    #  p['T_bot'] = T_bot
    #  pw = [0.070618]*2
    #  pw = [0.70618]*2
    pw = [P_initial]*2
    surf_temps = [T_0]*2
    #  bot_temps = [T_bot]*2
    #  outfile = 'air_static.boun'
    outfile = 'isothermal.boun'
    #  utils.write_boun(join(macro_dir,outfile),([0.,1e20],[0.,1e20]),[[pw,surf_temps],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]])
    utils.write_boun(join(macro_dir,outfile),[0.,1e20],[pw,surf_temps],"ti_linear",["pw","t"],zones=[-100])
    print('done.')



    #--------------------------------------------- 
    # FIGURE OUT TRAC PARAMETERS 
    #--------------------------------------------- 
    print('writing trac...')
    # Calculate concentration of the injected/produced methane
    c_inj = (p['methane_prod_rate']*vol_bulk) / abs(inj_rate) #[mol/kg 'vapor']
    p['c_inj_initial'] = (p['methane_prod_rate']*(1*1*1)) / abs(inj_rate)

    #--------------------------------------------- 
    # WRITE THE TRAC FILE 
    #--------------------------------------------- 
    # (Either use an initial source or continuous production rate)
    # INITIAL SOURCE (default)
    if p['methane_prod_rate'] == 0.0:
        sourcetype='initial source'
        with open(join(macro_dir,'trac.dat'),'w') as ft:
            ft.write('trac\n')
            ft.write('1.e-50     1.0      1.e-09   1.0\n')
            ft.write('0.       1.e20    1.0e6   0.0\n')
            ft.write('20      1.6     1.e-3  1.e-2   10\n')
            ft.write('1\n')
            ft.write('1  methane\n')
            #  ft.write('0  0 0 1 0 3.e-05  0 0 0\n\n') #no M-Q
            ft.write('0  0 0 1 1 3.e-05  0 0 0\n\n') #M-Q
            ft.write('1 0 0 1\n\n')
            ft.write('1  0 0   1.e-50\n')
            ft.write('{} 0 0   1.e0\n\n\n'.format(contamzonen))
            ft.write('stop')
    # CONTINUOUS METHANE PRODUCTION SOURCE RATE
    else:
        sourcetype='continuous source'
        # ADSORPTION "ON"
        if p['no_adsorption'] is False:
            with open(join(macro_dir,'trac.dat'),'w') as ft:
                ft.write('trac\n')
                ft.write('1.e-50     1.0      1.e-09   1.0\n')
                ft.write('0.       1.e20    1.0e6   0.0\n')
                ft.write('20      1.6     1.e-3  1.e-2   10\n')
                ft.write('2\n')  #liquid and solid (sorbed) methane
                #-- "LIQUID" METHANE -----
                ft.write('1  methane_l\n')
                #  ft.write('0  0 0 1 0 3.e-05  0 0 0\n\n') #no M-Q
                ft.write('0  0 0 1 1 3.e-05  0 0 0\n\n') #M-Q
                ft.write('1 0 0 1\n\n')
                if p['c_initial'] == 0.0:
                    ft.write('1  0 0   1.e-50\n\n')
                else:
                    ft.write('1  0 0   9.611e-5\n')
                    ft.write('{}  0 0   {}\n\n'.format(contamzonen,p['c_initial']))
                ft.write('{} 0 0   {}   0. 1e10\n\n'.format(contamzonen, c_inj))
                #-- SOLID/ADSORBED METHANE -----
                ft.write('0  methane_s\n')
                ft.write('1  0 0    1e-50\n\n\n')
                ft.write('stop')
            # Also write a trac w/o solid species for the noSorp part
            with open(join(macro_dir,'trac_noSorp.dat'),'w') as ft:
                ft.write('trac\n')
                ft.write('1.e-50     1.0      1.e-09   1.0\n')
                ft.write('0.       1.e20    1.0e6   0.0\n')
                ft.write('20      1.6     1.e-3  1.e-2   10\n')
                ft.write('1\n')
                ft.write('1  methane_l\n')
                #  ft.write('0  0 0 1 0 3.e-05  0 0 0\n\n') #no M-Q
                ft.write('0  0 0 1 1 3.e-05  0 0 0\n\n') #M-Q
                ft.write('1 0 0 1\n\n')
                if p['c_initial'] == 0.0:
                    ft.write('1  0 0   1.e-50\n\n')
                else:
                    ft.write('1  0 0   9.611e-5\n')
                    ft.write('{}  0 0   {}\n\n'.format(contamzonen,p['c_initial']))
                ft.write('{} 0 0   {}   0. 1e10\n\n'.format(contamzonen, c_inj))
                ft.write('stop')
        # ADSORPTION "OFF"
        else:
            with open(join(macro_dir,'trac.dat'),'w') as ft:
                ft.write('trac\n')
                ft.write('1.e-50     1.0      1.e-09   1.0\n')
                ft.write('0.       1.e20    1.0e6   0.0\n')
                ft.write('20      1.6     1.e-3  1.e-2   10\n')
                ft.write('1\n')
                ft.write('1  methane_l\n')
                #  ft.write('0  0 0 1 0 3.e-05  0 0 0\n\n') #no M-Q
                ft.write('0  0 0 1 1 3.e-05  0 0 0\n\n') #M-Q
                ft.write('1 0 0 1\n\n')
                if p['c_initial'] == 0.0:
                    ft.write('1  0 0   1.e-50\n\n')
                else:
                    ft.write('1  0 0   9.611e-5\n')
                    ft.write('{}  0 0   {}\n\n'.format(contamzonen,p['c_initial']))
                ft.write('{} 0 0   {}   0. 1e10\n\n'.format(contamzonen, c_inj))
                ft.write('stop')
            # Also make a copy with unique name
            shutil.copyfile(join(macro_dir,'trac.dat'), join(macro_dir, 'trac_noSorp.dat'))
    print('done.')

    #--------------------------------------------- 
    # WRITE THE RXN MACRO/FILE 
    #--------------------------------------------- 
    # Will have to write to a different rxn *.tpl file depending
    # on 
    print('writing rxn... ')
    # Methane adsorption parameters based on Gough2010
    sorp_exp_data_dir = '/lclscratch/jportiz/projects/gas_seepage/mars/runs/fehm_mars_tests/adsorption/Gough2010_tests/data/fig4'
    if not os.path.exists(sorp_exp_data_dir):
        sorp_exp_data_dir = RUNDIR+'/1d_heat_runs/heat_200m_synth/Gough2010_tests/data/fig4'

    #---- Sample parameters
    #  rho_dict = {114.9553571: 947.0039,
                #  119.3303571: 943.4958,
                #  121.2053571: 941.9659,
                #  123.5714286: 940.0130,
                #  128.6607143: 935.7279,
                #  135.0446429: 930.1898,
                #  142.1428571: 923.8194,
               #  }
    #  rho = 997.1827                          #[kg/m3]  "water" at 23°C and 0.08 MPa
    #  rho = rho_dict[T]                       #[kg/m3]  "water" at lookup temperature 
    #  phi = p['phi_m']                        #[-] just pick a value
    #  #  rho_rock =  m_sample/( V_sample * phi ) #[kg/m3] intrinsic density 
    #  # (Eventually, read density in and write to rock.dat, etc)
    #  rho_rock = 2900.                        #[kg/m3] intrinsic density 
    #  #  rho_bulk = m_sample/V_sample            #[kg/m3] bulk density (probably not used)
    M_w = molar_mass_ch4()/6.022e23         #[kg/molecule] molecular weight of CH4
    C_max = maxconc()
    deltaE_obs =  17991.2                    #[J/mol]
    #  deltaE_obs = -18068.79                   #[J/mol]
    rate = 1.e5                            #[1/hr]
    #-- Create rxn lookup table for T-dependent K_l
    sorp_exp_data = np.genfromtxt(join(sorp_exp_data_dir,'experimental.csv'),delimiter=',',skip_header=1)
    # List of lookup temps and distcoeffs
    temp_list = []
    kl_list = []
    for temp in sorp_exp_data[:,0]:
        temp_list.append(temp)           #[K]
        kl_list.append(K_eq(temp, molar_mass_ch4(), deltaE_obs))
    p['maxconc']    = maxconc()
    p['rate']       = rate
    #
    if p['no_adsorption'] is False:
        #-- If both active and constant regolith zones
        if all(rg is True for rg in regolith_type.values()):
            p['rxn_file'] = 'rxn_2types'
            print('  writing rxn file for both active and constant regolith adsorption...')
            # Calculate the constant distcoeff using the average regolith temperature
            temp = constT_average - T_upshift + 273.15
            constant_distcoeff = K_eq(temp, molar_mass_ch4(), deltaE_obs)
            p['constant_distcoeff'] = constant_distcoeff
        #-- If only an active regolith zone 
        elif regolith_type['active'] is True:
            p['rxn_file'] = 'rxn'
            print('  writing rxn file for just active regolith adsorption...')
        #-- If only a constant regolith zone
        elif  regolith_type['constant'] is True:
            p['rxn_file'] = 'rxn_constant'
            print('  writing rxn file for just constant regolith adsorption...')
            # Calculate the constant distcoeff using the average regolith temperature
            temp = constT_average - T_upshift + 273.15
            constant_distcoeff = K_eq(temp, molar_mass_ch4(), deltaE_obs)
            p['constant_distcoeff'] = constant_distcoeff
        #-- If there was some other option this missed...
        else:
            print('**********************')
            print(' WARNING: did not catch the type of rxn file to write!')
            print('**********************')

        print('  rxn template filename is {}.tpl'.format(p['rxn_file']))
        print('  rxn outputted filename is rxn.dat')
        #  pest_io.tpl_write(p,join(tpl_dir,'rxn.tpl'),join(macro_dir,'rxn.dat'))
        pest_io.tpl_write(p,join(tpl_dir,p['rxn_file']+'.tpl'),join(macro_dir,'rxn.dat'))
    print('done.')

    #--------------------------------------------- 
    # PLOT K_l vs T  (both experimental and projected for FEHM)
    #--------------------------------------------- 
    #---- Calculate K_l into the surface temperature range of the FEHM Mars sim
    T_exp_max = max(sorp_exp_data[:,0])            #[K]
    if p['synthetic'] == True:
        T_fehm_max = np.amax(td[:,-1])-T_upshift+273.15
        T_fehm_min = np.amin(td[:,-1])-T_upshift+273.15
    elif p['synthetic'] == False:
        T_fehm_max = np.amax(new_latetime_data[:,-1])-T_upshift+273.15
        T_fehm_min = np.amin(new_latetime_data[:,-1])-T_upshift+273.15
    T_range_fehm = np.linspace(T_fehm_min, T_fehm_max)
    Kl_range_fehm = [K_eq(temp, molar_mass_ch4(), deltaE_obs) for temp in T_range_fehm]
    #  #---- Plot experimental data points
    #  fig, ax1 = plt.subplots(1,figsize=(8,6.4))
    #  ax1.plot(sorp_exp_data[:,0], kl_list, color='k',ls='',marker='o',mec='k',mfc='k',ms=7.,zorder=100, label='experimental (Gough et al.,  2010)')
    #  #---- Plot exponential fit 
    #  exp_fit = np.polyfit(sorp_exp_data[:,0], np.log(kl_list),1)
    #  xfit = np.linspace(sorp_exp_data[:,0][0], sorp_exp_data[:,0][-1])
    #  yfit = np.exp(exp_fit[1]) * np.exp(exp_fit[0]*xfit)
    #  ax1.plot(xfit, yfit, ls='--',color='k')
    #  #---- Plot FEHM extrapolation 
    #  ax1.plot(T_range_fehm, Kl_range_fehm, ls='--', color='r',label='(eqn 13) FEHM Mars range')
    #  #Try using the exponential fit instead... (not sure which correct)
    #  #**** For now, let's go with Eqn 13 to get some results, but need to eventually figure out... ****
    #  Kl_exp = np.exp(exp_fit[1])*np.exp(exp_fit[0]*T_range_fehm)
    #  ax1.plot(T_range_fehm, Kl_exp, ls='--', color='k',label='(exponential fit) FEHM Mars range')
    #  #---- Axes Properties
    #  ax1.set_xlabel(r'$T$ [K]')
    #  ax1.set_ylabel(r'$K_{L}$')
    #  ax1.legend()
    #  ax1.set_yscale('log')
    #  plt.tight_layout()
    #  plt.savefig('temp_vs_Kl_plot.pdf')
    #  plt.close('all')
    #  #--------------------------------------------- 


    #--------------------------------------------- 
    # WRITE THE USERR SUBROUTINE FILE (CALLED BY RXN) 
    #--------------------------------------------- 
    # Only write userr input for rxn if there is thermally active
    # (i.e., time-varying) regolith adsorption
    if regolith_type['active'] is True:
        print('writing userr_data.dat... ')
        # Methane adsorption parameters based on Gough2010
        # (already defined earlier)
        # Augment the bc_dict dictionary with distribution coefficients for each zone
        #---------
        # TOP
        #---------
        kname = 'top'
        temps = bc_dict[kname]['vals'][1] - T_upshift + 273.15
        bc_dict[kname]['distcoeffs'] = K_eq(temps, molar_mass_ch4(), deltaE_obs)
        #---------
        # SUBSURFACE ACTIVE REGOLITH NODES 
        #   => z = -1m-regolith_base
        #---------
        counter = 0
        for y in T_depths[:-1]:
            counter+=1
            # Keys align with bc_dict
            kname = 'depth'+str(counter)
            #  temp_col = np.where(T_depths == y)[0][0]+1
            times_sub = bc_dict[kname]['times']
            temp_sub  = bc_dict[kname]['vals'] - T_upshift + 273.15
            #----------------------
            # Extend the timeseries cyclically since these are only 1 MY
            #----------------------
            dt = times_sub[1]-times_sub[0]
            timenew = times_sub
            tempnew = temp_sub
            # approx number of times to repeat cycle
            n_ts = math.ceil(p['sim_time']/times_sub[-1])
            # Add on to the time and temperature series cyclically
            for k in range(n_ts):
                #  if k==0: print(timenew[-1])
                tadd = times_sub+(timenew[-1]-timenew[0]+dt)
                #  if k==0: print(tadd[0])
                tempadd = temp_sub
                timenew = np.concatenate((timenew,tadd))
                tempnew = np.concatenate((tempnew,tempadd))
            # Cut down to be the same length as top boundary condition
            num_surf = len(bc_dict['top']['distcoeffs'])
            timenew = timenew[:num_surf]
            tempnew = tempnew[:num_surf]
            # Calculate new (extended) distcoeffs  and add to main BC directory
            bc_dict[kname]['distcoeffs'] = K_eq(tempnew, molar_mass_ch4(), deltaE_obs)
        # Write the dictionary entries to lists to be read into write_userr() function 
        distlist = []; zonelist=[]
        # Only want distcoeffs for top and active regolith zones
        #  keylist = list(bc_dict.keys()); keylist.remove('contam_zone')
        keylist = list(bc_dict.keys());
        [keylist.remove(item) for item in ['contam_zone','constant_regolith']]
        for key in keylist:
            #  print(key)
            distlist.append(bc_dict[key]['distcoeffs'])
            zonelist.append(bc_dict[key]['loc'])
        #---- Write the userr_data.dat file 
        utils.write_userr(join(macro_dir,'userr_data.dat'), bc_dict['top']['times'], distlist, 'distcoeffs',zonelist)
        print('done.')
    else:
        print('not writing a userr file because no active regolith adsorption present...')

    #--------------------------------------------- 
    # COPY AND/OR SYMLINK SOME IMPORTANT FILES 
    #--------------------------------------------- 
    #  symlink(join(tpl_dir,'perm.dat'),join(macro_dir,'perm.dat'),overwrite=True)
    #--Write the rock.dat file from tpl
    pest_io.tpl_write(p,join(tpl_dir,'rock.tpl'),join(macro_dir,'rock.dat'))
    #  symlink(join(tpl_dir,p['boun_file']),join(macro_dir,p['boun_file']),overwrite=True)
    #  symlink(join(tpl_dir,p['topzone_file']),join(macro_dir,os.path.split(p['topzone_file'])[1]),overwrite=True)
    #  symlink(join(tpl_dir,'geo-grad_air-static.boun'),join(macro_dir,'geo-grad_air-static.boun'),overwrite=True)
    #  symlink(join(tpl_dir,'eos_mars_co2.dat'),join(macro_dir,'eos_mars_co2.dat'),overwrite=True)
    #  #Copy eos file rather than symbolic link – symlinks can cause traceability issues
    #  shutil.copyfile(join('/lclscratch/jportiz/projects/gas_seepage/mars/runs/2d_fracture_network_runs/regolith_adsorption','eos_mars_co2.dat'),join(macro_dir,'eos_mars_co2.dat')) #newEOSmars
    #  shutil.copyfile(join('/lclscratch/jportiz/projects/gas_seepage/mars/runs/2d_fracture_network_runs/newEOSmars','eos_mars_co2.dat'),join(macro_dir,'eos_mars_co2.dat')) #newEOSmars
    #  symlink(join('/lclscratch/jportiz/projects/gas_seepage/mars/runs/2d_fracture_network_runs/newEOSmars','eos_mars_co2.dat'),join(macro_dir,'eos_mars_co2.dat'),overwrite=True) #newEOSmars

    #-------------------------------------------------- 
    ######### INITIALIZATION/SPINUP #################
    #  - NOTE: only initialize for noSorp simulations, since adsorption runs
    #  will restart from a noSorp_tracer simulation (for speed's sake)
    #  - if a ``noSorp_dir`` is provided manually in ``run.py`` for a noSorp
    #  simulation to restart from, have the code restart from that instead
    #-------------------------------------------------- 
    with open('v.pkl','wb') as f:
        pickle.dump([sourcetype,fehm_upshift,T_upshift],f)
        #  pickle.dump([fehm_upshift,T_upshift],f)

    #  if p['no_adsorption'] is True:
    # Don't do init0 if a ``noSorp_dir`` is specified manually in ``run.py``
    if p['noSorp_dir'] is not None:
        print()
        print('************************ WARNING *******************************')
        print('User specified {} as the noSorp_dir to restart from.'.format(p['noSorp_dir']))
        print()
        proceed = input('Do you want to proceed using {} as a starting point instead of the normal\ninit0 --> noSorp_tracer --> sorp/noSorp_restart progression? (y/n)'.format(p['noSorp_dir']))
        if proceed.lower() == 'y': pass  #skips to tracer restart simulation
        else:
            print('Ok. Going to do the normal init0 --> noSorp_tracer --> sorp/noSorp_restart progression instead.')


    # If init_dir not provided, create it:
    # Run pre-initialization first with uniform properties
    else:
        if init0_dir is None:
            print('Runing init0 simulation...')
            utils.run_fehm(p,'init0',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            print('    Done.')
            # Specify init0_dir
            #p['init0_dir'] = join(os.getcwd(),'init0')
            p['init0_dir'] = '../init0'
        else:
            p['init0_dir'] = init0_dir
        #  # Now do a 'spinup'  run
        #  # (For now, spinup may  be unnecessary since main sim is so long)
        #  if spinup_dir is None:
            #  utils.run_fehm(p,'spinup',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            #  p['spinup_dir'] = '../spinup'
        #  else:
            #  p['spinup_dir'] = spinup_dir

    #-------------------------------------------------- 
    #       TRANSPORT SIMULATIONS (EITHER noSorp OR sorp)
    #-------------------------------------------------- 

    # If a ``noSorp_dir`` was specified manually, do a _restart from it
    if p['noSorp_dir'] is not None:
        #  Make the noSorp_dir a full path back to main root dir (up 2 dir for this case)
        p['noSorp_dir'] = '../../{}'.format(p['noSorp_dir'])
        # Need init0_dir for *.storo file
        p['init0_dir'] = p['noSorp_dir']+'/init0'

        # NO ADSORPTION
        if p['no_adsorption'] is True:

            print('Running noSorp_restart_sim...')
            utils.run_fehm(p,'noSorp_restart',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            print('    Done.')
        # ADSORPTION
        else:
            utils.run_fehm(p,'sorp_tracer_restart',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)

    # Otherwise, do the Normal progression  ----------------
    # ---- _TRACER
    else:
        # Run noSorp_tracer regardless of adsorption mode
        # Do the noSorp_tracer sim before noSorp_restart or sorp_tracer_restart
        print('Running noSorp_tracer sim...')
        utils.run_fehm(p,'noSorp_tracer',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
        print('    Done.')
        # Need to assign the noSorp_dir variable (doesn't already exist)
        # Make the noSorp_dir a full path (only up 1 dir for normal progression)
        #  p['noSorp_dir'] = '../{}'.format(p['noSorp_dir'])  #doesn't work if noSorp_dir starts as None
        p['noSorp_dir'] = '..'  #the tpl is %noSorp_dir%/noSorp_tracer/run.fin

        # Now do the sorp_tracer_restart / noSorp_restart sim ----------
    # ---- _RESTART
        # NO ADSORPTION
        if p['no_adsorption'] is True:
            print('Running noSorp_restart_sim...')
            utils.run_fehm(p,'noSorp_restart',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            print('    Done.')
        # ADSORPTION
        else:
            print('Running sorp_tracer_restart sim...')
            utils.run_fehm(p,'sorp_tracer_restart',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            print('    Done.')

        #--- Eventually, set up a way to automatically run noSorp_restart simulation
        #    after this one concludes (remember to change the fracture aperture if necessary
        # NOTE: the below statement may not be true anymore
        #  print()
        #  print('************************ REMINDER ******************************')
        #  print('Because the ``no_adsorption`` parameter was set to "True" for')
        #  print('this simulation, this acts as the base run from which your ')
        #  print('other simulations will be restarted from (``noSorp_tracer``).')
        #  print()
        #  print('Currently, the followup noSorp run (e.g., ``noSorp_restart_b1mm``)')
        #  print('needs to be run manually (i.e., not using ``run.py``).')
        #  print()
        #  print('At this point, you can run your adsorption simulations using this')
        #  print('base (you can use ``run.py`` for relevant adsorption sims).')
        #  print('For example, you can use ``run.py`` in ``shallow50_base_case`` if')
        #  print('you have already run ``shallow50_noSorp_case``.')
        #  print('****************************************************************')
        #  #-------------------------------------------------- 

    #  #-------------------------------------------------- 
    #  ######### SORP_TRACER #################
    #  #-------------------------------------------------- 
    #  # - NOTE: All adsorption sims should restart from the noSorp_tracer simulation 
    #  # - It is up to the user to make sure that run has the right properties
    #  #-------------------------------------------------- 
    #  # ADSORPTION RXN ON (sorp)
    #  # - NOTE: All adsorption sims should restart from a noSorp_tracer simulation 
    #  else:
        #  # If a ``noSorp_dir`` was specified manually
        #  if p['noSorp_dir'] is not None:
            #  utils.run_fehm(p,'sorp_tracer_restart',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
        #  # Normal case 
        #  else:
            #  #  utils.run_fehm(p,'sorp_tracer',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            #  utils.run_fehm(p,'sorp_tracer_restart',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)


