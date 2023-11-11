import os,sys
from os.path import join
import numpy as np
from subprocess import call
sys.path.append('/project/gas_seepage/jportiz/scripts')
from tools import find_nearest, sci_notation
from nilson import *
sys.path.append('/project/gas_seepage/jportiz/scripts/gasmigra')
import utils
from mars_tools import solsec,days2sols,sols2days,mars_seasons
import tools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib import cbook
import scipy.integrate as integrate
from scipy import interpolate
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
#the below import allows you to outline text
import matplotlib.patheffects as path_effects
import pickle
import pandas as pd
import scipy.fft
from scipy.signal import blackman
from scipy import interpolate
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

#  def extract_Tsub(T_upshift=0., T_ref=-50.):
    #  '''
    #  Extract subsurface temperatures vs depth to then use as thermal boundary
    #  condition for future simulations (e.g., regolith adsorption).
#  
    #  :param T_upshift: float
        #  Value added to real Mars temperatures to allow FEHM to operate with a
        #  liquid water formulation (so that it is between 0-100 °C) [°C].
    #  :param T_ref: float
        #  Surface reference temperature used to calculate the EOS properties of
        #  CO2 for FEHM. Essentially the mean surface temperature at the location
        #  of interest [°C].
    #  '''
    #  cwd = os.getcwd()
    #  results_dir = join(cwd,'thermal_spinup')
    #  #-------------------------------------------------- 
    #  #           MESH DATA 
    #  #-------------------------------------------------- 
    #  # get name/path of meshfile used
    #  with open(join(results_dir,'fehmn.files'),'r') as ff:
        #  i=0
        #  for line in ff.readlines():
            #  if line.startswith('grida'):
                #  i+=1
                #  os.chdir(results_dir); wd = os.getcwd()
                #  gf = line.split(": ",1)[1].rstrip()
                #  os.chdir(os.path.split(gf)[0]); meshdir = os.getcwd()
                #  meshfile = join(meshdir,os.path.split(gf)[-1])
                #  os.chdir(cwd)
                #  if i==1: break
    #  grid = utils.grid_recarray(meshfile)
    #  fracnodes = np.where(grid['x']==0.) #left side
    #  fracdepths = grid['y'][fracnodes]   # ``   ``  
    #  mesh_dir1 = os.path.split(os.path.split(os.path.split(meshfile)[0])[0])[1]
    #  mesh_dir2 = os.path.split(os.path.split(meshfile)[0])[1]
    #  mesh_name = mesh_dir1+'/'+mesh_dir2
    #  mesh_name = mesh_name.replace('_','\_')
    #  #-------------------------------------------------- 
    #  #           FEHM  DATA 
    #  #-------------------------------------------------- 
    #  tf = os.path.join(results_dir, 'run_temp.his')
    #  #---- Temperatures 
    #  tdata = np.loadtxt(tf,skiprows=4)
    #  #  t_time = days2sols(tdata[:,0]) #convert days to sols
    #  t_cols = np.loadtxt(tf,skiprows=3,max_rows=1,dtype='str')[1::2]
    #  #  pdata = (pdata - fehm_upshift) #[Pa] remove the artificial pressure upshift from FEHM
    #  # Get depths of all nodes
    #  nodes =[]; ys = []
    #  for i in t_cols[1:]:
        #  nodes.append(int(i)-1)
        #  ys.append(grid[int(i)-1]['y'])
    #  #-------------------------------------------------- 
    #  #     TIME SERIES OF TEMPS AT EACH DEPTH 
    #  #-------------------------------------------------- 
    #  # Save columns: TIME  -200m ...  -2m  -1m  0m
    #  #               ----  -----      ---  ---  --- 
    #  for i in enumerate(ys):
        #  j = i[0]
        #  columns = np.columnstack( (tdata[:,0], t_data[) )

def plot_geothermal_gradient(geotherm, T_upshift=0., T_ref=-50.):
    '''
    Plot the theoretical geothermal gradient vs the output from FEHM.

    :param geotherm: float
        Prescribed geothermal gradient [°C/m].
    :param T_upshift: float
        Value added to real Mars temperatures to allow FEHM to operate with a
        liquid water formulation (so that it is between 0-100 °C) [°C].
    :param T_ref: float
        Surface reference temperature used to calculate the EOS properties of
        CO2 for FEHM. Essentially the mean surface temperature at the location
        of interest [°C].
    '''
    cwd = os.getcwd()
    results_dir = join(cwd,'geo_gradient')
    #-------------------------------------------------- 
    #           MESH DATA 
    #-------------------------------------------------- 
    # get name/path of meshfile used
    with open(join(results_dir,'fehmn.files'),'r') as ff:
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
    tf = os.path.join(results_dir, 'run_temp.his')
    #---- Temperatures 
    tdata = np.loadtxt(tf,skiprows=4)
    #  t_time = days2sols(tdata[:,0]) #convert days to sols
    t_cols = np.loadtxt(tf,skiprows=3,max_rows=1,dtype='str')[1::2]
    #  pdata = (pdata - fehm_upshift) #[Pa] remove the artificial pressure upshift from FEHM
    # Get depths of all nodes
    nodes =[]; ys = []
    for i in t_cols[1:]:
        nodes.append(int(i)-1)
        ys.append(grid[int(i)-1]['y'])
    #-------------------------------------------------- 
    #       CALCULATE THEORETICAL TEMP PROFILE 
    #-------------------------------------------------- 
    t_predict = abs(np.asarray(ys))*geotherm + T_ref
    #-------------------------------------------------- 
    #       PRINT RESULTS TO SCREEN 
    #-------------------------------------------------- 
    print("------------")
    print('Theoretical temperature gradient = {:.4f} °C/m'.format(geotherm))
    print('Model temperature gradient       = {:.4f} °C/m\n'.format(abs((max(tdata[-1,1:])-min(tdata[-1,1:]))/(ys[0]-ys[-1]))))
    #-------------------------------------------------- 
    #       PLOT GEOTHERMAL PROFILE 
    #-------------------------------------------------- 
    fig, ax = plt.subplots(1,figsize=(6,8))
    ax.plot(t_predict, abs(np.asarray(ys)), ls='-',color='k',label='theoretical')
    ax.plot(tdata[-1,1:]-T_upshift, abs(np.asarray(ys)), marker='o',mfc='None',ls='',color='darkorange',label='FEHM')
    ax.text(0.98,0.80,r'geothermal gradient: {}$^\circ$C/m'.format(geotherm),transform=ax.transAxes,ha='right')
    ax.legend()
    #---- PROPERTIES
    ax.invert_yaxis()
    ax.set_xlabel(r'Temperature [$^{\circ}$C]')
    ax.set_ylabel('Depth [m]')
    fig.savefig('output/geothermal_profile.pdf')
    plt.close('all')

def plot_timeseries(T_0,T_upshift):
    '''
    Plot timeseries of subsurface temperatures.

    Parameters
    ----------
    T_0 : float
        Surface temperature [°C].
    T_upshift : float
        Value added to real Mars temperatures to allow FEHM to operate with a
        liquid water formulation (so that it is between 0-100 °C) [°C].
    '''
    cwd = os.getcwd()
    results_dir = join(cwd,'thermal_spinup')
    #-------------------------------------------------- 
    #           MESH DATA 
    #-------------------------------------------------- 
    # get name/path of meshfile used
    with open(join(results_dir,'fehmn.files'),'r') as ff:
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
    tf = os.path.join(results_dir, 'run_temp.his')
    #---- Temperatures 
    tdata = np.loadtxt(tf,skiprows=4)
    tdata[:,0] = days2sols(tdata[:,0]) #convert days to sols
    t_cols = np.loadtxt(tf,skiprows=3,max_rows=1,dtype='str')[1::2]
    # Get depths of all nodes
    nodes =[]; ys = []
    for i in t_cols[1:]:
        nodes.append(int(i)-1)
        ys.append(grid[int(i)-1]['y'])

    #-------------------------------------------------- 
    #       PLOT SUBSURFACE TEMPERATURE TIME SERIES 
    #-------------------------------------------------- 
    #  depths = [0.0, 1.0, 5.0, 10.0, 30.0]
    #  nth=10 #plot every nth point
    #  fig, ax = plt.subplots(1)
    #  for depth in depths:
        #  color = next(ax._get_lines.prop_cycler)['color']
        #  #Plot FEHM 
        #  nodes_depth = np.where(grid['y']== -depth)[0] #nodes at desired depth
        #  for c in enumerate(t_cols):
            #  try:int(c[1])
            #  except: continue
            #  if int(c[1]) in nodes_depth:
                #  T_fehm = tdata[:,c[0]]
#  
        #  #FIX THIS EVENTUALLY !!
        #  ax.plot(tdata[:,0][::nth],T_fehm[::nth]-T_upshift, ls='-',color=color, label='{} m'.format(depth))
    depths = [0.0, 1.0, 5.0, 10.0, 30.0]
    #  depths = [0.0,0.5, 1.0, 5.0, 10.0, 30.0]
    nth=10 #plot every nth point
    fig, (ax1, ax2) = plt.subplots(2)
    for depth in depths:
        color = next(ax1._get_lines.prop_cycler)['color']
        #Plot FEHM 
        nodes_depth = np.where(grid['y']== -depth)[0] #nodes at desired depth
        for c in enumerate(t_cols):
            try:int(c[1])
            except: continue
            if int(c[1]) in nodes_depth:
                T_fehm = tdata[:,c[0]]

        #FIX THIS EVENTUALLY !!
        ax1.plot(tdata[:,0][::nth],T_fehm[::nth]-T_upshift, ls='-',color=color, label='{} m'.format(depth))
        #  ax2.plot(tdata[:,0][::nth],T_fehm[::nth]-T_upshift, ls='-',color=color, label='{} m'.format(depth))
        ax2.plot(tdata[:,0],T_fehm-T_upshift, ls='-',color=color, label='{} m'.format(depth))

        #  #---- Print some useful statistics
        #  print("Depth: {} m".format(depth))
        #  print("------------")
        #  print("max T_analytical = {:.4f} °C".format(max(T_s)))
        #  print("max T_FEHM       = {:.4f} °C".format(max(T_fehm)-T_upshift))
        #  print("min T_analytical = {:.4f} °C".format(min(T_s)))
        #  print("min T_FEHM       = {:.4f} °C".format(min(T_fehm)-T_upshift))
    #  ax.text(0.5,0.96,'Period = {} d'.format(period/86400.),transform=ax.transAxes,ha='center')
    #  ax.set_xlim(left=0., right=max(tdata[:,0][::nth]))
    #  ax.set_xlim(left=0., right=2714.)
    #  ax.set_xlim(left=0., right=668.)
    #  ax.set_xlim(left=tdata[-1,0]-668.,right=tdata[-1,0])
    #  ax.set_xlim(left=668.*74.,right=668.*75.)3
    ax1.set_xlim(left=668.*60.,right=668.*62.)
    #  ax2.set_xlim(left=668.*0.,right=668.*2.)
    ax2.set_xlim(left=668.*60.,right=668.*60.+5)
    ax1.legend(loc='lower left')
    ax2.set_xlabel('Time [sols]')
    ax1.set_ylabel(r'Temperature [$^\circ$ C]')
    ax2.set_ylabel(r'Temperature [$^\circ$ C]')
    fig.tight_layout()
    #  name = 'compare_subsurface_temps_'+signal+'.pdf'
    #  name = 'test3.pdf'
    name = 'output/subsurface_temp_timeseries.pdf'
    plt.savefig(name)
    plt.close('all')

def plot_ranges(T_0,T_upshift):
    '''
    Plot ranges of subsurface temperatures.

    Parameters
    ----------
    T_0 : float
        Surface temperature [°C].
    T_upshift : float
        Value added to real Mars temperatures to allow FEHM to operate with a
        liquid water formulation (so that it is between 0-100 °C) [°C].
    '''
    cwd = os.getcwd()
    results_dir = join(cwd,'thermal_spinup')
    #-------------------------------------------------- 
    #           MESH DATA 
    #-------------------------------------------------- 
    # get name/path of meshfile used
    with open(join(results_dir,'fehmn.files'),'r') as ff:
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

    #  surface_temp = 0.0
    #  geotherm = 0.0045 #[°C/m]

    #-------------------------------------------------- 
    #           FEHM  DATA 
    #-------------------------------------------------- 
    tf = os.path.join(results_dir, 'run_temp.his')
    #---- Temperatures 
    tdata = np.loadtxt(tf,skiprows=4)
    tdata[:,0] = days2sols(tdata[:,0]) #convert days to sols
    t_cols = np.loadtxt(tf,skiprows=3,max_rows=1,dtype='str')[1::2]
    # Get depths of all nodes
    nodes =[]; ys = []
    for i in t_cols[1:]:
        nodes.append(int(i)-1)
        ys.append(grid[int(i)-1]['y'])

    #-------------------------------------------------- 
    #       PLOT SUBSURFACE TEMPERATURE RANGES BY DEPTH 
    #           - similar to Fig. 4 in Jones2011
    #-------------------------------------------------- 
    #  depths = [0.0, 1.0, 5.0, 10.0, 30.0]
    #  depths = np.arange(0.0, max(abs(grid['y'])))
    depths = list(np.arange(0.0, 30.+1.))
    if max(abs(grid['y'])) not in depths: depths.append(max(abs(grid['y'])))
    T_ranges = {}
    nth=10 #plot every nth point
    for depth in depths:
        color = 'k'
        #  color = next(ax._get_lines.prop_cycler)['color']
        # Get the temperature range at each depth (put in dictionary)
        nodes_depth = np.where(grid['y']== -depth)[0] #nodes at desired depth
        for c in enumerate(t_cols):
            try:int(c[1])
            except: continue
            if int(c[1]) in nodes_depth:
                minT = min(tdata[:,c[0]])-T_upshift
                maxT = max(tdata[:,c[0]])-T_upshift
                T_ranges[depth] = [minT, maxT]
        #  ax.plot(T_ranges[depth], [depth]*2, ls='',marker='o',mfc='None',color=color, label='{} m'.format(depth)) #

    fig, ax = plt.subplots(1)
    # Plot the min and max as separate lines
    mins = []; maxs = []
    for i in T_ranges.values():mins.append(min(i))
    for i in T_ranges.values():maxs.append(max(i))
    # min(T) 
    ax.plot(mins, list(T_ranges.keys()), ls='-',color=color)
    # max(T) 
    ax.plot(maxs, list(T_ranges.keys()), ls='-',color=color)

    ax.set_yscale('symlog') #'symlog' allows 0's to plot 
    ax.set_ylim(0,depths[-1])
    ax.invert_yaxis()
    ax.set_xlabel(r'Temperature [$^\circ$ C]')
    ax.set_ylabel('Depth [m]')
    fig.tight_layout()
    name = 'output/subsurface_temp_ranges_log.pdf'
    plt.savefig(name)
    plt.close('all')




#----------------------------------------  
#TESTING
#----------------------------------------  
depths = [0.,0.5,1.]
period = 1.
T_0 = 0.
amplitude = 10.
thermal_cond = 2.5
rho_rock = 2900.
C_p = 800.

if __name__ == "__main__":

    #TESTING
    T_0 = 0.
    dT = 40.
    thermal_cond = 2.5
    rho_rock = 2900.
    C_p = 800.

    ys = [0.,0.1, 0.5, 1.0, 5, 10]
    #Plot temps at several depths
    fig, (ax1,ax2) = plt.subplots(2)
    # Diurnal Period
    period = 86400. #[s]
    delta_t = 100. #[s]
    ts = np.arange(0.,period+delta_t, delta_t)
    for y in ys:
        T_s = T_sub_analytical(y,ts,T_0,dT,period,thermal_cond,rho_rock,C_p)
        ax1.plot(ts/86400.,T_s,label='depth = {} m'.format(y))
    ax1.text(0.5,0.90,'Period = {} d'.format(period/86400.),transform=ax1.transAxes,ha='center')
    # Annual Period
    period = 86400.*365. #[s]
    delta_t = 100. #[s]
    ts = np.arange(0.,period+delta_t, delta_t)
    for y in ys:
        T_s = T_sub_analytical(y,ts,T_0,dT,period,thermal_cond,rho_rock,C_p)
        ax2.plot(ts/86400.,T_s,label='depth = {} m'.format(y))
    ax2.text(0.5,0.90,'Period = {} d'.format(period/86400.),transform=ax2.transAxes,ha='center')
    ax1.legend()
    ax1.set_xlabel('Time [d]')
    ax2.set_xlabel('Time [d]')
    ax1.set_ylabel(r'$\Delta T$ [$^\circ$ C]')
    ax2.set_ylabel(r'$\Delta T$ [$^\circ$ C]')
    fig.tight_layout()
    name = 'analytical_subsurface_temps.pdf'
    plt.savefig(name)
    plt.close('all')
    try: call(['evince', name])
    except: pass


