import os,sys
from os.path import join
sys.path.append('/project/gas_seepage/jportiz/scripts/')
from gasmigra import utils
import mars_tools
from mars_tools import sols2days, days2sols
from tools import symlink, sine_signal, cosine_signal, find_nearest
from mesh_tools import mesh_files
from glob import glob
import baro_decomp.convert_boun as convert_boun
import pandas as pd
#import ipdb
import pickle
import numpy as np
import numpy
from matk import pest_io

try:
    RUNDIR = os.environ["RUNDIR"]
except KeyError:
    print()
    print("Set the RUNDIR evironment variable first, e.g.:")
    print("    export RUNDIR='<FEHM-MARS_REPO_DIRECTORY>/example/' ")
    print()
try:
    FEHM_MARS_SRC = os.environ["FEHM_MARS_SRC"]
except KeyError:
    print()
    print("Set the FEHM_MARS_SRC evironment variable first, e.g.:")
    print("    export FEHM_MARS_SRC='<FEHM-MARS_REPO_DIRECTORY>/src/' ")
    print()


#  main_mesh_dir = '/project/gas_seepage/jportiz/mars/mesh/1d/'
#  tpl_dir = '/project/gas_seepage/jportiz/mars/runs/tpls/'
#  exe = '/home/jportiz/software/FEHM/src/xfehm_v3.3.1_mars'
main_mesh_dir = RUNDIR+'mesh/1d/'
tpl_dir = RUNDIR+'/runs/tpls/'
exe = FEHM_MARS_SRC+'/xfehm_v3.3.1_mars'

#DIRECTORY OF POTENTIAL MESH FILES 
mesh_dict = {
            1:  {'mesh_file':main_mesh_dir+'homogeneous_1x200m/grid.inp',
                },
            }

# Check if mesh name is too long for fehmn.files
meshname =  mesh_dict['depth200_fracDen035']['mesh_file']
print(meshname)
shortened=False
if len(meshname) > 100:
    shortened=True



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
    print_all = numpy.concatenate([gett,get])
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
        #  (used to be every 10th node).
        #  -``horz`` gives a horizontal line of nodes on the top of domain.
        #  -``both`` gives both (vertical then horizontal)
    #  '''
    #  # Make sure macro_dir exists
    #  if not os.path.isdir(macro_dir): os.mkdir(macro_dir)
    #  fh = open(os.path.join(macro_dir,'nodeline.dat'),'w')
    #  # Vertical line of nodes (x=0)
    #  #  get = grid[numpy.where(grid['x']==0)[0][0:-1:10]]
    #  get = grid[numpy.where(grid['x']==0)[0]]
    #  #  get = get[0::10] #<-- if you want only every 10th node
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

def perm_macro(grid, k_bg, macro_dir='.',zones=[],zone_perms=[]):
    '''
    Write perm macro file for FEHM (no gdkm).

    Parameters
    ----------
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
    fh = open(os.path.join(macro_dir,'perm.dat'),'w')
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

#  #TESTING
#  macro_dir = 'macros'
#  p = {'mesh':5,
     #  'sim_time':10.,
    #  }

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
    if 'df' not in p: p['df'] = 0.001 #[m]
    #  if 'methane_prod_rate' not in p: p['methane_prod_rate'] = 0. #initial tracer source (no source strength)
    #  # Tracer start time, used to shift starting time along baro pressure signal
    #  #if 'tracer_start' not in p: p['tracer_start'] = 0.0
    #  # Perm cutoff for gdkm nodes
    #  if 'k_gdkm_cutoff' not in p: p['k_gdkm_cutoff'] = 1.e-12
    #  if 'sim_time' not in p: p['sim_time'] = mars_tools.sols2days(2712.41) #sols-->days (extent of the MSL baro record)
    if 'phi_m' not in p: p['phi_m'] = 0.350 #default matrix porosity
    if 'phi_f' not in p: p['phi_f'] = 0.999 #default fracture porosity (upscaled if needed) 
    if 'k_m' not in p: p['k_m'] = 1.e-14 #default matrix perm
    if 'geothermal_gradient' not in p: p['geothermal_gradient'] = 0.0045 #default [°C/m]
    #  if 'rho_rock' not in p: p['rho_rock'] = 2900. #not used anywhere yet

    p['mesh_file'] = mesh_dict[p['mesh']]['mesh_file']
    #  p['stor_file'] = mesh_dict[p['mesh']]['stor_file']
    #  p['avs_file'] = mesh_dict[p['mesh']]['avs_file']
    #  p['matzone_file'] = mesh_dict[p['mesh']]['matzone_file']
    #  p['matzonn_file'] = mesh_dict[p['mesh']]['matzonn_file']
    #  p['outsidezone_file'] = mesh_dict[p['mesh']]['outsidezone_file']
    #  p['topzone_file'] = mesh_dict[p['mesh']]['topzone_file']
    #-------------------------------------------------- 
    # Throw a warning if exe does not have 'mars' in it
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
    #---- Write different perm.dat and nodeline.dat macros depending on mesh 
    if '2d_planar' in p['mesh_file']:
        zone_perms = [k_f]
        node_macro(g,macro_dir=macro_dir,dirxn='vert')
        #  node_macro(g,macro_dir=macro_dir,dirxn='vert')
    elif '1d' in p['mesh_file']:
        zones=[]        #no non-background zones
        zone_perms = []
        node_macro(g,macro_dir=macro_dir,dirxn='vert',spacing=1.0)
    perm_macro(g, k_bg=p['k_m'], macro_dir=macro_dir,zones=zones,zone_perms=zone_perms)

    #  # Write the fracture zone [-9]
    #  #-- All the leftmost nodes are the fracture
    #  fracnodes = g['node'][g['x']==0.0]
    #  utils.write_zone(join(macro_dir,'fracture.zone'),fracnodes,[-9])
    #  mesh_files.zone_to_zonn(join(macro_dir,'fracture.zone'))

    # Write topzone file [-100]
    with open(join(macro_dir,'top.zone'),'w+') as fz:
        fz.write('zone\n100\n')
        fz.write('-0.1 1.e4 1.e4 -0.1\n') #width coords
        fz.write('-0.1 -0.1  0.1  0.1\n\n') #up/down coords
        fz.write('stop')
    mesh_files.zone_to_zonn(join(macro_dir,'top.zone'))

    # Write bottom zone file [-200]
    botnodes = g['node'][g['y']==min(g['y'])]
    utils.write_zone(join(macro_dir,'bottom.zone'),botnodes,[200])
    mesh_files.zone_to_zonn(join(macro_dir,'bottom.zone'))

    #  rho_air = mars_tools.density_lookup(-4500.) #[kg/m^3] #check that this the same as in EOS macro

    #--------------------------------------------- 
    # WRITE THE EOS MACRO FOR CO2 "AIR" PROPERTIES
    #--------------------------------------------- 
    #---- (1) Get the pressures for the surface -------------------------------
    datadir='/project/gas_seepage/jportiz/mars/data/pressure/rems/processed_data'
    if mars_tools.days2sols(p['sim_time'])<=2712.41: #[sols] (extent of actual pressure record)
        pres_df = pd.read_csv(join(datadir,'pchip_press_elev_correction.csv'))
    else:
        #  pres_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_press_50.csv'))
        pres_df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_press_50_Ls.csv'))
    actual_press_soltime_start = pres_df['SOL_TIME'].iat[0]       #[sols]
    actual_press_Ls_start = pres_df['L_s'].iat[0]                 #[°]
    ts = np.array(pres_df['SOL_TIME']-pres_df['SOL_TIME'].iat[0]) #[sols] set first time to 'zero'
    ts = mars_tools.sols2days(ts) #[d] convert sols to days for FEHM boun file
    ps = np.array(pres_df['PRESS_ADJUSTED']/1e6) #[MPa]
    #  #--Shift the pressures up to 100x the initial pressure value for FEHM
    #  fehm_upshift = ps[0]*100. - ps[0] #MPa
    #--Shift the pressures up to 1000x the initial pressure value for FEHM 
    fehm_upshift = ps[0]*1000. - ps[0] #MPa
    ps = ps + fehm_upshift #this shifts up to FEHM range, but doesn't distort amplitudes
    print('writing the eos macro...')
    #  P_ref = 0.0006      # [MPa] equivalent to 600 Pa  (for some reason, setting P_ref to 0.0008 makes negative densities...)
    P_ref = 0.0007     # [MPa] equivalent to 700 Pa  (for some reason, setting P_ref to 0.0008 makes negative densities...)
    T_ref_eos = -50   # [°C]  the EOS T_ref can be -50, but T_ref in init should be mean surf. temp. @ Gale 
    #  P_upshift = P_ref * 100 - P_ref #shift pressures up to FEHM range without distorting amplitudes
    #  P_upshift = P_ref * 100 #shift pressures up to FEHM range without distorting amplitudes
    pupshift_factor = 1000 #shift pressure up to FEHM range without distorting amplitudes (for eos macro) 
    P_upshift = P_ref * pupshift_factor
    T_upshift = 95.    # [°C] trick FEHM to think it's above 0° (and keep T between 0-100°C)
    co2_props_file = '/project/gas_seepage/jportiz/mars/data/co2_properties/co2_isothermal_allTemps.csv'
    d_co2 = pd.read_csv(co2_props_file)
    # Calculate and set the eos parameters
    # EW1-3 -----------------------------
    # (EW1=P_ref; EW2=T_ref; EW3=rho_ref)
    #  EW1 = P_upshift
    # Make the approximation that P_ref is same as starting Pressure
    press_surf = ps[::4]
    EW1 = press_surf[0]
    EW2 = T_ref_eos+T_upshift
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
    EW7 = (H1-H0)/(P1-P0)
    # EW8 --------------------------------
    # (EW8=dE/dT)
    T0 = T_ref_eos+273
    index = min(d_co2.loc[(d_co2['Temperature (K)']>T0)&(round(d_co2['Pressure (MPa)'],5)==P_ref)].index)
    T1 = d_co2.loc[index]['Temperature (K)']
    H0 = EW6
    H1 = d_co2.loc[index]['Enthalpy (kJ/kg)']/1000 # [MJ/kg]
    EW8 = (H1-H0)/(T1-T0)
    # EW9 --------------------------------
    # (EW9=mu_ref)
    EW9 = float(d_co2.loc[(d_co2['Temperature (K)']==T_ref_eos+273)&(round(d_co2['Pressure (MPa)'],5)==P_ref)]['Viscosity (Pa*s)']) # [Pa s]
    # EW10 -------------------------------
    # (EW10=d(mu)/dP)
    P0 = P_ref
    index = min(d_co2.loc[(d_co2['Pressure (MPa)']>P_ref)&(d_co2['Temperature (K)']==T_ref_eos+273)].index)
    P1 = d_co2.loc[index]['Pressure (MPa)']
    mu0 = EW9
    mu1 = d_co2.loc[index]['Viscosity (Pa*s)'] # [MJ/kg]
    EW10 = (mu1-mu0)/(P1-P0)
    # EW11 -------------------------------
    # (EW11=d(mu)/dT)
    T0 = T_ref_eos+273
    index = min(d_co2.loc[(d_co2['Temperature (K)']>T0)&(round(d_co2['Pressure (MPa)'],5)==P_ref)].index)
    T1 = d_co2.loc[index]['Temperature (K)']
    mu0 = EW6
    mu1 = d_co2.loc[index]['Enthalpy (kJ/kg)']/1000 # [MJ/kg]
    EW11 = (mu1-mu0)/(T1-P0)
    # EV1-11 -----------------------------
    EV1=EW1; EV2=EW2; EV3=1.0; EV4=0; EV5=0;EV6=-0.020773
    EV7=0.827; EV8=11.398; EV9=1.372e-4; EV10=0.; EV11=0.
    #---------------------------------------------------- 
    # Write the EOS file
    eos_filename ='eos_mars_co2.dat'
    mars_tools.write_eos(join(macro_dir,eos_filename),EW1,EW2,EW3,EW4,EW5,EW6,EW7,EW8,EW9,EW10,EW11,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11)
    print('done.')

    #--------------------------------------------- 
    # Write the air_static.boun file (w/ geothermal gradient)
    #     - Used for initialization of geothermal gradient
    #--------------------------------------------- 
    T_ref = p['T_ref']    #[°C] mean surf. temp. @ Gale (different than T_ref_eos)
    T_0 = T_ref+T_upshift #50.#0.#-40.                    #[°C]
    p['T_0'] = T_0
    y_bot = min(g['y'])
    T_bot = T_0 + p['geothermal_gradient']*abs(y_bot)#calculate T at bottom of domain
    p['T_bot'] = T_bot
    #  pw = [0.070618]*2
    pw = [0.70618]*2
    surf_temps = [T_0]*2
    bot_temps = [T_bot]*2
    #  outfile = 'air_static.boun'
    outfile = 'geo-grad_air-static.boun'
    utils.write_boun(join(macro_dir,outfile),([0.,1e20],[0.,1e20]),[[pw,surf_temps],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]])

    #  #  #--------------------------------------------- 
    #  #  # Write the main sinusoidal Temp boun file to be used 
    #  #  #--------------------------------------------- 
    #  signals={'diurnal': 1.0,
             #  'monthly': 30.0,
             #  'annual': 365.,
            #  }
    #  #--TEMPERATURE
    #  #  period = 1.0                  #[d] 
    #  period = signals[p['signal']]   #[d]
    #  ampl = 10.                    #[°C]
#  
    #  #  dt = 0.010                     #[d]
    #  dt = period/100.                  #[d]
    #  ts = np.arange(0., period+dt, dt) #[d] 
    #  #  ts = np.arange(0., 1.+dt, dt) #[d] 
#  
#  
    #  #  temps = sine_signal(period,ts,T_0,ampl)
    #  temps = cosine_signal(period,ts,T_0,ampl)
    #  ps = [0.070618]*len(temps)
    #  #write to file
    #  #  outfile = 'synth_temp_diurnal.boun'
    #  #  # TESTING...
    #  #  macro_dir = '/lclscratch/jportiz/projects/gas_seepage/mars/runs/fehm_mars_tests/heat_flow/conduction/macros'        #TESTING
    #  #  out_file = outfile
    #  #  ts = np.arange(0., 1.+dt, dt) #[d] 
    #  #  t = [ts,[0.,1e20]]
    #  #  #  ys = [[temps,ps],[0.225,0.225]]
    #  #  ys = [[temps,ps],[T_bot,T_bot]]
    #  #  step = ['cy_linear','ti']
    #  #  varname = [["t","pw"],["t"]]
    #  #  zones = [[-100],[-200]]
#  
    #  # Boun with geothermal gradient
    #  outfile = 'synth_temp_'+p['signal']+'.boun'
    #  utils.write_boun(join(macro_dir,outfile),t=[ts,[0.,1e20]],ys=[[temps,ps],[T_bot,T_bot]],step=["cy_linear","ti_linear"],varname=[["t","pw"],["t"]],zones=[[-100],[-200]])
    #  p['boun_file'] = join('..',macro_dir,outfile)
    #  #  # Boun without geothermal gradient (don't think this is necessary here)
    #  #  outfile = 'synth_temp_'+p['signal']+'_no_geo.boun'
    #  #  utils.write_boun(join(macro_dir,outfile),t=ts,ys=[temps,ps],step="cy_linear",varname=["t","pw"],zones=-100)
    #  #  p['boun_file_no_geo'] = join('..',macro_dir,outfile)

    #--------------------------------------------- 
    # Write the boun file with Mars temperature fluctuations
    #     - include geothermal gradient
    #     - assign constant surface pressure
    #--------------------------------------------- 
    datadir='/project/gas_seepage/jportiz/mars/data/pressure/rems/processed_data'
    if p['synthetic']==False:
        df = pd.read_csv(join(datadir,'stitched_timeseries','stitched_amb_temp_80_Ls.csv'))
        ts = np.array(df['SOL_TIME']-df['SOL_TIME'].iat[0]) #[sols] set first time to 'zero'
        ts = mars_tools.sols2days(ts) #[d] convert sols to days for FEHM boun file
        surf_temps = np.array(df['AMB_TEMP']) + T_upshift # [°C]
        #  surf_temps = np.array(df['AMB_TEMP']) + 45 # [°C]  #TEST
        pw = [0.070618]*len(surf_temps) # [MPa]
        #  # Temps at bottom of domain from geothermal gradient are calculated earlier
        #  y_bot = min(g['y'])
        #  T_bot = T_0 + p['geothermal_gradient']*abs(y_bot) # calculate T at bottom of domain
        #  p['T_bot'] = T_bot
        #  bot_temps = [T_bot]*2
        #  topzonen = -100 #zone to apply BC changes to
        outfile = join('..',macro_dir,'mars_temps.boun')
        p['boun_file'] = outfile

        utils.write_boun(join(macro_dir,outfile),(ts[::4],[0.,1e20]),[[pw[::4],surf_temps[::4]],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]]) #hourly BC changes
    #  utils.write_boun(join(macro_dir,outfile),(ts,[0.,1e20]),[[pw,surf_temps],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]]) #15-min BC changes

    elif p['synthetic']==True:
        df = pd.read_csv(join(datadir,'pchip_amb_temp.csv'))
        # Load synthetic temperature components generated previously
        synth_components_dir = '/project/gas_seepage/jportiz/mars/analytical/fourier_analysis/generate_synthetic_data/output'
        with open(join(synth_components_dir, 'vars_tempsFullSignal.pkl'), 'rb') as f:
            periods,amps,shifts = pickle.load(f)

        # Make time array [sols] for 1 Mars Year
        #  dt = np.diff(df['SOL_TIME'])[0]             #[sols]
        #  dt = 0.01                                   #[sols]
        # Convert sols to days for FEHM
        dt = 0.04                                      #[days]
        #  ts = np.arange(0., sols2days(668.)+dt, dt)     #[sols]  #1-year cyclic
        ts = np.arange(0., p['sim_time']+dt, dt)     #[sols]     #full record
        mean_temp = p['T_ref']                         #[°C]

        # Generate synthetic temperatures
        surf_temps = mars_tools.synth_fourier_series(periods+amps+shifts, mean_temp, days2sols(ts), num_modulations=1) + T_upshift
        # Constant surface pressure
        pw = [0.070618]*len(surf_temps) # [MPa]

        # Write boun (cyclic time seems to have issues around 41,183.42 days)
        outfile = join('..',macro_dir,'mars_temps.boun')
        p['boun_file'] = outfile
        #  utils.write_boun(join(macro_dir,outfile),(ts,[0.,1e20]),[[pw,surf_temps],bot_temps],["cy_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]]) #roughly hourly BC changes
        utils.write_boun(join(macro_dir,outfile),(ts,[0.,1e20]),[[pw,surf_temps],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]]) #roughly hourly BC changes
        #  utils.write_boun(join(macro_dir,outfile),(ts,[0.,1e20]),[[pw,surf_temps],bot_temps],["ti_linear","ti"],[["pw","t"],["t"]],zones=[[-100],[-200]]) #roughly hourly BC changes

    #  # add a second model for air injection into contam zone (if needed) 
    #  contamzonen = -8
    #  fh = open(join(macro_dir,outfile),'r')
    #  fl = fh.readlines();fh.close()
    #  emptys=[] #emptylines in original boun file
    #  for i in enumerate(fl):
        #  if i[1]=='\n': emptys.append(i[0])
    #  fl.insert(emptys[0],'model 2\nti_linear\n2\n0.0    1e10\ndsw\n{}    {}\n'.format(inj_rate,inj_rate))
    #  fl.insert(emptys[1]+1,' {} 0 0 2\n'.format(int(contamzonen)))
    #  fh = open(join(macro_dir,outfile),'w')
    #  fl=''.join(fl)
    #  fh.write(fl);fh.close()

    #  #  ts = mars_tools.sols2days(ts) #[d] convert sols to days for FEHM boun file
    #  ps = np.array(df['PRESS_ADJUSTED']/1e6) #[MPa]
    #  temps = [0.05]*len(ps) #[deg C]
    #  #--shift the pressures up to 100x the initial pressure value for FEHM
    #  fehm_upshift = ps[0]*100. - ps[0] #MPa
    #  ps = ps + fehm_upshift #this shifts up to FEHM range, but doesn't distort amplitudes
    #  topzonen = -100 #zone to apply BC changes to
    #  outfile=p['boun_file']#'mars_press.boun'
    #  utils.write_boun(join(macro_dir,outfile),ts[::4],[ps[::4],temps[::4]],step='ti_linear',varname=['pw','t'],zones=topzonen) #hourly BC changes 







    #  #---- Write the trac file 
    #  # (Either use an initial source or continuous production rate)
    #  # INITIAL SOURCE (default)
    #  if p['methane_prod_rate'] == 0.0:
        #  sourcetype='initial source'
        #  with open(join(macro_dir,'trac.dat'),'w') as ft:
            #  ft.write('trac\n')
            #  ft.write('1.e-50     1.0      1.e-09   1.0\n')
            #  ft.write('0.       1.e20    1.0e6   0.0\n')
            #  ft.write('20      1.6     1.e-3  1.e-2   10\n')
            #  ft.write('1\n')
            #  ft.write('1  methane\n')
            #  #  ft.write('0  0 0 1 0 3.e-05  0 0 0\n\n') #no M-Q
            #  ft.write('0  0 0 1 1 3.e-05  0 0 0\n\n') #M-Q
            #  ft.write('1 0 0 1\n\n')
            #  ft.write('1  0 0   1.e-50\n')
            #  ft.write('{} 0 0   1.e0\n\n\n'.format(contamzonen))
            #  ft.write('stop')
    #  # CONTINUOUS METHANE PRODUCTION SOURCE RATE
    #  else:
        #  sourcetype='continuous source'
        #  with open(join(macro_dir,'trac.dat'),'w') as ft:
            #  ft.write('trac\n')
            #  ft.write('1.e-50     1.0      1.e-09   1.0\n')
            #  ft.write('0.       1.e20    1.0e6   0.0\n')
            #  ft.write('20      1.6     1.e-3  1.e-2   10\n')
            #  ft.write('1\n')
            #  ft.write('1  methane\n')
            #  #  ft.write('0  0 0 1 0 3.e-05  0 0 0\n\n') #no M-Q
            #  ft.write('0  0 0 1 1 3.e-05  0 0 0\n\n') #M-Q
            #  ft.write('1 0 0 1\n\n')
            #  if p['c_initial'] == 0.0:
                #  ft.write('1  0 0   1.e-50\n\n')
            #  else:
                #  #  ft.write('1  0 0   1.e-7\n')
                #  #  ft.write('1  0 0   1.e-50\n')
                #  ft.write('1  0 0   9.611e-5\n')
                #  ft.write('{}  0 0   {}\n\n'.format(contamzonen,p['c_initial']))
            #  #  ft.write('{} 0 0   1.e0\n\n\n'.format(contamzonen))
            #  ft.write('{} 0 0   {}   0. 1e10\n\n'.format(contamzonen, c_inj))
            #  ft.write('stop')

    #  symlink(join(tpl_dir,'perm.dat'),join(macro_dir,'perm.dat'),overwrite=True)
    #--Write the rock.dat file from tpl
    #  symlink(join(tpl_dir,'rock.dat'),join(macro_dir,'rock.dat'),overwrite=True)
    pest_io.tpl_write(p,join(tpl_dir,'rock.tpl'),join(macro_dir,'rock.dat'))

    #  symlink(join(tpl_dir,p['boun_file']),join(macro_dir,p['boun_file']),overwrite=True)l
    #  symlink(join(tpl_dir,p['topzone_file']),join(macro_dir,os.path.split(p['topzone_file'])[1]),overwrite=True)
    #  symlink(join(tpl_dir,'air_static.boun'),join(macro_dir,'air_static.boun'),overwrite=True)
    #  symlink(join(tpl_dir,'eos_mars_co2.dat'),join(macro_dir,'eos_mars_co2.dat'),overwrite=True)
    #  symlink(join(,'eos_mars_co2.dat'),join(macro_dir,'eos_mars_co2.dat'),overwrite=True)

    ######### INITIALIZATION/SPINUP #################
    # If init_dir not provided, create it:
    # Run a simulation to test a geothermal gradient (don't restart from it) 
    if init0_dir is None:
        if shortened==True:
            curr_sim_dir = join(sim_dir,'geo_gradient')
            common_prefix = os.path.commonprefix([meshname, os.getcwd()])
            short_mesh_name = os.path.relpath( meshname, start=curr_sim_dir)
            print('Shortened mesh_file name due to fehmn.files restriction:')
            print(short_mesh_name)
            p['mesh_file'] = short_mesh_name
        #  utils.run_fehm(p,'init0',exe=exe,tpl_dir=local_tpl_dir,macro_dir=macro_dir,verbose=verbose)
        utils.run_fehm(p,'geo_gradient',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
        # Specify init0_dir
        #p['init0_dir'] = join(os.getcwd(),'init0')
        #  p['init0_dir'] = '../init0'
        p['init0_dir'] = '../geo_gradient'
    else:
        p['init0_dir'] = init0_dir
    #  # Now change properties
    #  if init_dir is None:
        #  utils.run_fehm(p,'init',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
        #  # Specify init_dir
        #  #p['init_dir'] = join(os.getcwd(),'init')
        #  p['init_dir'] = '../init'
    #  else:
        #  p['init_dir'] = init_dir
    # Now do a 'thermal_spinup'  run
    if spinup_dir is None:
        if shortened==True:
            curr_sim_dir = join(sim_dir,'thermal_spinup')
            common_prefix = os.path.commonprefix([meshname, os.getcwd()])
            short_mesh_name = os.path.relpath( meshname, start=curr_sim_dir)
            print('Shortened mesh_file name due to fehmn.files restriction:')
            print(short_mesh_name)
            p['mesh_file'] = short_mesh_name

        utils.run_fehm(p,'thermal_spinup',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
        p['spinup_dir'] = '../thermal_spinup'
    else:
        p['spinup_dir'] = spinup_dir

    ### Create flow macro based on init run to use to set air-static max radial boundary condition (x=xmax)
    ##utils.flow_macro(g,p['init_dir'],macro_dir=macro_dir)

    #  # Find out what the steady-state concentration of the contaminated zone
    #  # should be (for continuous injection simulations)
    #  if 'continuous' in sourcetype:
        #  if conc_init_dir is None:
            #  utils.run_fehm(p,'conc_init',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
            #  p['conc_init_dir'] = '../conc_init'
        #  else:
            #  p['conc_init_dir'] = conc_init_dir

    ######### TRACER #################
    # Determine mass of air in cavity chimney to specify moles Xe/kg-air correctly
    # Run fehm single time step setting porosities to -1 everywhere except in cavity/chimney
    # and collect mass of air
    #if m_air_dir is None:
    #    utils.run_fehm(p,'m_air',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)
    #    m_air_dir = join(os.getcwd(),'m_air')
    #with open(join(m_air_dir,'run_glob.his')) as fh:
    #    lns = fh.readlines()
    #    m_a = float(lns[-1].split()[9])
    #    p['C'] = p['yield'] / m_a
    print('p = ')
    print(p)
    print('geotherm = ')
    print(p['geothermal_gradient'])
    with open('v.pkl','wb') as f:
        pickle.dump([p,p['geothermal_gradient'],T_ref,T_0,T_upshift],f)
        #  pickle.dump([p['geothermal_gradient'],period,ampl,T_0,T_upshift,p['signal']],f)
        #  #  pickle.dump([sourcetype,fehm_upshift],f)
        #  pickle.dump([fehm_upshift],f)
    #  utils.run_fehm(p,'tracer',exe=exe,tpl_dir=tpl_dir,macro_dir=macro_dir,verbose=verbose)


