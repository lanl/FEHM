
import os,sys
from os.path import join
import numpy
import numpy as np
from subprocess import call
sys.path.append(join(os.getcwd(),'scripts/'))
from tools import find_nearest, sci_notation
from nilson import *
sys.path.append(join(os.getcwd(),'scripts/stats_tools/'))
from stats_tools import stats_funcs
from stats_funcs import chisqg, redchisqg
#  sys.path.append('/project/gas_seepage/jportiz/scripts/gasmigra')
sys.path.append(join(os.getcwd(),'scripts/gasmigra/'))
import utils
import mars_tools
from mars_tools import solsec,days2sols,sols2days,mars_seasons, findPrincipleComponents, pbl_heights
import tools
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import tri
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
from matplotlib import cbook
import scipy.integrate as integrate
from scipy import interpolate, signal, stats
from scipy.interpolate import interp1d
from matplotlib import ticker
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
#the below import allows you to outline text
import matplotlib.patheffects as path_effects
import pickle
import pandas as pd
import math
from datetime import datetime
from plot import synth_seasons, read_mesh
from gasmigra import signalTools
from signalTools import bandPassFilter, lowPassFilter
import glob
from scipy import ndimage
from scipy.optimize import minimize, differential_evolution, NonlinearConstraint, LinearConstraint
import pickle
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from fluxComparisonAnalysis import getYlim
from mpl_toolkits.mplot3d import Axes3D
import argparse
import shutil
from matplotlib.animation import FuncAnimation, FFMpegWriter
from sklearn.linear_model import LinearRegression
#-----------------------------------------------------
#MATPLOTLIBRC PLOTTING PARAMETERS
# Load up sansmath so that math --> helvetica font
# Also need to tell tex to turn on sansmath package
plt.rcParams['text.latex.preamble'] = [
    r'\usepackage{helvet}',
    r'\usepackage{sansmath}',
    r'\sansmath']
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Helvetica'
#  plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.labelweight']=u'normal'
plt.rcParams['agg.path.chunksize'] = 10000  #agg has a hardcoded limit on points
#-----------------------------------------------------

def V_standard(T=223.15,P=600.):
    '''
    Calculate volume of 1 g molecular weight of gas at given P and T.
    Assumes n=1 [mol].
    Used in unit conversions to ppbv.
    '''
    R=8.314462e3   #[L*Pa/K/mol]
    v = R * T / P  #[L/mol]
    return v

def conc_gm3_to_ppbv(C_original, T=223.15, P=600., molec_weight=16.04):
    '''
    Convert concentration in g/m^3 to ppbv.
    Parameters
    ----------
    C_original : float
        Concentration of gas [g/m^3]
    T : float
        Temperature [K]
    P : float
        Atmospheric pressure [Pa]
    molec_weight : float
        Molecular weight of trace gas [g/mol]. Default is methane.
    '''
    ppbv = V_standard(T,P) * C_original * 1e6 / molec_weight
    return ppbv

def conc_ppbv_to_gm3(C_ppbv, T=223.15, P=600., molec_weight=16.04):
    '''
    Convert concentration in ppbv to  g/m^3.
    Parameters
    ----------
    C_ppbv : float
        Concentration of gas [ppbv]
    T : float
        Temperature [K]
    P : float
        Atmospheric pressure [Pa]
    molec_weight : float
        Molecular weight of trace gas [g/mol]. Default is methane.
    '''
    gm3 = ( C_ppbv * molec_weight) / V_standard(T,P) / 1e6
    return gm3

def calc_rmse(obs,pred):
    rmse = np.sqrt( np.sum( (obs-pred)**2 ) / len(obs) )
    return rmse

def read_fehm_results(outputDir, filename='fehm_results.csv'):
    '''
    Read in ``fehm_results.csv`` file that was created by ``plot.py``.
    '''
    print('Reading in {} from {}...'.format(filename, outputDir))
    dataFrame = pd.read_csv(join(outputDir, filename))
    print(    'Done.')
    return dataFrame

def interpolate_df_cols(dataFrame, t_interp, col_interp='SOL_TIME'):
    '''
    Interpolate each column in dataframe to the times/provided (t_interp) based
    on column provided (col_interp).
    '''
    print('Interpolating data to evenly spaced times ({})...'.format(col_interp))
    # Get the columns we need to interpolate
    #  cols_bool = [ c for c in dataFrame.columns != col_interp ]
    cols = [ c for c in dataFrame.columns.drop([col_interp, 'SEASON']) ]
    #  cols_bool = [ c for c in dataFrame.columns.any() not in [col_interp, 'SEASON'] ]
    #  cols = dataFrame.columns[cols_bool]
    # Interpolate each of those columns to the vals in t_interp
    dfi = pd.DataFrame(data=t_interp, columns=['SOL_TIME'])
    for col in cols:
        f0  = interp1d(dataFrame[col_interp], dataFrame[col], fill_value='extrapolate')
        ysi = f0(t_interp)
        # Add interpolated column tp dataFrame
        pos = np.shape(dfi)[1]
        dfi.insert(pos, col, ysi)
    # Now add the SEASON column (which couldn't be interpolated)
    season_i = [ mars_seasons(ls) for ls in dfi['L_s'] ]
    dfi.insert( np.shape(dfi)[1], 'SEASON', season_i )
    print('    Done.')
    return dfi

#  decimalHours = 15.005109 #DEBUG
#  decimalHours = 15.500000 #DEBUG
#  decimalHours = 15.250000 #DEBUG
def decimal2hhmmssTime(decimalHours):
    try: len(decimalHours)
    except TypeError: decimalHours = np.asarray([decimalHours])
    #  hours = [int(s) for s in (sols % 1) * 24 ]  # NOTE: these are "Mars" hours (24 'hrs' per sol)
    hours = [int(dH) for dH in decimalHours]     # NOTE: these are "Mars" hours (24 'hrs' per sol)
    minutes = [int(m) for m in (decimalHours%1 )%1*60 ]
    seconds = [int(sec) for sec in ((decimalHours%1 )%1*60%1)*60 ]
    hhmmss = ['{:02d}:{:02d}:{:02d}'.format(hours[n], minutes[n], seconds[n]) for n in range(len(hours))]
    return hhmmss

def calc_decimalTime(dataFrame, tname='SOL_TIME'):
    # Calculate the decimal times
    dec_time = dataFrame.copy()[tname]%1*24  #[hrs]
    # Insert into dataframe
    pos = dataFrame.columns.get_loc(tname)+1
    dataFrame.insert(pos, 'DECIMAL_TIME', dec_time)
    return dataFrame

def calc_localMarsTime(dataFrame):
    '''
    Calculate local time in Mars "hours" and add to dataFrame.
    '''
    sols    = dataFrame['SOL_TIME']
    solint = [int(sol) for sol in sols]
    hours = [int(s) for s in (sols % 1) * 24 ]  # NOTE: these are "Mars" hours (24 'hrs' per sol)
    minutes = [int(m) for m in (sols%1*24 )%1*60 ]
    seconds = [int(sec) for sec in ((sols%1*24 )%1*60%1)*60 ]
    hhmmss = ['{:02d}:{:02d}:{:02d}'.format(hours[n], minutes[n], seconds[n]) for n in range(len(hours))]
    # Add this to the a copy of the dataFrame (maybe)
    pos = np.shape(dataFrame)[1]
    dataFrame_new = dataFrame.copy()
    dataFrame_new.insert(pos, 'LOCAL_TIME', hhmmss)
    return dataFrame_new

def readAbundances(run_type='Enrichment',source='Webster2017'):
    '''
    Read literature atmospheric abundances into DataFrame from file.

    source :: str
        Flag indicating which literature values to read.
        Options are ``Webster2017`` (default)  and ``Webster2021``.
    '''
    if source == 'Webster2017':
        ab = pd.read_csv('/project/gas_seepage/jportiz/mars/data/methane/literature_data/methane_data_webster2017.csv', sep=',') # Drop all 'Unnamed' columns
    elif source == 'Webster2021':
        ab = pd.read_csv('/project/gas_seepage/jportiz/mars/data/methane/literature_data/methane_data_webster2021.csv', sep=',') # Drop all 'Unnamed' columns
    else: print('ERROR: source must be either ``Webster2017`` (default) or ``Webster2021``.')
    ab = ab.loc[:, ~ab.columns.str.contains('^Unnamed')]
    #Only want Enrichment method abundances
    em = ab[ab['RunType']==run_type]
    #  em.rename(columns={"Sol_time":"SOL_TIME"}, inplace=True)
    em = em.rename(columns={"Sol_time":"SOL_TIME","Sol":"SOL_INT"})
    #  # Calculate LMST as DECIMAL_TIME and add to dataFrame
    #  em = calc_decimalTime(em, tname='Sol_time')
    # Convert LMST from hh:mm to DECIMAL_TIME
    dct = []
    for t in em['LMST']:
        ts = t.split(':')
        hour = int(ts[0]); mm = int(ts[1])
        #  decimal_time = round((hour+mm/60)/24., 3)
        decimal_time = round((hour+mm/60), 3)
        dct.append(decimal_time)
    pos = em.columns.get_loc('LMST')+1
    em.insert(pos, 'DECIMAL_TIME', dct)
    return em

def interpolate_model2obsAbundance(sim_df, obs_df, time1=0, time2=1e10):
    '''
    Interpolate modeled abundances to the decimal times of the observed
    abundances. Gets the L_s from obs, then interpolates the Y value from the
    decimal_times on that given SOL_INT value.

    - Requires a choice of time1, time2 that makes the dataFrame subset at most 1
    Mars year long.
    '''
    #---- Subset the modeled values to the times specified by time1 and time2 (≤1 year)
    df_ss = sim_df.copy()[ sim_df['SOL_TIME'].between(time1,time2) ]
    #---- Get SOL_INT of model at the observed L_s times 
    f1 = interp1d(df_ss['L_s'], df_ss['SOL_TIME'], fill_value='extrapolate')
    ci_soltime = f1(obs_df['Ls(deg)']); ci_solint = [int(ci_st) for ci_st in ci_soltime]
    #---- Interpolate modeled C_ppbv to DECIMAL_TIME of observed SOL_INT
    ci_dectime = []
    for i,si in enumerate(ci_solint):
        #Subset for each measurobs_dfent
        df_si = df_ss.copy()[ df_ss['SOL_INT']==si ]
        #Interpolate model to obs time
        f0 = interp1d(df_si['DECIMAL_TIME'], df_si['C_ppbv'], fill_value='extrapolate')
        # Modeled C_ppbv at desired obs time 
        ci_dectime.append( float(f0(obs_df['DECIMAL_TIME'].iloc[i])) )
    #---- Assemble DataFrame
    result = pd.DataFrame()
    vals = [ ci_solint, np.asarray(obs_df['Ls(deg)']), np.asarray(obs_df['DECIMAL_TIME']), ci_dectime ]
    cols = [ 'SOL_INT', 'L_s', 'DECIMAL_TIME', 'C_ppbv' ]
    for i,v in enumerate(vals):
        line = pd.Series(v)
        result.insert(i, cols[i], line)
    return result


def calcAbundError(sim_df, obs_df, time1=0, time2=1e10, dof=2 ):
    '''
    Calculate errors (RMSE, Chi-square, and Reduced Chi-Square statistic)
    between simulated abundance and observed abundance.
    (Perhaps change eventually to use just the moving average 'C_mavgs'...)

    sim_df :: DataFrame
        Provide the dataframe with simulated C_ppbv values.
    obs_df :: DataFrame
        Provide the dataframe with observed values.
    time1 :: float
        Optional param to calculate error only on a subset of times ['SOL_TIME']
        between time1 and time2.
    time2 :: float
        (See time1).
    '''
    #-------------------------------------------------- 
    # Calculate Error b/w simulated Abundance and observation
    #   (only calculate error on subset data – final s.s. year)
    #   (specified by user)
    #-------------------------------------------------- 
    # (A) Observations
    xa = np.asarray(obs_df['Ls(deg)'])
    ya = np.asarray(obs_df['GlobalMeanCH4(ppbv)'])
    ya_sem = np.asarray(obs_df['Error(±1SEM).1'])
    # (B) Simulated 
    sim_interp_vals = interpolate_model2obsAbundance(sim_df, obs_df, time1, time2)
    ybi = np.asarray(sim_interp_vals['C_ppbv'])
    #-- Calculate RMSE
    #  error = calc_rmse(ya,ybi)
    rmse = calc_rmse(ya,ybi)
    #-- Calculate Chi-Square Statistic 
    chisq = chisqg(ydata=ya, ymod=ybi, sd=ya_sem)
    #-- Calculate Reduced Chi-Square Statistic 
    redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=ya_sem)
    return rmse, chisq, redchisq

#  def calcAbundError(sim_df, time1=0, time2=1e10, useMovingAvg=False):
    #  '''
    #  Calculate errors (RMSE, Chi-square, and Reduced Chi-Square statistic)
    #  between simulated abundance and observed abundance.
    #  (Perhaps change eventually to use just the moving average 'C_mavgs'...)
#  
    #  time1 :: float
        #  Optional param to calculate error only on a subset of times ['SOL_TIME']
        #  between time1 and time2.
    #  time2 :: float
        #  (See time1).
    #  useMovingAvg :: bool
        #  Flag to calculate error using the moving average of the Concentration ['C_mavgs'].
    #  '''
    #  copied_df = sim_df.copy()
    #  # Calculate the Error
    #  # Get Webster2018a abundances for comparison
    #  em = readAbundances(run_type='Enrichment')
    #  #-------------------------------------------------- 
    #  # Calculate Error b/w simulated Abundance and observation
    #  #   (only calculate error on subset data – final s.s. year)
    #  #   (specified by user)
    #  #-------------------------------------------------- 
    #  # Get final subset of data for error checking
    #  #  time2 = copied_df['SOL_TIME'].iloc[-1] - (668*1)
    #  #  time1 = time2 - (668*2)
    #  copied_df_ss = copied_df.copy()[ copied_df['SOL_TIME'].between(time1,time2) ]
    #  # (A) Observations
    #  xa = em['Ls(deg)']
    #  ya = em['GlobalMeanCH4(ppbv)']
    #  ya_sem = em['Error(±1SEM).1']
    #  # (B) Simulated 
    #  #-- Interpolate moving average to the observed Ls values 
    #  xb = copied_df_ss['L_s']
    #  if useMovingAvg == True:
        #  yb = copied_df_ss['C_mavgs']
    #  elif useMovingAvg == False:
        #  yb = copied_df_ss['C_ppbv']
    #  #  yb = C_mavgs   #moving average
    #  #  f0 = interp1d(xb,yb)
    #  f0 = interp1d(xb,yb,fill_value='extrapolate')
    #  ybi = f0(xa)
    #  #-- Calculate RMSE
    #  #  error = calc_rmse(ya,ybi)
    #  rmse = calc_rmse(ya,ybi)
    #  #-- Calculate Chi-Square Statistic 
    #  #  chisq = chisqg(ydata=ya, ymod=ybi, sd=None)
    #  chisq = chisqg(ydata=ya, ymod=ybi, sd=ya_sem)
    #  #-- Calculate Reduced Chi-Square Statistic 
    #  dof = 2 # no. of degrees of freedom
    #  #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=None)
    #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=ya_sem)
    #  return rmse, chisq, redchisq

def thomas(a,b,c,r):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(r) # number of equations
    ac, bc, cc, rc = map(np.array, (a, b, c, r)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        rc[it] = rc[it] - mc*rc[it-1]
    xc = bc
    xc[-1] = rc[-1]/bc[-1]
    for il in range(nf-2, -1, -1):
        xc[il] = (rc[il]-cc[il]*xc[il+1])/bc[il]
    return xc

def calc_moving_avg(dataFrame, xname='SOL_TIME', yname='C_ppbv', window_size=1.0):
    '''
    Parameters:
    ----------
    dataFrame : pandas DataFrame
        Has your data in it.
    xname : str
        Name of x or time column
    yname : str
        Name of y column
    window_size: float
        [sols] Size of time window to conduct moving average over.
        Should have same units as values in dataFrame[xname].
    '''
    # Calculate moving average in a given window  of size window_size
    dt = np.diff(dataFrame[xname])[0]
    win = int(round(window_size/dt))                #mv-win-sol window
    #  C_windows = C_ppbv.rolling(win)
    C_ppbv = dataFrame.copy()[yname]
    C_windows = C_ppbv.rolling(win, center=True)
    # Calculate average of each window
    C_mavgs = C_windows.mean()
    return C_mavgs


def diffusion_PBL(outdirname, dataFrame, D_collapsed, D_expanded, k_c, k_e, delx=1.0, delt=0.04, h_obs=1.0):
    '''
    Backward Euler for the one-dimensional diffusion equation.
    Fluxes feed from subsurface seepage model. Modeled diffusion into the PBL,
    which can vary in height (like a step function...). PBL can be in a
    "collapsed" state (quiescent, nighttime) or an "expanded" state (very deep,
    high mixing coefficient).

    Parameters:
    ----------
    D_collapsed : float
        Dispersion/diffusion coefficient [m2/s] for collapsed PBL state.
    D_expanded : float
        Dispersion/diffusion coefficient [m2/s] for expanded PBL state.
    k_c : float
        [1/s] 1st-order mass loss term for collapsed state. ( dC/dt = -kC )
        User should try to default to photochemical loss rate (~1.4e-10 [1/s]).
    k_e : float
        [1/s] 1st-order mass loss term for expanded state. ( dC/dt = -kC )
        Should always be ≥ k_c.
    delx : float
        [m] vertical grid spacing. Should do in increments of 0.5 or 1.0.
    delt : float
        [sols] time step size to use.
    h_obs : fload
        Observation node. Height of the SAM-TLS instrument [m]. This is where
        concentrations are being collected and calculated.
    '''
    #set up parameters
    #  D=50        #[m2/s] diffusivity
    #  k=1    #[1/t]  first order decay rate
    #  h = 10 #250.;   #[m]    PBL height
    #  h_obs = 1.0     #[m]  height of SAM-TLS (?)
    # Convert k loss rates from [1/s] to [1/sol]
    #  k = k_loss*solsec(668)  #[1/sol]   #<----- WRONG (?) (?) (?)
    k_c_sols = k_c*solsec(1)  #[1/sol] 
    k_e_sols = k_e*solsec(1)  #[1/sol] 
    k = k_c_sols              #initial k term
    #  delx=0.5 #0.005;#1;           #[m] grid spacing
    j_obs = int(h_obs/delx)       #idx of observation node

    # How long is simulation? 
    t_start = dataFrame['SOL_TIME'].iloc[0]
    t_end   = dataFrame['SOL_TIME'].iloc[-1]
    #  delt=0.001#5*delx^2; #made up a time-step  (make sure same to check units)
    duration = t_end-t_start  #[sols]
    nt=int(np.fix(duration/delt))   #number of time-steps to finish simulation 
    # Interpolate Fluxes and PBL_m for use later
    interp_obj_flux = interp1d(dataFrame['SOL_TIME'], dataFrame['FLUX']/1e6, fill_value='extrapolate')  #[kg/m2/s]
    interp_obj_pbl  = interp1d(dataFrame['SOL_TIME'], dataFrame['PBL_m'], fill_value='extrapolate')

    #-------------------------------------------------- 
    # Get initial PBL height 
    #-------------------------------------------------- 
    # Interpolate from dataFrame
    L_prev = int(interp_obj_pbl(t_start))
    # Get number of nodes (will change with PBL height in time)
    J=int(L_prev/delx + 1)             #number of nodes at new timestep 
    x_prev=np.linspace(0,L_prev,J)          #create x vector 

    C=np.zeros(J)       #initial condition 
    Cplot=[C]          #plot initial condition

    # Get max number of nodes there could be based on max PBL height
    Jmax = int(math.ceil(max(dataFrame['PBL_m'])/delx))
    #  Carr = np.zeros((nt,Jmax))  #<-- could initialize with negative ones, so all negative values at the end are from over the height of the PBL at that time
    Carr = np.ones((nt+1,Jmax))*-1  #<-- could initialize with negative ones, so all negative values at the end are from over the height of the PBL at that time
    C_obs = np.ones(nt+1)*-1
    Carr[0]=0.; C_obs[0]=0.0
    Larr = np.zeros(nt+1)
    Larr[0] = L_prev
    Fluxarr = np.zeros(nt+1)
    Fluxarr[0] = float(interp_obj_flux(t_start))

    ## loop over time
    #  ts=[0]
    #  ts=[]; real_ts=[]
    #  ts=np.zeros(nt+1); real_ts=np.zeros(nt+1)
    ts=np.zeros(nt+1); real_ts=np.zeros(nt+1); real_ts[0]=t_start
    #  print(f'nt = {nt}')
    for it in range(1,nt+1):
        time = it*delt;
        t_real = t_start + time
        #  print('--')
        #  print('it     = {}'.format(it))
        #  print('time   = {}'.format(time))
        #  print('t_real = {}'.format(t_real))
        #-------------------------------------------------- 
        # Get new FLUX and PBL height
        #-------------------------------------------------- 
        # Interpolate from dataFrame
        L    = int(interp_obj_pbl(t_real))
        # Convert flux to kg/m2/sol!
        flux = -1 * float(interp_obj_flux(t_real))*solsec()  #<--- MAKE NEGATIVE ??? ??? ???
        #  flux = -1 * float(interp_obj_flux(t_real))  #<--- MAKE NEGATIVE ??? ??? ???
        #  flux = float(interp_obj_flux(t_real))  #<--- MAKE NEGATIVE ??? ??? ???

        # Update Certain Parameters Based on PBL
        #-------------------------------------------------- 
        # Choose Diffusion Coefficient
        # ---- Convert Diffusion Coefficient from m2/s to m2/sol
        L_threshold = 300.  #[m] PBL height above which we use D_expanded
        if L < L_threshold:
            D = D_collapsed*solsec()  #[m2/sol]
            k = k_c_sols
            #  k = k_loss_sols
            #  D = D_collapsed
        else:
            D = D_expanded*solsec()   #[m2/sol]
            # TESTING !!! try increasing loss when PBL is expanded 
            k = k_e_sols
            #  k = k_loss_sols*1e5 # !!! !!!  <--------1e6 makes a difference (1e4 less-so)
            #  D = D_expanded
        # Get number of nodes (changes with PBL height)
        J=int(L/delx + 1)             #number of nodes at new timestep 
        mu=D*delt/delx**2             #mu shows up in coefficients
        x=np.linspace(0,L,J)          #create x vector 

        #-------------------------------------------------- 
        # create a,b,c vectors for Thomas algorithm - this could be outside the time
        #-------------------------------------------------- 
        #loop - for indexing reasons, we are using a(1)=0, c(J)=0; that way
        #the Thomas function is easier to write - these do not occur explicitly
        #in the Thomas algorithm but allow a(2) to line up with b(2), etc.
        #filling in and then overwriting the first and last row for B.C.
        # {B}
        b=1+(2*mu+k*delt)*np.ones(J)  #Interior Nodes - accounts for diffusion and decay
        b[0] =  -D/delx #specified flux on LEFT
        b[-1]=  1/delx   #no-flux on RIGHT   (x=h)
        #     b(1,1)= -1/delx;  #no-flux on LEFT

        # {A}
        a=-mu*np.ones(J)      #Interior Nodes
        a[0] = 0     #specified flux on LEFT   (not used by Thomas algorithm)
        a[-1]=-1/delx #no-flux on RIGHT

        # {C}
        c=-mu*np.ones(J)      #Interior Nodes
        #  c[0]  = 1/delx;  #specified flux on LEFT
        c[0]  = D/delx  #specified flux on LEFT
        c[-1] = 0.0      #specified flux on RIGHT (not used by Thomas algorithm)

        #-------------------------------------------------- 
        # Map Previous Concentration Profile to Current PBL Column for 'new' initial state
        #-------------------------------------------------- 
        # ALTERNATIVE WAY!
        # ----> Don't stretch conc profile going from collapsed to expanded
        # ----> Do compress   conc profile going from expanded  to collapsed
        C_prev = C  #previous concentration profile
        #  # Integrate concentration over profile to get total mass (mass conservation)
        #  m_prev = integrate.cumtrapz(y=C_prev, x=x_prev)[0]  #[kg]
        L_curr = L   #[m] current PBL height
        x_curr = x  # x arrray with current # of nodes
        #---- Going from EXPANDED to COLLAPSED
        if L_curr < L_prev:
            x_scaled = x_prev * L_curr/L_prev  #compress x array to current height (has prev # nodes)
            #Scale the Conc profile based on change in L
            C_scaled = C_prev * L_prev/L_curr
            #  C_scaled = C_prev #TESTING!!!
            # Interpolate scaled profile to x_curr
            f0 = interp1d(x_scaled, C_scaled)
            C_curr = f0(x_curr)
            x_scaled = x_prev * L_curr/L_prev
            # Make C_curr the length of updated PBL height, with all zeros 
            # (except for the lowest nodes that correspond to 
            C = C_curr
        #---- Going from COLLAPSED to EXPANDED 
        else:
            C_curr = np.zeros(len(x_curr))
            C_curr[:len(C_prev)] = C_prev
            C = C_curr


        #-------------------------------------------------- 
        # set up the right hand side vector - this should be in the time loop
        #-------------------------------------------------- 
        # (rhs is same as f vector)
        rhs=[[flux], list(C[1:-1]), [0.]];#f is the time-dependent flux on the left, 0 is no flux on the right
        rhs = [element for sublist in rhs for element in sublist]

        #call Thomas algorithm
        # Remove the non-existent values from a and c vectors (0.0)
        Cnext=thomas(a[1:],b,c[:-1],rhs)  #TESTINTG....
        #  Cnext=thomas(a,b,c,rhs)
        C=Cnext
        #---------------------
        # Store vals in array
        #---------------------
        Carr[it,:J] = C
        #  Carr[it,:len(C)+1] = C   #CHECK THIS.... make sure nothing gets cut off...
        C_obs[it]   = C[j_obs]  #C at observation node (SAM-TLS)
        Larr[it]    = L
        Fluxarr[it] = abs( flux/solsec() )
        #plot
        x=np.linspace(0,L,J)
        ts[it]      = time
        real_ts[it] = t_real
        x_prev=x  #set the x_prev for next timestep
        L_prev=L  #set the x_prev for next timestep

    # Save to dataFrame and CSV files
    #---- C_obs 
    # Calculate C_ppbv
    C_gm3 = C_obs*1000.  #convert kg/m3 to g/m3
    C_ppbv = conc_gm3_to_ppbv(C_gm3)
    # --
    data = np.column_stack( (real_ts, Fluxarr, Larr, C_obs, C_ppbv) )
    df = pd.DataFrame( data, columns=['SOL_TIME', 'FLUX_kg/m2/s', 'PBL_m', 'C_kgm3', 'C_ppbv'] )
    # -- Get other relevant columns by interpolating from input dataFrame
    cols = ['L_s', 'PRESS', 'TEMP']
    for ind,c in enumerate(cols):
        f0 = interp1d(dataFrame['SOL_TIME'], dataFrame[c])
        ci = f0(real_ts)
        df.insert(ind+1,c,ci)
    # -- Add decimal times
    df = calc_decimalTime(df, tname='SOL_TIME')
    #-- Insert integer sol
    sol_int = [ math.floor(s) for s in df['SOL_TIME'] ]
    pos = df.columns.get_loc('SOL_TIME')+1
    df.insert(pos, 'SOL_INT', sol_int)
    # -- SAVE
    outfilename = join(outdirname, 'conc_obs.csv')
    df.to_csv(outfilename, index=False)
    #  #---- Full Carr
    #  # >>>> Might use too much storage ... 1 year w/ delt=0.01 is 9.5GB...i <<<<
    # >>> (could just save last year...)
    #  outfilename = join(outdirname, 'conc_contours.csv')
    #  np.savetxt( outfilename,  Carr, delimiter=',' )
    #-- This is more efficient, uses less storage
    outfilename = join(outdirname, 'conc_contours.csv')
    np.save(outfilename, Carr)

    return df, Carr, C_obs, C_ppbv

def objFunction(params, dataFrame, t1_err, t2_err, constantPT=False):
    #  print(params[0])
    D_c      = 10**params[0]
    #  D_e      = 10**params[1]
    de_over_dc_factor = params[1]
    D_e      = D_c * de_over_dc_factor
    k_c      = 10**params[2]
    ke_over_kc_factor = 10**params[3]
    k_e      = k_c * ke_over_kc_factor
    #  k_e      = 10**params[2]
    print()
    print('D_c       = {:e}'.format(D_c))
    print('D_e       = {:e}'.format(D_e))
    print('k_c       = {:e}'.format(k_c))
    print('k_e       = {:e}'.format(k_e))

    # How many degrees of freedom?
    dof = len(params)

    # Run the diffusion model to calculate abundances
    # -----------------------------------------------
    #  df, c_all, c_obs, C_ppbv = diffusion_PBL(outputdir, dataFrame, D_collapsed, D_expanded, k_c, k_e, delx=delx, delt=delt, h_obs=h_obs)
    df, c_all, c_obs, C_ppbv = diffusion_PBL(outputdir, dataFrame, D_c, D_e, k_c, k_e, delx=delx, delt=delt, h_obs=h_obs)

    # -----------------------------------------------
    # Scale the FLUX so that mean of observation points is same as SAM-TLS background mean  (0.41ppbv)
    # -----------------------------------------------
    C_background = 0.41 #[ppbv]
    #---- Get Webster2018a abundances for comparison
    em = readAbundances(run_type='Enrichment',source='Webster2021')
    #  # Omit the outlier with C=20.5ppmv
    em = em.copy()[ em['GlobalMeanCH4(ppbv)']<20.]
    # Also omit the non-detections (they bring down the background avg) 
    em = em.copy()[ em['GlobalMeanCH4(ppbv)']>0.068]
    em = em.sort_values('Ls(deg)')
    #---- Subset the modeled values to the times specified (1 year)
    df_ss = df.copy()[ df['SOL_TIME'].between(t1_err,t2_err) ]
    #---- Interpolate sim concs at the time of observations
    sim_interp_vals = interpolate_model2obsAbundance(df_ss, em, t1_err, t2_err)
    #  #---- Calculate by what factor we need to reduce fluxes
    #  #  factor = np.mean(ci_dectime)/C_background
    #  factor = np.mean(sim_interp_vals['C_ppbv'])/C_background
    #  print('    (Reduction factor = {:.3f})'.format(factor))
    #  #---- NEW WAY !!! !!! !!! !!! !!! !!!
    #  #-------- Use a fixed value for average flux of 1e-13 kg/m2/s [V-M 2021]
    #  #---- Calculate by what factor we need to reduce fluxes
    #  flux_literature = 1e-13 #1e-16  #[kg/m2/s]
    #  # Get average annual flux from sim
    #  flux_avg = np.mean(df_ss['FLUX_kg/m2/s'])
    #  factor = flux_avg / flux_literature
    #  print('    (Reduction factor = {:.3f})'.format(factor))
    #  #---- Modify the fluxes and abundance
    #  flux_adj = df['FLUX_kg/m2/s'] / factor
    #  C_ppbv_adj = df['C_ppbv'] / factor
    #  #---- Add adjusted values into DataFrame
    #  try:
        #  pos1 = df.columns.get_loc('C_ppbv')
        #  df.insert(pos1+1, 'C_ppbv_ADJ',C_ppbv_adj)
        #  pos2 = df.columns.get_loc('FLUX_kg/m2/s')
        #  df.insert(pos2+1, 'FLUX_kg/m2/s_ADJ', flux_adj)
    #  except ValueError:
        #  df.pop('C_ppbv_ADJ'); df.pop('FLUX_kg/m2/s_ADJ')
        #  pos1 = df.columns.get_loc('C_ppbv')
        #  df.insert(pos1, 'C_ppbv_ADJ',C_ppbv_adj)
        #  pos2 = df.columns.get_loc('FLUX_kg/m2/s')
        #  df.insert(pos2, 'FLUX_kg/m2/s_ADJ', flux_adj)
    #  #  #---- Write this dataFrame to a CSV file
    #  #  df.to_csv(join(
    # NEW WAY
    # ---- factor is calculated outside of objFunction and passed to it

    # -----------------------------------------------
    # Calculate the Error 
    # -----------------------------------------------
    #---- Get Abundances again, but leave in the daytime non-detections for error calculation
    em = readAbundances(run_type='Enrichment',source='Webster2021')
    #  # Omit the outlier spike with C=20.5ppmv
    em = em.copy()[ em['GlobalMeanCH4(ppbv)']<20.]
    em = em.sort_values('Ls(deg)')
    #  # Replace C_ppbv with adjusted values for error calculation
    #  df_ss['C_ppbv'] = C_ppbv_adj
    #  rmse, chisq, redchisq = calcAbundError(df_ss, em, time1=t1_oneyear, time2=t2_oneyear, dof=dof)
    rmse, chisq, redchisq = calcAbundError(df_ss, em, time1=t1_err, time2=t2_err, dof=dof)
    print('RMSE               = {}'.format(rmse))
    print('Chi-Square         = {}'.format(chisq))
    print('Reduced Chi-Square = {}'.format(redchisq))
    # Choose one of these errors to tune the model
    #  error = rmse
    error = redchisq

    # -----------------------------------------------
    # Write pars and errors to file (manually for now) 
    # -----------------------------------------------
    restartdir = join(cwd, glob.glob('*restart*')[0])
    with open(join(restartdir, 'pbl_output/errors.csv'), 'a') as f:
    #  with open(join('/lclscratch/jportiz/projects/gas_seepage/mars/runs/2d_fracture_network_runs/regolith_adsorption_synth/sorp_depth200_b1mm_fracDen010/sorp_tracer_restart','new_pbl_test/errors.csv'), 'a') as f:
        outcols = [D_c,D_e,k_c,k_e,factor,rmse,redchisq]
        l='{},'*(len(outcols)-1)+'{}\n'
        f.write(l.format(*outcols))
        #  f.write(l.format([o for o in outcols]))
        #  f.write('{},{},{},{},{},{},{}\n'.format(D_c,D_e,k_c,k_e,factor,rmse,redchisq))
        #  f.write('{},{},{},{},{}\n'.format(D_c,D_e,k_e,factor,rmse,redchisq))
    return error


#--------------------------------------------------------------------- 
# Run plot script either on its own or from within ``run.py``
#--------------------------------------------------------------------- 
if __name__ == "__main__":
    print('Running pbl_diffusion.py as main...')

    #-------------------------------------------------- 
    # Define args for user to choose which part(s) of this script to run 
    #-------------------------------------------------- 
    descr ='Run the ``pbl_diffusion.py`` script. Running with no options defaults to running the whole script, which optimizes using differential_evolution to find the best parameters. Specifying an optional argument will skip the optimiztion and perform plotting using the provided parameters.'
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('-p', '--plotOnly', help='Full path to parameter file')
    #  parser.add_argument('-p', '--plotOnly',               action='store_true', help='path of parameter file')
    #  parser.add_argument('-p', '--plotOnly',               action='store_true', help='Do not optimize, just plot results for specific parameter combinations')
    args = parser.parse_args()
    #  print(args)                      #DEBUG

    if args.plotOnly:
        # [special case] special scenario name to add to output files if run using optional args
        special_id = os.path.split(args.plotOnly)[1][:-4]  #remove the '.csv' from path
        special_id_fn = '-'+special_id                     # for filenames
        print('Going to skip optimization and just do plotting for {}...'.format(special_id))
        #  input()
    else:
        # [default] special scenario name to add to output files if run using optional args
        special_id    = ''
        special_id_fn = ''
        print('Going to run the whole script...')
    #-------------------------------------------------- 


    cwd = os.getcwd()
    datadir = join(cwd, glob.glob('*restart*')[0], 'output')
    #  # Manual choice... 
    #  datadir = join(cwd, glob.glob('*tracer2*')[0], 'output')
    #  outputdir = datadir
    outputdir = join(cwd, glob.glob('*restart*')[0], 'pbl_output')  #TEMPORARY
    # Create outputdir if it doesn't exist
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    #-------------------------------------------------- 
    # Read in formatted DataFrame output by ``plot.py``
    #-------------------------------------------------- 
    df = read_fehm_results(datadir, filename='fehm_results.csv')

    #-------------------------------------------------- 
    # Interpolate to evenly spaced times 
    #-------------------------------------------------- 
    # Times for interpolation
    dt = 0.04  #[sols]
    min_t = 0.
    max_t = max(df['SOL_TIME'])
    ts = np.arange(min_t, max_t+dt, dt)
    dfi = interpolate_df_cols(df.copy(), ts, col_interp='SOL_TIME')


    #-------------------------------------------------- 
    # Subset the flux data to late time
    #   - because of some sim crashes, prob need 
    #     to keep calculations between 
    #     ~39100 and ~49000 sols
    #-------------------------------------------------- 
    # Several Years at Late Time
    # --------------------------
    #  t2 = dfi_ss['SOL_TIME'].iloc[-1] - (668*2)
    t2 = dfi['SOL_TIME'].iloc[-1] - (668*2)
    t1 = t2 - (668*5)
    #  print('t1 = {}'.format(t1))
    #  print('t2 = {}'.format(t2))
    dfi_ss = dfi.copy()[ dfi['SOL_TIME'].between(t1,t2) ]
    # A Single Year at Late time 
    # --------------------------
    t2_oneyear = t2
    t1_oneyear = t2_oneyear - (668*1)
    #  print('t1_oneyear = {}'.format(t1_oneyear))
    #  print('t2_oneyear = {}'.format(t2_oneyear))

    #-------------------------------------------------- 
    # Get mesh dimensions 
    #-------------------------------------------------- 
    # Go up 1 level from outputdir
    rundir = os.path.split(outputdir)[0]
    grid, mesh_depth, mesh_xwidth, mesh_zwidth, surfArea = read_mesh(rundir) #surfArea = km^2

    # 1-Minute PBL height per 30° bin (after Guzewich2017)
    #----------------------------------------------------- 
    pbl_data = pbl_heights(onlyMinima=False)
    all_Ls = pbl_data['L_s'].unique()
    # Get Ls bins of simulation results
    Ls_bin = []
    for ls in dfi_ss['L_s']:
        lbin = all_Ls[ all_Ls.searchsorted(ls) -1 ]
        Ls_bin.append( lbin )
    # Add Ls bins to dataFrame temporarily
    dfi_ss.insert(0,'Ls_bin',Ls_bin)
    # Also add decimal time
    #  dec_times = dfi_ss['SOL_TIME']%1
    dec_times = dfi_ss['SOL_TIME']%1*24
    #  dfi_ss.insert(1,'dec_times',dec_times)
    dfi_ss.insert(1,'DECIMAL_TIME',dec_times)
    # Insert integer sol
    sol_int = [ math.floor(s) for s in dfi_ss['SOL_TIME'] ]
    dfi_ss.insert(1,'SOL_INT', sol_int)
    #  # Interpolate PBLs to given times for each Ls bin 
    #  pbl_dfs = []
    #  for lsbin in all_Ls:
        #  print(lsbin)
        #  # Subset the dataframes to just current bin
        #  ss_p = pbl_data.copy()[ pbl_data['L_s']==lsbin ]
        #  ss_d = dfi_ss.copy()[ dfi_ss['Ls_bin']==lsbin ]
        #  #  f0 = interp1d( ss_p['LTST'], ss_p['PBL_m'] )
        #  # Interpolate PBL to times separately for each sol
        #  for s in np.unique(sol_int):
            #  ss_s = ss_d.copy()[ ss_d['SOL_INT']==s ]
            #  f0 = interp1d( ss_p['LMST'], ss_p['PBL_m'] )
            #  pbli = f0(ss_s['dec_times'])
            #  pbl_dfs.append(pbli)
        #  #  f0 = interp1d( ss_p['LMST'], ss_p['PBL_m'] )
        #  #  pbli = f0(ss_d['dec_times'])
        #  #  pbl_dfs.append(pbli)
    # Interpolate PBLS to dec_time first
    pbl_interp_dict = {}
    dec_times = dfi_ss.copy()[ dfi_ss['SOL_INT'] == dfi_ss['SOL_INT'].unique()[1] ]['DECIMAL_TIME']
    #  dec_times = dfi_ss.copy()[ dfi_ss['SOL_INT'] == dfi_ss['SOL_INT'].unique()[1] ]['dec_times']
    for lsbin in all_Ls:
        # Subset the PBL dataframe to just current Ls bin
        ss_p = pbl_data.copy()[ pbl_data['L_s']==lsbin ]
        # Interpolate to dec_times
        f0 = interp1d( ss_p['LMST'], ss_p['PBL_m'] )
        pbli = f0(dec_times)
        # Add these to dict
        pbl_interp_dict[lsbin] = [np.asarray(dec_times), pbli]
    # Now go by sol and add in the PBLs 
    #TESTING
    #  si = 45424
    #~~~~~~~~~
    all_si =dfi_ss['SOL_INT'].unique() #each sol_int
    saved_pbli = []
    for si in all_si:
        #  print(f'sol_int = {si}')  #DEBUG
        ss = dfi_ss.copy()[ dfi_ss['SOL_INT']==si ]
        # Determine which Ls bin it starts in
        lsbin = ss['Ls_bin'].iloc[0]
        # Interpolate the PBL heights to all times at this sol_int
        f0 = interp1d( pbl_interp_dict[lsbin][0], pbl_interp_dict[lsbin][1] )
        pbli = f0(ss['DECIMAL_TIME'])
        #  pbli = f0(ss['dec_times'])
        saved_pbli.append(pbli)
    dfi_ss.insert(np.shape(dfi_ss)[1], 'PBL_m', np.concatenate(saved_pbli) )




    C_0 = 0.41 #[ppbv] from [Klusman2022]  (0.45ppbv in  [Moores2019a])
    # Calculate air volume [m^3] of atmosphere 
    surface_area = mesh_xwidth*mesh_zwidth    #[m^2]
    # Get an approximate idea for the loss term  based on average influx into atmosphere
    # Assumes constant dt in dataframe
    dt = solsec(np.diff(dfi_ss['SOL_TIME'])[0])         #[s]


    #  #  #-------------------------------------------------- 
    #  #  # Plot PBL Heights for reference 
    #  #  #  (( TRY TO DISCRETIZE ))
    #  #  #--------------------------------------------------
    #  #only check Lsbin 150
    #  ss120 = dfi_ss.copy()[ dfi_ss['Ls_bin']==90 ]
    #  ss120 = dfi_ss.copy()[ dfi_ss['Ls_bin']==120 ]
    #  ss150 = dfi_ss.copy()[ dfi_ss['Ls_bin']==150 ]
    #  ss240 = dfi_ss.copy()[ dfi_ss['Ls_bin']==240 ]
    #  #  targ = ss240.copy()[ ss240['SOL_INT']==45575 ]
    #  targ = ss240.copy()[ (ss240['SOL_INT']>=45575) & (ss240['SOL_INT']<45585) ]
    #  x = targ['dec_times']
    #  y = targ['PBL_m']
#  
    #  # Discretize
    #  # Interpolate the PBL
    #  f0 = interp1d(x,y)
    #  dx = np.diff(x)[0]/2
    #  xi = np.arange(x.iloc[0], x.iloc[-1], dx)
    #  #  yi = f0(xi)
    #  #  ynew = [y.iloc[0]]
    #  ynew = []
    #  #  for i,xp in enumerate(x[:-1]):
        #  #  ynew[i] = y[
#  
#  
#  
    #  fig,ax1 = plt.subplots(1, figsize=(12,4) )
    #  #  for si in ss120['SOL_INT'].unique():
        #  #  x = ss120['dec_times']
        #  #  y = ss120['PBL_m']
        #  #  ax1.plot(x,y,c='g')
    #  #  for si in ss150['SOL_INT'].unique():
        #  #  x = ss150['dec_times']
        #  #  y = ss150['PBL_m']
        #  #  ax1.plot(x,y,c='b')
    #  #  for si in ss240['SOL_INT'].unique():
        #  #  x = ss240['dec_times']
        #  #  y = ss240['PBL_m']
        #  #  ax1.plot(x,y,c='r')
    #  targ = ss240.copy()[ ss240['SOL_INT']==45575 ]
    #  t = targ['dec_times']
    #  L = targ['PBL_m']
    #  ax1.plot(t,L,c='r')
    #  ax1.xaxis.set_minor_locator(MultipleLocator(1))
    #  ax1.xaxis.set_major_locator(MultipleLocator(2))
    #  #  ax1.set_xlim(7,10)
    #  plt.savefig(join(outputdir, 'pbls.pdf'))
    #  plt.close('all')



    #================================================== 
    #  OPTIMIZE 
    #================================================== 

    # TESTING ------
    outdirname = outputdir
    start =  t1
    end = t2
    dataFrame = dfi_ss.copy()[ dfi_ss['SOL_TIME'].between(start,end) ]
    #---- Desired Average Flux
    flux_literature = 2e-16#5.e-14  #1e-13  [kg/m2/s]  [V-M 2021]
    #  print('Desired average flux: {:.2e} kg/m2/s'.format(flux_literature))
    flux_avg = np.mean(dataFrame.copy()[ dataFrame['SOL_TIME'].between(t1_oneyear, t2_oneyear) ]['FLUX']/1e6)
    # Calculate factor by which to multiply by current fluxes to get desired flux
    factor = flux_avg / flux_literature  #value of factor will be saved in objFunction()
    #  print('Factor              : {:.4f}'.format(factor))
    flux_orig = dataFrame.copy()['FLUX'] #[kg/km2/s]
    flux_adj  = flux_orig / factor
    dataFrame['FLUX'] = flux_adj  #replace with adjusted values
    #  print('New     average flux: {:.2e} kg/m2/s'.format(np.mean(flux_adj/1e6)) )
    #----
    #  #  dataFrame['FLUX']=dataFrame['FLUX']/factor
    #  D_collapsed = 0.035#3.5  #[m2/s]
    #  D_expanded = 5#500.  #[m2/s]
    #  #  k_loss = 1/( 304*(86400*365.) ) #[1/s]
    #  #  k_c    = 1/( 304*(86400*365.) ) #[1/s]
    k_c_photochemical    = 1/( 304*(86400*365.) ) #[1/s]
    #  k_e    = k_c*1e5
    delx=1.0  #[m]
    delt=0.04 #[sols]
    h_obs=1.0
    #---------------


    #---------------
    # Set Parameter Bounds 
    #---------------
    bnds_Dc = ( -4., 1.2 )          # [0] log scale
    bnds_De = (1, 1000 )            # [1] factor by which D_e ≥ D_c 
    bnds_kc = (np.log10(k_c_photochemical), -1)   # [2] log scale
    bnds_ke = (0, 6)    #[3] LOG factor by which k_e ≥ k_c
    #  #---- Nonlinear Constraint (log(k_e) be no bigger than -1)
    #  # https://towardsdatascience.com/introduction-to-optimization-constraints-with-scipy-7abd44f6de25
    #  # Add constraint that the value of log(k_e) be no bigger than -1
    #  #  constraint_fun = lambda x: x[3] <= -1  # log(k_e) ≤ -1
    #  #  constraint_fun = lambda x: 10**(x[3]*x[2])
    #  #  constraint_fun = lambda x: (10**x[3])*(10**x[2])
    #  constraint_fun = lambda x: x[3]+x[2]    #log space
    #  #  nlc = NonlinearConstraint(constraint_fun, -np.inf, 0.0)
    #  #  nlc = NonlinearConstraint(constraint_fun, -np.inf, 0.0)
    #  #  nlc = NonlinearConstraint(constraint_fun, -np.inf, -1.)
    #---- Linear Constraint (log(k_e) be no bigger than -1)
    # https://towardsdatascience.com/introduction-to-optimization-constraints-with-scipy-7abd44f6de25
    #  lc = LinearConstraint( np.eye(4), [-np.inf]*4, [np.inf, np.inf, np.inf, -1] )
    A = np.eye(4)
    A[-1] = [0.,0.,1,1]
    lc = LinearConstraint( A, [-np.inf]*4, [np.inf, np.inf, np.inf, -1] )

    #---- Assemble Bounds Tuple
    #  bnds = ( bnds_Dc, bnds_De, bnds_ke)        # 3 PARAMS
    bnds = ( bnds_Dc, bnds_De, bnds_kc, bnds_ke)  # 4 PARAMS

    # ================================================== 
    # Only optimize if ``plotOnly`` arg was not specified
    # ================================================== 
    if not args.plotOnly:

        #---------------
        # Save values and Errors to a Text File
        #---------------
        start = time.perf_counter() #[s]
        outfilename = join(outputdir, 'errors.csv')
        errfilename = outfilename
        with open(outfilename, 'w') as f:
            f.write('D_c,D_e,k_c,k_e,factor,rmse,redchisq\n')
            #  f.write('D_c,D_e,k_e,factor,rmse,redchisq\n')
            #  f.write('D_c,D_e,k_e,rmse,redchisq\n')

        #---------------
        # DIFFERENTIAL EVOLUTION
        #---------------
        print('Performing differential evolution...')
        #---- 3 Params, 0 constraints
        #  results = differential_evolution(objFunction, args=(dataFrame, t1_oneyear, t2_oneyear), bounds=bnds,tol=0.01)# tol=1e-4) #tol=1e-5)
        #---- 4 Params, 1 constraint
        num_workers = 8#4
        results = differential_evolution(objFunction, args=(dataFrame, t1_oneyear, t2_oneyear, factor), bounds=bnds, constraints=lc, tol=0.01,popsize=2,polish=False,workers=num_workers) #popsize1-->16min; popsize15-->7hrs
        print('    Done.')

        #---------------
        # GRADIENT METHOD  (MUCH MUCH FASTER!)
        #---------------
        #  results = minimize(objFunction, args=(dfi_ss.copy(), surface_area), x0=x0, bounds=bnds, method='SLSQP', tol=1e-5)
        #  results = minimize(objFunction, args=(dfi_ss.copy(), t1_oneyear, t2_oneyear, surface_area), x0=x0, bounds=bnds, method='SLSQP', tol=1e-5)

        #---- Display how long the optimization took
        end = time.perf_counter()  #[s]
        elapsed = (end-start)/60.  #[min]
        print()
        print("Elapsed time = {:.1f} min".format(elapsed))



        #--------------------------------------------------- 
        #  PLOT THE ERROR SPACE
        #--------------------------------------------------- 
        #---- Read in Data
        err_df = pd.read_csv(errfilename)
        err_name = 'redchisq'
        err_label_dic = {'redchisq': r'$\chi_{{\nu}}^2$',
                         'rmse': 'RMSE'
                        }
        #-- Calculate parameter ratios and insert into dataframe (then save)
        de_dc = err_df['D_e'] / err_df['D_c']
        ke_kc = err_df['k_e'] / err_df['k_c']
        try:
            err_df.insert(2, 'De/Dc',de_dc)
            err_df.insert(5, 'ke/kc',ke_kc)
        except ValueError:
            err_df.pop('De/Dc'); err_df.pop('ke/kc')
            err_df.insert(2, 'De/Dc',de_dc)
            err_df.insert(5, 'ke/kc',ke_kc)
        err_df.to_csv(errfilename,index=False)

        # Subset the errors to some # of best runs
        #  #---- Error less than XX% of total error range
        #  minerr = min(err_df[err_name])
        #  err_factor = 0.05#0.01          # % difference in error
        #  maxerr = minerr + (err_factor*minerr)
        #  err_df_ss = err_df.copy()[err_df[err_name] <= maxerr]
        #---- Error less than some absolute delta value
        err_factor = 0.5                # absolute delta error 
        minerr = min(err_df[err_name])
        maxerr = minerr + err_factor
        err_df_ss = err_df.copy()[err_df[err_name] <= maxerr]
        #----
        #---- Best error and parameter set
        best_err = min(err_df[err_name])
        best_pars = err_df.copy()[ err_df[err_name]==best_err ]
        print()
        print('Best Pars:')
        print('----------')
        print(best_pars)
        # Calculate relevant ratios
        de_dc = float(best_pars['D_e'] / best_pars['D_c'])
        ke_kc = float(best_pars['k_e'] / best_pars['k_c'])
        de_dc_col = err_df_ss['D_e'] / err_df_ss['D_c']
        ke_kc_col = err_df_ss['k_e'] / err_df_ss['k_c']
        #---- Relevant Statistical Regressions 
        #-- De vs Dc
        xx = np.asarray( np.log10(err_df_ss['D_c']) ).reshape(-1,1)
        yy = np.asarray( np.log10(err_df_ss['D_e']) ).reshape(-1,1)
        lm = LinearRegression()
        lm.fit(xx,yy)
        D_y_intercept = 10**lm.intercept_[0]        #y-intercept (in linear space)
        D_slope = 10**lm.coef_[0][0]                #slope       (in linear space) 
        D_r_squared = 1 - (1-lm.score(xx, yy))*(len(yy)-1)/(len(yy)-xx.shape[1]-1)
        #-- ke vs kc
        xx = np.asarray( err_df_ss['k_c'] ).reshape(-1,1)
        yy = np.asarray( err_df_ss['k_e'] ).reshape(-1,1)
        lm = LinearRegression()
        lm.fit(xx,yy)
        k_y_intercept = lm.intercept_[0]        #y-intercept (in linear space)
        k_slope = lm.coef_[0][0]                #slope       (in linear space) 
        k_r_squared = 1 - (1-lm.score(xx, yy))*(len(yy)-1)/(len(yy)-xx.shape[1]-1)
        #  #-- log(ke/kc) / log(De/Dc) 
        #  xx = np.asarray( err_df_ss['De/Dc'] ).reshape(-1,1)
        #  yy = np.asarray( err_df_ss['ke/kc'] ).reshape(-1,1)
        #  lm = LinearRegression()
        #  lm.fit(xx,yy)
        #  Dk_y_intercept = lm.intercept_[0]        #y-intercept (in linear space)
        #  Dk_slope = lm.coef_[0][0]                #slope       (in linear space) 
        #  Dk_r_squared = 1 - (1-lm.score(xx, yy))*(len(yy)-1)/(len(yy)-xx.shape[1]-1)



        print('best D_e / D_c Ratio              : {:.2f}'.format(de_dc))
        print('best k_e / k_c Ratio              : {:.2f}'.format(ke_kc))
        print('best k_c / k_photochemical        : {:.2f}'.format( float(best_pars['k_c'])/k_c_photochemical) )
        print('----')
        print('Range of D_e/D_c             : {:.2f} - {:.2f}'.format( float(min(de_dc_col)), float(max(de_dc_col)) ) )
        print('Average  D_e/D_c             : {:.2f}         '.format( np.mean(de_dc_col) ) )
        print('Slope    D_e/D_c             : {:.2f}         '.format( D_slope) )
        print('R^2      D_e/D_c             : {:.3f}         '.format( D_r_squared) )
        print('Range of k_e/k_c             : {:.2f} - {:.2f}'.format( float(min(ke_kc_col)), float(max(ke_kc_col)) ) )
        print('--')
        print('Average  k_e/k_c             : {:.2f}         '.format( np.mean(ke_kc_col) ) )
        print('Slope    k_e/k_c             : {:.2f}         '.format( k_slope) )
        print('R^2      k_e/k_c             : {:.3f}         '.format( k_r_squared) )
        print('Range of k_c/k_photochemical : {:.2f} - {:.2f}'.format(float(min(err_df_ss['k_c']))/k_c_photochemical, float(max(err_df_ss['k_c']))/k_c_photochemical))
        print()



        #-----------------------------------------------------
        #  3D PLOT (D_c, D_e, k_e/k_e, and error) 
        #-----------------------------------------------------
        fig = plt.figure(figsize=(9,7))
        ax1 = fig.add_subplot(111, projection='3d')
        cname = 'viridis' #viridis, plasma, inferno, magma 
        cm = plt.cm.get_cmap(cname)
        cm = cm.reversed()  # Make best results brighter (lower error=brighter)
        #  linnorm = matplotlib.colors.Normalize(vmin=c.min(), vmax=c.max())
        lognorm = matplotlib.colors.LogNorm()  #normalize colors to Log scale (some v. big errors)
        #  cm = plt.cm.get_cmap('viridis')
        #---- LOG SCALE
        x = err_df['D_c'];y=err_df['D_e'];z=err_df['k_e']/err_df['k_c']
        x = np.log10(x); y = np.log10(y); z = np.log10(z)
        e = ax1.scatter(x,y,z,c=err_df[err_name],vmin=min(err_df[err_name]),cmap=cm,zorder=0,s=40,norm=lognorm)
        #  e = ax1.scatter(x,y,z,c=err_df[err_name],vmin=min(err_df[err_name]),cmap=cm,zorder=0,s=40)

        ax1.set_xlabel(r'log$(D_c)$')
        ax1.set_ylabel(r'log$(D_e)$')
        ax1.set_zlabel(r'log$(k_e / k_c)$')
        #  ax1.view_init(elev=24., azim=128)
        #  ax1.set_title(titles[gas]+r', $S_w$ = '+str(sat)+'\% (tuff)')
        #
        cb = fig.colorbar(e)
        #  cb.ax.invert_yaxis()
        #  cb.ax.set_ylabel(r"Error")
        cb.ax.set_ylabel('Error, '+err_label_dic[err_name],rotation=-90,va='bottom')
        #  fig.savefig('parameterSpace_3d_{}_{}.pdf'.format(gas,satname))
        plt.tight_layout()
        fig.savefig(join(outputdir,'parameterSpace_3d.pdf'))
        plt.close('all')

        #-----------------------------------------------------
        #  3D PLOT OF BEST SUBSET  (same axes as above) 
        #-----------------------------------------------------
        fig = plt.figure(figsize=(9,7))
        ax1 = fig.add_subplot(111, projection='3d')
        cname = 'viridis' #viridis, plasma, inferno, magma 
        cm = plt.cm.get_cmap(cname)
        cm = cm.reversed()  # Make best results brighter (lower error=brighter)
        #  cm = plt.cm.get_cmap('viridis')
        #---- LOG SCALE
        #  x = err_df_ss['D_c'];y=err_df_ss['D_e'];z=err_df_ss['k_e']
        # Try making same as non-subset plot using subset indices
        x = x[err_df_ss.index];y=y[err_df_ss.index];z=z[err_df_ss.index]
        e = ax1.scatter(x,y,z,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm,zorder=0,s=40)

        ax1.set_xlabel(r'log$(D_c)$')
        ax1.set_ylabel(r'log$(D_e)$')
        #  ax1.set_zlabel(r'log$(k_e)$')
        ax1.set_zlabel(r'log$(k_e / k_c)$')
        #  ax1.view_init(elev=24., azim=128)
        #  ax1.set_title(titles[gas]+r', $S_w$ = '+str(sat)+'\% (tuff)')
        #  ax1.set_title(r'Subset Parameter Space (best {:.1f}\%)'.format(err_factor*100)) #PERCENT
        ax1.set_title(r'Subset Parameter Space (w/in {:.1f} of $min(${}$)$ )'.format(err_factor, err_label_dic[err_name]))  #ABSOLUTE DELTA ERROR
        #
        cb = fig.colorbar(e)
        #  cb.ax.invert_yaxis()
        #  cb.ax.set_ylabel(r"Error")
        cb.ax.set_ylabel('Error, '+err_label_dic[err_name],rotation=-90,va='bottom')
        plt.tight_layout()
        fig.savefig(join(outputdir,'parameterSpace_3d-subset.pdf'))
        plt.close('all')

        #-----------------------------------------------------
        #  2D PLOTS  (SUBSET)
        #-----------------------------------------------------
        fig,ax = plt.subplots(3,figsize=(5,12))
        #  fig,ax = plt.subplots(2,figsize=(5,8))
        alpha = 1.0#0.4  #how opaque
        # [0] D_e vs D_c  -------------------------------------------
        #---- LOG SCALE
        x = err_df_ss['D_c'];y=err_df_ss['D_e']
        # (LOG SCALE?)
        x = np.log10(x); y = np.log10(y)
        #  sc = ax[0].scatter(x,y,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm, alpha=alpha, norm=lognorm)
        sc = ax[0].scatter(x,y,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm, alpha=alpha)
        #--
        ax[0].set_xlabel(r'log$(D_c)$')
        ax[0].set_ylabel(r'log$(D_e)$')

        # [1] k_e vs k_c  -------------------------------------------
        #---- LOG SCALE
        x = err_df_ss['k_c'];y=err_df_ss['k_e']
        # (LOG SCALE?)
        #  x = np.log10(x); y = np.log10(y)
        #  sc = ax[1].scatter(x,y,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm, alpha=alpha, norm=lognorm)
        sc = ax[1].scatter(x,y,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm, alpha=alpha)
        #--
        ax[1].set_xlabel(r'$k_c$')
        ax[1].set_ylabel(r'$k_e$')
        #  ax[1].set_xlabel(r'log$(k_c)$')
        #  ax[1].set_ylabel(r'log$(k_e)$')

        # [2] (D_e/D_c) vs (k_e/k_c) --------------------------------
        x = err_df_ss['De/Dc']
        y = err_df_ss['ke/kc']
        # (LOG SCALE?)
        x = np.log10(x)
        y = np.log10(y)

        #  sc = ax[2].scatter(x,y,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm, alpha=alpha, norm=lognorm)
        sc = ax[2].scatter(x,y,c=err_df_ss[err_name],vmin=min(err_df_ss[err_name]),cmap=cm, alpha=alpha)
        ax[2].set_xlabel(r'log$(D_e$/$D_c$)')
        #  ax[2].set_ylabel(r'$k_e$/$k_c$')
        ax[2].set_ylabel(r'log$(k_e$/$k_c)$')

        #  ax[2].set_ylabel(r'log$(D_e$/$D_c)$')
        #  ax[2].set_ylim(bottom=1.0) # D_e/D_c ratio is at minimum=1.0 (D_e ≥ D_c) 

        #---- Properties
        fig.subplots_adjust(right=0.75)
        fig.subplots_adjust(left=.2)
        #  ax[0].set_title(titles[gas]+r', $S_w$ = '+str(sat)+'\% (tuff)')
        cbar_ax = fig.add_axes([0.8, 0.15, 0.04, 0.7])
        cb = fig.colorbar(sc,cax=cbar_ax)
        #-- Turn these off if log colorbar
        #  tick_locator = ticker.MaxNLocator(nbins=10)
        #  cb.locator = tick_locator
        #  cb.update_ticks()
        #--
        cb.ax.set_ylabel('Error, '+err_label_dic[err_name],rotation=-90,va='bottom')
        #  plt.tight_layout()
        fig.savefig(join(outputdir,'parameterSpace_2d-subset.pdf'))
        plt.close('all')

        #-----------------------------------------------------
        #  2D ERROR CONTOUR PLOTS  (SUBSET)
        #-----------------------------------------------------
        fig,ax = plt.subplots(3,figsize=(5,12))
        alpha = 0.4  #how opaque
        nlevels = 25
        # [0] D_e vs D_c  -------------------------------------------
        x = err_df_ss['D_c'];y=err_df_ss['D_e']
        # (LOG SCALE?)
        x = np.log10(x); y = np.log10(y)
        z = err_df_ss[err_name]

        triang = tri.Triangulation(x,y)
        sc = ax[0].tricontourf(triang, z, nlevels, vmin=min(err_df_ss[err_name]),cmap=cm)
        #  sc = ax[0].tricontourf(triang, z, nlevels, vmin=min(err_df_ss[err_name]),cmap=cm,norm=lognorm)
        #----------
        # Try to remove weird lines between contours
        # < https://discourse.matplotlib.org/t/unwanted-lines-between-contourf-contour-levels/12539/2 >
        #     make sure the ax.contourf call is the same as above (can't use "cs"...)
        for c in ax[0].tricontourf(triang, z,nlevels, vmin=min(err_df_ss[err_name]),cmap=cm).collections:
            c.set_linewidth(0.1)
        #----------
        #--
        ax[0].set_xlabel(r'log$(D_c)$')
        ax[0].set_ylabel(r'log$(D_e)$')

        # [1] k_e vs k_c  -------------------------------------------
        #---- LOG SCALE
        x = err_df_ss['k_c'];y=err_df_ss['k_e']
        # (LOG SCALE?)
        #  x = np.log10(x); y = np.log10(y)
        z = err_df_ss[err_name]

        triang = tri.Triangulation(x,y)
        sc = ax[1].tricontourf(triang, z, nlevels, vmin=min(err_df_ss[err_name]),cmap=cm)
        #  sc = ax[1].tricontourf(triang, z, nlevels, vmin=min(err_df_ss[err_name]),cmap=cm,norm=lognorm)
        #----------
        # Try to remove weird lines between contours
        # < https://discourse.matplotlib.org/t/unwanted-lines-between-contourf-contour-levels/12539/2 >
        #     make sure the ax.contourf call is the same as above (can't use "cs"...)
        for c in ax[1].tricontourf(triang, z,nlevels, vmin=min(err_df_ss[err_name]),cmap=cm).collections:
            c.set_linewidth(0.1)
        #----------
        #--
        ax[1].set_xlabel(r'$k_c$')
        ax[1].set_ylabel(r'$k_e$')
        #  ax[1].set_xlabel(r'log$(k_c)$')
        #  ax[1].set_ylabel(r'log$(k_e)$')

        # [2] (D_e/D_c) vs (k_e/k_c) ----------------------------------
        x = err_df_ss['De/Dc']
        y = err_df_ss['ke/kc']
        # (LOG SCALE?)
        x = np.log10(x); y = np.log10(y)
        z = err_df_ss[err_name]

        triang = tri.Triangulation(x,y)
        sc = ax[2].tricontourf(triang, z, nlevels, vmin=min(err_df_ss[err_name]),cmap=cm)
        #  sc = ax[2].tricontourf(triang, z, nlevels, vmin=min(err_df_ss[err_name]),locator=ticker.LogLocator(),cmap=cm,norm=lognorm)
        #----------
        # Try to remove weird lines between contours
        # < https://discourse.matplotlib.org/t/unwanted-lines-between-contourf-contour-levels/12539/2 >
        #     make sure the ax.contourf call is the same as above (can't use "cs"...)
        for c in ax[2].tricontourf(triang, z,nlevels, vmin=min(err_df_ss[err_name]),cmap=cm).collections:
            c.set_linewidth(0.1)
        #----------
        #  ax[2].set_xlabel(r'$D_e$/$D_c$')
        #  ax[2].set_ylabel(r'$k_e$/$k_c$')
        ax[2].set_xlabel(r'log$(D_e$/$D_c$)')
        ax[2].set_ylabel(r'log$(k_e$/$k_c)$')
        #  ax[2].set_ylim(bottom=1.0) # D_e/D_c ratio is at minimum=1.0 (D_e ≥ D_c) 

        #---- Properties
        fig.subplots_adjust(right=0.75)
        fig.subplots_adjust(left=.2)
        #  ax[0].set_title(titles[gas]+r', $S_w$ = '+str(sat)+'\% (tuff)')
        cbar_ax = fig.add_axes([0.8, 0.15, 0.04, 0.7])
        cb = fig.colorbar(sc,cax=cbar_ax)
        #  tick_locator = ticker.MaxNLocator(nbins=10)
        #  cb.locator = tick_locator
        #  cb.update_ticks()
        cb.ax.set_ylabel('Error, '+err_label_dic[err_name],rotation=-90,va='bottom')
        #  plt.tight_layout()
        fig.savefig(join(outputdir,'parameterSpace_2d-contour-subset.pdf'))
        plt.close('all')
        #--------------------------------------------------- 
        # END ERROR SPACE PLOTS


    # ==================================================
    # BEGIN NORMAL (TIME SERIES) PLOTS
    # ==================================================

    # -----------------------------------------------
    # Run the model again with the best Params to get FULL Time Series
    # -----------------------------------------------


    # NORMAL CASE ~~~~~~~~~~~~~~~~~~~~
    #  if args.plotOnly is False:
    if not args.plotOnly:
        #---- Read in best_pars from Saved Data File (errors.csv)
        err_df = pd.read_csv(errfilename)
        # Best error and parameter set (best_pars)
        best_err = min(err_df[err_name])
        best_pars = err_df.copy()[ err_df[err_name]==best_err ]
        pars = best_pars
        #  #---- Run the diffusion model using best_pars
        #  c_df, c_all, c_obs, C_ppbv = diffusion_PBL(outputdir, dataFrame, D_collapsed=float(best_pars['D_c']), D_expanded=float(best_pars['D_e']), k_c=float(best_pars['k_c']), k_e=float(best_pars['k_e']), delx=delx, delt=delt, h_obs=h_obs)
    # RUN WITH SPECIFIED PARS ~~~~~~~~
    #  elif args.plotOnly is True:
    elif args.plotOnly:
        # Read in the pars from file provided by argument
        scenario_pars = pd.read_csv(args.plotOnly)
        pars = scenario_pars
        #  pars = pars[ pars['redchisq']==min(pars['redchisq'])]
        #---- Modify the outputdir to a separate sub-directory (e.g. ``scenarioA``)
        # Make sure path hasn't already been added
        if outputdir.count('lclscratch')==1:
            outputdir = join(outputdir, special_id)
        else:
            input('WARNING: scenario output path as already been added')
        # Create outputdir if it doesn't exist
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        # Copy the scenario file into scenario<> dir for posterity
        shutil.copy( args.plotOnly, join(outputdir,os.path.split(args.plotOnly)[-1]) )
        #  c_df, c_all, c_obs, C_ppbv = diffusion_PBL(outputdir, dataFrame, D_collapsed=float(scenario_pars['D_c']), D_expanded=float(scenario_pars['D_e']), k_c=float(scenario_pars['k_c']), k_e=float(scenario_pars['k_e']), delx=delx, delt=delt, h_obs=h_obs)

    #---- Run the diffusion model using the specified pars 
    c_df, c_all, c_obs, C_ppbv = diffusion_PBL(outputdir, dataFrame, D_collapsed=float(pars['D_c']), D_expanded=float(pars['D_e']), k_c=float(pars['k_c']), k_e=float(pars['k_e']), delx=delx, delt=delt, h_obs=h_obs)

    # Write periodically to a summary file with useful diagnostics
    summaryFile = join(outputdir,'summaryFile.md')
    # Create file here. Append to it throughout the rest of the script
    with open(summaryFile, 'w') as f:
        if not args.plotOnly:
            case='Best Case'
            desc='Case with best overall error.'
        else:
            case=special_id
        f.write('SUMMARY of {}\n'.format(case))
        f.write('==================================================\n\n')
        f.write('Parameters Used:\n')
        f.write('--------------------------------------------------\n')
        f.write(pars.to_string())
        f.write('\n\n')

    #  # -----------------------------------------------
    #  # Scale the FLUX so that mean of observation points is same as SAM-TLS background mean  (0.41ppbv)
    #  # -----------------------------------------------
    #  C_background = 0.41 #[ppbv]
    #  #  #---- Get Webster2018a abundances for comparison
    #  #  em = readAbundances(run_type='Enrichment',source='Webster2021')
    #  #  #  # Omit the outlier with C=20.5ppmv
    #  #  em = em.copy()[ em['GlobalMeanCH4(ppbv)']<20.]
    #  #  # Also omit the non-detections (they bring down the background avg) 
    #  #  em = em.copy()[ em['GlobalMeanCH4(ppbv)']>0.068]
    #  #  em = em.sort_values('Ls(deg)')
    #  #  #---- Subset the modeled values to the times specified (1 year)
    #  #  df_ss = c_df.copy()[ c_df['SOL_TIME'].between(t1_oneyear,t2_oneyear) ]
    #  #  #---- Interpolate sim concs at the time of observations
    #  #  sim_interp_vals = interpolate_model2obsAbundance(df_ss, em, t1_oneyear, t2_oneyear)
    #  #  #---- Calculate by what factor we need to reduce fluxes
    #  #  #  factor = np.mean(ci_dectime)/C_background
    #  #  factor = np.mean(sim_interp_vals['C_ppbv'])/C_background
    #  #  print('    (Reduction factor = {:.3f})'.format(factor))
    #  #  #---- Modify the fluxes and abundance
    #  #  C_ppbv_orig = c_df['C_ppbv'].copy()
    #  #  flux_orig   = c_df['FLUX_kg/m2/s'].copy()
    #  #  flux_adj = flux_orig / factor
    #  #  C_ppbv_adj = C_ppbv_orig / factor
    #  #  #  C_ppbv_adj = c_df['C_ppbv'] / c_df['FLUX_kg/m2/s'] / factor
    #  #  #  C_ppbv_adj = c_df['C_ppbv'] / factor
    #  #  #---- Add adjusted values into DataFrame
    #  #  try:
        #  #  pos1 = c_df.columns.get_loc('C_ppbv')
        #  #  c_df.insert(pos1+1, 'C_ppbv_ADJ',C_ppbv_adj)
        #  #  pos2 = c_df.columns.get_loc('FLUX_kg/m2/s')
        #  #  c_df.insert(pos2+1, 'FLUX_kg/m2/s_ADJ', flux_adj)
    #  #  except ValueError:
        #  #  c_df.pop('C_ppbv_ADJ'); c_df.pop('FLUX_kg/m2/s_ADJ')
        #  #  pos1 = c_df.columns.get_loc('C_ppbv')
        #  #  c_df.insert(pos1, 'C_ppbv_ADJ',C_ppbv_adj)
        #  #  pos2 = c_df.columns.get_loc('FLUX_kg/m2/s')
        #  #  c_df.insert(pos2, 'FLUX_kg/m2/s_ADJ', flux_adj)
    #  #  #---- For now, just make C_ppbv and FLUX equal to their adjusted values
    #  #  c_df['C_ppbv'] = c_df['C_ppbv_ADJ']
    #  #  c_df['FLUX_kg/m2/s'] = c_df['FLUX_kg/m2/s_ADJ']
    #  #---- Modify the fluxes and abundance
    #  C_ppbv_orig = c_df['C_ppbv'].copy()
    #  flux_orig   = c_df['FLUX_kg/m2/s'].copy()
    #  flux_adj = flux_orig / float(best_pars['factor'])
    #  C_ppbv_adj = C_ppbv_orig / float(best_pars['factor'])
    #  #  C_ppbv_adj = c_df['C_ppbv'] / c_df['FLUX_kg/m2/s'] / factor
    #  #  C_ppbv_adj = c_df['C_ppbv'] / factor
    #  #---- Add adjusted values into DataFrame
    #  try:
        #  pos1 = c_df.columns.get_loc('C_ppbv')
        #  c_df.insert(pos1+1, 'C_ppbv_ADJ',C_ppbv_adj)
        #  pos2 = c_df.columns.get_loc('FLUX_kg/m2/s')
        #  c_df.insert(pos2+1, 'FLUX_kg/m2/s_ADJ', flux_adj)
    #  except ValueError:
        #  c_df.pop('C_ppbv_ADJ'); c_df.pop('FLUX_kg/m2/s_ADJ')
        #  pos1 = c_df.columns.get_loc('C_ppbv')
        #  c_df.insert(pos1, 'C_ppbv_ADJ',C_ppbv_adj)
        #  pos2 = c_df.columns.get_loc('FLUX_kg/m2/s')
        #  c_df.insert(pos2, 'FLUX_kg/m2/s_ADJ', flux_adj)
    #  #---- For now, just make C_ppbv and FLUX equal to their adjusted values
    #  c_df['C_ppbv'] = c_df['C_ppbv_ADJ']
    #  c_df['FLUX_kg/m2/s'] = c_df['FLUX_kg/m2/s_ADJ']




    # 1-sol
    C_mavgs_1 = calc_moving_avg(c_df, window_size=1.0)
    # 2-sol
    C_mavgs_2 = calc_moving_avg(c_df, window_size=2.0)
    # 10-sol
    C_mavgs_10 = calc_moving_avg(c_df, window_size=10.0)
    # 668-sol
    C_mavgs_668 = calc_moving_avg(c_df, window_size=668.)
    #  # Calculate moving average in a given window  of size mv_win
    #  mv_win = 1. #[sol] window for moving avg
    #  dt = np.diff(c_df['SOL_TIME'])[0]
    #  win = int(round(mv_win/dt))                #mv-win-sol window
    #  #  C_windows = C_ppbv.rolling(win)
    #  C_ppbv = c_df.copy()['C_ppbv']
    #  C_windows = C_ppbv.rolling(win, center=True)
    #  # Calculate average of each window
    #  C_mavgs = C_windows.mean()
    #  #---- Also calculate for a window of 1 Mars year
    #  mv_win_annual = 668. #[sol]
    #  win_annual = int(round(mv_win_annual/dt))                #mv-win-sol window
    #  C_windows_annual = C_ppbv.rolling(win_annual, center=True)
    #  # Calculate average of each window
    #  C_mavgs_annual = C_windows_annual.mean()
    #  print( 'Mean Conc = {:.4} ppbv'.format(np.mean(C_mavgs)) )
    #-------------------------------------------------- 
    # Calculate Night-Time Average Abundance (for plot) 
    #-------------------------------------------------- 
    # Define "night-time" as between 0.0 and 2.0 decimal time (want same-day)
    nightDecTime = (0., 2.0)
    # Get night abundances for each SOL_INT 
    sol_int = c_df['SOL_INT']
    meanNightConc_arr = np.zeros( (len(np.unique(sol_int)), 3) )
    array = np.zeros( (len(np.unique(sol_int)), 3) )
    for i,s in enumerate(np.unique(sol_int)):
        # Subset only the current SOL_INT for Night-Time hours
        ss_si = c_df.copy()[ (c_df['SOL_INT']==s) & (c_df['DECIMAL_TIME'].between(nightDecTime[0], nightDecTime[1])) ]
        # Calculate the mean of night-time Concentrations
        meanNightConc = np.mean(ss_si['C_ppbv'])
        #  meanNightConc = np.mean(ss_si['PBL_m'])
        # Store in array [ cols= SOL_INT, L_s, mean(C) ]
        meanNightConc_arr[i,0]  = int(s)
        meanNightConc_arr[i,1]  = np.mean(ss_si['L_s'])
        meanNightConc_arr[i,-1] = meanNightConc
    # Create a DataFrame with the night averages in it
    c_night_df = pd.DataFrame(meanNightConc_arr, columns=['SOL_INT', 'L_s', 'C_ppbv'])
    c_night_df['SOL_INT'] = c_night_df['SOL_INT'].astype('int')
    #  c_night_df


    # Get Webster2018a abundances for comparison
    em = readAbundances(run_type='Enrichment',source='Webster2021')
    #  # Omit the outlier with C=20.5ppmv
    em = em.copy()[ em['GlobalMeanCH4(ppbv)']<20.]
    #  # Also omit the non-detections (will be dealt with later)
    #  em = em.copy()[ em['GlobalMeanCH4(ppbv)']>0.068]
    em = em.sort_values('Ls(deg)')

    #-------------------------------------------------- 
    # Calculate Error b/w simulated Abundance and observation
    #-------------------------------------------------- 
    dof = np.shape(bnds)[0]
    rmse, chisq, redchisq = calcAbundError(c_df.copy(), em, time1=t1_oneyear, time2=t2_oneyear,dof=dof)
    rmse_text     = 'RMSE = {:.2f}'.format(rmse)
    chisq_text    = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(chisq)
    redchisq_text = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)
    #  redchisq_text = r'Reduced $\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)
    print('RMSE               = {}'.format(rmse))
    print('Chi-Square         = {}'.format(chisq))
    print('Reduced Chi-Square = {}'.format(redchisq))


    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    # MULTI-YEAR TEST PLOT
    fig, ax1 = plt.subplots(1, sharex=True, figsize=(12,5))
    #  # [1] Modeled Abundance ---------------------------------------------
    ax1.plot(c_df['SOL_TIME'], c_df['C_ppbv'], color='k')
    plt.savefig(join(outputdir, 'TEST_multiYear_atm_abundance_plot.pdf'))
    plt.close('all')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

    #-------------------------------------------------- 
    # Plot atmospheric abundance (atm_abundance_and_flux_plot.pdf) 
    #-------------------------------------------------- 
    # Get a copy of c_df that is just one year (call it 'df_copy')
    df_copy = c_df.copy()[c_df['SOL_TIME'].between(t1_oneyear, t2_oneyear)]
    # Sort by L_s
    df_copy = df_copy.sort_values(by='L_s')
    #  df_copy = c_df.copy()
    fig, (ax1,ax2) = plt.subplots(2, sharex=False, figsize=(8,10))
    #  # [1] Modeled Abundance ---------------------------------------------
    ax1.plot(df_copy['L_s'], df_copy['C_ppbv'], color='lightgrey',label='simulated',zorder=0)
    #  # Moving average
    #  ax1.plot(df_copy['L_s'], C_mavgs_1[df_copy['C_ppbv'].index], color='k', lw=3, label='sim. moving avg. 1')
    #  #  ax1.plot(df_copy['L_s'], C_mavgs_2[df_copy['C_ppbv'].index], color='b', lw=3, label='sim. moving avg. 2')
    #  #  ax1.plot(df_copy['L_s'], C_mavgs_10[df_copy['C_ppbv'].index], color='g', lw=3, label='sim. moving avg. 10')
    #  TRY PLOTTING NIGHT-TIME AVERAGE CONCENTRATION ~~~~~~~~~~~~~~~~~~~~~~
    night_df = c_night_df.copy()[ c_night_df['SOL_INT'].between(int(min(df_copy['SOL_TIME'])), int(max(df_copy['SOL_TIME']))) ]
    night_df = night_df.sort_values(by='L_s')
    ax1.plot(night_df['L_s'], night_df['C_ppbv'], color='black', lw=3, label='night-time avg.')
    #  ax1.plot(night_df['L_s'], night_df['C_ppbv'], color='purple', lw=3, label='night-time avg.')
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # [ALSO] Webster2018a Abundances -----------------------------------
    ax1.errorbar(em['Ls(deg)'], em['GlobalMeanCH4(ppbv)'], yerr=em['Error(±1SEM).1'],marker='o',ms=10,capsize=5,capthick=1,ls='',color='red',mec='k',label=r'Webster \textit{et al.} (2021)',zorder=50)
    # [2] Flux -----------------------------------
    ax2.plot(df_copy['L_s'], df_copy['FLUX_kg/m2/s'], color='k')
    #  # Adjusted Flux
    #  ax2.plot(df_copy['L_s'], df_copy['FLUX_ADJ'], color='k')
    # Properties
    #[1]
    ax1.legend()
    #  max_y = max( df_copy['C_ppbv'] )
    #  max_y = max( max(df_copy['C_ppbv']), max(em['GlobalMeanCH4(ppbv)']) )
    max_y = max( max(df_copy['C_ppbv']), max(em['GlobalMeanCH4(ppbv)']+em['Error(±1SEM).1']) )
    ax1.set_ylim(0,max_y*1.04) #TESTING
    #  ax1.set_ylim(0,1.7) #TESTING
    ax1.set_xlabel(r'\textbf{Solar Longitude, L$_s$ [°]}')
    ax1.set_xlim(0., 360.)
    #  ax1.set_xlim(160., 165.)
    ax1.set_ylim(bottom=0.)
    ax1.xaxis.set_minor_locator(MultipleLocator(30)) #15-degree minor ticks 
    ax1.xaxis.set_major_locator(MultipleLocator(60)) #30-degree major ticks 
    ax1.set_ylabel(r'\textbf{CH$_4$ Abundance [ppbv]}')
    #[2]
    ax2.set_xlim(0., 360.)
    ax2.set_ylim(bottom=0.)
    ax2.xaxis.set_minor_locator(MultipleLocator(30)) #15-degree minor ticks 
    ax2.xaxis.set_major_locator(MultipleLocator(60)) #30-degree major ticks 
    ax2.set_xlabel(r'\textbf{Solar Longitude, L$_s$ [°]}')
    ax2.set_ylabel(r'\textbf{Surface Flux [kg/m$^2$/s]}')
    #-- Print the mean abundance on plot
    meanCppbv = df_copy['C_ppbv'].mean()
    # What we really want is the mean abundance at time of night-time observations...
    obs_night = em.copy()[(em['DECIMAL_TIME'].between(0,4)) | (em['DECIMAL_TIME']>19.)]
    meanCobs_night = np.mean(obs_night['GlobalMeanCH4(ppbv)'])
    # (A) Observations
    xa = np.asarray(obs_night['Ls(deg)'])
    ya = np.asarray(obs_night['GlobalMeanCH4(ppbv)'])
    # (B) Simulated 
    sim_interp_vals = interpolate_model2obsAbundance(df_copy, obs_night)
    ybi = np.asarray(sim_interp_vals['C_ppbv']); ybi[ybi<0] = 0.
    meanCsim_night = np.mean(ybi)

    y0 = 0.95
    #  ax1.text(0.1,0.95, r'$C_{{mean}}$ = {:.2f} ppbv'.format(meanCppbv), transform=ax1.transAxes)
    ax1.text(0.1,0.95, r'$C_{{mean}}$ = {:.2f} ppbv'.format(meanCsim_night), transform=ax1.transAxes)
    vshift=0.05
    #  ax1.text(0.1,y0-vshift,   rmse_text,     transform=ax1.transAxes)
    #  ax1.text(0.1,y0-vshift*2, redchisq_text, transform=ax1.transAxes)
    ax1.text(0.1,y0-vshift, redchisq_text, transform=ax1.transAxes)
    #  plt.tight_layout()
    plt.savefig(join(outputdir, 'atm_abundance_and_flux_plot.pdf'))
    plt.close('all')

    #-------------------------------------------------- 
    # ADD TO summaryFile 
    #-------------------------------------------------- 
    with open(summaryFile, 'a') as f:
        f.write('Overall Error (Annual):\n')
        f.write('--------------------------------------------------\n')
        f.write('Reduced Chi-Square : {:.4f}\n'.format(redchisq))
        f.write('\n')


    #================================================== 
    # Plot all abundances against LMST 
    #================================================== 
    #-------------------------------------------------- 
    # Read in Webster2021 abundances  
    #-------------------------------------------------- 
    em = readAbundances(run_type='Enrichment', source='Webster2021')
    # Remove 20pppv measurement
    em = em.copy()[ em['GlobalMeanCH4(ppbv)']<20 ]
    #-------------------------------------------------- 
    # Read Simulated Abundances 
    #-------------------------------------------------- 
    df_all = df_copy.copy()#[ df_copy['L_s'].between(90,180) ]

    #-------------------------------------------------- 
    # Plot day-night Measured vs Simulated Atm Abundances for All Obs
    #-------------------------------------------------- 
    print('Plotting day-night differences in atm abundance for all observations...')
    #  gridspec = {'width_ratios': [1, 1, 0.1, 0.1]}  #colorbar axes shoudld be small
    gridspec = {'width_ratios': [1, 0.1]}  #colorbar axes should be small
    # [Subplot Map]
    #--------------
    # [0] [1]
    # [2] [3] 
    #--------------
    # 0=top panel,    1=colorbar, 
    # 2=bottom panel, 3=blank, 4=second axis on [2]
    pblcount = 0  #counter for plotting only 1 PBL height (for non-detects)
    fig, ax = plt.subplots(2,2, sharex=False,figsize=(7,10), gridspec_kw=gridspec)
    axs = ax.ravel()
    ax4 = axs[2].twinx()
    # Assign colors increasing with L_s(°)
    cmapName='plasma' #uniform
    cmap=matplotlib.cm.get_cmap(cmapName)
    # Custom colorbar delineating L_s of each point
    cbar_ticks = 60
    c = np.arange(0, 360+cbar_ticks, cbar_ticks)
    norm = matplotlib.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cmap.set_array([])
    # Measured Abundances
    #Save the Abundance values interpolated to the observation times 
    obs_abund_dict = {}
    sim_abund_dict = {}
    for idx, row in em.iterrows():
        targ_ls = row['Ls(deg)']
        targ_dec_time = row['DECIMAL_TIME']

        # [0] ABUNDANCES -----------------------------------------------
        circle_size=5#10
        star_size=10#18
        #---- Measured Abundances
        x = em['DECIMAL_TIME'].loc[idx]
        y = em['GlobalMeanCH4(ppbv)'].loc[idx]
        yerr = em['Error(±1SEM).1'].loc[idx]
        axs[0].errorbar(x,y, yerr=yerr,marker='o',ms=circle_size,capsize=5,ls='',color=cmap.to_rgba(targ_ls),mec='k')

        #---- Simulated Abundances
        # Find closest simulated L_s to target obs L_s
        fuzz = 0.02
        a = df_all.copy()[ df_all['L_s'].between(targ_ls-fuzz, targ_ls+fuzz) ]
        if len(a)>1: a = a.iloc[-1]
        #  if len(a)>1: a = a.iloc[0]
        # Get the SOL_TIME of that day and get that whole day
        targ_sol = a['SOL_TIME']
        b = df_all.copy()[ df_all['SOL_TIME'].between(int(targ_sol), int(targ_sol)+1.0) ]
        #make midnight next day equal to hour 24
        if b['SOL_TIME'].iloc[-1] == int(targ_sol)+1.0:
            b['DECIMAL_TIME'].iloc[-1] = 24.00
        #  # Interpolate Abundance to time of observation
        # Numpy Interpolate seems to work better for this 
        x = targ_dec_time
        y = np.interp(x, b['DECIMAL_TIME'], b['C_ppbv'])
        axs[0].errorbar(x,y, yerr=0., marker='*',ms=star_size,capsize=0,ls='',color=cmap.to_rgba(targ_ls),mec='k')
        #  # [2] FLUXES AND PBL HEIGHT ------------------------------------
        #  # (for now, only do for the non-detects...)
        #  if em['GlobalMeanCH4(ppbv)'].loc[idx] < 0.07:
            #  #---- FLUXES
            #  axs[2].plot(b['DECIMAL_TIME'], b['FLUX_kg/m2/s'], ls='-', color=cmap.to_rgba(targ_ls), label=r'flux (L$_s$ = {}°)'.format(targ_ls))
            #  #  axs[2].plot(b['DECIMAL_TIME'], b['FLUX_ADJ'], ls='-', color=cmap.to_rgba(targ_ls), label=r'flux (L$_s$ = {}°)'.format(targ_ls))
            #  #  axs[2].semilogy(b['DECIMAL_TIME'], b['FLUX_ADJ'], ls='-', color=cmap.to_rgba(targ_ls), label=r'flux (L$_s$ = {}°)'.format(targ_ls))         #LOG SCALE FLUXES (?)
            #  #---- PBL HEIGHT (only need one, should be same)
            #  if pblcount==0:
                #  ax4.plot(b['DECIMAL_TIME'], b['PBL_m'], ls=':', color='k', label='PBL')
                #  #  ax4.plot(b['DECIMAL_TIME'], b['PBL_m'], ls=':', color=cmap.to_rgba(targ_ls), label='PBL')
            #  #---- PRESSURE (don't need a scale, only plot one)
                #  #  ax10 = ax4.twinx()
                #  ax10 = axs[2].twinx()
                #  ax10.plot(b['DECIMAL_TIME'], b['PRESS'], ls='-', color='blue', label='pressure')
                #  # Note: Don't plot the yticks or labels for pressure (can mention in caption)
                #  plt.setp(ax10.get_yticklabels(), visible=False)
                #  ax10.tick_params(axis='y', which='both', length=0)
                #  pblcount+=1
        # Save values to dictionary so they're easy to access if needed
        obs_abund_dict[targ_ls] = {'DECIMAL_TIME': targ_dec_time,
                                   'SOL'         : int(targ_sol),
                                   'C_ppbv'      : em['GlobalMeanCH4(ppbv)'].loc[idx],
                                   #  'PBL_m'       : np.interp(x,b['DECIMAL_TIME'], b['PBL_m']),
                                   'Error(±1SEM).1'  : em['Error(±1SEM).1'].loc[idx],
                                  }
        sim_abund_dict[targ_ls] = {'DECIMAL_TIME': targ_dec_time,
                                   'SOL'         : int(targ_sol),
                                   'C_ppbv'      : y,
                                   'PBL_m'       : np.interp(x,b['DECIMAL_TIME'], b['PBL_m']),
                                  }
    # Calculate the error b/w sim and obs
    rmse, chisq, redchisq = calcAbundError(df_all, em,  time1=t1_oneyear, time2=t2_oneyear, dof=dof)
    # Write error to file
    with open(join(outputdir,'allDayNight_error'), 'w') as f:
        f.write('RMSE              = {:.4f}\n'.format(rmse))
        f.write('Chi-square        = {:.4f}\n'.format(chisq))
        f.write('Reduced Chi-square = {:.4f}\n'.format(redchisq))
    # Create text fields if desired
    rmse_text     = 'RMSE = {:.2f}'.format(rmse)
    chisq_text    = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(chisq)
    redchisq_text = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)
    #  redchisq_text = r'Reduced $\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)

    #  # CALCULATE ERROR, OMITTING Ls=103.48° [idx=0]
    #  #-- Calculate RMSE
    #  rmse = calc_rmse(obs[1:], sim[1:])
    #  #-- Calculate Chi-Square Statistic 
    #  obs_sem = np.asarray([obs_abund_dict[c]['Error(±1SEM).1'] for c in obs_abund_dict])[1:]
    #  #  chisq = chisqg(ydata=ya, ymod=ybi, sd=None)
    #  chisq = chisqg(ydata=obs[1:], ymod=sim[1:], sd=obs_sem)
    #  #-- Calculate Reduced Chi-Square Statistic 
    #  dof = np.shape(bnds)[0]#2 # no. of degrees of freedom
    #  #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=None)
    #  redchisq = redchisqg(ydata=obs[1:], ymod=sim[1:], deg=dof, sd=obs_sem)
    #  #  # Write error to file
    #  #  with open(join(outputdir,'summerDayNight_error_omitLs_103'), 'w') as f:
        #  #  f.write('RMSE              = {:.4f}\n'.format(rmse))
        #  #  f.write('Chi-square        = {:.4f}\n'.format(chisq))
        #  #  f.write('Reduced Chi-square = {:.4f}\n'.format(redchisq))

    # [1] Colorbar
    cax = axs[1]
    cb1 = plt.colorbar(cmap, ticks=c, cax=cax)
    cb1.set_label(r'\textbf{Solar Longitude, L$_s$ [°]}', rotation=-90, va='bottom')
    # Make [3] invisible
    axs[3].set_visible(False)
    #  #-- Properties
    # [0]
    axs[0].set_title('All Day-Night Variations')
    vshift=0.05; y0 = 0.95
    axs[0].text(0.97, y0,rmse_text, transform=axs[0].transAxes, ha='right')
    axs[0].text(0.97, y0-vshift,redchisq_text, transform=axs[0].transAxes, ha='right')
    # Legend
    # Only need legend for one marker each
    axs[0].lines[0].set_label('observed')
    axs[0].lines[3].set_label('simulated')
    #  for i in range(10):
        #  axs[0].lines[i].set_label(i)
    axs[0].legend(loc='upper center')
    leg = axs[0].get_legend()
    # Force legend colors to be black (doesn't seem to work...)
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('black')
    axs[0].set_ylabel(r'\textbf{CH$_4$ Abundance [ppbv]}')
    axs[0].set_xlabel(r'\textbf{LMST}')
    #  axs[0].set_ylim(0., max(max_y) )
    #  axs[0].set_ylim(bottom=0.)
    #  ax1.legend()
    pad=0.5
    axs[0].set_xlim(0-pad, 24+pad)
    axs[0].xaxis.set_major_locator(MultipleLocator(2)) #2-hour major ticks 
    axs[0].xaxis.set_minor_locator(MultipleLocator(1)) #1-hour minor ticks 
    # [2]
    axs[2].set_ylim(bottom=0)
    axs[2].set_xlabel(r'\textbf{LMST}')
    axs[2].set_ylabel(r'\textbf{Surface Flux [kg/m$^2$/s]}')
    axs[2].set_xlim(0-pad, 24+pad)
    axs[2].xaxis.set_major_locator(MultipleLocator(2)) #2-hour major ticks 
    axs[2].xaxis.set_minor_locator(MultipleLocator(1)) #1-hour minor ticks 
    #  axs[2].set_yscale('semilog')  #LOG SCALE FLUXES (?)
    # [ax4]
    ax4.set_ylim(bottom=0)
    ax4.set_ylabel(r'\textbf{PBL height [m]}', rotation=-90, va='bottom')
    h2, l2 = axs[2].get_legend_handles_labels()
    h4, l4 = ax4.get_legend_handles_labels()
    #  h10, l10 = ax10.get_legend_handles_labels()
    #  axs[2].legend(h2+h4+h10, l2+l4+l10, loc='best')

    #------------------------------------------------------------  
    # Now add full simulated abundance time series to top plot 
    #------------------------------------------------------------  
    max_y = []  #track max abund for ylimits
    for idx, row in em.iterrows():
        targ_ls = row['Ls(deg)']
        targ_dec_time = row['DECIMAL_TIME']
        #  fuzz = 0.01
        a = df_all.copy()[ df_all['L_s'].between(targ_ls-fuzz, targ_ls+fuzz) ]
        if len(a)>1: a = a.iloc[0]
        # Get the SOL_TIME of that day and get that whole day
        targ_sol = a['SOL_TIME']
        b = df_all.copy()[ df_all['SOL_TIME'].between(int(targ_sol), int(targ_sol)+1.0) ]
        #make midnight next day equal to hour 24
        if b['SOL_TIME'].iloc[-1] == int(targ_sol)+1.0:
            b['DECIMAL_TIME'].iloc[-1] = 24.00
        # ------- LINES ----------
        # Plot whole simulated abundance for the affirmative detections (solid line)
        if em['GlobalMeanCH4(ppbv)'].loc[idx] >= 0.07:
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='-', color=cmap.to_rgba(targ_ls))
            #  max_y.append( max(b['C_ppbv']) )
        # Also plot the whole simulated abundance prediction for the non-detects  (dashed line)
        if em['GlobalMeanCH4(ppbv)'].loc[idx] < 0.07:
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='--', color=cmap.to_rgba(targ_ls))
        max_y.append( max(b['C_ppbv']) )
    #-- Properties
    axs[0].set_ylim(0., max(max_y)*1.04 )
    plt.savefig(join(outputdir,'all_abundance_dayNight_differences.pdf'))
    plt.close('all')
    print('    Done.')


    #================================================== 
    # Take a closer look at day-night differences for N. Summer (Ls = 90-180°)
    #================================================== 
    #-------------------------------------------------- 
    # Read in Webster2021 abundances  
    #-------------------------------------------------- 
    em = readAbundances(run_type='Enrichment', source='Webster2021')
    em_summer = em.copy()[ em['Ls(deg)'].between(90,180) ]
    #  # Calculate the decimal times
    #  times = em_summer.copy()['LMST']
    #  dec_time = np.asarray([float(t.split(':')[0])+float(t.split(':')[1])/60 for t in times])  #[hrs]
    #  # Insert into dataframe
    #  pos = em_summer.columns.get_loc('LMST')+1
    #  em_summer.insert(pos, 'DECIMAL_TIME', dec_time)

    #-------------------------------------------------- 
    # Read Simulated Abundances 
    #-------------------------------------------------- 
    df_summer = df_copy.copy()[ df_copy['L_s'].between(90,180) ]
    #  # Calculate the decimal times
    #  dec_time = df_summer.copy()['SOL_TIME']%1*24  #[hrs]
    #  # Insert into dataframe
    #  pos = df_summer.columns.get_loc('SOL_TIME')+1
    #  df_summer.insert(pos, 'DECIMAL_TIME', dec_time)

    #-------------------------------------------------- 
    # Plot day-night Measured vs Simulated Atm Abundances for N. Summer 
    #-------------------------------------------------- 
    print('Plotting day-night differences in atm abundance for N. Summer...')
    #  gridspec = {'width_ratios': [1, 1, 0.1, 0.1]}  #colorbar axes shoudld be small
    gridspec = {'width_ratios': [1, 0.1]}  #colorbar axes should be small
    # [Subplot Map]
    #--------------
    # [0] [1]
    # [2] [3] 
    #--------------
    # 0=top panel,    1=colorbar, 
    # 2=bottom panel, 3=blank, 4=second axis on [2]
    pblcount = 0  #counter for plotting only 1 PBL height (for non-detects)
    fig, ax = plt.subplots(2,2, sharex=False,figsize=(7,10), gridspec_kw=gridspec)
    axs = ax.ravel()
    #  ax4 = ax2.twinx()
    ax4 = axs[2].twinx()
    # Assign colors increasing with L_s(°)
    #  cmapName='viridis'  #sequential
    #  cmapName='cool' #sequential
    cmapName='plasma' #uniform
    #  cmapName='twilight' #cyclic
    #  cmapName='twilight_shifted' #cyclic
    cmap=matplotlib.cm.get_cmap(cmapName)
    #  colors = iter(cmap(np.linspace(0,1,len(em_summer))))

    # Custom colorbar delineating L_s of each point
    cbar_ticks = 10
    #  c = np.arange(90, 180+cbar_ticks, cbar_ticks)
    c = np.arange(100, 160+cbar_ticks, cbar_ticks)
    norm = matplotlib.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cmap.set_array([])

    # Measured Abundances
    #Save the Abundance values interpolated to the observation times 
    obs_abund_dict = {}
    sim_abund_dict = {}
    for idx, row in em_summer.iterrows():
        #  c = next(colors)
        targ_ls = row['Ls(deg)']
        targ_dec_time = row['DECIMAL_TIME']
        #  print('Ls = {}'.format(targ_ls))
        #  print('targ time = {}'.format(targ_dec_time))

        # [0] ABUNDANCES -----------------------------------------------
        #---- Measured Abundances
        x = em_summer['DECIMAL_TIME'].loc[idx]
        y = em_summer['GlobalMeanCH4(ppbv)'].loc[idx]
        yerr = em_summer['Error(±1SEM).1'].loc[idx]
        #  ax1.errorbar(x,y, yerr=yerr,marker='o',ms=10,capsize=5,ls='',color='red',mec='k')
        axs[0].errorbar(x,y, yerr=yerr,marker='o',ms=10,capsize=5,ls='',color=cmap.to_rgba(targ_ls),mec='k')

        #---- Simulated Abundances
        # Find closest simulated L_s to target obs L_s
        #  a = find_nearest(df_summer['L_s'], targ_ls)
        fuzz = 0.01
        a = df_summer.copy()[ df_summer['L_s'].between(targ_ls-fuzz, targ_ls+fuzz) ]
        if len(a)>1: a = a.iloc[0]
        # Get the SOL_TIME of that day and get that whole day
        targ_sol = a['SOL_TIME']
        b = df_summer.copy()[ df_summer['SOL_TIME'].between(int(targ_sol), int(targ_sol)+1.0) ]
        #make midnight next day equal to hour 24
        if b['SOL_TIME'].iloc[-1] == int(targ_sol)+1.0:
            b['DECIMAL_TIME'].iloc[-1] = 24.00
        #  # Interpolate Abundance to time of observation
        # Numpy Interpolate seems to work better for this 
        x = targ_dec_time
        #  x = targ_dec_time + 1.0  #TEST    TEST     TEST     TEST     TEST     TEST     TEST      TEST
        y = np.interp(x, b['DECIMAL_TIME'], b['C_ppbv'])
        axs[0].errorbar(x,y, yerr=0., marker='*',ms=18,capsize=0,ls='',color=cmap.to_rgba(targ_ls),mec='k')
        # [2] FLUXES AND PBL HEIGHT ------------------------------------

        # (for now, only do for the non-detects...)
        if em_summer['GlobalMeanCH4(ppbv)'].loc[idx] < 0.07:
            #---- FLUXES
            # Dashed line for non-detects
            axs[2].plot(b['DECIMAL_TIME'], b['FLUX_kg/m2/s'], ls='--', color=cmap.to_rgba(targ_ls))
            #  axs[2].plot(b['DECIMAL_TIME'], b['FLUX_kg/m2/s'], ls='-', color=cmap.to_rgba(targ_ls), label=r'flux (L$_s$ = {}°)'.format(targ_ls))
            #---- PBL HEIGHT (only need one, should be same)
            if pblcount==0:
                ax4.plot(b['DECIMAL_TIME'], b['PBL_m'], ls=':', color='k', label='PBL')
                #  ax4.plot(b['DECIMAL_TIME'], b['PBL_m'], ls=':', color=cmap.to_rgba(targ_ls), label='PBL')
            #---- PRESSURE (don't need a scale, only plot one)
                #  ax10 = ax4.twinx()
                ax10 = axs[2].twinx()
                ax10.plot(b['DECIMAL_TIME'], b['PRESS'], ls='-', color='blue', label='pressure')
                # Note: Don't plot the yticks or labels for pressure (can mention in caption)
                plt.setp(ax10.get_yticklabels(), visible=False)
                ax10.tick_params(axis='y', which='both', length=0)
                pblcount+=1
        else:
            #---- FLUXES
            # Solid line for positive detects
            axs[2].plot(b['DECIMAL_TIME'], b['FLUX_kg/m2/s'], ls='-', color=cmap.to_rgba(targ_ls))
        # Save values to dictionary so they're easy to access if needed
        obs_abund_dict[targ_ls] = {'DECIMAL_TIME': targ_dec_time,
                                   'SOL'         : int(targ_sol),
                                   'C_ppbv'      : em_summer['GlobalMeanCH4(ppbv)'].loc[idx],
                                   #  'PBL_m'       : np.interp(x,b['DECIMAL_TIME'], b['PBL_m']),
                                   'Error(±1SEM).1'  : em_summer['Error(±1SEM).1'].loc[idx],
                                  }
        sim_abund_dict[targ_ls] = {'DECIMAL_TIME': targ_dec_time,
                                   'SOL'         : int(targ_sol),
                                   'C_ppbv'      : y,
                                   'PBL_m'       : np.interp(x,b['DECIMAL_TIME'], b['PBL_m']),
                                  }
        #  print()
        #  fig.colorbar(cmap, ticks=c)
    #  # --------------------------------------
    #  # Custom colorbar
    #  #  divider = make_axes_locatable(plt.gca())
    #  #  divider = make_axes_locatable(ax1)
    #  ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    #  #  #  cb1 = matplotlib.colorbar.ColorbarBase(ax_cb, cmap=cmap, orientation='vertical')
    #  # Set the colormap and norm to correspond to the data for which
    #  # the colorbar will be used.
    #  norm = matplotlib.colors.Normalize(vmin=90, vmax=180)
    #  #  cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
    #  cb1 = matplotlib.colorbar.ColorbarBase(ax_cb, cmap=cmap, norm=norm, orientation='vertical')
    #  cb1.set_label(r'L$_s$ [°]')
    #  # --------------------------------------
    # Calculate the error b/w sim and obs
    obs = np.asarray([obs_abund_dict[c]['C_ppbv'] for c in obs_abund_dict])
    sim = np.asarray([sim_abund_dict[c]['C_ppbv'] for c in sim_abund_dict])
    #-- Calculate Errors
    #  dof = 4 # no. of degrees of freedom
    rmse, chisq, redchisq = calcAbundError(df_summer, em_summer, dof=dof)
    #  #-- Calculate RMSE
    #  rmse = calc_rmse(obs, sim)
    #  #-- Calculate Chi-Square Statistic 
    #  obs_sem = np.asarray([obs_abund_dict[c]['Error(±1SEM).1'] for c in obs_abund_dict])
    #  #  chisq = chisqg(ydata=ya, ymod=ybi, sd=None)
    #  chisq = chisqg(ydata=obs, ymod=sim, sd=obs_sem)
    #  #-- Calculate Reduced Chi-Square Statistic 
    #  dof = np.shape(bnds)[0]#2 # no. of degrees of freedom
    #  #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=None)
    #  redchisq = redchisqg(ydata=obs, ymod=sim, deg=dof, sd=obs_sem)
    # Write error to file
    with open(join(outputdir,'summerDayNight_error'), 'w') as f:
        f.write('RMSE              = {:.4f}\n'.format(rmse))
        f.write('Chi-square        = {:.4f}\n'.format(chisq))
        f.write('Reduced Chi-square = {:.4f}\n'.format(redchisq))
    # Create text fields if desired
    rmse_text     = 'RMSE = {:.2f}'.format(rmse)
    chisq_text    = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(chisq)
    redchisq_text = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)
    #  redchisq_text = r'Reduced $\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)

    #  # CALCULATE ERROR, OMITTING Ls=103.48° [idx=0]
    #  #-- Calculate RMSE
    #  rmse = calc_rmse(obs[1:], sim[1:])
    #  #-- Calculate Chi-Square Statistic 
    #  obs_sem = np.asarray([obs_abund_dict[c]['Error(±1SEM).1'] for c in obs_abund_dict])[1:]
    #  #  chisq = chisqg(ydata=ya, ymod=ybi, sd=None)
    #  chisq = chisqg(ydata=obs[1:], ymod=sim[1:], sd=obs_sem)
    #  #-- Calculate Reduced Chi-Square Statistic 
    #  dof = np.shape(bnds)[0]#2 # no. of degrees of freedom
    #  #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=None)
    #  redchisq = redchisqg(ydata=obs[1:], ymod=sim[1:], deg=dof, sd=obs_sem)
    #  # Write error to file
    #  with open(join(outputdir,'summerDayNight_error_omitLs_103'), 'w') as f:
        #  f.write('RMSE              = {:.4f}\n'.format(rmse))
        #  f.write('Chi-square        = {:.4f}\n'.format(chisq))
        #  f.write('Reduced Chi-square = {:.4f}\n'.format(redchisq))
#  
    # [1] Colorbar
    cax = axs[1]
    cb1 = plt.colorbar(cmap, ticks=c, cax=cax)
    cb1.set_label(r'\textbf{Solar Longitude, L$_s$ [°]}', rotation=-90, va='bottom')
    # Make [3] invisible
    axs[3].set_visible(False)
    #  #-- Properties
    # [0]
    axs[0].set_title('N. Summer Day-Night Variations')
    vshift=0.05; y0 = 0.95
    #  axs[0].text(0.97, y0,rmse_text, transform=axs[0].transAxes, ha='right')
    axs[0].text(0.97, y0,redchisq_text, transform=axs[0].transAxes, ha='right')
    #  axs[0].text(0.97, y0-vshift,redchisq_text, transform=axs[0].transAxes, ha='right')
    # Legend
    # Only need legend for one marker each
    axs[0].lines[0].set_label('observed')
    axs[0].lines[3].set_label('simulated')
    #  for i in range(10):
        #  axs[0].lines[i].set_label(i)
    axs[0].legend(loc='upper center')
    leg = axs[0].get_legend()
    # Force legend colors to be black (doesn't seem to work...)
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('black')
    axs[0].set_ylabel(r'\textbf{CH$_4$ Abundance [ppbv]}')
    axs[0].set_xlabel(r'\textbf{LMST}')
    #  axs[0].set_ylim(0., max(max_y) )
    #  axs[0].set_ylim(bottom=0.)
    #  ax1.legend()
    pad=0.5
    axs[0].set_xlim(0-pad, 24+pad)
    axs[0].xaxis.set_major_locator(MultipleLocator(2)) #2-hour major ticks 
    axs[0].xaxis.set_minor_locator(MultipleLocator(1)) #1-hour minor ticks 
    # [2]
    axs[2].set_ylim(bottom=0)
    axs[2].set_xlabel(r'\textbf{LMST}')
    axs[2].set_ylabel(r'\textbf{Surface Flux [kg/m$^2$/s]}')
    axs[2].set_xlim(0-pad, 24+pad)
    axs[2].xaxis.set_major_locator(MultipleLocator(2)) #2-hour major ticks 
    axs[2].xaxis.set_minor_locator(MultipleLocator(1)) #1-hour minor ticks 
    #  axs[2].set_yscale('semilog')  #LOG SCALE FLUXES (?)
    # [ax4]
    ax4.set_ylim(bottom=0)
    ax4.set_ylabel(r'\textbf{PBL height [m]}', rotation=-90, va='bottom')
    h2, l2 = axs[2].get_legend_handles_labels()
    h4, l4 = ax4.get_legend_handles_labels()
    h10, l10 = ax10.get_legend_handles_labels()
    axs[2].legend(h2+h4+h10, l2+l4+l10, loc='best')
    #  axs[2].legend(h2+h4, l2+l4, loc='best')
    #  #  plt.tight_layout()
    #  plt.subplots_adjust(wspace=0.05)
    #  plt.savefig(join(outputdir,'summer_abundance_dayNight_differences_noLines.pdf'))
    #  #  plt.close('all')

    #------------------------------------------------------------  
    # Now add full simulated abundance time series to top plot 
    #------------------------------------------------------------  
    max_y = []  #track max abund for ylimits
    for idx, row in em_summer.iterrows():
        targ_ls = row['Ls(deg)']
        targ_dec_time = row['DECIMAL_TIME']
        fuzz = 0.01
        a = df_summer.copy()[ df_summer['L_s'].between(targ_ls-fuzz, targ_ls+fuzz) ]
        if len(a)>1: a = a.iloc[0]
        # Get the SOL_TIME of that day and get that whole day
        targ_sol = a['SOL_TIME']
        b = df_summer.copy()[ df_summer['SOL_TIME'].between(int(targ_sol), int(targ_sol)+1.0) ]
        #make midnight next day equal to hour 24
        if b['SOL_TIME'].iloc[-1] == int(targ_sol)+1.0:
            b['DECIMAL_TIME'].iloc[-1] = 24.00
        # ------- LINES ----------
        # Plot whole simulated abundance for the affirmative detections (solid line)
        if em_summer['GlobalMeanCH4(ppbv)'].loc[idx] >= 0.07:
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='-', color=cmap.to_rgba(targ_ls))
            #  max_y.append( max(b['C_ppbv']) )
        # Also plot the whole simulated abundance prediction for the non-detects  (dashed line)
        if em_summer['GlobalMeanCH4(ppbv)'].loc[idx] < 0.07:
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='--', color=cmap.to_rgba(targ_ls))
        my = max( max(df_copy['C_ppbv']), max(em_summer['GlobalMeanCH4(ppbv)']+em_summer['Error(±1SEM).1']) )
        max_y.append( my )
        #  max_y.append( max(b['C_ppbv']) )
    #-- Properties
    axs[0].set_ylim(0., max(max_y)*1.04 )
    plt.savefig(join(outputdir,'summer_abundance_dayNight_differences.pdf'))
    plt.close('all')
    print('    Done.')

    #-------------------------------------------------- 
    # ADD TO summaryFile 
    #-------------------------------------------------- 
    with open(summaryFile, 'a') as f:
        f.write('Day-Night Error (N. Summer):\n')
        f.write('--------------------------------------------------\n')
        f.write('Reduced Chi-Square : {:.4f}\n'.format(redchisq))
        f.write('\n')

    #================================================== 
    # Take a closer look at day-night differences for N. Spring (Ls = 41.5°)
    #================================================== 
    #-------------------------------------------------- 
    # Read in Webster2021 abundances  
    #-------------------------------------------------- 
    em = readAbundances(run_type='Enrichment', source='Webster2021')
    em_spring = em.copy()[ em['Ls(deg)'].between(0,90) ]
    # (TRY THIS) Leave out 20ppbv point
    em_spring = em_spring.copy()[ em_spring['GlobalMeanCH4(ppbv)']<20. ]
    #  # Calculate the decimal times
    #  times = em_spring.copy()['LMST']
    #  dec_time = np.asarray([float(t.split(':')[0])+float(t.split(':')[1])/60 for t in times])  #[hrs]
    #  # Insert into dataframe
    #  pos = em_spring.columns.get_loc('LMST')+1
    #  em_spring.insert(pos, 'DECIMAL_TIME', dec_time)
        #  # (for now, only do for the non-detects...)

    #-------------------------------------------------- 
    # Read Simulated Abundances 
    #-------------------------------------------------- 
    df_spring = df_copy.copy()[ df_copy['L_s'].between(0,90) ]
    #  # Calculate the decimal times
    #  dec_time = df_spring.copy()['SOL_TIME']%1*24  #[hrs]
    #  # Insert into dataframe
    #  pos = df_spring.columns.get_loc('SOL_TIME')+1
    #  df_spring.insert(pos, 'DECIMAL_TIME', dec_time)

    #-------------------------------------------------- 
    # Plot day-night Measured vs Simulated Atm Abundances for N. Spring 
    #-------------------------------------------------- 
    print('Plotting day-night differences in atm abundance for N. Spring...')
    #  gridspec = {'width_ratios': [1, 1, 0.1, 0.1]}  #colorbar axes shoudld be small
    gridspec = {'width_ratios': [1, 0.1]}  #colorbar axes should be small
    # [Subplot Map]
    #--------------
    # [0] [1]
    # [2] [3] 
    #--------------
    # 0=top panel,    1=colorbar, 
    # 2=bottom panel, 3=blank, 4=second axis on [2]
    pblcount = 0  #counter for plotting only 1 PBL height (for non-detects)
    fig, ax = plt.subplots(2,2, sharex=False,figsize=(7,10), gridspec_kw=gridspec)
    axs = ax.ravel()
    #  ax4 = ax2.twinx()
    ax4 = axs[2].twinx()
    # Assign colors increasing with L_s(°)
    cmapName='viridis'  #sequential
    #  cmapName='cool' #sequential
    #  cmapName='plasma' #uniform
    #  cmapName='twilight' #cyclic
    #  cmapName='twilight_shifted' #cyclic
    cmap=matplotlib.cm.get_cmap(cmapName)
    #  colors = iter(cmap(np.linspace(0,1,len(em_spring))))

    # Custom colorbar delineating L_s of each point
    cbar_ticks = 10
    #  c = np.arange(100, 160+cbar_ticks, cbar_ticks)
    c = np.arange(30, 70+cbar_ticks, cbar_ticks)
    norm = matplotlib.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cmap.set_array([])

    # Measured Abundances
    #Save the Abundance values interpolated to the observation times 
    obs_abund_dict = {}
    sim_abund_dict = {}
    for idx, row in em_spring.iterrows():
        #  c = next(colors)
        targ_ls = row['Ls(deg)']
        targ_dec_time = row['DECIMAL_TIME']
        #  print('Ls = {}'.format(targ_ls))
        #  print('targ time = {}'.format(targ_dec_time))

        # [0] ABUNDANCES -----------------------------------------------
        #---- Measured Abundances
        x = em_spring['DECIMAL_TIME'].loc[idx]
        y = em_spring['GlobalMeanCH4(ppbv)'].loc[idx]
        yerr = em_spring['Error(±1SEM).1'].loc[idx]
        #  ax1.errorbar(x,y, yerr=yerr,marker='o',ms=10,capsize=5,ls='',color='red',mec='k')
        axs[0].errorbar(x,y, yerr=yerr,marker='o',ms=10,capsize=5,ls='',color=cmap.to_rgba(targ_ls),mec='k')

        #---- Simulated Abundances
        # Find closest simulated L_s to target obs L_s
        #  a = find_nearest(df_spring['L_s'], targ_ls)
        fuzz = 0.01
        a = df_spring.copy()[ df_spring['L_s'].between(targ_ls-fuzz, targ_ls+fuzz) ]
        if len(a)>1: a = a.iloc[0]
        # Get the SOL_TIME of that day and get that whole day
        targ_sol = a['SOL_TIME']
        b = df_spring.copy()[ df_spring['SOL_TIME'].between(int(targ_sol), int(targ_sol)+1.0) ]
        #make midnight next day equal to hour 24
        if b['SOL_TIME'].iloc[-1] == int(targ_sol)+1.0:
            b['DECIMAL_TIME'].iloc[-1] = 24.00
        #  # Interpolate Abundance to time of observation
        # Numpy Interpolate seems to work better for this 
        x = targ_dec_time
        y = np.interp(x, b['DECIMAL_TIME'], b['C_ppbv'])
        #  print('y = {}'.format(y))
        axs[0].errorbar(x,y, yerr=0., marker='*',ms=18,capsize=0,ls='',color=cmap.to_rgba(targ_ls),mec='k')

        # Plot whole simulated abundance for the affirmative detections (solid line)
        if em_spring['GlobalMeanCH4(ppbv)'].loc[idx] >= 0.07:
            #  ax1.plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='-', color=cmap.to_rgba(targ_ls))
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='-', color=cmap.to_rgba(targ_ls))
        # Also plot the whole simulated abundance prediction for the non-detects  (dashed line)
        if em_spring['GlobalMeanCH4(ppbv)'].loc[idx] < 0.07:
            #  ax1.plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='--', color=cmap.to_rgba(targ_ls))
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='--', color=cmap.to_rgba(targ_ls))

        # [2] FLUXES AND PBL HEIGHT ------------------------------------
        # (for now, only do for the big 20ppbv spike at Ls=41.5)...
        #  if em_spring['GlobalMeanCH4(ppbv)'].loc[idx] > 20.0:
        if em_spring['GlobalMeanCH4(ppbv)'].loc[idx] > 0.0:
            #---- FLUXES
            #  axs[2].plot(b['DECIMAL_TIME'], b['FLUX_ADJ'], ls='-', color=cmap.to_rgba(targ_ls), label=r'flux (L$_s$ = {}°)'.format(targ_ls))
            #  axs[2].plot(b['DECIMAL_TIME'], b['FLUX_kg/m2/s_ADJ'], ls='-', color=cmap.to_rgba(targ_ls))
            axs[2].plot(b['DECIMAL_TIME'], b['FLUX_kg/m2/s'], ls='-', color=cmap.to_rgba(targ_ls))
            #  axs[2].plot(b['DECIMAL_TIME'], b['FLUX_ADJ'], ls='-', color=cmap.to_rgba(targ_ls))
            #  axs[2].semilogy(b['DECIMAL_TIME'], b['FLUX_ADJ'], ls='-', color=cmap.to_rgba(targ_ls), label=r'flux (L$_s$ = {}°)'.format(targ_ls))         #LOG SCALE FLUXES (?)
            #---- PBL HEIGHT (only need one, should be same)
            if pblcount==0:
                ax4.plot(b['DECIMAL_TIME'], b['PBL_m'], ls=':', color='k', label='PBL')
                #  ax4.plot(b['DECIMAL_TIME'], b['PBL_m'], ls=':', color=cmap.to_rgba(targ_ls), label='PBL')
            #---- PRESSURE (don't need a scale, only plot one)
                #  ax10 = ax4.twinx()
                ax10 = axs[2].twinx()
                ax10.plot(b['DECIMAL_TIME'], b['PRESS'], ls='-', color='blue', label='pressure')
                # Note: Don't plot the yticks or labels for pressure (can mention in caption)
                plt.setp(ax10.get_yticklabels(), visible=False)
                ax10.tick_params(axis='y', which='both', length=0)
                pblcount+=1
            #  #TEST


        # Save values to dictionary so they're easy to access if needed
        obs_abund_dict[targ_ls] = {'DECIMAL_TIME': targ_dec_time,
                                   'SOL'         : int(targ_sol),
                                   'C_ppbv'      : em_spring['GlobalMeanCH4(ppbv)'].loc[idx],
                                   #  'PBL_m'       : np.interp(x,b['DECIMAL_TIME'], b['PBL_m']),
                                   'Error(±1SEM).1'  : em_spring['Error(±1SEM).1'].loc[idx],
                                  }
        sim_abund_dict[targ_ls] = {'DECIMAL_TIME': targ_dec_time,
                                   'SOL'         : int(targ_sol),
                                   'C_ppbv'      : y,
                                   'PBL_m'       : np.interp(x,b['DECIMAL_TIME'], b['PBL_m']),
                                  }

        #  fig.colorbar(cmap, ticks=c)
    #  # --------------------------------------
    #  # Custom colorbar
    #  #  divider = make_axes_locatable(plt.gca())
    #  #  divider = make_axes_locatable(ax1)
    #  ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    #  #  #  cb1 = matplotlib.colorbar.ColorbarBase(ax_cb, cmap=cmap, orientation='vertical')
    #  # Set the colormap and norm to correspond to the data for which
    #  # the colorbar will be used.
    #  norm = matplotlib.colors.Normalize(vmin=90, vmax=180)
    #  #  cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
    #  cb1 = matplotlib.colorbar.ColorbarBase(ax_cb, cmap=cmap, norm=norm, orientation='vertical')
    #  cb1.set_label(r'L$_s$ [°]')
    #  # --------------------------------------
    # Calculate the error b/w sim and obs
    obs = np.asarray([obs_abund_dict[c]['C_ppbv'] for c in obs_abund_dict])
    sim = np.asarray([sim_abund_dict[c]['C_ppbv'] for c in sim_abund_dict])
    #-- Calculate RMSE
    rmse = calc_rmse(obs, sim)
    #-- Calculate Chi-Square Statistic 
    obs_sem = np.asarray([obs_abund_dict[c]['Error(±1SEM).1'] for c in obs_abund_dict])
    #  chisq = chisqg(ydata=ya, ymod=ybi, sd=None)
    chisq = chisqg(ydata=obs, ymod=sim, sd=obs_sem)
    #-- Calculate Reduced Chi-Square Statistic 
    dof = np.shape(bnds)[0]#2 # no. of degrees of freedom
    #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=None)
    redchisq = redchisqg(ydata=obs, ymod=sim, deg=dof, sd=obs_sem)
    # Write error to file
    with open(join(outputdir,'springDayNight_error'), 'w') as f:
        f.write('RMSE              = {:.4f}\n'.format(rmse))
        f.write('Chi-square        = {:.4f}\n'.format(chisq))
        f.write('Reduced Chi-square = {:.4f}\n'.format(redchisq))
    # Create text fields if desired
    rmse_text     = 'RMSE = {:.2f}'.format(rmse)
    chisq_text    = r'$\chi_{{\nu}}^2$ = {:.2f}'.format(chisq)
    redchisq_text = r'Reduced $\chi_{{\nu}}^2$ = {:.2f}'.format(redchisq)

    # CALCULATE ERROR, OMITTING Ls=41.5° [idx=3]
    idx=3  #omit index
    #-- Calculate RMSE
    #  rmse = calc_rmse(obs[1:], sim[1:])
    rmse = calc_rmse(np.delete(obs,idx,0), np.delete(sim,idx,0))
    #-- Calculate Chi-Square Statistic 
    #  obs_sem = np.asarray([obs_abund_dict[c]['Error(±1SEM).1'] for c in obs_abund_dict])[1:]
    obs_sem = np.asarray([obs_abund_dict[c]['Error(±1SEM).1'] for c in obs_abund_dict])
    obs_sem = np.delete(obs_sem, idx, 0)
    #  chisq = chisqg(ydata=obs[1:], ymod=sim[1:], sd=obs_sem)
    chisq = chisqg(ydata=np.delete(obs,idx,0), ymod=np.delete(sim,idx,0), sd=obs_sem)
    #-- Calculate Reduced Chi-Square Statistic 
    dof = np.shape(bnds)[0]#2 # no. of degrees of freedom
    #  redchisq = redchisqg(ydata=ya, ymod=ybi, deg=dof, sd=None)
    redchisq_omit = redchisqg(ydata=np.delete(obs,idx,0), ymod=np.delete(sim,idx,0), deg=dof,sd=obs_sem)
    # Write error to file
    with open(join(outputdir,'springDayNight_error_omitLs_41'), 'w') as f:
        f.write('RMSE              = {:.4f}\n'.format(rmse))
        f.write('Chi-square        = {:.4f}\n'.format(chisq))
        f.write('Reduced Chi-square = {:.4f}\n'.format(redchisq_omit))
        redchisq_text_omit = r'Reduced $\chi_{{\nu}}^2$ = {:.2f} ({:.2f})'.format(redchisq, redchisq_omit)


    # [1] Colorbar
    cax = axs[1]
    cb1 = plt.colorbar(cmap, ticks=c, cax=cax)
    cb1.set_label(r'\textbf{Solar Longitude, L$_s$ [°]}', rotation=-90, va='bottom')
    # Make [3] invisible
    axs[3].set_visible(False)

    #  #-- Properties
    # [0]
    #  axs[0].set_yscale('log')
    axs[0].set_title('N. Spring Day-Night Variations')
    vshift=0.05; y0 = 0.95
    #  axs[0].text(0.97, y0,rmse_text, transform=axs[0].transAxes, ha='right')
    axs[0].text(0.97, y0,redchisq_text_omit, transform=axs[0].transAxes, ha='right')
    #  axs[0].text(0.97, y0-vshift,redchisq_text_omit, transform=axs[0].transAxes, ha='right')
    #  axs[0].text(0.97, y0-vshift,redchisq_text, transform=axs[0].transAxes, ha='right')
    # Legend
    # Only need legend for one marker each
    axs[0].lines[0].set_label('observed')
    axs[0].lines[3].set_label('simulated')
    #  for i in range(10):
        #  axs[0].lines[i].set_label(i)
    #  axs[0].legend(loc='upper center')
    axs[0].legend(loc='center')
    leg = axs[0].get_legend()
    # Force legend colors to be black (doesn't seem to work...)
    leg.legendHandles[0].set_color('black')
    leg.legendHandles[1].set_color('black')
    axs[0].set_ylabel(r'\textbf{CH$_4$ Abundance [ppbv]}')
    axs[0].set_xlabel(r'\textbf{LMST}')
    axs[0].set_ylim(bottom=0.)
    #  ax1.legend()
    pad=0.5
    axs[0].set_xlim(0-pad, 24+pad)
    axs[0].xaxis.set_major_locator(MultipleLocator(2)) #2-hour major ticks 
    axs[0].xaxis.set_minor_locator(MultipleLocator(1)) #1-hour minor ticks 
    # [2]
    axs[2].set_ylim(bottom=0)
    axs[2].set_xlabel(r'\textbf{LMST}')
    axs[2].set_ylabel(r'\textbf{Surface Flux [kg/km$^2$/s]}')
    axs[2].set_xlim(0-pad, 24+pad)
    axs[2].xaxis.set_major_locator(MultipleLocator(2)) #2-hour major ticks 
    axs[2].xaxis.set_minor_locator(MultipleLocator(1)) #1-hour minor ticks 
    #  axs[2].set_yscale('semilog')  #LOG SCALE FLUXES (?)
    # [ax4]
    ax4.set_ylim(bottom=0)
    ax4.set_ylabel(r'\textbf{PBL height [m]}', rotation=-90, va='bottom')
    h2, l2 = axs[2].get_legend_handles_labels()
    h4, l4 = ax4.get_legend_handles_labels()
    h10, l10 = ax10.get_legend_handles_labels()
    axs[2].legend(h2+h4+h10, l2+l4+l10, loc='best')
    #  axs[2].legend(h2+h4, l2+l4, loc='best')

    #  #  plt.tight_layout()
    #  plt.subplots_adjust(wspace=0.05)
    #  plt.savefig(join(outputdir,'spring_abundance_dayNight_differences_noLines.pdf'))
    #  #  plt.close('all')

    #------------------------------------------------------------  
    # Now add full simulated abundance time series to top plot 
    #------------------------------------------------------------  
    for idx, row in em_spring.iterrows():
        targ_ls = row['Ls(deg)']
        targ_dec_time = row['DECIMAL_TIME']
        fuzz = 0.01
        a = df_spring.copy()[ df_spring['L_s'].between(targ_ls-fuzz, targ_ls+fuzz) ]
        if len(a)>1: a = a.iloc[0]
        # Get the SOL_TIME of that day and get that whole day
        targ_sol = a['SOL_TIME']
        b = df_spring.copy()[ df_spring['SOL_TIME'].between(int(targ_sol), int(targ_sol)+1.0) ]
        #make midnight next day equal to hour 24
        if b['SOL_TIME'].iloc[-1] == int(targ_sol)+1.0:
            b['DECIMAL_TIME'].iloc[-1] = 24.00
        # ------- LINES ----------
        # Plot whole simulated abundance for the affirmative detections (solid line)
        if em_spring['GlobalMeanCH4(ppbv)'].loc[idx] >= 0.07:
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='-', color=cmap.to_rgba(targ_ls))
        # Also plot the whole simulated abundance prediction for the non-detects  (dashed line)
        if em_spring['GlobalMeanCH4(ppbv)'].loc[idx] < 0.07:
            axs[0].plot(b['DECIMAL_TIME'], b['C_ppbv'], ls='--', color=cmap.to_rgba(targ_ls))
    plt.savefig(join(outputdir,'spring_abundance_dayNight_differences.pdf'))
    plt.close('all')
    print('    Done.')




    #--------------------------------------------------------------------------------  
    # CONTOUR PLOT OF CONC PROFILES IN TIME 
    #--------------------------------------------------------------------------------  
    # 
    # Make this a time slice at approx the time (Ls) that will be sampled...  !!! !!! !!!
    #
    Ls_target = 130. #[°] approximate L_s to get a couple sols of contours for 
    # Get the index of the corresponding L_s in data (last year)
    ctarg = find_nearest( np.asarray(df_copy['L_s']), Ls_target )
    ltarg = df_copy.copy()[ df_copy['L_s']==ctarg ]
    li = ltarg.index[0]
    si = int(ltarg['SOL_INT'])
    # Get dataFrame of all times on that sol and the following (2 sols total)
    cont_df = df_copy.copy()[ df_copy['SOL_INT'].between(si,si+1) ]
    # Using those indices, grab correct contours in time
    targ_idx = np.asarray(cont_df.index)


    #-----------------------------------------------------
    # Generate meshgrid
    #-----------------------------------------------------
    # X: time array
    # Y: vertical height
    # Z: concentration value

    x = cont_df['SOL_TIME']
    y = np.linspace(0.,max(c_df['PBL_m']),np.shape(c_all)[-1]) #vertical coordinate 
    #----  Make all negative ones a smaller negative 
    neg_value = -1e-20
    c_all_copy = c_all[targ_idx]
    c_all_copy[c_all_copy==-1.] = neg_value #set to  v. small neg. value so contours aren't messed up
    #---
    z = c_all_copy
    # Convert kg/m3 to ppbv
    z_ppbv = np.zeros(np.shape(z))
    for i,profile in enumerate(z):
        profile_ppbv = conc_gm3_to_ppbv(profile*1000) #first convert kg/m3 to g/m3
        z_ppbv[i] = profile_ppbv
    z_ppbv[ z_ppbv<0. ] = neg_value
    #----
    X,Y = np.meshgrid(x,y)
    Z   = z_ppbv.T

    #-----------------------------------------------------
    # Create Contour Plot
    #-----------------------------------------------------
    fig, ax1 = plt.subplots(1, figsize=(10,5))
    #-- Contour Properties
    cmap_name = matplotlib.cm.inferno
    cmap_name.set_under(color='grey')
    vmin = 0.; #vmax=0.005
    nlevels=50#100#50
    #  levels = np.linspace(vmin, vmax, nlevels+1)
    #  cmap_name = 'jet'
    #  cmap_name = 'inferno'
    #--
    cs = ax1.contourf(X,Y,Z, nlevels, cmap=cmap_name, extend='both',vmin=vmin)
    # Try to remove weird lines between contours
    # < https://discourse.matplotlib.org/t/unwanted-lines-between-contourf-contour-levels/12539/2 >
    #---- make sure the ax1.contourf call is the same as above (can't use "cs"...)
    for c in ax1.contourf(X,Y,Z,nlevels, cmap=cmap_name, extend='both',vmin=vmin).collections:
        #  c.set_linewidth(0.15)
        c.set_linewidth(0.1)
    #-- Colorbar
    cbar = fig.colorbar(cs)
    cbar.set_label(r'\textbf{CH$_4$ Abundance [ppbv]}', rotation=-90, va='bottom')
    #  cbar.ax.set_ylabel(r'CH$_4$ Abundance [ppbv]', rotation=-90, va='bottom')
    cs.set_clim(cs.get_clim()[0], cs.get_clim()[1])
    #  cs.set_clim(vmin,vmax)
    #  cs.cmap.set_under('grey')
    #--
    #-- Properties
    ax1.set_title(r'$L_s =$ {:.0f}$^{{\circ}}$'.format(Ls_target))
    ax1.set_ylim(0, max(cont_df['PBL_m'])*1.04)
    ax1.set_xlabel(r'\textbf{Time [sols]}')
    ax1.set_ylabel(r'\textbf{Vertical Coordinate [m]}')
    #  # XLIMITS in SOL_TIME
    #  xlimits_sols = ax1.get_xlim()
    #  # XLIMITS in DECIMAL_TIME
    #  hmin=cont_df['DECIMAL_TIME'][ cont_df['SOL_TIME']==xlimits_sols[0] ]
    #  hmax=cont_df['DECIMAL_TIME'][ cont_df['SOL_TIME']==xlimits_sols[1] ]
    #  xlimits_lmst = (hmin, hmax)
    #  #  ax1.set_ylim(0, 500)
    #  # # Manually set xticklabels as Hourly time
    #  # Get every 2 hours
    #  nh = 2  #ticks every 2 hours
    #  #  for i in x: if
    #  ax1.set_xticks( x )
    #  ax1.set_xticklabels( cont_df['DECIMAL_TIME'] )
    #  #  xtl = list( cont_df['DECIMAL_TIME'] )

    #--
    plt.savefig(join(outputdir, 'contour_plot.pdf'))
    plt.close('all')


    #-----------------------------------------------------
    # Create Contour ANIMATION  
    #-----------------------------------------------------
    #  fig, ax1 = plt.subplots(1, figsize=(10,5))
    fig, ax1 = plt.subplots(1, figsize=(8,10))
    #-- Contour Properties
    cmap_name = matplotlib.cm.inferno
    cmap_name.set_under(color='grey')
    vmin = 0.; #vmax=0.005
    vmax = max(cont_df['C_ppbv'])
    nlevels=100#50#20#50#100#50
    levels=np.linspace(0., max(cont_df['C_ppbv']), nlevels)
    #-- Colorbar
    cs.set_clim(vmin,vmax)
    cbar = fig.colorbar(cs)
    cbar.set_label(r'\textbf{CH$_4$ Abundance [ppbv]}', rotation=-90, va='bottom')
    #  cs.cmap.set_under('grey')
    def animate(i):
        #  # Try limiting to just a few frames at first
        #  if i<10:
        cs = ax1.contourf(X[:,i:i+2],Y[:,i:i+2],Z[:,i:i+2], levels, cmap=cmap_name, extend='min',vmin=vmin,vmax=vmax)
        # Try to remove weird lines between contours
        for c in ax1.contourf(X[:,i:i+2],Y[:,i:i+2],Z[:,i:i+2],levels, cmap=cmap_name, extend='min',vmin=vmin,vmax=vmax).collections:
            c.set_linewidth(0.1)
        #  #-- Properties
        #  ax1.set_title(r'$L_s =$ {:.0f}$^{{\circ}}$'.format(Ls_target))
        lmst = cont_df['DECIMAL_TIME'].iloc[i]
        ax1.set_title(r'$t=$ {:.2f} (LMST)'.format(lmst),y=-0.05)
        #  ax1.text(0.5,-0.05, r'\textbf{$t=$ {:.2f} (LMST)}'.format(lmst), transform=ax1.transAxes, ha='center')
        ax1.set_ylim(0, max(cont_df['PBL_m'])*1.04)
        # Make x axis v small so it's essentially just one time slice
        ax1.set_xlim(X[0][i], X[0][i]+0.0001)
        # Turn off x ticks and labels
        ax1.ticklabel_format(useOffset=False)
        ax1.tick_params(axis='x',which='both',
                        bottom=False, top=False, labelbottom=False)
        #  ax1.set_xlabel(r'\textbf{Time [sols]}')
        ax1.set_ylabel(r'\textbf{Vertical Coordinate [m]}')


    # Run the animation
    ani = FuncAnimation(fig, animate, frames=np.shape(X)[1]-1, interval=10, repeat=False) # full animation (minus 2 frames) 
    #  ani = FuncAnimation(fig, animate, frames=np.shape(X)[1]-45, interval=10, repeat=False) # only about 10 frames...
    # Save as mp4
    f = join(outputdir, 'contour_animation.mp4')
    writervideo = FFMpegWriter(fps=2)
    #  writervideo = FFMpegWriter(fps=60)
    ani.save(f, writer=writervideo)
    plt.close('all')
    #--------------------------------------------------- 
    # END NORMAL (TIME SERIES PLOTS) 


    #-----------------
    # MISC CALCULATONS
    #-----------------
    # Average annual flux
    #  flux_avg = np.mean(df_copy.copy()['FLUX_kg/m2/s_ADJ'])
    flux_avg = np.mean(df_copy.copy()['FLUX_kg/m2/s'])
    # How much time each year is in collapsed/expanded state?
    pbls = c_df[ c_df['SOL_TIME'].between(t1_oneyear,t2_oneyear) ]['PBL_m']
    n_c = len(pbls[ pbls<300.  ])  #number of entries collapsed state
    n_e = len(pbls[ pbls>=300. ])  #number of entries expanded  state
    # Amount of time in each
    c_sec = n_c * dt #[s]
    e_sec = n_e * dt #[s]


    #  for i in range(1,12):
        #  red = 0.905
        #  #  red = 15
        #  print('dof='+str(i))
        #  print('pdf   = {}'.format(stats.chi2.pdf(red,i)))
        #  print('1-pdf = {}'.format(1-stats.chi2.pdf(red,i)))
        #  print('cdf   = {}'.format(stats.chi2.cdf(red,i)))
        #  print('1-cdf = {}'.format(1-stats.chi2.cdf(red,i)))
        #  print('ppf   = {}'.format(stats.chi2.ppf(red,i)))
        #  print('----')
