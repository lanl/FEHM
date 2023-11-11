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
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
#the below import allows you to outline text
import matplotlib.patheffects as path_effects
import pickle
import pandas as pd
import glob

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
plt.rcParams['agg.path.chunksize'] = 10000  #agg has a hardcoded limit on points
#-----------------------------------------------------


def draw_brace(ax, xspan, yy, text):
    """
    Draws an annotated brace on the axes spanning a range.
    Note: xspan should be a tuple.
    https://stackoverflow.com/questions/18386210/annotating-ranges-of-data-in-matplotlib
    """
    xmin, xmax = xspan
    xspan = xmax - xmin
    ax_xmin, ax_xmax = ax.get_xlim()
    xax_span = ax_xmax - ax_xmin
    ymin, ymax = ax.get_ylim()
    yspan = ymax - ymin
    resolution = int(xspan/xax_span*100)*2+1 # guaranteed uneven
    beta = 300./xax_span # the higher this is, the smaller the radius
    x = np.linspace(xmin, xmax, resolution)
    x_half = x[:int(resolution/2)+1]
    y_half_brace = (1/(1.+np.exp(-beta*(x_half-x_half[0])))
                    + 1/(1.+np.exp(-beta*(x_half-x_half[-1]))))
    y = np.concatenate((y_half_brace, y_half_brace[-2::-1]))
    y = yy + (.05*y - .01)*yspan # adjust vertical position
    ax.autoscale(False)
    #  ax.plot(x, y, color='black', lw=1)
    ax.plot(x, y, color='black', lw=1, clip_on=False)
    ax.text((xmax+xmin)/2., yy+.07*yspan, text, ha='center', va='bottom')

def local_max(xvals, yvals, x_range):
    '''
    Find the index, xvalue (often time), and local max
    given a span of x values to look through.
    Outputs the index, xvalue, and yvalue of the local maxima.

    Parameters
    ----------
    xvals : ndarray
        All x-values of the data.
    yvals : ndarray
        All y-values of the data.
    x_range : list
        Range of x values over which to find the local maxima.
    '''
    #Get the range indices
    ix1 = np.where(xvals == find_nearest(xvals,x_range[0]))[0][0]
    ix2 = np.where(xvals == find_nearest(xvals,x_range[1]))[0][0]
    #Get index of max value
    i_lmax = np.argmax(yvals[ix1:ix2])
    #Get the x value of the local max
    if isinstance(xvals, pd.Series):  x_lmax = xvals[ix1:ix2].iloc[i_lmax]
    else: x_lmax = xvals[ix1:ix2][i_lmax]
    #Get the y value of the local max
    if isinstance(xvals, pd.Series): y_lmax = yvals[ix1:ix2].iloc[i_lmax]
    else: y_lmax = yvals[ix1:ix2][i_lmax]
    #Now get the global index of the max (since the full dataset is larger)
    global_ilmax = np.where(xvals==x_lmax)[0][0]
    return global_ilmax, x_lmax, y_lmax

def timeAveragedFlux(all_times, all_fluxes, time_range):
    ''' Calculate time-integrated flux over ``time_range``). '''
    # Get time range indices
    stime = time_range[0]
    etime = time_range[1]
    sti = np.where(all_times == find_nearest(all_times, stime))[0][0] #start time index
    eti = np.where(all_times == find_nearest(all_times, etime))[0][0] #end   time index
    # Calculate cumulative flux in the time range
    cum_flux = integrate.cumtrapz(all_fluxes[sti:eti], x=all_times[sti:eti]*solsec())[-1] #[kg/km2]
    # Get time-averaged flux by dividing by number of seconds in time range
    avg_flux = cum_flux / ((etime-stime) * solsec())    #[kg/km2/s]
    return avg_flux

def read_mesh(run_dir, fehmn_files_name='fehmn.files'):
    ''' Read mesh used from fehmn.files and extract key quantities. '''
    fehmn_files_path = join(run_dir, fehmn_files_name)
    cwd = os.getcwd()
    with open(fehmn_files_path,'r') as ff:
        i=0
        for line in ff.readlines():
            if line.startswith('grida'):
                i+=1
                os.chdir(run_dir); wd = os.getcwd()
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
    #---- Calculate mesh depth
    mesh_depth = max(grid['y'])-min(grid['y']) #[m]
    #---- Calculate mesh width for flux calculations
    mesh_xwidth = max(grid['x']) - min(grid['x'])           #[m]
    mesh_zwidth = 1.                                        #[m]
    #---- Calculate mesh "surface" area in km^2 for flux calcs
    surface_area = mesh_xwidth/1000. * mesh_zwidth/1000.    #[km^2]
    return grid, mesh_depth, mesh_xwidth, mesh_zwidth, surface_area

def getIndex(arr, val):
    ''' Get index of closest value in array. '''
    if isinstance(val, str):
        ind = np.where(arr==val)[0][0]
    else:
        ind =np.where(arr==find_nearest(arr, val))[0][0]
    return ind

def synth_seasons(justNames=False):
    '''
    Read in the synthetic mission sols + solar longitude (L_s) data.

    Create a dictionary with the time ranges [sols] of the bounds of each
    season (useful when dataset is very long).
    '''
    if justNames==False:
        synthdir='/project/gas_seepage/jportiz/mars/analytical/fourier_analysis/generate_synthetic_data/output/'
        synthfile='synthetic_SOL_TIME_Ls.csv'
        df = pd.read_csv(join(synthdir,synthfile), index_col=None)
        # Add SEASON name to data frame
        sss = []
        for l in df['L_s']: sss.append(mars_seasons(l))
        df.insert(df.shape[1],column='SEASON',value=sss)
    # If justNames==True, skip reading df to save time
    else: df = None
    # Definitions of seasons
    season_dict = {'N. Spring': [0,90],
                   'N. Summer': [90,180],
                   'N. Autumn': [180,270],
                   'N. Winter': [270,360]
                  }
    #Season colors (for plotting)
    #https://imagecolorpicker.com/
    season_colors = {'N. Spring': '#c1db75',  #nice, yellow green
                     'N. Summer': '#ffde7a',  #gold/yellow "Grandis" #ffd17a 
                     'N. Autumn': '#ff8f00',  #"Pizzaz" '#fc7c00'a nice tangerine
                     'N. Winter': '#a4bce3',  #a smoky "Regent St." blue
                    }
    return season_dict, season_colors, df

def add_seasonsToPlot(axis, tLimits=[47000, 47668], num_cycles=1, extraSeason=None, axis_label_fontsize=14):
    '''
    Adds seasons as highlighted areas to provided plot axis.
    NOTE: only works on time series plots. Very simple to just code up the seasons on a solar Longitude (L_s) plot.

    Parameters
    ----------
    axis : AxesSubPlot object
        Matplotlib axis to add seasons to.
    tLimits : list or tuple
        min and max time you want to generate season ranges list. For
        simplicity, can just be your plot/steady-state starting time to the
        end.
    num_cycles : int
        Number of seasonal cycles you want to plot.
    extraSeason : str
        In case the first season plotted only has a bit showing, you can
        specify an extra season to plot in addition to those in ``num_cycles``.
        Must be one of ['N. Spring', 'N. Summer', 'N. Autumn', 'N. Winter'].
    axis_label_fontsize : int or float
        Font size of season name label at top middle of each season.
    '''
    # Call the helper function above
    season_dict, season_colors, df_seasons = synth_seasons()
    # Get all date ranges (t1,t2)  for each season in plot window
    season_ranges = {'N. Spring': [],
                     'N. Summer': [],
                     'N. Autumn': [],
                     'N. Winter': [],
                     }
    season = df_seasons['SEASON'].iloc[getIndex(df_seasons['SOL_TIME'], tLimits[0])] #get first season (of shortened results)
    t1 = tLimits[0]              #get first sol_time (of shortened results)
    t2 = 0.
    df = df_seasons.copy()
    # Shorten dataframe
    df = df[df.SOL_TIME>=t1]
    c = 0
    while t2 <= tLimits[1]:
        df = df[df.SOL_TIME>t2]         #shorten the dataframe each iteration 
        t1 = df.SOL_TIME.iloc[0]               #update t1
        try: i = df[df.SEASON != season].index[0]   #first index not equal to t1 season
        except IndexError: break               #end loop if no more seasons 
        t2 = df.SOL_TIME.loc[i-1]              #last time w/ same season as t1
        season_ranges[season].append([t1,t2])  #add this range to the dictionary
        c+=1
        season = df.SEASON.loc[i]  #now get the "next" season, restart loop

    # Now Plot all the season ranges in the season_ranges dict
    for season in season_ranges:
        for r in season_ranges[season][:num_cycles]:
            axis.axvspan(r[0],r[1], color=season_colors[season],alpha=0.75,zorder=0)
    # also plot the the extraSeason season if chosen
    if extraSeason is not None:
        r = season_ranges[extraSeason][:num_cycles+1][-1]
        axis.axvspan(r[0],r[1], color=season_colors[season],alpha=0.75,zorder=0)
    # Add text saying "autumn", "spring", etc to the first cycle
    for season in season_ranges:
        for r in season_ranges[season][:num_cycles]:
            # Skip the 1st extraSeason'd season label and only label it the second time
            if season != extraSeason:
                mdp = (r[1]+r[0])/2 #midpoint of seasonal range
                #  axis.text(mdp, 0.93*max_flux, season.replace(' ','\n'), ha='center',size=axis_label_fontsize)
                axis.text(mdp, 0.93*axis.get_ylim()[-1], season.replace(' ','\n'), ha='center',size=axis_label_fontsize)
    if extraSeason is not None:
        r = season_ranges[extraSeason][:num_cycles+1][-1]
        mdp = (r[1]+r[0])/2 #midpoint of seasonal range
        #  axis.text(mdp, 0.93*max_flux, season.replace(' ','\n'), ha='center',size=axis_label_fontsize)
        axis.text(mdp, 0.93*axis.get_ylim()[-1], season.replace(' ','\n'), ha='center',size=axis_label_fontsize)

    return season_ranges


def make_plots(fehm_upshift, T_upshift):

    cwd = os.getcwd()
    #  results_dir = join(cwd,'tracer')
    #  results_dir = join(cwd) #for "/continue_originalb "
    if len(glob.glob('*restart*'))>1:
        print('******************************************')
        print('                Warning ')
        print('******************************************')
        print('More than one *restart*/ directory detected...')
        print(glob.glob('*restart*'))
        print()
        results_base = input('Please enter name of results directory for plotting:')
    else:
        results_base = glob.glob('*restart*')[0]
    results_dir = join(cwd, results_base) #for the "restart" tracer sim in this dir
    output_dir = join(results_dir,'output')
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print('Going to put all output in {}...'.format(join(results_base,'output')))
    print('Plotting results from "{}" directory...'.format(results_base))
    #-------------------------------------------------- `
    #           MESH DATA 
    #-------------------------------------------------- 
    #  fehmnFiles = join(results_dir,'fehmn.files')
    # Extract domain metrics from mesh
    grid, depth, mesh_xwidth, mesh_zwidth, surface_area = read_mesh(results_dir)

    #  rho_w = dd[1] #[kg/m3]  (use the density in FEHM)
    topnode = np.where((grid['x']==0) & (grid['y']==0.))[0][0] + 1 #any top node
    botnode = np.where((grid['x']==0) & (grid['y']==min(grid['y'])))[0][0] + 1 #bot node

    #-------------------------------------------------- 
    #           FEHM  DATA 
    #-------------------------------------------------- 
    pf = os.path.join(results_dir, 'run_presWAT.his')
    df = os.path.join(results_dir, 'run_denWAT.his')
    #  cf = os.path.join(results_dir, 'run_methane_l.trc')
    gf = os.path.join(results_dir, 'run_glob.his')
    zf = os.path.join(results_dir, 'runmethane_l.cflx')
    #  zf = os.path.join(results_dir, 'runmethane.cflx')
    #---- Pressures
    #  pdata = np.loadtxt(pf,skiprows=4)
    pdata = np.genfromtxt(pf,skip_header=4, invalid_raise=False) #skips bad rows
    p_time = days2sols(pdata[:,0]) #convert days to sols
    p_cols = np.loadtxt(pf,skiprows=3,max_rows=1,dtype='str')[1::2]
    pdata = (pdata - fehm_upshift) #[Pa] remove the artificial pressure upshift from FEHM
    # Put the top pressures in a pandas DataFrame
    press_df = pd.DataFrame({'SOL_TIME':p_time, 'PRESS':pdata[:,getIndex(p_cols, str(topnode))]*1e6})
    #---- Temperatures (from a different simulation)
    temps_df = pd.read_csv( 'surface_temperatures.csv', sep=',', names=['SOL_TIME', 'TEMP'], index_col=None )
    # Correct the temperatures and convert to Kelvin
    temps_df['TEMP'] = temps_df['TEMP'] - T_upshift + 273.15
    #  #---- Density                                    (not needed)
    #  dd = np.genfromtxt(df,skip_header=4,max_rows=1)
    #  #---- Concentrations                             (not needed)
    #  cdata = np.loadtxt(cf,skiprows=4)
    #  c_time = days2sols(cdata[:,0]) #convert days to sols
    #  c_cols = np.loadtxt(cf,skiprows=3,max_rows=1,dtype='str')[1::2]
    #  #---- Global balances                            (not needed)
    #  gdata = np.loadtxt(gf,skiprows=4)
    #  g_time = days2sols(gdata[:,0]) #convert days to sols
    #  g_mass = gdata[:,1] #[kg] delta mass (doesn't seem quite right)
    #---- Zone solute fluxes
    with open(zf,'r') as zh:
        z_cols = zh.readline().strip().split()[4:]
    z_cols.insert(0,'Time')
    #  zd = np.loadtxt(zf,skiprows=1,dtype=[(i,'float') for i in z_cols]) #[kg/d]
    zd = np.genfromtxt(zf,skip_header=1,dtype=[(i,'float') for i in z_cols], invalid_raise=False) #[kg/d]
    z_time = days2sols(zd['Time'])                          #[sols]
    zdata = zd.view(np.recarray)                            #[kg/d]  FEHM "flux"
    zflux = zdata.Sink*days2sols(1)/solsec()/surface_area   #[kg/km^2/s] flux 

    #-------------------------------------------------- 
    #       SAVE PROCESSED FEHM DATA TO CSV
    #-------------------------------------------------- 
    # Put the fluxes in a pandas DataFrame
    flux_df = pd.DataFrame({'SOL_TIME':z_time, 'FLUX':zflux})
    # Add L_s to DataFrame
    sd, sc, ls_df = synth_seasons()
    # Interpolate to the times in flux_df
    f0 = interp1d(ls_df['SOL_TIME'], ls_df['L_s'])
    flux_df.insert( np.shape(flux_df)[1]-1, 'L_s', f0(flux_df['SOL_TIME']) )
    # Add season name to DataFrame
    flux_df.insert( np.shape(flux_df)[1], 'SEASON', [mars_seasons(l) for l in flux_df['L_s']] )
    #---- Create a FULL DataFrame with all the FEHM data in it 
    full_df = flux_df.copy()
    # Interpolate the pressures to the flux times 
    f0 = interp1d(press_df['SOL_TIME'], press_df['PRESS'], fill_value='extrapolate')
    full_df.insert( np.shape(full_df)[1], 'PRESS', f0(full_df['SOL_TIME']) )
    # Interpolate the temperatures to the flux times 
    f0 = interp1d(temps_df['SOL_TIME'], temps_df['TEMP'])
    full_df.insert( np.shape(full_df)[1], 'TEMP', f0(full_df['SOL_TIME']) )
    #---- Save full DataFrame to CSV 
    filename = 'fehm_results.csv'
    full_df.to_csv( join(results_dir, 'output', filename), sep=',', index=None)

    #  #---- TEST (read the csv)
    #  a = pd.read_csv(join(results_dir,'output',filename))

    #-------------------------------------------------- 
    #           CALCULATE TIME-AVERAGED FLUX 
    #-------------------------------------------------- 
    #  #---- Print this info to a file
    #  tem = sys.stdout
    #  name = 'fluxCalcs'
    #  sys.stdout = f = open(name, 'w')
    #---- Print flux info to a file
    print('Writing ``fluxCalcs`` file...')
    name = 'fluxCalcs'
    fff = open(join(output_dir,name),'w')

    cum_sink = integrate.cumtrapz(flux_df['FLUX']*surface_area,x=flux_df['SOL_TIME']*solsec()) #[kg]
    #---- Choose the zoomed-in timespan for the zoomed publication plot
    targ_zstart = 47000             #[sols] target start of zoomed-in time
    targ_zend = targ_zstart+2700    #[sols] target end   of zoomed-in time
    #---- Find indices of data in zoomed-in window
    pstart_i = getIndex(p_time, targ_zstart)
    pend_i   = getIndex(p_time, targ_zend)
    zstart_i = getIndex(z_time, targ_zstart)
    zend_i   = getIndex(z_time, targ_zend)

    #---- Save pressure and flux data from these windows as .csv
    presscols = np.column_stack( (p_time[pstart_i:pend_i], pdata[pstart_i:pend_i,np.where(p_cols==str(topnode))[0][0]]*1e6) )
    fluxcols = np.column_stack( (z_time[zstart_i:zend_i], zflux[zstart_i:zend_i]) )
    np.savetxt(join(output_dir,'ssPressures.csv'), presscols, delimiter=',')
    np.savetxt(join(output_dir,'ssFluxes.csv'), fluxcols, delimiter=',')


    stime2 = 45000. #start
    etime2 = 47000. #end
    sti2 = getIndex(z_time, stime2)                            #start time index
    eti2 = getIndex(z_time, etime2)                            #end time index
    stime3 = 45250. #start
    etime3 = 45450. #end
    sti3 = getIndex(z_time, stime3)                            #start time index
    eti3 = getIndex(z_time, etime3)                            #end time index
    stime4 = 45550. #start
    etime4 = 45850. #end
    sti4 = getIndex(z_time, stime4)                            #start time index
    eti4 = getIndex(z_time, etime4)                            #end time index
    stime5 = 47200. #start         #NEW
    etime5 = 47600.#47500. #end    #NEW
    sti5 = getIndex(z_time, stime5)                            #start time index
    eti5 = getIndex(z_time, etime5)                            #end time index
    stime6 = 47600.#47500. #start  #NEW
    etime6 = 47900.#47800. #end    #NEW
    sti6 = getIndex(z_time, stime6)                            #start time index
    eti6 = getIndex(z_time, etime6)                            #end time index

    flux_zoomed  = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[targ_zstart, targ_zend])
    flux_zoomed2 = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime2, etime2])
    flux_zoomed3 = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime3, etime3])
    flux_zoomed4 = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime4, etime4])
    flux_zoomed5 = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime5, etime5])
    flux_zoomed6 = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime6, etime6])


    #---- Calculate time-averaged flux for different window sizes around the first "big" spike and "little" spike
    #Indices, times, and fluxes of the big and little methane spikes
    #  b_spike = local_max(z_time, zflux, x_range=[stime5, etime5]) #big    spike
    #  l_spike = local_max(z_time, zflux, x_range=[stime6, etime6]) #little spike
    b_spike = local_max(flux_df['SOL_TIME'], flux_df['FLUX'], x_range=[stime5, etime5]) #big    spike
    l_spike = local_max(flux_df['SOL_TIME'], flux_df['FLUX'], x_range=[stime6, etime6]) #little spike
    #Windows
    win_marsYear  = 668          #[sols] 
    win_long      = 334.8        #[sols] (Mumma et al, 2009)
    win_med       = 116.8        #[sols] (Mumma et al, 2009)
    win_short     = 58.4         #[sols] (Mumma et al, 2009)
    win_5         = 5.           #[sols] (Giuranna et al, 2009) [max release duration]
    win_vshort    = 0.72         #[sols] (Giuranna et al, 2019) [17.5 hrs; from SI]

    #Store results in dictionay {window: [(start, end), flux]}
    bigspike_dict = {}
    littlespike_dict = {}
    for win in [win_marsYear,win_long,win_med,win_short,win_5,win_vshort]:
        #Do big spike first
        stime = b_spike[1]-win/2
        etime = b_spike[1]+win/2
        #  flx = timeAveragedFlux(z_time, zflux, time_range=[stime,etime])
        flx = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime,etime])
        bigspike_dict[win] = [ (stime, etime), flx ]
        fff.write('"Big spike" flux (window={:.1f})     = {:.4g} kg/s/km^2 (b/w sols {:.0f} - {:.0f})\n'.format(win,flx,stime,etime))
        #Now do  little spike
        stime = l_spike[1]-win/2
        etime = l_spike[1]+win/2
        #  flx = timeAveragedFlux(z_time, zflux, time_range=[stime,etime])
        flx = timeAveragedFlux(flux_df['SOL_TIME'], flux_df['FLUX'], time_range=[stime,etime])
        littlespike_dict[win] = [ (stime, etime), flx ]
        fff.write('"Little spike" flux (window={:.1f})  = {:.4g} kg/s/km^2 (b/w sols {:.0f} - {:.0f})\n'.format(win,flx,stime,etime))
    #--- Print key results in same format as the article table (for easy reading)
    factor = 1e-9  #[kg/km2/s] print all results relative to factor
    fff.write('\n')
    fff.write('==========================================================\n')
    fff.write('                 Time-averaged flux        \n')
    fff.write('                  10^-9 [kg/km^2/s]        \n')
    fff.write('----------------------------------------------------------\n')
    fff.write('    Northern Summer    ||     Northern Winter       \n')
    fff.write('Window duration [sols] || Window duration [sols]\n')
    fff.write(' 58       117     334  ||  58       117     334 \n')
    fff.write('----------------------------------------------------------\n')
    fff.write('{:.4f}  {:.4f}  {:.4f} || {:.4f}  {:.4f}  {:.4f}\n'.format(bigspike_dict[58.4][-1]/factor, bigspike_dict[116.8][-1]/factor, bigspike_dict[334.8][-1]/factor, littlespike_dict[58.4][-1]/factor, littlespike_dict[116.8][-1]/factor, littlespike_dict[334.8][-1]/factor))
    fff.write('==========================================================\n')
    fff.write('\n\n')
    # Also print out calculation of surface area required to supply yearly methane budget (Korablev2019; 4.0 kg/sol)
    fff.write('----------\n')
    fff.write('AREA REQUIRED TO SUPPLY ANNUAL METHANE BUDGET:\n')
    req_area = 1/ (bigspike_dict[win_marsYear][1]*88775.24/4.0) #[km2]
    fff.write('Required area                 = {:.0f} km^2\n'.format(req_area))
    fff.write('% comparison to Gale crater   = {:.2%}\n'.format(req_area/18626.))
    fff.write('% comparison to total Mars SA = {:%}\n'.format(req_area/144806235.))
    fff.close()
    #  os.system('cat '+name)
    os.system('cat '+join(output_dir,name))
    print('    Done.')



    #-------------------------------------------------- 
    #           PLOT
    #-------------------------------------------------- 
    ms = 5 #markersize
    lw = 1.5  #linewidth
    axls = 14.     # axis label size

    fig = plt.figure(figsize=(12,9))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    # PRESSURE
    ax2.plot(p_time[::10],pdata[:,np.where(p_cols==str(topnode))[0][0]][::10]*1e6,ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=0)
    p_lims = ax2.get_ylim(); p_range = p_lims[-1]-p_lims[0]
    ax2.set_ylabel('Barometric pressure [Pa]',rotation=-90,va='bottom')
    ax2.set_ylim(p_lims[0]-.5*p_range, p_lims[1]+0.5*p_range)
    # MASS FLOW RATE
    ax1.plot(z_time,zdata.Sink*days2sols(1)/solsec(), label='sink',c='k',ls='-')
    #  ax1.set_yscale('log')
    #  ax1.set_ylim(bottom=1e-10,top=1e-4)
    #  ax1.set_title('Baro pumping, initial source\nmesh: {}'.format(mesh_name))
    #  ax1.set_title('Baro pumping, {}\nmesh: {}'.format(sourcetype,mesh_name))
    ax1.set_xlabel('Time [sols]')
    ax1.set_ylabel(r'CH$_4$ mass flow rate [kg/s]')
    ax1.xaxis.set_minor_locator(MultipleLocator(100))
    fig.tight_layout()
    ax1.set_xlim(right=60000)
    ax2.set_xlim(right=60000)
    #  ax1.legend()
    #  plt.figtext(0.02, 0.95, r'Total CH$_4$ seepage: {:0.2g} kg'.format(cum_sink[-1]),fontsize=axls,transform=ax1.transAxes)
    plt.figtext(0.02, 0.95, r'Total CH$_4$ seepage: '+tools.sci_notation(cum_sink[-1])+ ' kg',fontsize=axls,transform=ax1.transAxes)
    fig.savefig(join(output_dir,'non-log_massflowrate_plot.pdf'))
    plt.close(fig)

    #---- ZOOMED-IN MASS FLOW RATE PLOT
    ax2 = ax1.twinx()
    # PRESSURE
    ax2.plot(p_time[::10],pdata[:,np.where(p_cols==str(topnode))[0][0]][::10]*1e6,ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=0)
    p_lims = ax2.get_ylim(); p_range = p_lims[-1]-p_lims[0]
    ax2.set_ylabel('Barometric pressure [Pa]',rotation=-90,va='bottom')
    ax2.set_ylim(p_lims[0]-.5*p_range, p_lims[1]+0.5*p_range)
    ax2.set_xlim(left=22000.,right=60000.)
    # MASS FLOW RATE
    ax1.plot(z_time,zdata.Sink*days2sols(1)/solsec(), label='sink',c='k',ls='-')
    #  ax1.set_yscale('log')
    #  ax1.set_ylim(bottom=1e-10,top=1e-4)
    #  ax1.set_title('Baro pumping, initial source\nmesh: {}'.format(mesh_name))
    #  ax1.set_title('Baro pumping, {}\nmesh: {}'.format(sourcetype,mesh_name))
    ax1.set_xlabel('Time [sols]')
    ax1.set_ylabel(r'CH$_4$ mass flow rate [kg/s]')
    ax1.set_xlim(left=22000.,right=52000.)
    #  ax1.xaxis.set_minor_locator(MultipleLocator(100))
    # Get max flux in this window
    #  ax1.set_ylim(bottom=0.0,top=0.2e-7)
    t1 = ax1.get_xlim()[0]; t2 = ax1.get_xlim()[1]
    maxval= max(np.asarray(zdata.Sink*days2sols(1))[getIndex(z_time,t1):getIndex(z_time,t2)]/solsec())
    ax1.set_ylim(0., maxval)
    fig.tight_layout()
    #  ax1.legend()
    #  plt.figtext(0.02, 0.95, r'Total CH$_4$ seepage: {:0.2g} kg'.format(cum_sink[-1]),fontsize=axls,transform=ax1.transAxes)
    plt.figtext(0.02, 0.95, r'Total CH$_4$ seepage: '+tools.sci_notation(cum_sink[-1])+ ' kg',fontsize=axls,transform=ax1.transAxes)
    fig.savefig(join(output_dir,'zoomed_non-log_massflowrate_plot.pdf'))
    plt.close(fig)


    #---- CUMULATIVE MASS FLOW RATE PLOT
    fig = plt.figure(figsize=(12,9))
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    # PRESSURE
    ax2.plot(p_time[::10],pdata[:,np.where(p_cols==str(topnode))[0][0]][::10]*1e6,ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=0)
    p_lims = ax2.get_ylim(); p_range = p_lims[-1]-p_lims[0]
    ax2.set_ylabel('Barometric pressure [Pa]',rotation=-90,va='bottom')
    ax2.set_ylim(p_lims[0]-.5*p_range, p_lims[1]+0.5*p_range)

    # CUMULATIVE MASS FLOW OUT
    cum_sink = integrate.cumtrapz(zdata.Sink*days2sols(1)/solsec(),x=z_time*solsec())
    ax1.plot(z_time[:-1], cum_sink, label='cum_sink',c='k',ls='-')

    #  ax1.set_title('Baro pumping, initial source\nmesh: {}'.format(mesh_name))
    #  ax1.set_title('Baro pumping, {}\nmesh: {}'.format(sourcetype,mesh_name))
    ax1.set_xlabel('Time [sols]')
    ax1.set_ylabel(r'Cumulative CH$_4$ mass flow [kg]')
    ax1.xaxis.set_minor_locator(MultipleLocator(100))
    ax1.set_xlim(right=60000)
    ax2.set_xlim(right=60000)
    fig.tight_layout()
    fig.savefig(join(output_dir,'cumulative_massflow_plot.pdf'))
    plt.close(fig)

    #--- PLOT FOR PUBLICATIONS ---
    fig = plt.figure(figsize=(12,9))
    ax1 = fig.add_subplot(111)
    axls = 28.     # axis label size
    cum_sink = integrate.cumtrapz(zdata.Sink*days2sols(1)/solsec(),x=z_time*solsec())
    ax2 = ax1.twinx()
    # PRESSURE
    ax2.plot(p_time[::10],pdata[:,np.where(p_cols==str(topnode))[0][0]][::10]*1e6,ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=0)
    p_lims = ax2.get_ylim(); p_range = p_lims[-1]-p_lims[0]
    ax2.set_ylabel('Barometric pressure [Pa]',rotation=-90,va='bottom',fontsize=axls)
    ax2.yaxis.label.set_color('darkgrey')
    ax2.tick_params(axis='y', colors='darkgrey')
    ax2.set_ylim(p_lims[0]-.5*p_range, p_lims[1]+0.5*p_range)
    # MASS FLOW RATE
    ax1.plot(z_time,zdata.Sink*days2sols(1)/solsec(), label='sink',c='k',ls='-')
    #  ax1.set_title('Baro pumping, initial source\nmesh: {}'.format(mesh_name))
    #  ax1.set_title('Baro pumping, {}\nmesh: {}'.format(sourcetype,mesh_name))
    #  ax1.set_xticklabels(fontsize=24.)
    #  ax1.set_yticklabels(fontsize=24.)
    #  ax2.set_yticklabels(fontsize=24.)
    #  plt.rcParams['xtick.labelsize']=24.
    #  plt.rcParams['ytick.labelsize']=24.
    ax1.tick_params(axis="x", labelsize=24)
    ax1.tick_params(axis="y", labelsize=24)
    ax2.tick_params(axis="y", labelsize=24)
    ax1.set_xlabel('Time [sols]',fontsize=axls)
    ax1.set_ylabel(r'CH$_4$ mass flow rate [kg/s]',fontsize=26.)
    ax1.xaxis.set_minor_locator(MultipleLocator(100))
    ax1.set_xlim(right=60000)
    ax2.set_xlim(right=60000)
    fig.tight_layout()
    #  ax1.legend()
    #  plt.figtext(0.02, 0.95, r'Total CH$_4$ seepage: {:0.2g} kg'.format(cum_sink[-1]),fontsize=axls,transform=ax1.transAxes)
    #  plt.figtext(0.02, 0.95, r'Total CH$_4$ seepage: '+tools.sci_notation(cum_sink[-1])+ ' kg',fontsize=axls,transform=ax1.transAxes)
    fig.savefig(join(output_dir,'pub_non-log_massflowrate_plot.pdf'))
    plt.close(fig)

    #--- ZOOM FLUX PLOT FOR PUBLICATIONS ---
    fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(111)
    axls = 28.     # axis label size
    cum_sink = integrate.cumtrapz(flux_df['FLUX']*surface_area,x=flux_df['SOL_TIME']*solsec())
    ax2 = ax1.twinx()
    # PRESSURE
    #  ax2.plot(p_time[::10],pdata[:,np.where(p_cols==str(topnode))[0][0]][::10]*1e6,ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=0)
    ax2.plot(full_df['SOL_TIME'][::10],full_df['PRESS'][::10],ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=0)
    p_lims = ax2.get_ylim(); p_range = p_lims[-1]-p_lims[0]
    ax2.set_ylabel('Barometric pressure [Pa]',rotation=-90,va='bottom',fontsize=axls)
    ax2.yaxis.label.set_color('darkgrey')
    ax2.tick_params(axis='y', colors='darkgrey')
    ax2.set_ylim(p_lims[0]-.5*p_range, p_lims[1]+0.5*p_range)
    # MASS FLUX 
    max_flux = max(full_df['FLUX'][zstart_i:zend_i])
    min_flux = min(full_df['FLUX'][zstart_i:zend_i])
    #  ax1.plot(z_time,zdata.Sink*days2sols(1)/solsec()/surface_area, label='sink',c='k',ls='-')
    ax1.plot(full_df['SOL_TIME'],full_df['FLUX'], label='sink',c='k',ls='-')
    ax1.tick_params(axis="x", labelsize=24)
    ax1.tick_params(axis="y", labelsize=24)
    ax2.tick_params(axis="y", labelsize=24)
    ax1.set_xlabel('Time [sols]',fontsize=axls)
    ax1.set_ylabel(r'Surface CH$_4$ flux [kg/km$^2$/s]', fontsize=axls)
    # SETTINGS
    ax1.xaxis.set_minor_locator(MultipleLocator(100))
    ylim_flux = max_flux*1.1
    ax1.set_ylim(min_flux,ylim_flux)
    #  ax1.set_ylim(min_flux,max_flux*1.1)
    #  ax1.set_ylim(min_flux,8.e-7) #same ylim for all the run_14d* plots
    ax1.set_xlim(targ_zstart,targ_zend) #x limits (steady-state)
    ax2.set_xlim(targ_zstart,targ_zend) #x limits ( ''      '' )
    #Annotate with curly braces stating time-average flux for certain time ranges
    an1 = r"{}".format(sci_notation(flux_zoomed5,2)) + "\n" + "kg/km$^2$/s"
    an2 = r"{}".format(sci_notation(flux_zoomed6,2)) + "\n" + "kg/km$^3$/s"
    an_ypos = ylim_flux #eventually may change this position
    #  an_ypos = max_flux #eventually may change this position
    xshift = 11. #space between otherwise-touching braces [sols]
    draw_brace(ax2, (stime5,etime5), an_ypos, an1)
    draw_brace(ax2, (stime6+xshift,etime6+xshift), an_ypos, an2)
    fig.tight_layout()
    fig.savefig(join(output_dir,'pub_zoom_non-log_flux_plot.pdf'))
    plt.close(fig)

    #--- ZOOM FLUX PLOT WITH INSET FOR PUBLICATIONS ---
    fig = plt.figure(figsize=(13,6))
    ax1 = fig.add_subplot(111)
    axls = 24#28.     # axis label size
    cum_sink = integrate.cumtrapz(flux_df['FLUX']*surface_area,x=flux_df['SOL_TIME']*solsec())
    ax2 = ax1.twinx()
    # PRESSURE
    ax2.plot(full_df['SOL_TIME'][::10],full_df['PRESS'][::10],ls='-',lw=lw,color='lightgrey',alpha=0.7,zorder=2)
    p_lims = ax2.get_ylim(); p_range = p_lims[-1]-p_lims[0]
    ax2.set_ylabel(r'\textbf{Barometric pressure [Pa]}',rotation=-90,va='bottom',fontsize=axls-6)
    ax2.yaxis.label.set_color('darkgrey')
    ax2.tick_params(axis='y', colors='darkgrey')
    ax2.set_ylim(p_lims[0]-.5*p_range, p_lims[1]+0.5*p_range)
    # MASS FLUX 
    max_flux = max(full_df['FLUX'][zstart_i:zend_i])
    min_flux = min(full_df['FLUX'][zstart_i:zend_i])
    ax1.plot(z_time,zdata.Sink*days2sols(1)/solsec()/surface_area, label='sink',c='k',ls='-',zorder=1)
    # PROPERTIES
    ax1.tick_params(axis="x", labelsize=axls-4)
    ax1.tick_params(axis="y", labelsize=axls-4)
    ax2.tick_params(axis="y", labelsize=axls-4)
    ax1.set_xlabel(r'\textbf{Time [sols]}',fontsize=axls-6)
    ax1.set_ylabel(r'\textbf{Surface CH$_4$ flux [kg/km$^2$/s]}', fontsize=axls-6)
    # SETTINGS
    ax1.xaxis.set_minor_locator(MultipleLocator(100))
    ax1.set_ylim(min_flux, ylim_flux)
    #  ax1.set_ylim(min_flux,max_flux*1.1)
    #  ax1.set_ylim(min_flux,8.e-7) #same ylim for all the run_14d* plots
    ax1.set_xlim(targ_zstart,targ_zend)
    ax2.set_xlim(targ_zstart,targ_zend)
    #---- ANNOTATE WITH CURLY BRACES STATING TIME-AVERAGE FLUX FOR CERTAIN TIME RANGES
    win = 334.8     #[sols] flux window size to plot (334.8, 116.8, 58.4)
    an1 = r"{}".format(sci_notation(bigspike_dict[win][1],2)) + "\n" + "kg/km$^2$/s"
    an2 = r"{}".format(sci_notation(littlespike_dict[win][1],2)) + "\n" + "kg/km$^2$/s"
    an_ypos = ylim_flux #eventually may change this position
    #  an_ypos = max_flux #eventually may change this position
    xshift = 10. #space between otherwise-touching braces [sols]
    sb = bigspike_dict[win][0][0]; eb = bigspike_dict[win][0][1]
    sl = littlespike_dict[win][0][0]; el = littlespike_dict[win][0][1]
    draw_brace(ax1, (sb,eb), an_ypos, an1)
    draw_brace(ax1, (sl,el), an_ypos, an2)
    fig.tight_layout()
    #---- DO AN INSET PLOT ZOOMING IN TO DIURNAL FLUCTUATIONS --------
    # Create a set of inset Axes: these should fill the bounding box allocated to
    # them.
    ins = local_max(full_df['SOL_TIME'], full_df['FLUX'], [stime6,etime6]) #find local max flux info within this time range
    ins_win = 10 #[sols]  window size of inset
    stime_ins = ins[1]-ins_win/2 #start time
    etime_ins = ins[1]+ins_win/2 #end   time
    sti_ins = getIndex(full_df['SOL_TIME'],stime_ins) #start time index
    eti_ins = getIndex(full_df['SOL_TIME'],etime_ins) #end time index
    axins = plt.axes([0,0,1,1])
    # Manually set the position and relative size of the inset axes within ax1
    ip = InsetPosition(ax1, [0.7,0.4,0.25,0.55]) #leftedge,bottomedge,width,height
    axins.set_axes_locator(ip)
    # The data: only display for the selected span 
    #  axins.plot(z_time[sti_ins:eti_ins], zflux[sti_ins:eti_ins], c='k',zorder=3)
    axins.plot(full_df['SOL_TIME'][sti_ins:eti_ins], full_df['FLUX'][sti_ins:eti_ins], c='k',zorder=3)
    # Mark the region corresponding to the inset axes on ax1 and draw lines
    # in grey linking the two axes.
    mark_inset(ax1, axins, loc1=2, loc2=4, fc='none', ec='0.5',zorder=5)#ec='0.5',fc='none'
    #------------------------------------------------------------------
    fig.savefig(join(output_dir,'pub_zoom_non-log_flux_wInset_plot.pdf'))

    #------------------------------------------------------------------
    #   ADD ON HIGHLIGHTED ZONES SHOWING "SEASONS"
    #------------------------------------------------------------------
    add_seasonsToPlot(ax1, tLimits=[targ_zstart, targ_zend], num_cycles=1, extraSeason='N. Winter', axis_label_fontsize=axls/2-1)
    fig.savefig(join(output_dir,'pub_zoom_non-log_flux_wInsetAndSeasons_plot.pdf'))
    #  fig.savefig(join(output_dir,'pub_zoom_non-log_flux_wInsetAndSeasons_plot.png'),format='png',transparent=True,dpi=600) #PNG plot for AGU poster
    plt.close(fig)
    plt.close('all')


    #--- L_S FLUX PLOT COMPARING TO ABUNDANCES IN WEBSTER_2018a (Fig 1b)
    #Read in the Webster2018a methane abundances
    #  cols = ['RunType','Sol_time','Ls(deg)','InSituCH4 (ppbv)','Error(±1SEM)','GlobalPressureMultiplierFp','GlobalMeanCH4(ppbv)','Error(±1SEM)']
    ab = pd.read_csv('/project/gas_seepage/jportiz/mars/data/methane/literature_data/methane_data_webster2017.csv', sep=',')
    #Only want Enrichment method abundances
    em = ab[ab['RunType']=='Enrichment']
    fig = plt.figure(figsize=(8,6))
    #  fig = plt.figure(figsize=(12,6))
    ax1 = fig.add_subplot(111)
    axls = 24#28.     # axis label size
    ax2 = ax1.twinx()
    # METHANE ABUNDANCES (WEBSTER_2018a) 
    ax2.errorbar(em['Ls(deg)'], em['GlobalMeanCH4(ppbv)'], yerr=em['Error(±1SEM).1'],marker='o',ms=10,capsize=5,ls='',color='red',mec='k',label=r'Webster \textit{et al.} (2018)')
    ax2.set_ylabel(r'\textbf{CH$_4$ Abundance [ppbv]}',rotation=-90,va='bottom',fontsize=axls-6,color='red')
    ax2.yaxis.label.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='black'), path_effects.Normal()]) #outline ylabel text in black
    #  ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red')
    # MASS FLUX 
    #  # Interpolate time [sols] to L_s [deg]    
    #  f0 = interpolate.interp1d(sdf['SOL_TIME'], sdf['L_s'])
    #  lsi = f0(z_time) #interpolated L_s for all FEHM flux data
    #  max_flux = max(zdata.Sink[zstart_i:zend_i]*days2sols(1)/solsec()/surface_area)
    #  min_flux = min(zdata.Sink[zstart_i:zend_i]*days2sols(1)/solsec()/surface_area)
    #  ax1.plot(lsi[zstart_i:zend_i],zdata.Sink[zstart_i:zend_i]*days2sols(1)/solsec()/surface_area, label='sink',c='k',ls='-',zorder=1)
    # Interpolate time [sols] to L_s [deg]  #<-- may not need anymore 
    max_flux = max(full_df['FLUX'][zstart_i:zend_i])
    min_flux = min(full_df['FLUX'][zstart_i:zend_i])
    ax1.plot(full_df['L_s'][zstart_i:zend_i], full_df['FLUX'][zstart_i:zend_i], label='sink',c='k',ls='-',zorder=1)
    #  ax1.plot(lsi[zstart_i:zend_i],zdata.Sink[zstart_i:zend_i]*days2sols(1)/solsec()/surface_area, label='sink',c='k',ls='-',zorder=1)
    ax1.tick_params(axis="x", labelsize=axls-4)
    ax1.tick_params(axis="y", labelsize=axls-4)
    ax2.tick_params(axis="y", labelsize=axls-4)
    ax1.set_xlabel(r'\textbf{Solar Longitude, $L_s$ [$^\circ$]}',fontsize=axls-6)
    ax1.set_ylabel(r'\textbf{Surface CH$_4$ flux [kg/km$^2$/s]}', fontsize=axls-6)
    # SETTINGS
    ax1.xaxis.set_minor_locator(MultipleLocator(20)) #20-degree minor ticks 
    ax1.xaxis.set_major_locator(MultipleLocator(40)) #20-degree minor ticks 
    ax1.set_ylim(min_flux,ylim_flux)
    #  ax1.set_ylim(min_flux,max_flux*1.1)
    ax1.set_xlim(0, 360)
    ax2.set_xlim(0, 360)
    #  ax2.set_ylim(0.0, 0.8) #[ppbv]
    ax2.legend(loc='upper center', facecolor='white',edgecolor='k',fancybox=True,shadow=True,bbox_to_anchor=(0.5,1.12),fontsize=axls-10)
    # Plot the season ranges with colors zones 
    sdict, sc, season_df = synth_seasons()
    for season in sdict:
        ax1.axvspan(sdict[season][0],sdict[season][1], color=sc[season],alpha=0.75,zorder=0)
    # Add text saying "autumn", "spring"
    for season in sdict:
        mdp = (sdict[season][1]+sdict[season][0])/2 #midpoint of seasonal range
        ax1.text(mdp, ax1.get_ylim()[-1]*0.95, season, ha='center',fontsize=axls-4)

    fig.tight_layout()
    fig.savefig(join(output_dir,'pub_zoom_non-log_flux_wAbundancesWebster2018a_plot.pdf'))
    fig.savefig(join(output_dir,'pub_zoom_non-log_flux_wAbundancesWebster2018a_plot.png'),format='png',transparent=True, dpi=600) #PNG for transparent background... AGU poster
    plt.close(fig)
    plt.close('all')

#--------------------------------------------------------------------- 
# Run plot script either on its own or from within ``run.py``
if __name__ == "__main__":
    print('NOTE: likely better results plotting from run.py')
    #load certain vars from master_model.py run
    with open('v.pkl','rb') as f:
        #  sourcetype,fehm_upshift = pickle.load(f)
        sourcetype,fehm_upshift,T_upshift = pickle.load(f)
    #  T_upshift=95.
    make_plots(fehm_upshift=fehm_upshift, T_upshift=T_upshift)
    #  make_plots(sourcetype='unknown source type',fehm_upshift=0.0)

