'''
Collection of useful tools for use in Mars analysis.
'''
import os,sys
import numpy
import numpy as np
sys.path.append('/project/gas_seepage/jportiz/scripts')
from gasmigra import signalTools
from scipy import ndimage
import pandas as pd
from scipy.interpolate import interp1d

def g_earth():
    '''Return Earth's gravitational acceleration. [m^2/s]'''
    return 9.807

def g_mars():
    '''Return Mars's gravitational acceleration. [m^2/s]'''
    return 3.711

def solsec(num_sols=1):
    '''
    Return number of seconds in a martian sol.
    Optional argument can provide the number of sols (default=1).
    '''
    return 88775.24409*num_sols

def days2sols(days):
    '''Convert Earth days to martian sols.'''
    return days*86400./solsec()

def sols2days(sols):
    '''Convert martian sols to Earth days.'''
    return sols*solsec()/86400.

def mars_radius():
    '''Return radius in meters.'''
    return 3389439.4

def molar_mass_ch4():
    '''Return g/mol CH_4'''
    return 16.04

'''
Mars altitude-air density lookup table.
Viking 1 Lander site, daytime temps
http://www.braeunig.us/space/atmmars.htm
'''
table = numpy.array([
    [-7000,249,1.31e+3,2.74e-2],
    [-6000,248,1.20e+3,2.52e-2],
    [-5000,247,1.10e+3,2.31e-2],
    [-4000,246,1.00e+3,2.12e-2],
    [-3000,245,9.16e+2,1.94e-2],
    [-2000,244,8.37e+2,1.78e-2],
    [-1000,243,7.65e+2,1.64e-2],
    [0,242,6.99e+2,1.50e-2],
    [1000,241,6.39e+2,1.38e-2],
    [2000,240,5.84e+2,1.27e-2],
    [3000,239,5.34e+2,1.16e-2],
    [4000,238,4.88e+2,1.07e-2],
    [5000,237,4.46e+2,9.79e-3],
    [6000,236,4.07e+2,8.98e-3],
    [7000,234,3.72e+2,8.26e-3],
    [8000,232,3.40e+2,7.64e-3],
    [9000,230,3.11e+2,7.05e-3],
    [10000,228,2.84e+2,6.50e-3],
    [15000,216,1.81e+2,4.36e-3],
    [20000,205,1.16e+2,2.93e-3],
    [25000,194,7.37e+1,1.97e-3],
    [30000,183,4.70e+1,1.34e-3],
    [40000,161,1.91e+1,6.18e-4]
                 ])

def density_lookup(altitude):
    '''
    Interpolate approximate atmospheric density given the altitude/elevation.
    (e.g. Curiosity Rover landing site = -4500 m)
    '''
    altitudes=table[:,0] #[m]
    rhos = table[:,-1]   #[kg/m3]
    return numpy.interp(altitude,altitudes,rhos)

def mars_seasons(L_s):
    '''
    Northern hemisphere seasons by solar longitude (L_s).
    http://www-mars.lmd.jussieu.fr/mars/time/solar_longitude.html
    https://www.planetary.org/articles/mars-calendar
    '''
    #  summer = 45.  <= L_s <= 135. #these are wrong
    #  autumn = 135. <= L_s <= 225.
    #  winter = 225. <= L_s <= 315.
    #  if not any([summer,autumn,winter]): spring = True
    #  spring = 315. <= L_s <= 360. and <=45.
    spring = 0.   <= L_s <= 90.
    summer = 90.  <= L_s <= 180.
    autumn = 180. <= L_s <= 270.
    winter = 270. <= L_s <= 360.
    #  if not any([spring,summer,autumn]): winter = True
    # Return the season (1=Summer, 2=Autumn, 3=Winter, 4=Spring)
    if spring   is True: s = 'N. Spring'
    elif summer   is True: s = 'N. Summer'
    elif autumn is True: s = 'N. Autumn'
    elif winter is True: s = 'N. Winter'

    return s

def pbl_heights(lookup_Ls='all',onlyMinima=True,timeType='LMST'):
    '''
    Planetary Boundary Layer (PBL) heights broken up into Ls 30° bins.
    Can specify either the daily minimum heights for each 30° bin, or
    1-min (LTST) heights for each 30° bin.

    Reference: Newman2017
               Used by Moores2019a (daily minimums)
               Used by Guzewich20 (hourly

    Parameters
    ----------
    lookup_Ls : float or ndarray (if onlyMinima=False)
        if onlyMinima=True:
            float, L_s [°] to lookup the PBL height for.
        if onlyMinima=False:
            ndarray, of form [ L_s, decimal_time ]
            Decimal time can be calculated like:  df['SOL_TIME'] % 1
        if empty, function returns all the pbl_heights as an array.
    onlyMinima : bool
        True (default) gives the minimum PBL heights in 30° bins (Moores2019a).
        False gives time-dependent PBL heights in 30° bins (Guzewich2017).
    timeType : str
        Flag saying which local time type to use ('LMST' (default) or 'LTST').
    '''

    if onlyMinima==True:
        # Daily Minimums 
        #-------------------------
        #PBL heights, based on Newman et al 2017, Grid B, daily Minimums (Ls, PBLh in meters)
        pbl = [ [0, 248.3385],
                [30.0000, 207.2998],
                [60.0000, 216.8055],
                [90.0000, 232.3460],
                [120.0000, 227.8850],
                [150.0000, 246.0847],
                [180.0000, 245.5930],
                [210.0000, 247.8800],
                [240.0000, 256.5392],
                [270.0000, 254.9117],
                [300.0000, 256.1173],
                [330.0000, 250.3968],
                [360.0000, 248.3385],
              ]
        if lookup_Ls=='all': pass
        else:
            # Lookup a single PBL height at the L_s provided
            #  for row in range(len(pbl)):
            for row in enumerate(pbl):
                ls = row[1][0]; i = row[0]
                # None of my data have L_s == 360.0, so this shouldn't break
                if (lookup_Ls >= ls) and (lookup_Ls<pbl[i+1][0]):
                    pbl_individual = pbl[i][1]
                    break
            pbl=pbl_individual

    elif onlyMinima==False:
        # 1-min PBL Heights 
        #-------------------------
        #PBL heights after Guzewich2017
        path = '/project/gas_seepage/jportiz/scripts/mars/pbl_heights/Guzewich2017'
        #  datafile = 'pbl_heights_Guzewich2017.csv'
        datafile = 'pbl_heights_Guzewich2017_processed.csv'
        df = pd.read_csv(os.path.join(path,datafile))
        # Return all PBL heights in data
        if lookup_Ls=='all': pbl=df
        # Return single PBL height at a specific L_s and decimal time
        else:
            #  print('Can only lookup single PBL height if ``onlyMinima=True``.')
            if len(lookup_Ls)==1:
                print('If ``onlyMinima=True``, then ``lookup_Ls`` must be a list or tuple of [L_s, decimal_time].')
            else:
                # Lookup the PBL heights per time at the L_s provided 
                all_Ls = df['L_s'].unique()
                # Get Ls bins 
                Ls_bin = all_Ls[ all_Ls.searchsorted(lookup_Ls[0]) - 1 ]
                #  Ls_bin_next = Ls_bin + 30
                # Get all the PBL heights per time for desired bin 
                pbl_times = df[ df['L_s']==Ls_bin ]
                # Interpolate the PBL based on provided time 
                if timeType=='LMST':
                    f0 = interp1d(pbl_times['LMST'], pbl_times['PBL_m'])
                elif timeType == 'LTST':
                    f0 = interp1d(pbl_times['LTST'], pbl_times['PBL_m'])
                else: print('timeType must be either LMST or LTST.')
                pbl_interp = f0(lookup_Ls[1]) #[m]
                pbl = pbl_interp
    else: print('onlyMinima must be a boolean.')

    return pbl

def findPrincipleComponents(press,time,comps = ['s','d','e','y']):
#  def findPrincipleComponents(press,time,comps = ['s','d','w','m','y']):
    '''
    decomposes a pressure signal into subdiurnal, diurnal, seasonal, and Mars-yearly components.

    :param press: barometric pressures
    :type press: lst(float) or ndarray(float)
    :param time: times (sols)
    :type time: lst(float) or ndarray(float)
    :returns: lst(float) -- a list of signal periods (sols)
    :returns:lst(float) -- a list of signal amplitudes
    :returns: float -- mean pressure
    '''

    #get frequency spectrum    
    #  periodFFT,magFFT, phaseFFT = signalTools.simpleFFT(press,time)
    periodFFT,magFFT, phaseFFT = signalTools.simpleFFT(press,time, useFilter=True)

    #mean pressure (DCoffset)
    avgMag = numpy.mean(press)  

    #list of compnent periods and amps to build
    amp = list()
    period = list()    
    phase_shifts = list()

    #semidiurnal
    if 's' in comps:
        #  i,j = indexBetween(periodFFT,0.4,0.6)
        i,j = signalTools.indexBetween(periodFFT,0.4,0.6)
        period.append(0.5)
        amp.append(sum(magFFT[i:j]))   
        phase_shifts.append(phaseFFT[i]) #JPO TESTING

    #diurnal
    if 'd' in comps:
        i,j = signalTools.indexBetween(periodFFT,0.8,1.2)
        period.append(1)
        amp.append(sum(magFFT[i:j]))
        phase_shifts.append(phaseFFT[i]) #JPO TESTING

    #  #Weather Fronts
    #  if 'w' in comps:    
        #  i,j = indexBetween(periodFFT,3,10)
        #  period.append(periodFFT[int(ndimage.measurements.center_of_mass(np.array(magFFT[i:j]))[0])+i])
        #  amp.append(sum(magFFT[i:j]))
#  
    #  #Weather Fronts 2
    #  if 'm' in comps:     
        #  i,j = indexBetween(periodFFT,10,30)
        #  period.append(periodFFT[int(ndimage.measurements.center_of_mass(np.array(magFFT[i:j]))[0])+i])
        #  amp.append(sum(magFFT[i:j])*2)

    #Seasonal
    if 'e' in comps:
        #  i,j = signalTools.indexBetween(periodFFT,220,550)
        i,j = signalTools.indexBetween(periodFFT,330,342)
        period.append(periodFFT[int(ndimage.measurements.center_of_mass(np.array(magFFT[i:j]))[0])+i])
        amp.append(sum(magFFT[i:j])*2)
        #  phase_shifts.append(phaseFFT[i]) #JPO TESTING (<--- THIS IS WRONG!!!)
        phase_shifts.append(phaseFFT[i]) #JPO TESTING (<--- THIS IS WRONG!!!)

    #Year
    if 'y' in comps:
        i,j = signalTools.indexBetween(periodFFT,660,680)
        period.append(periodFFT[int(ndimage.measurements.center_of_mass(np.array(magFFT[i:j]))[0])+i])
        amp.append(sum(magFFT[i:j]))
        phase_shifts.append(phaseFFT[i]) #JPO TESTING

    #Convert reletive mag to absolute
    a_rms = np.sqrt(np.mean(np.square(press-avgMag)))
    a_tot = a_rms*np.sqrt(2)
    norm = sum(amp)
    for i,a in enumerate(amp):
        amp[i] = a/norm*a_tot

    return period,amp,phase_shifts,avgMag

def createSyntheticSignal(period,amp,mean,time):
    '''
    combines pressure signal components into a single synthetic signal

    :param period: list of periods for signals
    :type period: lst(float) or ndarray(float)
    :param amp: lsit of amplitudes for signals
    :type amp: lst(float) or ndarray(float)
    :param mean: the mean pressure (acts as DC offset)
    :type mean: float
    :param time: list of times 
    :type time: lst(float) or ndarray(float)
    :returns: lst(float) -- the pressure of the synthetic signal at passed times
    '''

    synthSignal = np.full_like(time,mean)

    if type(period) is not list:
        period = [period]
        amp = [amp]

    for p,a in zip(period,amp):
        for i,t in enumerate(time):
            synthSignal[i] += a*np.sin(6.283*t/p)
    return synthSignal

#  def createSyntheticSignal(period,amp,mean,time):
    #  '''
    #  combines pressure signal components into a single synthetic signal.
    #  (Add ability to specify phase of each component.)
#  
    #  :param period: list of periods for signals
    #  :type period: lst(float) or ndarray(float)
    #  :param amp: lsit of amplitudes for signals
    #  :type amp: lst(float) or ndarray(float)
    #  :param mean: the mean pressure (acts as DC offset)
    #  :type mean: float
    #  :param time: list of times
    #  :type time: lst(float) or ndarray(float)
    #  :returns: lst(float) -- the pressure of the synthetic signal at passed times
    #  '''
#  
    #  synthSignal = np.full_like(time,mean)
#  
    #  if type(period) is not list:
        #  period = [period]
        #  amp = [amp]
#  
    #  #  for p,a in zip(period,amp):
        #  #  for i,t in enumerate(time):
            #  #  synthSignal[i] += a*np.sin(6.283*t/p)
    #  for p,a,phi in zip(period,amp,pshifts):
        #  for i,t in enumerate(time):
            #  synthSignal[i] += a*np.sin(6.283*t/p)
    #  return synthSignal

def synth_fourier_series(pars, mean_pa, ts, num_modulations=0):
    '''
    [[[[ ONE OR TWO MODULATIONS ]]]]
    Attempt to generate a Fourier series with a (2) seasonally modulated diurnal
    amplitude.
    Detects the diurnal component, and modifies it with one of the other
    components (use the last and 2nd-to-last  period/amplitude/phase shift in
    pars because all the ones at beginning are fixed values (either from prior
    calibration, or fixed as 0.25-, 0.5-, and 1.0-sol period components).

    Same args as convert_boun.fourier_series().
    '''
    y_synth = mean_pa
    ps,ampls,pshifts = numpy.reshape(pars,(3,-1))
    #  lastEntry = ps[-1] #get the last period in pars
    if   num_modulations == 0: lastEntry = -999   #bogus value
    elif num_modulations == 1: lastEntry = ps[-1] #get the last period in pars
    elif num_modulations == 2: lastEntry = ps[-2] #get the 2nd-to-last period in pars 
    else: print('ERROR: currently only supports 1 or 2 diurnal modulations.')
    counter = -1
    for p,a,phi in zip(ps,ampls,pshifts):
        counter+=1
        w = 2*numpy.pi/p # Frequency
        # Modify diurnal signal with the last component in pars
        if p==1.0:
            # Save the diurnal parameters 
            w_d  = 2*numpy.pi/p
            a_d  = a
            ph_d = phi
            continue  #skip adding diurnal until the end
            #  pass  #skip adding diurnal until the end
        else:
            #  if p == lastEntry:
            #  if (p == lastEntry) & (1 in ps):
            if (counter == len(ps)-1) & (1 in ps):
                # 1 modulation of diurnal component
                if   num_modulations == 1:
                    # Get the seasonal parameters to modulate with
                    w_se  = 2*numpy.pi/p
                    a_se  = a
                    ph_se = phi
                    #-- Add signal
                    y_synth -= ( a_d + a_se* np.sin(w_se*ts + ph_se) ) * np.sin(w_d*ts + ph_d)
                # 2 modulations of diurnal component
                elif num_modulations ==2:
                    # Get the seasonal parameters to modulate with
                    #-- 1st seasonal modulation
                    w_se1  = 2*numpy.pi/p
                    a_se1  = a
                    ph_se1 = phi
                    #-- 2nd seasonal modulation
                    w_se2  = 2*numpy.pi/ps[-1]
                    a_se2  = ampls[-1]
                    ph_se2 = pshifts[-1]
                    #-- Add signal
                    y_synth -= ( a_d + a_se1* np.sin(w_se1*ts + ph_se1) + a_se2* np.sin(w_se2*ts + ph_se2) ) * np.sin(w_d*ts + ph_d)
                break
            else:
                y_synth -= a*(numpy.sin(w*ts+phi))
    return y_synth



#testing
outfile = 'eos_mars_co2_test.dat'
EW1 = 0.06
EW2 = 0.05
EW3 = 0.011637
EW4 = -4.26e-5
EW5 = 19.379
EW6 = 0.88
EW7 = 0.
EW8 = 4.2e-3
EW9 = 1.18278e-4
EW10 = 0.
EW11 = 2.5418e-8
####
EV1 = 0.06
EV2 = 0.05
EV3 = 1.00
EV4 = 0
EV5 = 0
EV6 = -0.020773
EV7 = 0.827
EV8 = -11.398
EV9 = 1.372e-4
EV10 = 0.
EV11 = 0.

def write_eos(outfile,EW1,EW2,EW3,EW4,EW5,EW6,EW7,EW8,EW9,EW10,EW11,EV1,EV2,EV3,EV4,EV5,EV6,EV7,EV8,EV9,EV10,EV11):
    '''
    Write an FEHM ``eos`` macro file.
    For now, assumes you will be doing a single-phase "water"
    simulation.

    Parameters
    ----------
    outfile : str
        Output file name.
    EW1 : float
        Liquid reference pressure [MPa].
    EW2 : float
        Liquid reference temperature [°C].
    EW3 : float
        Liquid reference density [kg/m^3].
    EW4 : float
        Derivative of liquid density with respect to pressure at reference
        conditions.
    EW5 : float
        Derivative of liquid density with respect to temperature at reference
        conditions.
    EW6 : float
        Liquid reference enthalpy [MJ/kg].
    EW7 : float
        Derivative of liquid enthalpy with respect to pressure at reference
        conditions.
    EW8 : float
        Derivative of liquid enthalpy with respect to temperature at reference
        conditions.
    EW9 : float
        Liquid reference viscosity [Pa s].
    EW10 : float
        Derivative of liquid viscosity with respect to pressure at reference
        conditions.
    EW11 : float
        Derivative of liquid viscosity with respect to temperature at reference
        conditions.
    EV1-11 : float
        Same as EW1-11, but for Vapor.
    '''

    f = open(outfile, 'w')
    f.write('eos\n')
    # GROUP 1
    f.write('3 0 0\n')
    # GROUP 2 (LIQUID)
    f.write('{} '.format(EW1).ljust(10," "))
    f.write('{} '.format(EW2).ljust(10," "))
    f.write('{} '.format(EW3).ljust(10," "))
    f.write('{} '.format(EW4).ljust(10," "))
    f.write('{} '.format(EW5).ljust(10," "))
    f.write('{} '.format(EW6).ljust(10," "))
    f.write('{} '.format(EW7).ljust(10," "))
    f.write('{} '.format(EW8).ljust(10," "))
    f.write('{} '.format(EW9).ljust(10," "))
    f.write('{} '.format(EW10).ljust(10," "))
    f.write('{} '.format(EW11).ljust(10," "))
    f.write('\n')
    # GROUP 3 (VAPOR)
    f.write('{} '.format(EV1).ljust(10," "))
    f.write('{} '.format(EV2).ljust(10," "))
    f.write('{} '.format(EV3).ljust(10," "))
    f.write('{} '.format(EV4).ljust(10," "))
    f.write('{} '.format(EV5).ljust(10," "))
    f.write('{} '.format(EV6).ljust(10," "))
    f.write('{} '.format(EV7).ljust(10," "))
    f.write('{} '.format(EV8).ljust(10," "))
    f.write('{} '.format(EV9).ljust(10," "))
    f.write('{} '.format(EV10).ljust(10," "))
    f.write('{} '.format(EV11).ljust(10," "))
    f.write('\nstop')
    f.close()
