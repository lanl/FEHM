import numpy
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.signal import blackman
#  import cPickle as pickle
import pickle 
from os.path import join
import os

def baro_fft(pas,ts,c,ax=None):
    '''
    Perform spectral analysis of barometric data.
    It is assumed that data is continuous with constant delta time
    between pressure observations.

    :param pas: barometric pressures
    :type pas: lst(float) or ndarray(float)
    :param ts: times
    :type ts: lst(float) or ndarray(float)
    :param ax: Axis for plotting
    :type ax: matplotlib axis
    '''

    # Create blackman window filter
    w = blackman(len(pas))
    # fft with blackman
    yf = fft(pas)
    # fft with blackman
    yfw = fft(pas*w)
    # Find delta times, assume all are constant
    T = numpy.diff(ts)[0]
    # Length of fft array
    N = yf.shape[0]
    # If number of pressures is odd
    if len(pas) % 2:
        xf = numpy.linspace(0.0, 1.0/(2.0*T), (N-1)//2)
        xf[1:] = 1/xf[1:]
        x0 = numpy.zeros_like(xf)
        #axarr2[fi].semilogx(xf[1:], x0[1:], (2.0/N * numpy.abs(yf[1:(N-1)/2])),'.')
        #axarr2[fi].vlines(xf[1:], x0[1:], (2.0/N * numpy.abs(yf[1:(N-1)/2])))
        if ax is not None:
            ax.vlines(xf[1:], x0[1:], (2.0/N * numpy.abs(yfw[1:(N-1)//2])),color=c)
        pf = numpy.column_stack([xf[1:],(2.0/N * numpy.abs(yfw[1:(N-1)//2]))])
    # Else, number of pressures is even
    else:
        xf = numpy.linspace(0.0, 1.0/(2.0*T), N//2)
        xf[1:] = 1/xf[1:]
        x0 = numpy.zeros_like(xf)
        #axarr2[fi].semilogx(xf[1:], x0[1:], (2.0/N * numpy.abs(yf[1:N/2])),'.')
        #axarr2[fi].vlines(xf[1:], x0[1:], (2.0/N * numpy.abs(yf[1:N/2])))
        if ax is not None:
            ax.vlines(xf[1:], x0[1:], (2.0/N * numpy.abs(yfw[1:N//2])),color=c)
        pf = numpy.column_stack([xf[1:], (2.0/N * numpy.abs(yfw[1:N//2]))])

    return pf

def fourier_series(pars,mean_pa,ts):
    y_synth = mean_pa
    ps,ampls,pshifts = numpy.reshape(pars,(3,-1))
    for p,a,phi in zip(ps,ampls,pshifts):
        w = 2*numpy.pi/p # Frequency
        y_synth -= a*(numpy.sin(w*ts+phi))
    return y_synth

def write_boun(out_file, t, y, step="ti_linear"):
    fout = open(out_file,'w')
    #fout.write('# Mean pressure: '+str(numpy.mean(y))+'\n')
    fout.write('boun\nmodel 1\n'+step+'\n')
    fout.write(str(len(t))+'\n')
    cnt = 0
    for i in range(len(t)):
        fout.write(' %12g' % t[i])
        cnt += 1
        if cnt == 5 or i == len(t)-1:
            fout.write('\n')
            cnt = 0

    fout.write(varname+'\n')
    cnt = 0
    for i in range(len(y)):
        fout.write(' %12g' % y[i])
        cnt += 1
        if cnt == 5 or i == len(y)-1:
            fout.write('\n')
            cnt = 0

    fout.write('\n')
    fout.write(' -6 0 0 1\n')
    fout.write(' -106 0 0 1\n')
    fout.write('\nstop\n')
    fout.close()

def write_pflotran_boun_tpl(out_file, ts, ys):
    with open(out_file,'w') as fout:
        fout.write('ptf #\n')
        fout.write('TIME_UNITS d\n')
        fout.write('DATA_UNITS Pa\n')
        for t,y in zip(ts,ys):
            fout.write(' %12g #%12g + pressure_offset#\n' % (t,y*1e6))

if __name__=='__main__':

    files = [
        #'inputdata/boun_pres',
        #'inputdata/boun_pres3',
        'inputdata/boun_pres_ancjan',
        #'inputdata/boun_pres_ancjun',
        #'inputdata/boun_pres_denvaltjan',
        #'inputdata/boun_pres_denvaltjun',
        'inputdata/boun_pres_denvjan',
        #'inputdata/boun_pres_denvjun',
        'inputdata/boun_pres_honjan',
        #'inputdata/boun_pres_honjun',
        #'inputdata/boun_pres_kdec',
        #'inputdata/boun_pres_kjan',
        #'inputdata/boun_pres_kjun',
        #'inputdata/boun_pres_kleen',
        #'inputdata/boun_pres_kmar',
        #'inputdata/boun_pres_ksep'
        ]

    f, axarr = plt.subplots(len(files), sharex=True, figsize=(25,16))
    f2, axarr2 = plt.subplots(len(files), sharex=True, figsize=(15,6))
    f3, axarr3 = plt.subplots(1, figsize=(15,4))
    f4, axarr4 = plt.subplots(1, figsize=(15,4))

    # Create dictionary of details
    baro_dict = {}

    cs = ['r','g','k','m','b']
    for fi,amys_file in enumerate(files):

        fh = open(amys_file,'r')
        nn = int(fh.readline().strip())
        ts = numpy.fromfile(fh,count=nn,sep='\n')
        varname = fh.readline().strip()
        pas = numpy.fromfile(fh,count=nn,sep='\n')
        fh.close()

        out_file = amys_file.split('/')[1]+'.dat'
        pout_file = amys_file.split('/')[1]+'.tpl'
        write_boun(out_file, ts, pas,step='ti_linear')
        write_pflotran_boun_tpl(pout_file, ts, pas)

        mean_pas = numpy.mean(pas)
        axarr[fi].plot(ts, pas, label='Measured')
        axarr[fi].axhline(y=mean_pas,c='k',linestyle='--',label='Mean = '+str(mean_pas))

        if fi == 0:
            ts_ind = numpy.where(ts<365)[0][-1]
            axarr4.plot(ts[0:ts_ind], pas[0:ts_ind],label='Measured',c='gray')
            axarr4.text(ts[ts_ind], pas[ts_ind],'Measured',color='gray')
            baro_dict[101] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 32.625, # Determined manually
                              'description': 'Anchorage measured',
                             }

            # Create base case, constant mean pressure
            ts_base = [0,1.e20]
            pas_base = [mean_pas,mean_pas]
            out_file = amys_file.split('/')[1]+'_base.dat'
            pout_file = amys_file.split('/')[1]+'_base.tpl'
            write_boun(out_file, ts_base, pas_base, step='ti_linear')
            write_pflotran_boun_tpl(pout_file, ts_base, pas_base)
            baro_dict[100] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 0,
                              'description': 'Anchorage synthetic constant pressure',
                             }

            periods = [7.305]
            ampls = [0.001]*1
            pshifts = [0.]*1
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.07305)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            ts_ind = numpy.where(ts<134)[0][-1]
            axarr3.plot(ts[0:ts_ind], y_synth[0:ts_ind], label='Weekly',c=cs[0])
            axarr3.text(ts[ts_ind],y_synth[ts_ind],'Weekly',color=cs[0])
            out_file = amys_file.split('/')[1]+'_weekly.dat'
            pout_file = amys_file.split('/')[1]+'_weekly.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[102] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 4*periods[0], # Multiply by 4 so that spinup is longer
                              'description': 'Anchorage synthetic weekly',
                             }
            periods = [182.625]
            ampls = [0.001]*1
            pshifts = [0.]*1
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.07305)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            ts_ind = numpy.where(ts<134)[0][-1]
            axarr3.plot(ts[0:ts_ind], y_synth[0:ts_ind], label='Bi-annual',c=cs[1])
            axarr3.text(ts[ts_ind],y_synth[ts_ind],'Bi-annual',color=cs[1])
            out_file = amys_file.split('/')[1]+'_biannual.dat'
            pout_file = amys_file.split('/')[1]+'_biannual.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[103] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0],
                              'description': 'Anchorage synthetic bi-annual',
                             }
            periods = [30.4375]
            ampls = [0.001]*1
            pshifts = [0.]*1
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.01521875)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            ts_ind = numpy.where(ts<134)[0][-1]
            axarr3.plot(ts[0:ts_ind], y_synth[0:ts_ind], label='Monthly',c=cs[2])
            axarr3.text(ts[ts_ind],y_synth[ts_ind],'Monthly',color=cs[2])
            out_file = amys_file.split('/')[1]+'_monthly.dat'
            pout_file = amys_file.split('/')[1]+'_monthly.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[104] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0],
                              'description': 'Anchorage synthetic monthly',
                             }
            periods = [1]
            ampls = [0.002*0.03]
            pshifts = [0.]
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.05)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            ts_ind = numpy.where(ts<134)[0][-1]
            axarr3.plot(ts[0:ts_ind], y_synth[0:ts_ind], label='Daily',c=cs[3])
            axarr3.text(ts[ts_ind],y_synth[ts_ind],'Daily',color=cs[3],va='top')
            out_file = amys_file.split('/')[1]+'_daily.dat'
            pout_file = amys_file.split('/')[1]+'_daily.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[105] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0]*30,
                              'description': 'Anchorage synthetic daily',
                             }
            periods = [0.5]
            ampls = [0.002*0.025]
            pshifts = [0.]
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.05)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            ts_ind = numpy.where(ts<134)[0][-1]
            axarr3.plot(ts[0:ts_ind], y_synth[0:ts_ind], label='Bi-Daily',c=cs[4])
            axarr3.text(ts[ts_ind],y_synth[ts_ind],'Bi-daily',color=cs[4])
            out_file = amys_file.split('/')[1]+'_bidaily.dat'
            pout_file = amys_file.split('/')[1]+'_bidaily.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[106] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0]*60,
                              'description': 'Anchorage synthetic bi-daily',
                             }
            periods = [182.625,30.4375,7.305,1.,0.5]
            ampls = [0.001]*3+[0.001*0.1]+[0.001*0.1]
            pshifts = [0.]*5
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.07305)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            ts_ind = numpy.where(ts<365)[0][-1]
            axarr[fi].plot(ts, y_synth, label='Combined')
            axarr4.plot(ts[0:ts_ind], y_synth[0:ts_ind], label='Combined',c='r')
            axarr4.text(ts[ts_ind], y_synth[ts_ind],'Synthetic',color='r')
            out_file = amys_file.split('/')[1]+'_combined.dat'
            pout_file = amys_file.split('/')[1]+'_combined.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[107] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0],
                              'description': 'Anchorage synthetic combined',
                             }
            periods = [7.305,1.,0.5]
            ampls = [0.001]*1+[0.001*0.03]+[0.001*0.025]
            pshifts = [0.]*3
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,periods[0]+0.001,0.07305)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            axarr[fi].plot(ts, y_synth, label='Combined 02')
            out_file = amys_file.split('/')[1]+'_combined02.dat'
            pout_file = amys_file.split('/')[1]+'_combined02.tpl'
            write_boun(out_file, ts_boun, y_boun,step='cy_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[108] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 4*periods[0],
                              'description': 'Anchorage synthetic combined (weekly, daily, bi-daily)',
                             }
            # Reset periods to total combined
            periods = [182.625,30.4375,7.305,1.,0.5]
        elif fi == 1:
            baro_dict[201] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 33.5, # Determined manually
                              'description': 'Denver measured',
                             }

            # Create base case, constant mean pressure
            ts_base = [0,1.e20]
            pas_base = [mean_pas,mean_pas]
            out_file = amys_file.split('/')[1]+'_base.dat'
            pout_file = amys_file.split('/')[1]+'_base.tpl'
            write_boun(out_file, ts_base, pas_base, step='ti_linear')
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[200] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 0,
                              'description': 'Denver constant pressure',
                             }

            periods = [365.25,22,5.6,1.,0.5]
            ampls = [0.0005,0.0005]+[0.0004]+[0.0001]*2
            pshifts = [0.]*5
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,730.5+0.001,0.07305)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            axarr[fi].plot(ts, y_synth, label='Synthetic')
            out_file = amys_file.split('/')[1]+'_synth.dat'
            pout_file = amys_file.split('/')[1]+'_synth.tpl'
            write_boun(out_file, ts_boun, y_boun)
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[202] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0],
                              'description': 'Denver synthetic combined',
                             }
        elif fi == 2:
            baro_dict[301] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 30.75, # Determined manually
                              'description': 'Honolulu measured',
                             }

            # Create base case, constant mean pressure
            ts_base = [0,1.e20]
            pas_base = [mean_pas,mean_pas]
            out_file = amys_file.split('/')[1]+'_base.dat'
            pout_file = amys_file.split('/')[1]+'_base.tpl'
            write_boun(out_file, ts_base, pas_base, step='ti_linear')
            write_pflotran_boun_tpl(pout_file, ts_base, pas_base)
            baro_dict[300] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': 0,
                              'description': 'Honolulu constant pressure',
                             }

            periods = [365.25,22,8,1.,0.5,0.33]
            ampls = [0.00015]*3+[0.00005]*3
            pshifts = [-365.25/4,0,0,0,0,0]
            y_synth = fourier_series(periods+ampls+pshifts,mean_pas,ts)
            ts_boun = numpy.arange(0.,730.5+0.001,0.07305)
            y_boun = fourier_series(periods+ampls+pshifts,mean_pas,ts_boun)
            axarr[fi].plot(ts, y_synth, label='Synthetic')
            out_file = amys_file.split('/')[1]+'_synth.dat'
            pout_file = amys_file.split('/')[1]+'_synth.tpl'
            write_boun(out_file, ts_boun, y_boun)
            write_pflotran_boun_tpl(pout_file, ts_boun, y_boun)
            baro_dict[302] = {'file': join(os.getcwd(),out_file),
                              'p_mean': mean_pas,
                              'spinup_period': periods[0],
                              'description': 'Honolulu synthetic combined',
                             }


        axarr[fi].set_title(amys_file)
        axarr[fi].set_xlabel("Time [d]")
        axarr[fi].set_ylabel("Pressure [MPa]")
        axarr[fi].legend()

        # Perform fft and plot amplitude vs period
        pf = baro_fft(pas,ts,cs[fi],axarr2[fi])

        # Add synthetic periods as dashed vertical lines to plot
        for p in periods:
            axarr2[fi].axvline(x=p,linestyle='--')
        #axarr2[fi].set_title(amys_file)
        axarr2[fi].set_ylim(6e-7,100)
        axarr2[fi].set_xscale("log")
        axarr2[fi].set_yscale("log")

        # Dump period, amplitude to file
        numpy.savetxt(amys_file.split('/')[1]+'.fft',pf,header='period_[days] amplitude_[MPa]')

    axarr2[2].set_xlabel("Period [d]")
    axarr2[1].set_ylabel("Amplitude [MPa]")
    pickle.dump(baro_dict,open('baro_dict.pkl','wb'))
    f2.savefig('fft.png')


    axarr3.set_xlabel("Time [d]")
    axarr3.set_ylabel("Pressure [MPa]")
    axarr3.set_xlim(0,142)
    f3.tight_layout()
    f3.savefig('Anchorage_pressure_components.png')

    axarr4.set_xlabel("Time [d]")
    axarr4.set_ylabel("Pressure [MPa]")
    axarr4.set_xlim(0,390)
    f4.tight_layout()
    f4.savefig('Anchorage_meas_combined.png')
    plt.show(block=True)


