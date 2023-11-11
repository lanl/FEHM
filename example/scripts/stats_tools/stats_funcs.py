'''
Collection of functions for statistical analysis, goodness-of-fit, etc.
'''
import numpy
import numpy as np


def calc_rmse(obs,pred):
    '''
    Calculate the Root Mean Squared Error.
    '''
    rmse = np.sqrt( np.sum( (obs-pred)**2 ) / len(obs) )
    return rmse

def chisqg(ydata,ymod,sd=None):
    '''
    Returns the chi-square error statistic as the sum of squared errors between
    Ydata(i) and Ymodel(i). If individual standard deviations (array sd) are supplied,
    then the chi-square error statistic is computed as the sum of squared errors
    divided by the standard deviations.     Inspired on the IDL procedure linfit.pro.
    See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.

    x,y,sd assumed to be Numpy arrays. a,b scalars.
    Returns the float chisq with the chi-square statistic.

    Rodrigo Nemmen
    http://goo.gl/8S1Oo
    '''
    # Chi-square statistic (Bevington, eq. 6.9)  
    #  if sd==None:
    if sd is None:
         chisq=numpy.sum((ydata-ymod)**2)
    else:
         chisq=numpy.sum( ((ydata-ymod)/sd)**2 )
    return chisq

def redchisqg(ydata,ymod,deg=2,sd=None):
    '''
    Returns the reduced chi-square error statistic for an arbitrary model,   
    chisq/nu, where nu is the number of degrees of freedom. If individual   
    standard deviations (array sd) are supplied, then the chi-square error   
    statistic is computed as the sum of squared errors divided by the standard
    deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.

    ydata,ymod,sd assumed to be Numpy arrays. deg integer.

    Usage:
    >>> chisq=redchisqg(ydata,ymod,n,sd)
    where  
     ydata : data  
     ymod : model evaluated at the same x points as ydata  
     n : number of free parameters in the model  
     sd : uncertainties in ydata  

    Rodrigo Nemmen
    http://goo.gl/8S1Oo
    '''
    # Chi-square statistic  
    #  if sd==None:
    if sd is None:
        chisq=numpy.sum((ydata-ymod)**2)
    else:
        chisq=numpy.sum( ((ydata-ymod)/sd)**2 )
    # Number of degrees of freedom assuming 2 free parameters  
    nu=ydata.size-1-deg
    return chisq/nu

