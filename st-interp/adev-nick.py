#############################################################################
######This is the adev code:
######Adev.py
#############################################################################
 
import numpy as np
from decimal import Decimal as d
import math
from scipy import stats
from scipy.optimize import curve_fit
 
 
def prepareData(freq, t, tFormat, tau, dx, theEnd):
    # Remove NaN values from data
    if type(freq[0]) == d:
        t = [d(i) for i in t]
        arr = [ [i,j] for i,j in zip(t,freq) if not (i.is_nan() or j.is_nan()) ]
        t = [ i[0] for i in arr ]
        freq = [ i[1] for i in arr ]
    else:
        arr = [ [i,j] for i,j in zip(t,freq) if not (math.isnan(i) or math.isnan(j)) ]
        t = [ i[0] for i in arr ]
        freq = [ i[1] for i in arr ]
       
    
    if tFormat == 'mjd':
        t = [float(i - t[0])*24*3600 for i in t]
    else:
        t = [float(i - t[0]) for i in t]
   
    
    # find mean time different between data points (dt)
    dtArray = np.diff(t)
    if len(set(dtArray)) != 1:
        dtArray = dtArray[abs(dtArray - np.mean(dtArray)) < 2 * np.std(dtArray)]
    if type(freq[0]) == d:
        dt = d(np.mean(dtArray))
    else:
        dt = np.mean(dtArray)
    
    # remove frequency mean so it's just frequency differences
    m = np.mean(freq)
    freq = [f - m for f in freq]
   
    # integrate frequency data to get phase
    # (it's faster to calculate adev with phase data)
    s = [freq[0]]
    phase = [s.append(s[0]+n) or s.pop(0) for n in freq]
    phase = [p*dt for p in phase]
   
    M = len(freq) - 1
    end = (M+1)/theEnd
   
    # specify averaging time values
    if tau == "full":
        # full range
        R = [i for i in range(1,int(end)) if i < end]
    elif tau == "exp2":
        # exponential 2 range
        x = range(1,1000)
        s = [1]
        R = [s.append(s[0]*2) or s.pop(0) for n in x if s[0] <= end]
    elif tau == "exp10":
    # exponential 10 range
        R = [i*10**exp for exp in range(0, 10) for i in range(1, 10) if i*10**exp <= end]
    elif tau == 'other':
        # dx = 1 to 2, for more points but faster than tau="full"
        x = range(1,int(end))
        s = [1]
        R = [s.append(s[0]*dx) or s.pop(0) for n in x if s[0] <= end]
        R = [int(i) for i in R]
       
    
    return phase, M, R, dt
 
 
### mandatory arguments ###
# freq: list of frequencies
# t: list of times corresponding to frequencies
 
### optional arguments ###
# overlapping: (True or False), specifies whether to use overlapping adev or not.
#              Overlapping adev is smoother than normal adev.
# tau: ("full", "exp2", "exp10" or "other")
# dx: if tau="other", this specifies density of adev points
 
### output ###
# returns [t_adev, x_adev, error]
# t_adev: averages times, specified by "tau" input parameter
# x_adev: calculated adev for each t_adev
# error: 1-sigma error for each x_adev
 
 
def adev(freq, t, overlapping = False, tFormat="mjd", tau = "full", dx=2):
   
    phase, M, R, dt = prepareData(freq, t, tFormat, tau, dx, 2) # use 1/2 of the dataset
#     phase, M, R, dt = prepareData(freq, t, tFormat, tau, dx, 3) # use 1/3 of the dataset
   
    
    x = []
    error = []
   
    #########################################
    # this is where the adev is calculated ##
    #########################################
   
    # overlapping allan deviation (smoother than normal adev)
    if overlapping == True:
        for n in R:
            # three terms in summation
            x1 = phase[2*n:M+1]
            x2 = phase[n:M-n+1]
            x3 = phase[0:M-2*n+1]
           
            # calculate adev from summation terms
            sigma = np.sqrt( sum((a - 2*b + c)**2 for a,b,c in zip(x1,x2,x3)) \
                           / (2*(M + 1 - 2*n)*(n*dt)**2) )
           
            # append to list of adevs
            x.append( sigma )
           
            # add error bar - 1sigma confidence intervals (assuming white frequency noise data)
           
            # degrees of freedom for white FM
            df = (3*(M-1)/(2*n) - 2*(M-2)/M)* (4*n**2)/(4*n**2 + 5)
            # degrees of freedom for white PM
#             df = (M+1)*(M-2*n)/(2*(M-n))
               
            if type(phase[0]) == d:
                error.append( d(0.87) * sigma/d((M/n)**0.5) )
            else:
                errorHigh = sigma * (( df/stats.chi2.ppf(0.16 , df) )**0.5 - 1)
                errorLow = sigma * (1 - ( df/stats.chi2.ppf(0.84 , df) )**0.5 )
                error.append([errorLow, errorHigh])
#                 error.append( 0.87 * sigma/(M/n)**0.5 )
    
    # normal allan deviation
    else:
        for n in R:
            # take every nth phase data point
            ph = phase[::n]
            L = len(ph)
           
            # three terms in summation
            x1 = ph[2:L]
            x2 = ph[1:L-1]
            x3 = ph[0:L-2]
           
            # calculate adev from summation terms
            sigma = np.sqrt( sum((a - 2*b + c)**2 for a,b,c in zip(x1,x2,x3)) \
                            / (2*(L -2)*(n*dt)**2) )
           
            # append to list of adevs
            x.append( sigma )
           
            # add error bar (assuming white noise data)
            if type(phase[0]) == d:
                error.append( d(0.87) * sigma/d(L**0.5) )
            else:
                error.append( 0.87 * sigma/L**0.5 )
 
       
    
    return [float(i*dt) for i in R], x, np.array(error).T
 
 
 
def modAdev(freq, t, overlapping = False, tFormat="mjd", tau = "full", dx=2):
   
    phase, M, R, dt = prepareData(freq, t, tFormat, tau, dx, 3)
   
    
    x = []
    error = []
   
    #########################################
    # this is where the adev is calculated ##
    #########################################
   
    for n in R:
        # three terms in summation
        x1 = phase[2*n : 3*n]
        x2 = phase[n : 2*n]
        x3 = phase[0 : n]
       
        # calculate adev from summation terms
        a1 =  sum((a - 2*b + c) for a,b,c in zip(x1,x2,x3)) 
        
        a = []
        a.append(a1)
        for i in range(0, M-3*n+1):
            delta = phase[3*n+i] - 3*phase[2*n+i] + 3*phase[n+i] - phase[i]
            a.append( a[-1] + delta )
 
        sigma = sum( [j**2 for j in a] )
        x.append( np.sqrt( sigma / (2*(M + 2 - 3*n)* n**4 * dt**2) ) )
       
        
        
        # add error bar (assuming white noise data)
        if type(phase[0]) == d:
            error.append( d(0.87) * sigma/d((M/n)**0.5) )
        else:
            error.append( 0.87 * sigma/(M/n)**0.5 )
    
        
    
    return float([i*dt for i in R]), x, error
 
 
 
 
 
 
 
def getWhiteNoiseFit(t, adev0, weights, tCutoff=10):
    # Define the nonlinear function
    def whitenoise_func(x, a):
        return a/np.sqrt(x)
   
    tfit = [i for i in t if i > tCutoff]
    adevfit = [j for i,j in zip(t,adev0) if i > tCutoff]
    weightsfit = [j for i,j in zip(t,weights) if i > tCutoff]
   
    # Unweighted fit the nonlinear model
    # popt, pcov = curve_fit(whitenoise_func, tfit, adevfit)
   
    # Weighted fit the nonlinear model
    popt, pcov = curve_fit(whitenoise_func, tfit, adevfit, sigma=weightsfit, absolute_sigma=True)
   
    whiteNoiseFit = whitenoise_func(t, *popt)
   
    return whiteNoiseFit


####################################################
###### This is how I would call the adev function:
####################################################
 
# get the allan deviation points
tAlYb,adevAlYb,AlYberror = adev.adev(AlYbOffset, comb.mjd, overlapping=True, tau='other',dx=1.5)
# fit a white noise line excluding the points below 100 seconds
whiteNoiseFitAlYb = adev.getWhiteNoiseFit(tAlYb, adevAlYb, AlYberror[1], tCutoff=100)
 
 