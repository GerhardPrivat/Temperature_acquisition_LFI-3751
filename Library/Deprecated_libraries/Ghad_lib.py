#!/usr/bin/env mpython_q
import numpy as np
import sys
from numpy import NaN, Inf, arange, isscalar, array, asarray
import ntpath

def getdiskparam(disknum=111):
    """
        dir=0 "ordinary" - big radius, dir=1 "extraordinary" - small radius
    """
    if disknum==103:
        from WGM_lib import nLNO1
        R0_mm = 1.6
        r0_mm = 1./6.
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==111:
        from WGM_lib import nLNO1
        R0_mm = 1.0
        r0_mm = 1./6.
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==106:
        from WGM_lib import nLNO1
        R0_mm = 1.27
#        R0_mm = 0.1626
        r0_mm = R0_mm
#        r0_mm = 0.97
        MgO = 5.4 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]

    if disknum==666:
        from WGM_lib import nLNO1
        R0_mm = 1.0
        r0_mm = R0_mm/6.
        r0_mm
        MgO = 5.7 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==966:
        from WGM_lib import nLNO1
#for JmO long paper
        R0_mm = 1.5
        r0_mm = .6
        MgO = 5.0
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==9191:
        from WGM_lib import nLNO1
#for JmO long paper
        R0_mm = 2.4978
        r0_mm = .581000
        r0_mm = .65
        MgO = 5.0 #Mg conc
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
        
    if disknum==919:
        from WGM_lib import nLNO1
#        R0_mm = 2.4601  MgO = 5.3 #Mg conc
#MgO = 5.4 #Mg conc  R0_mm = 2.500
#        R0_mm = 2.499
#        R0_mm = 2.489
#        R0_mm = 2.486

        R0_mm = 2.4978
#        R0_mm = 2.
#        r0_mm = .581000
        r0_mm = .65
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]

    if disknum==966:
        from WGM_lib import nLNO1
        R0_mm = 2.5
        r0_mm = .5
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]

    if disknum==977:
#        xcut disk
        from WGM_lib import nLNO1
        R0_mm = 2.5
        r0_mm = .5
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
#        tilted extension coefficients
        lamdir0 = [stretchx* (13.31  *10**(-6)+3.828 *10**(-6))/2., stretchx2*(0.03059*10**(-6)+stretchx2*0.00814*10**(-6))] 
        lamdir1 = [stretchx*  13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        
        
    if disknum==9193:
        from WGM_lib import nLNO1

        R0_mm = 2.4978

        r0_mm = .1
        r0_mm = .581000
        r0_mm = 1.5
        r0_mm = R0_mm
        r0_mm = 4.
        
        MgO = 5.0 #Mg conc
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        lamdir1 = lamdir0
        
    if disknum==1701:
        from WGM_lib import nMgF2
        R0_mm = 0.366
        r0_mm = R0_mm / 6.        
        MgO = 5.0 #Mg conc
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = nMgF2
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        lamdir1 = lamdir0

    if disknum==2601:
        from WGM_lib import nLNO1
        R0_mm = 2.42
        r0_mm = R0_mm/5.
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
                
    if disknum==2604:
        from WGM_lib import nLNO1
        R0_mm = 0.75
        r0_mm = R0_mm/6.1
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = nLNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    return R0_mm, r0_mm, material, MgO, MgO_th, lamdir0, lamdir1
    
        
#def scaler(x_in,disknum=919,x0=0,slope=1):
#    if disknum == 919:
#        slope = 1.2
#        offset = 20
#    return x_in*slope + offset
    
def scaler(x_in,disknum=919,x0=0,slope=1):
    if disknum == 106: 
        #worked for Cs long letter
        x0 = 0
        slope = 1.
        offset = 0

    if disknum == 9192: 
        #worked for Cs long letter
        x0 = 141.91
        slope = 1.45
        offset = 0.001
        
    if disknum == 9191: 
        #worked for Cs long letter
        x0 = 100
        slope = 1.
        offset = 18.2
        
    if disknum == 9193: 
        #worked for Cs short letter
        x0 = 100
        slope = 1.22
        offset = -3
        offset = 11.5
        
    if disknum == 919: 
        #worked for Cs short letter
        x0 = 100
#        x0 = 200
        slope = 1.22
        offset = -3
        offset = 11.5
        
    if disknum == 966: 
        #worked for Cs short letter
        x0 = 100
        slope = 1.22
        offset = -3
        offset = 11.5

    x_out = (x_in - x0)*slope + x0 + offset
    return x_out
    
def ntest(lda,T,pol=0,MgO=0,MgO_th=5.0,E = 0): return 1.5

def path_leaf(path):
    head, tail = ntpath.split(path)
    tail =  ('.').join(tail.split('.')[:-1])
    return tail or ntpath.basename(head)
    
def tiltcorrect(ydata,windowlen=100,degree = 1):
#    import poly_iter
    x_ground = np.arange(len(ydata))
    x_ground_fit = smooth(x_ground,window_len=int(len(x_ground)/windowlen),window='hanning')
    ydata_fit    = smooth(ydata   ,window_len=int(len(ydata)   /windowlen),window='hanning')
    coeff = np.polyfit(x_ground_fit,ydata_fit, degree)
    
#    if degree == 1:
#        ytilt = tuple(x*coeff[0] + coeff[1] for x in x_ground)
#    if degree == 2:
#        ytilt = tuple(x*x*coeff[0] + x*coeff[1] + coeff[2]  for x in x_ground)
#    if degree == 3:
#        ytilt = tuple(x*x*x*coeff[0] + x*x*coeff[1] + x*coeff[2] + coeff[3]  for x in x_ground)
    ytilt = np.array([], dtype='double')
    for xx in x_ground:
        ytilt = np.append(ytilt,poly_horner(np.flipud(coeff),xx))
    return  ydata - ytilt 

def poly_horner(A, x):
    p = A[-1]
    i = len(A) - 2
    while i >= 0:
        p = p * x + A[i]
        i -= 1
    return p
    
def poly_naive(x, coeff):
    result = coeff[0]
    for i in range(1, len(coeff)):
        result = result + coeff[i] * x**i
    return result
    
def poly_iter(A, x):
    p = 0
    xn = 1
    for a in A:
        p += xn * a
        xn *= x
    return p

def nonorthprism(alpha,beta,d,angtype='deg'): 
#	calculates distances in a triangle with given baseline and its two adjacent angles
	import math
	if angtype=='deg':
		alpha = alpha *math.pi / 180
		beta = beta *math.pi / 180
	h = d / (1 / math.tan(alpha) + 1 / math.tan(beta))	
	d1 = h / math.tan(alpha)
	d2 = h / math.tan(beta)
	return h, d1, d2

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    http://wiki.scipy.org/Cookbook/SignalSmooth
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2-1):-(window_len/2)]
    


def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    https://gist.github.com/endolith/250860
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)

def WGM_norm_freq(L,q = 1,p = 0,Roverr = 1,n=2.22,pol=1):
      """
      This iteratively solves the normalized WGM dispersion equation for a spheroid at the target
      orbital number L and returs the normalized frequency (GHz)
      Roverr stands for ellipticity, R, r are big and small radia in mm
      pol=0 "ordinary", pol=1 "extraordinary"
	lambda =     freq_term  * (    simple dispersion        + (      geom_term                       ) -            pol term        + small contribution
	f      = c / (2 pi Rsemi n) * ( l + a(q) * (l/2)^(1/3)  + (2* p * (Rsemi-rsemi) + Rsemi)/(2*rsemi) - n^(-1+2*pol)/sqrt(n^2 -1  ) + 3*a(q)^2 * (l/2)^(-1/3)/20 )
      """
      import math
      #from WGM_lib import nLNO1
      from scipy.special import ai_zeros
      q -= 1
      AiRoots =  -ai_zeros(40)[0] # >0
      geom_term = p*(math.sqrt(Roverr)-1) + math.sqrt(Roverr)/2.0
      freq_term = 1
      pol_term = n**(-1+2*pol)/math.sqrt(n**2-1)
      
      f = freq_term * (L + AiRoots[q]*(L/2.0)**(1.0/3.0) + geom_term - pol_term + 3*AiRoots[q]**2*(L/2.0)**(-1.0/3.0)/20)
      
      return f

def autocorr(x):
    result = np.correlate(x, x, mode='full')
#    result = np.convolve(x, x, mode='same')
#    result = np.convolve(x, x, mode='full')
    return result[result.size/2.:]
#    return result

def nu2lambda(x):
    #for example: lambda in nm to frequency in GHz
    cms = 299792458 # speed of light
    return cms/x
    
def lambda2nu(x):
    #for example: lambda in nm to frequency in GHz
    cms = 299792458 # speed of light
    return cms/x
    
def lorentz(x=range(1,100),x0=0,gamma=1,amp=1,y0=1):
    return amp / (1.0 + 4*((x - x0) / gamma)**2 ) + y0
    
def gauss(x, *p):
    A, mu, sigma, y0 = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y0
    
def lowpass(cutoff, fs, order=5):
    import numpy as np
    from scipy.signal import butter, lfilter, freqz
    import matplotlib.pyplot as plt
    
    def butter_lowpass(cutoff, fs, order=5):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        return b, a
    def butter_lowpass_filter(data, cutoff, fs, order=5):
        b, a = butter_lowpass(cutoff, fs, order=order)
        y = lfilter(b, a, data)
        return y

    # Filter requirements.
    order = 6
    fs = 30.0       # sample rate, Hz
    cutoff = 3.667  # desired cutoff frequency of the filter, Hz
    
    # Get the filter coefficients so we can check its frequency response.
    b, a = butter_lowpass(cutoff, fs, order)
    
    # Plot the frequency response.
    w, h = freqz(b, a, worN=8000)
    plt.subplot(2, 1, 1)
    plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
    plt.plot(cutoff, 0.5*np.sqrt(2), 'ko')
    plt.axvline(cutoff, color='k')
    plt.xlim(0, 0.5*fs)
    plt.title("Lowpass Filter Frequency Response")
    plt.xlabel('Frequency [Hz]')
    plt.grid()

def read_data(filename,delim = ' '):
# Funktion zum Einlesen der selbsterzeugten kalibrationsdatein
    import csv
    rawdata=[]
	
    reader = csv.reader(open(filename, 'r'), delimiter=delim)#, quotechar='|')

    for line in reader:
		#print len(line)	
          rawdata.append(line)
	#ardata = np.array(rawdata, dtype='float64')
    return rawdata#ardata
    
def read_agilent_DSOx2024a(expdatfolder,expfilenr=0,expfiletype='csv',grdsubtract=False):
#This function read a measured frequeny spectrum and prepares it a bit  
#
#@author: Gerhard Schunk  
    from Ghad_lib import read_data, smooth, path_leaf
    import numpy as np
    import csv
    import glob
    
    headerlines = 4
    x1 = np.array([], dtype='float64')
    y1 = np.array([], dtype='float64')
    y2 = np.array([], dtype='float64')
    y3 = np.array([], dtype='float64')
    y4 = np.array([], dtype='float64')
    
    if expfiletype=='csv':
        expfilepaths = glob.glob(expdatfolder + '//*.csv')
        expfilename =  path_leaf(expfilepaths[expfilenr])   
        reader = csv.reader(open(expfilepaths[expfilenr], 'r'), delimiter=',')#, quotechar='|')
        
        
        line_count = 0
        for line in reader:
            if headerlines > line_count:
                line_count = line_count + 1
                continue
            x1 = np.append(x1,float(line[0]))
            y1 = np.append(y1,float(line[1]))
            y2 = np.append(y2,float(line[2]))
            y3 = np.append(y3,float(line[3]))
            y4 = np.append(y4,float(line[4]))
            
    
    print 'exp data name: ', expfilename, '\n'

    if grdsubtract:
        groundsmoothlen = int(round(len(yfreq)/10))
        yfreq_ground = smooth(yfreq,groundsmoothlen,window='hanning')
        yfreq_ground = yfreq_ground[0:len(yfreq)]
        yfreq_filter = yfreq - yfreq_ground
        for ite in range(yfreq_filter.size):
            if yfreq_filter[ite]<0:
                yfreq_filter[ite] = 0
                
    return  x1, y1, y2, y3, y4, expfilename              
#    return xfreq_GHz, yfreq, expfilename
                
                
                
