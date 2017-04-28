#!/usr/bin/env mpython_q
def factorial_002(n,exact=0):
    """
    The factorial function, n! = special.gamma(n+1).

    If exact is 0, then floating point precision is used, otherwise
    exact long integer is computed.

    - Array argument accepted only for exact=0 case.
    - If n<0, the return value is 0.

    Parameters
    ----------
    n : int or array_like of ints
        Calculate ``n!``.  Arrays are only supported with `exact` set
        to False.  If ``n < 0``, the return value is 0.
    exact : bool, optional
        The result can be approximated rapidly using the gamma-formula
        above.  If `exact` is set to True, calculate the
        answer exactly using integer arithmetic. Default is False.

    Returns
    -------
    nf : float or int
        Factorial of `n`, as an integer or a float depending on `exact`.

    Examples
    --------
    >>> arr = np.array([3,4,5])
    >>> sc.factorial(arr, exact=False)
    array([   6.,   24.,  120.])
    >>> sc.factorial(5, exact=True)
    120L

    """
    import numpy as np
    if exact:
        if n < 0:
            return 0L
        val = 1L
        for k in xrange(1,n+1):
            val *= k
        return val
    else:
        from scipy import special
        n = np.asarray(n)
        sv = special.errprint(0)
        vals = special.gamma(n+1)
        sv = special.errprint(sv)
        return np.where(n>=0,vals,0)

def alert(text):
	print '---------------------------'
	print text
	print '---------------------------\n'

def show_func():
      print 'ntest(lda,T,pol=0,MgO=0,MgO_th=5.0,n = nLNO1,E = 0)'
      print 'nLNO1(lda,T,pol,MgO=0,MgO_th=5.0,n = nLNO1,E = 0)'
      print 'nBBO(lda,T,pol=0,MgO=0,MgO_th=5.0,n = nLNO1,E=0)'
      print 'nBBO1(lda,T,pol=0,MgO=0,MgO_th=5.0,n = nLNO1,E=0)'
      print 'FSR_simple(R,lda,T,pol,L,q = 1,MgO=0,MgO_th=5.0,n = nLNO1,E = 0)'
      print 'WGM_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,n = nLNO1,E = 0)'
      print 'ndiamond(lda,T=20,pol=0)'
      print 'coupling_Q_factor(refr_diamond,refr_disk,radius,wavelength)'

def ntest(lda,T,pol=0,MgO=0,MgO_th=5.0,E = 0): return 1.5

def nLNO1(lda,T,pol,MgO=0,MgO_th=5.0,E = 0):
      """
      This is Sellmeyer eq. for Lithium Niobate from "Influence of the defect structure on the refractive indices
      of undoped and Mg-doped lithium niobate", U.Schlarb and K.Betzler, Phys. Rev. B 50, 751 (1994).
      Temperature in Centigrades, wavelength in microns. pol=0 "ordinary", pol=1 "extraordinary" MgO
      =0 MgO concentration in % MgO_th the threshold concentration. External electric field E in V/cm
      """
      import math
      if pol: # extraordinary
            w0 = 218.203
            mu = 6.4047E-6
            A0 = 3.9466E-5
            A_NbLi = 23.727E-8
            A_Mg = 7.6243E-8
            A_IR = 3.0998E-8
            A_UV = 2.6613
            rr = 33E-10 # r_33 cm/V
      else: # ordinary
            w0 = 223.219
            mu = 1.1082E-6
            A0 = 4.5312E-5
            A_NbLi = -1.4464E-8
            A_Mg = -7.3548E-8
            A_IR = 3.6340E-8
            A_UV = 2.6613
            rr = 11E-10 # r_13 cm/V
      B = A0+A_Mg*MgO
      if MgO < MgO_th: B+=(MgO_th-MgO)*A_NbLi
      def f(T):
             return (T + 273)**2 + 4.0238E5*(1/math.tanh(261.6/(T + 273)) - 1)
      def n(lda,T):
            return math.sqrt(B/((w0+(f(T)-f(24.5))*mu)**(-2)-lda**(-2))-A_IR*lda**2+A_UV)
      return n(lda,T)+0.5*rr*E*n(lda,T)**3

def nLNO1_disp(lda,T,pol,MgO=0,MgO_th=5.0):
      """
      dispersion dn/dlda (in nm) for Lithium Niobate dispersion nLNO1
      """
      import math
      if pol: # extraordinary
            w0 = 218.203
            mu = 6.4047E-6
            A0 = 3.9466E-5
            A_NbLi = 23.727E-8
            A_Mg = 7.6243E-8
            A_IR = 3.0998E-8
            A_UV = 2.6613
      else: # ordinary
            w0 = 223.219
            mu = 1.1082E-6
            A0 = 4.5312E-5
            A_NbLi = -1.4464E-8
            A_Mg = -7.3548E-8
            A_IR = 3.6340E-8
            A_UV = 2.6613
      B = A0+A_Mg*MgO
      if MgO < MgO_th: B+=(MgO_th-MgO)*A_NbLi
      def f(T):
             return (T + 273)**2 + 4.0238E5*(1/math.tanh(261.6/(T + 273)) - 1)
      def n(lda,T):
            return math.sqrt(B/((w0+(f(T)-f(24.5))*mu)**(-2)-lda**(-2))-A_IR*lda**2+A_UV)
      return -(B*lda**(-3)*((w0+(f(T)-f(24.5))*mu)**(-2)-lda**(-2))**(-2)+A_IR*lda)/n(lda,T)

def nBBO(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
      """
      This is Sellmeyer eq. for BBO from http://refractiveindex.info/
      (Handbook of Optics, 3rd edition, Vol. 4. McGraw-Hill 2009);
      T is the angle (in degrees) between the optical axis and polarization:
      T = 0 extraordinary, T = 90 ordinary.
      """
      import math 
      a = 2.7405,2.3730 # ordinary, extraordinary
      b = 0.0184,0.0128
      c = 0.0179,0.0156
      d = 0.0155,0.0044
      lda = lda/1000.0 # nm -> microns
      no = math.sqrt(a[0]-d[0]*lda*lda+b[0]/(lda*lda-c[0]))
      ne = math.sqrt(a[1]-d[1]*lda*lda+b[1]/(lda*lda-c[1]))
      o = math.sin(math.pi*T/180)/no
      e = math.cos(math.pi*T/180)/ne
      return (o*o+e*e)**(-0.5)

def nBBO1(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
      """
      This is Sellmeyer eq. for BBO from Optics Communications 184 (2000) 485-491;
      T is the angle (in degrees) between the optical axis and polarization:
      T = 0 extraordinary, T = 90 ordinary.
      """
      import math 
      a = 2.7359,2.3753 # ordinary, extraordinary
      b = 0.01878,0.01224
      c = 0.01822,0.01667
      d = 0.01471,0.01627
      x = 0.0006081,0.0005716
      y = 0.00006740,0.00006305
      lda = lda/1000.0 # nm -> microns
      no = math.sqrt(a[0]-d[0]*lda*lda+b[0]/(lda*lda-c[0])+x[0]*lda**4-y[0]*lda**6)
      ne = math.sqrt(a[1]-d[1]*lda*lda+b[1]/(lda*lda-c[1])+x[1]*lda**4-y[1]*lda**6)
      o = math.sin(math.pi*T/180)/no
      e = math.cos(math.pi*T/180)/ne
      return (o*o+e*e)**(-0.5)

def nKBBF1(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
      """
      This is Sellmeyer eq. form IEEE JOURNAL OF QUANTUM ELECTRONICS, VOL. 44, NO. 7, JULY 2008 ;
      T is the angle (in degrees) between the optical axis and polarization:
      T = 0 extraordinary, T = 90 ordinary.
      """
      import math 
      a = 1.1713,0.9316 # ordinary, extraordinary
      b = 0.00733,0.00675
      c = 0.01022,0.00169
      lda = lda/1000.0 # nm -> microns
      no = math.sqrt(1+ (a[0]*lda*lda)/(lda*lda - b[0]) - c[0]*lda*lda)
      ne = math.sqrt(1+ (a[1]*lda*lda)/(lda*lda - b[1]) - c[1]*lda*lda)
      o = math.sin(math.pi*T/180)/no
      e = math.cos(math.pi*T/180)/ne
      return (o*o+e*e)**(-0.5)

def ndiamond(lda):
      """
      This is Sellmeyer eq. for diamond from from F.Peter, Z Phys 15, 358 (1923)
      """
      import math 
      a = 0.3306*lda*lda/(lda*lda-175*175)
      b = 4.3356*lda*lda/(lda*lda-106*106)
      return math.sqrt(1+a+b)

def FSR_simple(R,lda,T,pol,L,q = 1,MgO=0,MgO_th=5.0, n = nLNO1,E = 0):
      import math
      #from WGM_lib import nLNO1
      from scipy.special import ai_zeros
      q -= 1
      AiRoots = -ai_zeros(40)[0]
      a = 29.9792/(2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E))
      b = 1+0.5*AiRoots[q]*((L/2.0+0.5)**(1.0/3.0)-(L/2.0-0.5)**(1.0/3.0))
      return a*b

def Get_frac_L(R,r,lda,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0, n = nLNO1,E = 0):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target wavelength lda
      and returs the fractional orbital number L
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in mm
      """
      import math
      #from WGM_lib import nLNO1
      from scipy.special import ai_zeros
      q -= 1
      r = math.sqrt(r*R) # redefine the rim radius to spheroid semi-axis
      AiRoots = -ai_zeros(40)[0]
      freq_term = 2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E)*1E7/lda
      geom_term = (2*p*(R-r)+R)/2.0/r
      pol_term = n(lda,T,pol,MgO,MgO_th,E)**(1-2*pol)/math.sqrt(n(lda,T,pol,MgO,MgO_th,E)**2-1)
      L = freq_term
      dL = 1
      while math.fabs(dL) >= 1E-9:
            L1 = freq_term - AiRoots[q]*(L/2.0)**(1.0/3.0) - geom_term + pol_term\
                 - 3*AiRoots[q]**2*(L/2.0)**(-1.0/3.0)/20
            dL = L1-L
            L = L1
      return L

def WGM_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,n = nLNO1,E = 0):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target
      orbital number L and returs the wavelength (nm) and frequency (GHz)
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in mm
      """
      import math
      #from WGM_lib import nLNO1
      from scipy.special import ai_zeros
      q -= 1
      r = math.sqrt(r*R) # redefine the rim radius to spheroid semi-axis
      AiRoots =  -ai_zeros(40)[0] # >0
      geom_term = (2*p*(R-r)+R)/2.0/r
      def nm2GHz(lda): return 2.99792E8/lda
      def wl1(lda):
            freq_term = 2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E)*1E7
            pol_term = n(lda,T,pol,MgO,MgO_th,E)**(1-2*pol)/math.sqrt(n(lda,T,pol,MgO,MgO_th,E)**2-1)
            return freq_term/(L+AiRoots[q]*(L/2.0)**(1.0/3.0)\
                        +geom_term-pol_term+3*AiRoots[q]**2*(L/2.0)**(-1.0/3.0)/20)
      lda0 = 2*math.pi*R*n(1000,T,pol,MgO,MgO_th,E)*1E7/L # initial wavelength guess based on n(1000 nm)
      lda = wl1(lda0)
      #print lda0, lda
      df = 1
      while math.fabs(df) > 1E-6: #1 kHz accuracy
            lda1 = wl1(lda)
            df = nm2GHz(lda1)-nm2GHz(lda)
            lda = lda1
            #print df,lda1,n(lda,T,pol,MgO,MgO_th,E)
      return lda, nm2GHz(lda)

def Gaunt(L1,m1,L2,m2,L3,m3):
      """
      Gaunt formula for Legendre polynomials overlap. Large factorials are treated in Stirling's approximation.
      """
      import math
      def fa(n):
            if n<=10: fac = float(math.factorial(n))
            else: fac = math.sqrt(2*math.pi*n)*(n/math.e)**n
            return fac
      params = [(L1,m1),(L2,m2),(L3,m3)]
      params = sorted(params,key=lambda list: list[1])
      l,u = params[-1]
      prms = params[:-1]
      prms = sorted(prms,key=lambda list: list[0])
      n,w = prms[0]
      m,v = prms[1]
      if not v+w == u: return 0
      else:
            s = (l+m+n)/2.0
            if not s-round(s) == 0 : return 0
            elif m+n <l : return 0
            elif m-n >l : return 0
            else:
                  p = int(max(0, n-m-u))
                  q = int(min(n+m-u,l-u,n-w))
                  a = (-1)**(s-m-w)*math.sqrt((2*l+1)*(2*m+1)*(2*n+1)*fa(l-u)*fa(n-w)/2.0/fa(m-v))/4/math.pi/math.pi
                  b = a*math.sqrt(math.sqrt(2*math.pi*(m+v)*(n+w)/(l+u))/(s-m)/(s-n))/2.0/math.pi
                  c = b*math.exp(0.5*(m+n-l+u-v-w))*(s-l)**(l-s)
                  d = c*((n+w)/float(l+n-m))**(s-m)*((n+w)/2.0)**((w+m-l)/2.0)
                  t0 = d*2**((w-n)/2.0)*((m+v)/float(l+m-n))**((l+m-n)/2.0)*(m+v)**((v-l+n)/2.0)*(l+m-n)**(l-u)
                  sum  = 0
                  for t in range (p,q+1):
                        f1 = (-1)**t*fa(m+n-u-t)/fa(t)/fa(l-u-t)/fa(n-w-t)
                        f2 = f1*((l+u+t)/float(m-n+u+t))**(0.5+t)
                        f3 = f2*((l+u+t)/float(l+u))**(0.5*l+0.5*u)
                        f4 = f3*((l+u+t)/float(l+m+n))**((l+m+n)/2.0)*(l+u+t)**((u-m-n)/2.0)
                        sum+= f4*((l+m-n)/float(u+m-n+t))**(u+m-n)
                  return t0*sum

def Gaunt_low(L1,m1,L2,m2,L3,m3):
      """
      Gaunt formula for Legendre polynomials overlap computed exactly. Good only for low orders.
      """
      import math
      def fa(n):
            if n<=10: fac = float(math.factorial(n))
            else: fac = math.sqrt(2*math.pi*n)*(n/math.e)**n
            return fac
      def fa_ratio(n1, n2): #calculates n1!/n2!
            p = 1.0
            for j in range (min(n1, n2)+1, max(n1, n2)+1):
                  p = p*j
            if n1 >= n2: return p
            else: return 1/p
      params = [(L1,m1),(L2,m2),(L3,m3)]
      params = sorted(params,key=lambda list: list[1])
      l,u = params[-1]
      prms = params[:-1]
      prms = sorted(prms,key=lambda list: list[0])
      n,w = prms[0]
      m,v = prms[1]
      if not v+w == u: return 0
      else:
            s = (l+m+n)/2.0
            if not s-round(s) == 0 : return 0
            elif m+n <l : return 0
            elif m-n >l : return 0
            else:
                  s = int(s)
                  p = int(max(0, n-m-u))
                  q = int(min(n+m-u,l-u,n-w))
                  #print l,m,n,s
                  #print fa(l-u),fa(n-w), fa(m-v),fa(s-l)
                  a1 = (-1)**(s-m-w)*math.sqrt((2*l+1)*(2*m+1)*(2*n+1)*fa(l-u)*fa(n-w)/fa(m-v)/math.pi)/fa(s-l)/float(2*s+1)/4.0
                  a2 = math.sqrt(fa(m+v)*fa(n+w)/fa(l+u))*fa(s)/fa(s-m)/fa(s-n)
                  sum  = 0
                  for t in range (p,q+1):
                        b1 = (-1)**t*fa(m+n-u-t)/fa(t)/fa(l-u-t)/fa(n-w-t)    
                        b2 = fa_ratio(l+u+t,2*s)*fa_ratio(2*s-2*n,m-n+u+t)
                        sum+= b1*b2
                  #print a1, a2, sum
                  return a1*a2*sum

def RadOvlp(R,L1,q1,L2,q2,L3,q3):
      """
      Calculates radial overlap of three Airy functions
      """
      import math
      from scipy.special import airy, ai_zeros
      from scipy.integrate import quad
      AiRoots =  -ai_zeros(40)[0]
      q1, q2, q3 = q1-1, q2-1, q3-1
      def nk(L,q): return L*(1+AiRoots[q]*(2*L**2)**(-1.0/3.0))
      def AiArg(x,L,q): return -AiRoots[q]*(L-nk(L+0.5,q)*x)/(L-nk(L+0.5,q))
      def RadFunc(x,L,q): return airy(AiArg(x,L,q))[0]/math.sqrt(x)
      L0 = min(L1,L2,L3)
      q0 = max(q1,q2,q3)+1
      Rmin = 1-5*math.sqrt(q0)*L0**(-2.0/3.0)
      N1 = R**3*quad(lambda x:  x*x*RadFunc(x,L1,q1)**2, Rmin, 1)[0]
      N2 = R**3*quad(lambda x:  x*x*RadFunc(x,L2,q2)**2, Rmin, 1)[0]
      N3 = R**3*quad(lambda x:  x*x*RadFunc(x,L3,q3)**2, Rmin, 1)[0]
      #print N1, N2, N3
      def Core(x): return RadFunc(x,L1,q1)*RadFunc(x,L2,q2)*RadFunc(x,L3,q3)*x*x
      return R**3*quad(lambda x: Core(x), Rmin, 1)[0]/math.sqrt(N1*N2*N3)

def AngOvlp(L1,p1,L2,p2,L3,p3,limit=3):
      """
      Calculates angular overlap of three spherical functions. Default limit is sufficient for up to p = 10
      """
      import numpy
      from scipy.special import hermite
      from scipy.integrate import quad
      #from factorial import factorial
      #def HG(p,L,x): return (2**L*math.factorial(p))**-0.5*hermite(p)(x*math.sqrt(L))*(L/math.pi)**0.25*math.exp(-L*x**2/2.0)
      #def HG(p,L,x): return hermite(p)(x*math.sqrt(L))*math.exp(-L*x**2/2.0) # removed a large factor that would drop out in normalization anyway
      # michls try
      #def HG(p,L,x): return hermite(p)(x*numpy.sqrt(L))*numpy.exp(-L*x**2/2.0)/(numpy.sqrt(factorial(p)))
      def HG(p,L,x): return hermite(p)(x*numpy.sqrt(L))*numpy.exp(-L*x**2/2.0)/(numpy.sqrt(factorial_002(p)))
      def HGnorm(p,L,x): return HG(p,L,x)*(4*numpy.pi*quad(lambda x: HG(p,L,x)**2, 0, ((p+1)*numpy.pi/L/2)**0.5)[0])**(-0.5)
      return 4*numpy.pi*quad(lambda x: HGnorm(p1,L1,x)*HGnorm(p2,L2,x)*HGnorm(p3,L3,x), 0, ((max(p1,p2,p3)+1)*numpy.pi/min(L1,L2,L3)/2)**0.5)[0]

def T_phasematch(L1,p1,q1,L2,p2,q2,L3,p3,q3,R,r,T0=70,MgO=0,MgO_th=5.0,n = nLNO1,E = 0,Tstep=0.1):
      '''
      Finds 1/wl2 + 1/wl2 = 1/wl3, L1+L2 -> L3 phase matching temperature down to the accuracy deltaf (GHz)
      based on the initial guess T0 with initial search step Tstep
      Requires importing Sellmeyer equations. This program is unaware of the orbital selection rules!
      '''
      import math
      from WGM_lib import WGM_freq
      def Df(T): # freqency detuning in GHz
            lda_1, freq_1 = WGM_freq(R,r,L1,T,0,q1,p1,MgO,MgO_th,E,n)
            lda_2, freq_2 = WGM_freq(R,r,L2,T,0,q2,p2,MgO,MgO_th,E,n)
            lda_3, freq_3 = WGM_freq(R,r,L3,T,1,q3,p3,MgO,MgO_th,E,n)
            return freq_3-freq_1-freq_2
      T = T0
      df = Df(T)
      while math.fabs(df)>1E-6: #1 kHz accuracy
            #print "T = ",T," df = ", df, 
            df1 = Df(T+Tstep)
            ddfdt = (df1-df)/Tstep
            dT = df/ddfdt
            T -= dT
            Tstep = math.fabs(dT/2.0)   
            df = Df(T)
            #print " dT = ", dT, " T -> ",T, " df -> ", df
            if T<-200 or T> 500: break
      return T

def T_phasematch_blk(wl1,wl2, T0=70,MgO=0,MgO_th=5.0,E = 0,Tstep=0.1, n = nLNO1):
      '''
      Finds 1/wl2 + 1/wl2 = 1/wl3 (wl in nanometers), collinear bulk phase matching temperature 
      based on the initial guess T0 with initial search step Tstep
      Requires importing Sellmeyer equations. 
      '''
      import math
      wl3 = 1/(1.0/wl1+1.0/wl2)
      def Dk(T):
            k1 = 2*math.pi*n(wl1,T,0,MgO,MgO_th,E)*1e7/wl1 # cm^-1
            k2 = 2*math.pi*n(wl2,T,0,MgO,MgO_th,E)*1e7/wl2
            k3 = 2*math.pi*n(wl3,T,1,MgO,MgO_th,E)*1e7/wl3
            return k3-k2-k1
      T = T0
      dk = Dk(T)
      while math.fabs(dk)>0.1: #10 cm beat length
            dk1 = Dk(T+Tstep)
            ddkdt = (dk1-dk)/Tstep
            dT = dk/ddkdt
            T -= dT
            Tstep = math.fabs(dT/2)   
            dk = Dk(T)
            #print " dT = ", dT, " T -> ",T, " dk -> ", dk
            if T<-200 or T> 500:
                  print "No phase matching!"
                  break
      return T

def SH_phasematch_blk(angle, wl0=1000,step=1.01, n = nBBO):
      '''
      Finds the fundamental wavelength (in nm) for frequency doubling o+o->e at a given angle between the optical axis and the e.
      Requires importing Sellmeyer equations. 
      '''
      import math
      def Dk(wl):
            k1 = 2*math.pi*n(wl,90)*1e7/wl # ordinary cm^-1
            k2 = 4*math.pi*n(wl/2.0,angle)*1e7/wl # extraordinary with angle
            return k2-2*k1
      wl = wl0
      dk = Dk(wl)
      while math.fabs(dk)>0.1: #10 cm beat length
            dk1 = Dk(wl*step)
            ddkdt = (dk1-dk)/step
            dwl = dk/ddkdt
            wl = wl/dwl
            #step = min(step,math.fabs(dwl/2) ) 
            dk = Dk(wl)
            print " dwl = ", dwl, " wl -> ",wl, " dk -> ", dk
      return wl

def Phasematch(wlP,p1,q1,p2,q2,p3,q3,R,r,T,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      '''
      
      NOT FULLY FUNCTIONAL!!      
                  
      Finds wl1 and wl2 such that 1/wl2 + 1/wl2 = 1/wlP and  L1-p1+L2-p2 = L3-p3 at a given temperature T.
      Requires importing Sellmeyer equations. R is evaluated at T. This program is unaware of the other orbital
      selection rules! L1, L2, L3 are allowed to be non-integer.
      '''
      import math
      from WGM_lib import Get_frac_L
      def Dm(wl): # m3-m1-m2, orbital detuning
            L3 = Get_frac_L(R,r,wlP,T,1,q3,p3,MgO,MgO_th)
            L1 = Get_frac_L(R,r,wl,T,0,q1,p1,MgO,MgO_th)
            wl2 = 1.0/(1.0/wlP-1.0/wl)
            L2 = Get_frac_L(R,r,wl2,T,0,q2,p2,MgO,MgO_th)
            return L2-L3+p3+L1-p1-p2
      wl = 2*wlP #try at degeneracy
      dm = Dm(wl)
      print "dm = ", dm
      ddm = 0
      wlstep = 1.0 # chosing an appropriate wavelength step so the m variation is not too large or too small
      while math.fabs(ddm)<1:
            wlstep = 2*wlstep
            ddm = Dm(wl+wlstep)-dm
      print "wlstep = ", wlstep, " initial slope = ", ddm/wlstep
      while math.fabs(dm)>1E-6: # equiv *FSR accuracy
            dm1 = Dm(wl+wlstep)
            print "dm1 = ", dm1
            slope = (dm1-dm)/wlstep
            print "slope = ", slope, "nm^-1"
            dwl = dm/slope
            wl -= dwl
            print "new wl = ", wl
            wlstep = min(math.fabs(dwl/2), wlstep)   
            dm = Dm(wl)
            print " d wl = ", dwl, " wl -> ",wl, " dm -> ", dm
            if math.fabs(dm)> Get_frac_L(R,r,2*wlP,T,0,q1,p1,MgO,MgO_th): break
      L1 = Get_frac_L(R,r,wl,T,0,q1,p1,MgO,MgO_th)
      wl2 = 1.0/(1.0/wlP-1.0/wl)
      L2 = Get_frac_L(R,r,wl2,T,0,q2,p2,MgO,MgO_th)
      return L1, L2, wl, wl2

def n_eff_FSR(wl,R,T,FSR,pol,MgO=0,MgO_th=5.0,n = nLNO1,E = 0,n_disp=nLNO1_disp):
      '''
      Finds effective refraction index (material+geometrical) based on the FSR measuremert (in GHz).
      R in cm, wl in nm. Also tries to guess the appropriate q and L. The first order approximation.
      '''
      import math
      from scipy.special import ai_zeros
      import bisect
      nr = n(wl,T,pol,MgO,MgO_th,E)
      wl_deriv = n_disp(wl,T,pol,MgO,MgO_th)
      denom = 2*math.pi*3*R*FSR*(nr-wl_deriv*wl)/29.9792-2
      neff = nr/denom
      AiRoots = -ai_zeros(40)[0]
      L = 2*math.pi*R*neff*1e7/wl
      a_q = (nr/neff-1)*(2*L*L)**(1.0/3.0)
      j = bisect.bisect(AiRoots, a_q)
      #q = 0
      #if a_q>AiRoots[0]/2.0: q = 1
      if j: q = j+(a_q-AiRoots[j-1])/(AiRoots[j]-AiRoots[j-1])
      else: q = 1+(a_q-AiRoots[0])/(AiRoots[1]-AiRoots[0])
      return neff,q,L

def n_eff(wl,R,r,T,pol,q,p,MgO=0,MgO_th=5.0,n = nLNO1,E = 0):
      '''
      Finds effective refraction index (material+geometrical) for a given WGM.
      R,r in cm, wl in nm. 
      '''
      import math
      from WGM_lib import Get_frac_L
      L = Get_frac_L(R,r,wl,T,pol,q,p,MgO,MgO_th,E,n)
      neff = L*wl*1e-7/(2*math.pi*R)
      return neff,L

def Qc(wl,R,d,pol,ns,nc):
      '''
      Coupling Q for a spherical resonator with radius R (cm) for wavelength wl (nm) and gap d (nm)
      From Gorodetsky's book, pg 246, for the "o" and with a factor n**4 for the "e"
      '''
      import math
      R = R*1e7 # cm -> nm
      k = 2*math.pi/wl
      root = math.sqrt(2*math.pi**5*ns/(nc*nc-ns*ns))
      a = root*(ns*ns-1)*(R/wl)**1.5
      b = math.exp(2*d*k*math.sqrt(ns*ns-1))
      qq = 1
      if pol: qq = ns*ns*ns*ns
      return a*b*qq

def coupling_Q_factor(radius,lda,distance,T=60,pol=0,MgO=5.0,MgO_th=5.0,nprism=ndiamond,ndisk=nLNO1):
      '''
      This function caculates the coupling bandwidth depended on the:
            - distance between prisma and disk
            - the wavelength
            - thecoupling prism material
      '''      
      import numpy as np
      #define the refractive indices
      refr_disk=ndisk(lda,T,pol,MgO,MgO_th,E=0)
      refr_coupler = nprism(lda,T,pol)
      
      lda=lda*10**(-9)
      distance=distance*10**(-9)
      
      Q_load = np.sqrt(2.0)*np.pi**(2.5)*refr_disk**(0.5)*(refr_disk-1) /\
      np.sqrt(refr_coupler**2-refr_disk**2) * (radius/lda)**(1.5)* \
      np.exp(2.0*distance*2.0*np.pi/lda*np.sqrt((refr_disk**2)-1.0))
      return (Q_load)
      