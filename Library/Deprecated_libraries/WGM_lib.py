#!/usr/bin/env mpython_q
def ntest(lda,T,pol=0,MgO=0,MgO_th=5.0,E = 0): return 1.5

def thermexpLiNbO3V3(R0,T,dir,lamdir0,lamdir1):
	"""
	This calculates thermal expansion of LiNbo3 congruent material. Depends on temperature and crystall axis
	Same notation for orientation as for refractive index pol (pol = 1 # ordinar xy axis (0) or extraordinary (1) z axis)
	R0 is at 0 degree celsius.	
	see http://books.google.de/books?id=dukadDN0iL8C&pg=PA73&lpg=PA73&dq=lithium+niobate+thermal+expansion+coefficient&source=bl&ots=wiiCzA56bZ&sig=JEVfX-gI5lpWZqJsjEuuU8ZUKGM&hl=en&sa=X&ei=kGsgUsjHL4_EsgaQ04CIDw&ved=0CFMQ6AEwBDgK#v=onepage&q=lithium%20niobate%20thermal%20expansion%20coefficient&f=false	
	plot at C:\Users\gschunk\Documents\Libary\python	
      """
#	T = T - 25 #zero expansion is at 0 degree
	if dir == 0:
		Rt = R0 *(1 + lamdir0[0]*T + .5*(lamdir0[0]**2 + lamdir0[1])*T**2 + 0.333*(-3.901*10**(-11) + lamdir0[0]*lamdir0[1])*T**3)
	if dir == 1:
		Rt = R0 *(1 + lamdir1[0]*T + .5*(lamdir1[0]**2 + lamdir1[1])*T**2 + 0.333*(-3.175*10**(-11) + lamdir1[0]*lamdir1[1])*T**3)
	return Rt 
 
def thermexpLiNbO3V2(R0,T,dir,lamdir0=[13.33863*10**(-6),.5*0.03078*10**(-6),-0.333333*4.81345*10**(-11)],lamdir1=[3.91101*10**(-6),.5*0.00868*10**(-6),-0.333333*5.401055*10**(-11)]):
	"""
	This calculates thermal expansion of LiNbo3 congruent material. Depends on temperature and crystall axis
	Same notation for orientation as for refractive index pol (pol = 1 # ordinary (0) or extraordinary (1))
	R0 is at 0 degree celsius.	
      expansion along xy axis = pol =0	
	see http://books.google.de/books?id=dukadDN0iL8C&pg=PA73&lpg=PA73&dq=lithium+niobate+thermal+expansion+coefficient&source=bl&ots=wiiCzA56bZ&sig=JEVfX-gI5lpWZqJsjEuuU8ZUKGM&hl=en&sa=X&ei=kGsgUsjHL4_EsgaQ04CIDw&ved=0CFMQ6AEwBDgK#v=onepage&q=lithium%20niobate%20thermal%20expansion%20coefficient&f=false	
	plot at C:\Users\gschunk\Documents\Libary\python	
      coeff: lamdir0=[x,x**2,x**3]
	"""
	if dir == 0:
		Rt = R0 + R0* (lamdir0[0]*T + lamdir0[1]*T**2 + + lamdir0[2]*T**3  )
	if dir == 1:
		Rt = R0 + R0* (lamdir1[0]*T + lamdir1[1]*T**2 + + lamdir1[2]*T**3  )
	return Rt 
 
def thermexpLiNbO3(R0,T,dir):
	"""
	This calculates thermal expansion of LiNbo3 congruent material. Depends on temperature and crystall axis
	Same notation for orientation as for refractive index pol (pol = 1 # ordinary (0) or extraordinary (1))
	R0 is at 0 degree celsius.			
	see http://books.google.de/books?id=dukadDN0iL8C&pg=PA73&lpg=PA73&dq=lithium+niobate+thermal+expansion+coefficient&source=bl&ots=wiiCzA56bZ&sig=JEVfX-gI5lpWZqJsjEuuU8ZUKGM&hl=en&sa=X&ei=kGsgUsjHL4_EsgaQ04CIDw&ved=0CFMQ6AEwBDgK#v=onepage&q=lithium%20niobate%20thermal%20expansion%20coefficient&f=false	
	plot at C:\Users\gschunk\Documents\Libary\python	
	"""
	if dir == 0:
#		expansion along xy axis = pol =0, Temp in degree celsius
#		alpha = 1.67
#		alpha = .8
		#R0*(1 + alpha* 10**-5*(T-20))
#		Rt = R0*(13.33863 + 0.03078*T  - 4.81345*10**(-5) * T**2 + 2.82544*10**(-8)*T**3)
		#delta a T delta t -> a  = integral delta a delta t dt = A + a (Rt)
#		Rt = R0* (1 + (13.33863 + 0.03078*T  - 4.81345*10**(-5) * T**2 + 2.82544*10**(-8)*T**3 ) * T * 10**(-6)  ) 
		Rt = R0 + R0* ( (13.33863*T + .5 * 0.03078*T**2  - .333333* 4.81345*10**(-5) * T**3 + .25 *2.82544*10**(-8)*T**4) * 10**(-6) ) 
	if dir == 1:
#		expansion along z axis = pol =1
#		alpha = .4
#		alpha = 15
#		Rt = R0* (1 + (3.91101 + 0.00868*T  -5.401055*10**(-5)  * T**2 + 6.89656*10**(-8)*T**3) * T * 10**(-6)  )  
		Rt = R0 + R0* ( (3.91101*T + .5 * 0.00868*T**2  - .333333* 5.401055*10**(-5) * T**3 + .25 *6.89656*10**(-8)*T**4) * 10**(-6) ) 
	return Rt 
	
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

def nLNO1xcut(lda,T,pol,MgO=0,MgO_th=5.0,E = 0):
	#this uses the seelmeier equation for lithium niobate in an x-cut configuration
	from WGM_lib import nLNO1
	if pol == 0:
		neff = (nLNO1(lda,T,0,MgO=0,MgO_th=5.0,E = 0) + nLNO1(lda,T,1,MgO=0,MgO_th=5.0,E = 0)) / 2.
#		neff = nLNO1(lda,T,1,MgO,MgO_th,E)
	if pol == 1:
		neff = nLNO1(lda,T,0,MgO,MgO_th,E)
	return neff
	
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

def nMgF2(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
      """
      Numeric phase matching analysis for 3 dr harmonic generation e, e, e \
      -> o in xy -cut Subscript[MgF, 2]
      In this version we fix R, r and vary the wavelength
      only lambda in nm is used
      pol=0 "ordinary", pol=1 "extraordinary" 
      """
      import math
      a = 0.48755108, 0.41344023 # (* ordinary,extraordinary*)
      b = 0.39875031, 0.50497499
      c = 2.3120353, 2.4904862
      d = 43.38408, 36.84262
      u = 94.61442, 90.76162
      v = 23793.604, 23771.995
      n = math.sqrt(1 + (a[pol]*lda**2)/(lda**2 - d[pol]**2) + (b[pol]*lda**2)/(lda**2 - u[pol]**2) + (c[pol]*lda**2)/(lda**2 - v[pol]**2))
      return n
      
def nT_MgF2(lda,A,T,*args):
      """
      This uses the n_MgF2 at 19C and dn/dT data from Corning to calculate temperature-dependent n
      """
      import math
      import lagrange_int
      from WGM_lib1 import n_MgF2
      no19, ne19 = n_MgF2(lda,90), n_MgF2(lda,0)
      wavelengths = [457.9,632.8,1150.0,3390.0]
      temperatures = range(-180,220,20)
      dne_dT = [[ 0.177,  0.165,  0.144,  0.15 ],[ 0.168,  0.154,  0.133,  0.14 ],[ 0.159,  0.144,  0.122,  0.13 ],
       [ 0.151,  0.133,  0.11 ,  0.12 ],[ 0.142,  0.122,  0.099,  0.12 ],[ 0.133,  0.112,  0.088,  0.11 ],[ 0.124,  0.101,  0.077,  0.1  ],
       [ 0.116,  0.09 ,  0.065,  0.09 ],[ 0.107,  0.08 ,  0.054,  0.08 ],[ 0.098,  0.069,  0.043,  0.07 ],[ 0.089,  0.058,  0.032,  0.06 ],
       [ 0.081,  0.048,  0.02 ,  0.05 ],[ 0.072,  0.037,  0.009,  0.04 ],[ 0.063,  0.027, -0.002,  0.03 ],[ 0.054,  0.016, -0.013,  0.02 ],
       [ 0.046,  0.005, -0.025,  0.01 ],[ 0.037, -0.005, -0.036,  0.   ],[ 0.028, -0.016, -0.047, -0.01 ],[ 0.019, -0.027, -0.059, -0.02 ],
       [ 0.011, -0.037, -0.07 , -0.03 ]]
      dno_dT = [[ 0.244,  0.223,  0.199,  0.2  ],[ 0.234,  0.212,  0.188,  0.2  ],[ 0.225,  0.201,  0.177,  0.19 ],
       [ 0.215,  0.19 ,  0.166,  0.18 ],[ 0.205,  0.179,  0.155,  0.17 ],[ 0.196,  0.168,  0.144,  0.16 ],
       [ 0.186,  0.157,  0.132,  0.15 ],[ 0.176,  0.146,  0.121,  0.14 ],[ 0.166,  0.135,  0.11 ,  0.13 ],
       [ 0.157,  0.124,  0.099,  0.12 ],[ 0.147,  0.112,  0.088,  0.11 ],[ 0.137,  0.101,  0.077,  0.1  ],
       [ 0.128,  0.09 ,  0.066,  0.1  ],[ 0.118,  0.079,  0.054,  0.09 ],[ 0.108,  0.068,  0.043,  0.08 ],
       [ 0.099,  0.057,  0.032,  0.07 ],[ 0.089,  0.046,  0.021,  0.06 ],[ 0.079,  0.035,  0.01 ,  0.05 ],
       [ 0.07 ,  0.024, -0.001,  0.04 ],       [ 0.06 ,  0.013, -0.013,  0.03 ]]
      # interpolating to given temperature wavelength lda
      dne_dT1,dno_dT1 = [],[]
      for x in dne_dT:
            xa,xintderiv,xerror_codex = lagrange_int(lda,4,wavelengths,x,3,0)
            dne_dT1.append(xa)
      for x in dno_dT:
            xa,xintderiv,xerror_codex = lagrange_int(lda,4,wavelengths,x,3,0)
            dno_dT1.append(xa)
      # computing crude n(T) from derivatives
      nos,nes = [], []
      for ind in range(20):
            a,b = min(ind,10),max(ind,10)
            s = 0
            if ind>10: s = 1
            if ind<10: s = -1
            nos.append(no19+1e-5*20*s*sum(dno_dT1[a:b]))
            nes.append(ne19+1e-5*20*s*sum(dne_dT1[a:b]))      
      # interpolating n to given temperature T
      ne,xintderiv,xerror_codex = lagrange_int(T,20,temperatures,nes,3,0)
      no,xintderiv,xerror_codex = lagrange_int(T,20,temperatures,nos,3,0)
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
      
def ndiamond(lda):
      """
      This is Sellmeyer eq. for diamond from from F.Peter, Z Phys 15, 358 (1923)
      """
      import math 
      a = 0.3306*lda*lda/(lda*lda-175*175)
      b = 4.3356*lda*lda/(lda*lda-106*106)
      return math.sqrt(1+a+b)

def n_fused_silica(lda):
      """
      This is Sellmeyer eq. for fused silica http://refractiveindex.info/?group=GLASSES&material=F_SILICA
      """
      import math 
      a = 0.6961663*lda*lda/(lda*lda-68.4043*68.4043)
      b = 0.4079426*lda*lda/(lda*lda-116.2414*116.2414)
      c = 0.8974794*lda*lda/(lda*lda-9896.161*9896.161)
      return math.sqrt(1+a+b+c)
						
def FSR_bulk(R,lda,T,pol,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      import math
      #from WGM_lib import nLNO1
      FSR = 29.9792/(2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E))
      return FSR
						
def FSR_simple(R,lda,T,pol,L,q = 1,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      import math
      #from WGM_lib import nLNO1
      from scipy.special import ai_zeros
      q -= 1
      AiRoots = -ai_zeros(40)[0]
      a = 29.9792/(2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E))
      b = 1+0.5*AiRoots[q]*((L/2.0+0.5)**(1.0/3.0)-(L/2.0-0.5)**(1.0/3.0))
      return a*b

def Get_frac_L(R,r,lda,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target wavelength lda
      and returs the fractional orbital number L
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in mm
      """
      import math
      #from WGM_lib import nLNO1
	#there are five terms contributing to L
	# frequency, geom and pol term stay constant, the rest converge to the correct value
	
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
#
def Find_modeinsweep(fsweep,Lstart,R,r,lda,T,pol,q,p,MgO,MgO_th=5.0,E = 0, n = nLNO1):
      """
      This iteratively reduces L until mode is in sweep range
      """
      import math
      wl, f = WGM_freq(R,r,Lstart,T,pol,q,p,MgO,MgO_th,E,n)
      FSR = WGM_FSR(R,r,Lstart,T,pol,q,p,MgO,MgO_th,E,n)
#      print 'n, FSR',n, FSR
      fdiff = f- fsweep
      f1 = f
      L1 = Lstart
#      x = 1
      while math.fabs(fdiff) >= FSR:
#      while math.fabs(4) >= x:
            L1 -= 1
            wl1, f1 = WGM_freq(R,r,L1,T,pol,q,p,MgO,MgO_th,E,n)
            fdiff = f1- fsweep
#            print 'f diff', fdiff, ' GHz'
#            x += 1
      return f1, L1
	
def Get_Rfromf(fmeas,Lmeas,R,r,lda,T,pol,q,p,MgO,MgO_th=5.0,E = 0, n = nLNO1):
      """
      This iteratively changes Radius until mode with given mode numbers fits exactly to measured frequency
      """
      import math
      R1 = R

      fmism = 1
      while math.fabs(fmism) >= 0.000001:
           wl1, f1 = WGM_freq(R1,r,Lmeas,T,pol,q,p,MgO,MgO_th,E,n) 
           fmism = fmeas-f1      
           deltaR = fmism/f1 * R1 / 2.									
#           print 'R loop is', R1, 'delta R is', deltaR
           R1 = R1 - deltaR					
      return R1, fmism
						
def WGM_FSR(R,r,L,T,pol,q,p ,MgO,MgO_th=5.0,E = 0, n = nLNO1):
     wlmode_minus, frmode_minus = WGM_freq(R,r,L-1.,T,pol,q,p,MgO,MgO_th,E,n) 
#     wlmode, frmode = WGM_freq(R,r,L,T,pol,q,p,MgO,MgO_th,E, n) 
     wlmode_plus, frmode_plus   = WGM_freq(R,r,L+1.,T,pol,q,p,MgO,MgO_th,E,n) 
     FSR = (frmode_plus - frmode_minus ) / 2.
     return FSR
    
def WGMxcut_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target
      orbital number L and returs the wavelength (nm) and frequency (GHz)
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in cm (Dmitry had mm here)
	lambda = frequ + geom_term - pol term
	lambda = 2*pi*R*n / ( l + a(q) * (l/2)^(1/3)  + (2* p * (R-r) + R)/(2*r) - n^(1-2*pol)/sqrt(n^2 -1  ) +3*a(q)^2 * (l/2)^(-1/3)/20 )
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
            if pol ==0:
                n_eff = (n(lda,T,0,MgO,MgO_th,E)+n(lda,T,1,MgO,MgO_th,E))/2.
            else:
                n_eff = n(lda,T,0,MgO,MgO_th,E)
            freq_term = 2*math.pi*R*n_eff*1E7
            pol_term = n_eff**(-1+2*pol)/math.sqrt(n_eff**2-1)
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
      
def WGM_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target
      orbital number L and returs the wavelength (nm) and frequency (GHz)
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in cm (Dmitry had mm here)
	lambda = frequ + geom_term - pol term
	lambda = 2*pi*R*n / ( l + a(q) * (l/2)^(1/3)  + (2* p * (R-r) + R)/(2*r) - n^(1-2*pol)/sqrt(n^2 -1  ) +3*a(q)^2 * (l/2)^(-1/3)/20 )
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
            pol_term = n(lda,T,pol,MgO,MgO_th,E)**(-1+2*pol)/math.sqrt(n(lda,T,pol,MgO,MgO_th,E)**2-1)
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

def WGM_freq_m(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = nLNO1):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target
      orbital number L and returs the wavelength (nm) and frequency (GHz)
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in cm (Dmitry had mm here)
	lambda = frequ + geom_term - pol term
	lambda = 2*pi*R*n / ( l + a(q) * (l/2)^(1/3)  + (2* p * (R-r) + R)/(2*r) - n^(1-2*pol)/sqrt(n^2 -1  ) +3*a(q)^2 * (l/2)^(-1/3)/20 )
      [1] I. Breunig, B. Sturman, F. Sedlmeir, H. G. L. Schwefel, and K. Buse, Optics Express 21, 30683 (2013).
      [1] Y. Demchenko and M. Gorodetsky, Journal of the Optical Society of America B 30, 3056 (2013).
      """
      import math
      #from WGM_lib import nLNO1
      from scipy.special import ai_zeros
      q -= 1
      r = math.sqrt(r*R) # redefine the rim radius to spheroid semi-axis
      AiRoots =  -ai_zeros(40)[0] # >0
      geom_term = (2*p*+1)*R/2.0/r
      m = L-p
      def nm2GHz(lda): return 2.99792E8/lda
      def wl1(lda):
            freq_term = 2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E)*1E7
            pol_term = n(lda,T,pol,MgO,MgO_th,E)**(-1+2*pol)/math.sqrt(n(lda,T,pol,MgO,MgO_th,E)**2-1)
            
            return freq_term/(m+AiRoots[q]*(m/2.0)**(1.0/3.0)\
                        +geom_term-pol_term+3*AiRoots[q]**2*(m/2.0)**(-1.0/3.0)/20)
                        
      lda0 = 2*math.pi*R*n(1000,T,pol,MgO,MgO_th,E)*1E7/m # initial wavelength guess based on n(1000 nm)
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
      import math
      from scipy.special import hermite
      from scipy.integrate import quad
      #def HG(p,L,x): return (2**L*math.factorial(p))**-0.5*hermite(p)(x*math.sqrt(L))*(L/math.pi)**0.25*math.exp(-L*x**2/2.0)
      def HG(p,L,x): return hermite(p)(x*math.sqrt(L))*math.exp(-L*x**2/2.0) # removed a large factor that would drop out in normalization anyway
      def HGnorm(p,L,x): return HG(p,L,x)*(4*math.pi*quad(lambda x: HG(p,L,x)**2, 0, limit*((p+1)*math.pi/L/2)**0.5)[0])**(-0.5)
      return 4*math.pi*quad(lambda x: HGnorm(p1,L1,x)*HGnorm(p2,L2,x)*HGnorm(p3,L3,x), 0, limit*((max(p1,p2,p3)+1)*math.pi/min(L1,L2,L3)/2)**0.5)[0]

def T_phasematch(L1,p1,q1,L2,p2,q2,L3,p3,q3,R,r,T0=70,MgO=0,MgO_th=5.0,E = 0,Tstep=0.1, n = nLNO1):
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

def n_eff_FSR(wl,R,T,FSR,pol,MgO=0,MgO_th=5.0,E = 0,n = nLNO1,n_disp=nLNO1_disp):
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

def n_eff(wl,R,r,T,pol,q,p,MgO=0,MgO_th=5.0,E = 0,n = nLNO1):
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


      