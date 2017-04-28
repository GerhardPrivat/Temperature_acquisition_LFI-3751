# This is supposed to be a combined WGM librabry file containing all functions from MJF, GS, DS, and others, respectively



### --------------------------------------
### ----- REFRACTIVE INDEX FUNCTIONS -----
### --------------------------------------

def n_BBO(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
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
	
def nT_BBO(lda,A,T,*args):
      """
      This is Sellmeyer eq. for BBO from Optics Communications 184 (2000) 485-491;
      With temperature coefficients from Appl.Phys.A 52, 359-368 (1991) 
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
      no = no*(1-(T-20)*16.6e-6)
      ne = ne*(1-(T-20)*9.3e-6)
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
	
def n_BBO1(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
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
	
def n_CaF2(lda,*args):
      """
      This is Sellmeyer eq. for fluorite from http://refractiveindex.info/?group=CRYSTALS&material=CaF2
      """
      import math 
      a = 0.5675888*lda*lda/(lda*lda-50.263605*50.263605)
      b = 0.4710914*lda*lda/(lda*lda-100.3909*100.3909)
      c = 3.8484723*lda*lda/(lda*lda-34649.040 *34649.040 )
      return math.sqrt(1+a+b+c)

def nT_CaF2(lda,T,*args):
      """
      This uses the n_CaF2 at 19C and dn/dT data from Corning to calculate temperature-dependent n
      """
      import math
      from lagrange_int import *
      from WGM_lib1 import n_CaF2
      n19 = n_CaF2(lda)
      wavelengths = [460.0,630.0,1150.0,3390.0]
      temperatures = range(-180,220,20)
      dn_dT = [[ -0.39, -0.40, -0.41, -0.40 ],[ -0.53, -0.54, -0.56, -0.52],[ -0.64, -0.66, -0.68, -0.63],
       [-0.74, -0.77, -0.78, -0.73],[ -0.83, -0.85, -0.87, -0.82 ],[-0.90, -0.93, -0.95, -0.89 ],
       [ -0.95, -0.99, -1.01, -0.95],[-1.00, -1.03, -1.06, -1.00],[ -1.04, -1.07, -1.10, -1.05],
       [-1.07, -1.10, -1.13, -1.09 ],[ -1.10, -1.13, -1.15, -1.12],[ -1.12, -1.15, -1.18, -1.14],
       [-1.14, -1.17, -1.20, -1.17],[ -1.16, -1.19, -1.22, -1.19],[-1.18, -1.21, -1.24, -1.21],
       [-1.20, -1.23, -1.26, -1.23],[-1.22, -1.26, -1.29, -1.25],[-1.26, -1.30, -1.32, -1.27 ],
       [-1.29, -1.34, -1.36, -1.30], [-1.34, -1.40, -1.41, -1.34 ]]
      # interpolating to given temperature wavelength lda
      dn_dT1 = []
      for x in dn_dT:
            xa,xintderiv,xerror_codex = lagrange_int(lda,4,wavelengths,x,3,0)
            dn_dT1.append(xa)
      # computing crude n(T) from derivatives
      ns = []
      for ind in range(20):
            a,b = min(ind,10),max(ind,10)
            s = 0
            if ind>10: s = 1
            if ind<10: s = -1
            ns.append(n19+1e-5*20*s*sum(dn_dT1[a:b]))    
      # interpolating n to given temperature T
      n,xintderiv,xerror_codex = lagrange_int(T,20,temperatures,ns,3,0)
      return n
	
def n_const(lda,A,*args):
    import math
    if A==0: return math.sqrt(4.66893170463187) #e
    else: return math.sqrt(5.00135666372516) #o
	
def n_fused_silica(lda,*args):
      """
      This is Sellmeyer eq. for fused silica http://refractiveindex.info/?group=GLASSES&material=F_SILICA
      """
      import math 
      a = 0.6961663*lda*lda/(lda*lda-68.4043*68.4043)
      b = 0.4079426*lda*lda/(lda*lda-116.2414*116.2414)
      c = 0.8974794*lda*lda/(lda*lda-9896.161*9896.161)
      return math.sqrt(1+a+b+c)
	
def n_diamond(lda,*args):
      """
      This is Sellmeyer eq. for diamond from from F.Peter, Z Phys 15, 358 (1923)
      """
      import math 
      a = 0.3306*lda*lda/(lda*lda-175*175)
      b = 4.3356*lda*lda/(lda*lda-106*106)
      return math.sqrt(1+a+b)
	
def n_KBBF1(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
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
	
def n_LNO1(lda,T,pol,MgO=0,MgO_th=5.0,E = 0):
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
	
def n_LNO1_disp(lda,T,pol,MgO=0,MgO_th=5.0):
      """
      dispersion dn/dlda (in nm) for Lithium Niobate dispersion n_LNO1
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
	
def n_LNO1xcut(lda,T,pol,MgO=0,MgO_th=5.0,E = 0):
	#this uses the seelmeier equation for lithium niobate in an x-cut configuration
	from WGM_lib import n_LNO1
	if pol == 0:
		neff = (n_LNO1(lda,T,0,MgO=0,MgO_th=5.0,E = 0) + n_LNO1(lda,T,1,MgO=0,MgO_th=5.0,E = 0)) / 2.
		#neff = n_LNO1(lda,T,1,MgO,MgO_th,E)
	if pol == 1:
		neff = n_LNO1(lda,T,0,MgO,MgO_th,E)
	return neff

def n_LTaO(lda,A,T,*args):
      """
      This is Sellmeyer eq. for Lithium Tantalate from J. Appl. Phys.80 (11), 1 December 1996
      "Temperature-dependent dispersion relation of ferroelectric lithium tantalate",
      Kazi Sarwar Abedin and Hiromasa Ito
      Temperature in Centigrades, wavelength in nm.
      T is the temperature in C, A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      MgO=0 MgO concentration in % MgO_th the threshold concentration. External electric field E in V/cm
      """
      import math
      A1 =4.5122, 4.5299
      A2 =0.0847522, 0.0844313
      A3 =0.19876, 0.20344
      A4 =-0.0239046, -0.0237909
      B1 =-9.66449e-9, 1.72995e-7
      B2 = 8.815e-8, -4.7733e-7
      B3 = 4.25637e-8, -8.31467e-8
      T0=25.0 # deg. C
      lda = lda*1.0e-3 # mn->microns
      def f(T): return (T-T0)*(T+T0+546)
      def g(lda,T,A2,A3,B1,B2): return (A2+B1*f(T))/(lda*lda-(A3+B2*f(T))**2)
      no = math.sqrt(A1[0]+(g(lda,T,A2[0],A3[0],B1[0],B2[0]))+B3[0]*f(T)+A4[0]*lda*lda)
      ne = math.sqrt(A1[1]+(g(lda,T,A2[1],A3[1],B1[1],B2[1]))+B3[1]*f(T)+A4[1]*lda*lda)
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
	
def n_MgF2(lda,T,pol=0,MgO=0,MgO_th=5.0,E=0):
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
	  #x = 94.61442,90.76162
      #y = 23793.604,23771.995 
      #no = math.sqrt(1+a[0]*lda*lda/(lda*lda-d[0]*d[0])+b[0]*lda*lda/(lda*lda-x[0]*x[0])+c[0]*lda*lda/(lda*lda-y[0]*y[0]))
      #ne = math.sqrt(1+a[1]*lda*lda/(lda*lda-d[1]*d[1])+b[1]*lda*lda/(lda*lda-x[1]*x[1])+c[1]*lda*lda/(lda*lda-y[1]*y[1]))
      #o = math.sin(math.pi*A/180)/no
      #e = math.cos(math.pi*A/180)/ne
      #return (o*o+e*e)**(-0.5)
	  return n
	  
def nT_MgF2(lda,A,T,*args):
      """
      This uses the n_MgF2 at 19C and dn/dT data from Corning to calculate temperature-dependent n
      """
      import math
      from lagrange_int import *
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
	
def n_quartz(lda,A,*args):
      """
      This is Sellmeyer eq. for alpha-quartz from Handbook of Optics, 3rd edition, Vol. 4. McGraw-Hill 2009;
      T is the angle (in degrees) between the optical axis and polarization:
      T = 0 extraordinary, T = 90 ordinary.
      """
      import math 
      a_ord = [0.663044,0.517852,0.175912,0.565380,1.675299]
      b_ord = [60.0,106.0,119.0,8844.0,20742.0]
      a_extr = [0.665721,0.503511,0.214792,0.539173,1.807613]
      b_extr = [60.0,106.0,119.0,8792.0,19700.0]
      def n(a,b):
            ss = 1
            for i in range(5):
                  qq = a[i]*lda*lda/(lda*lda-b[i]*b[i])
                  ss+=qq
            return math.sqrt(ss)
      no = n(a_ord,b_ord)
      ne = n(a_extr,b_extr)
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
	
def n_rutile1(lda,A,*args):
      """
      This is Sellmeyer eq. for TiO2 (rutile) from
      OPTICAL MATERIALS EXPRESS Vol. 2, No. 12, 1797 (2012) Borne et al.
      A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      """
      import math 
      a_ord = [3.2089,3.4e-5,1.227e-5,3.2545e-8]
      a_extr = [2.9713,5.1891e-5,1.228e-5,4.295e-8]
      no = math.sqrt(a_ord[0]+a_ord[1]/(a_ord[2]-lda**(-2))-a_ord[3]*lda*lda)
      ne = math.sqrt(a_extr[0]+a_extr[1]/(a_extr[2]-lda**(-2))-a_extr[3]*lda*lda)
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
	
def n_rutile2(lda,A,T,*args):
      """
      This is Sellmeyer eq. for TiO2 (rutile) from Rams, Tejeda, and Cabrera, 
      J. Appl. Phys., Vol. 82, No. 3, 1 August 1997; T in centigrades.
      A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      """
      import math 
      a_ord = [3.3427,3.119e-5,4.6173e-8,1.0407e-10]
      a_extr = [3.0993,5.0407e-5,6.344e-9,1.5444e-10]
      zzo = (291.04-0.016144*T+4.5057e-5*T*T)**(-2)-lda**(-2)
      zze = (284.19-0.01385*T+2.9614e-5*T*T)**(-2)-lda**(-2)
      no = math.sqrt(a_ord[0]+a_ord[1]/zzo-(a_ord[2]+a_ord[3]*T)*lda*lda)
      ne = math.sqrt(a_extr[0]+a_extr[1]/zze+(a_extr[2]-a_extr[3]*T)*lda*lda)
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
	
def n_test(lda,T,pol=0,MgO=0,MgO_th=5.0,E = 0):
	  return 1.5
	
def n_YVO4(lda,A,*args):
      """
      See the Laser components data sheet.
      A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      """
      import math 
      a = 3.77834, 4.59905 # ordinary, extraordinary
      b = 0.069736, 0.110534
      c = 0.04724, 0.04813
      d = 0.0108133, 0.0122676
      lda = lda/1000.0 # nm -> microns
      no = math.sqrt(a[0]-d[0]*lda*lda+b[0]/(lda*lda-c[0]))
      ne = math.sqrt(a[1]-d[1]*lda*lda+b[1]/(lda*lda-c[1]))
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
	


### --------------------------------------
### ------- EXPANSION COEFFICIENTS -------
### --------------------------------------

def alpha_CaF2(T): # Thermal expansion coefficient of fluorite from Corning data sheet
    import math
    from lagrange_int import *
    temperatures = range(-180,220,20)
    alphas=[6.7,9.1,11.1,12.8,14.1,15.2,16.2,17.0,17.7,18.3,18.7,19.1,19.4,19.7,20.0,20.4,20.8,21.3,21.7,22.2]
    alpha,xintderiv,xerror_codex = lagrange_int(190,20,temperatures,alphas,3,0)
    return 1e-6*alpha

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



### --------------------------------------
### ------- MATHEMATICAL FUNCTIONS -------
### --------------------------------------

def AngOvlp(L1,p1,L2,p2,L3,p3,limit=3):
      """
      Calculates angular overlap of three spherical functions. Default limit is sufficient for up to p = 10
      """
      import numpy
	  import math
      from scipy.special import hermite
      from scipy.integrate import quad
      #from factorial import factorial
      #def HG(p,L,x): return (2**L*math.factorial(p))**-0.5*hermite(p)(x*math.sqrt(L))*(L/math.pi)**0.25*math.exp(-L*x**2/2.0)
      #def HG(p,L,x): return hermite(p)(x*math.sqrt(L))*math.exp(-L*x**2/2.0) # removed a large factor that would drop out in normalization anyway
      # michls try
	  #def HG(p,L,x): return hermite(p)(x*math.sqrt(L))*math.exp(-L*x**2/2.0) # removed a large factor that would drop out in normalization anyway
      #def HGnorm(p,L,x): return HG(p,L,x)*(4*math.pi*quad(lambda x: HG(p,L,x)**2, 0, limit*((p+1)*math.pi/L/2)**0.5)[0])**(-0.5)
      #return 4*math.pi*quad(lambda x: HGnorm(p1,L1,x)*HGnorm(p2,L2,x)*HGnorm(p3,L3,x), 0, limit*((max(p1,p2,p3)+1)*math.pi/min(L1,L2,L3)/2)**0.5)[0]
      #def HG(p,L,x): return hermite(p)(x*numpy.sqrt(L))*numpy.exp(-L*x**2/2.0)/(numpy.sqrt(factorial(p)))
      def HG(p,L,x): return hermite(p)(x*numpy.sqrt(L))*numpy.exp(-L*x**2/2.0)/(numpy.sqrt(factorial_002(p)))
      def HGnorm(p,L,x): return HG(p,L,x)*(4*numpy.pi*quad(lambda x: HG(p,L,x)**2, 0, ((p+1)*numpy.pi/L/2)**0.5)[0])**(-0.5)
      return 4*numpy.pi*quad(lambda x: HGnorm(p1,L1,x)*HGnorm(p2,L2,x)*HGnorm(p3,L3,x), 0, ((max(p1,p2,p3)+1)*numpy.pi/min(L1,L2,L3)/2)**0.5)[0]

def AngOvlp_new(R_ratio,m1,p1,m2,p2,m3,p3,limit=3):
      """
      Calculates angular overlap of three normalized spherical functions  using Breuing et al. Opt. Expr. 21 2013.
      R_ratio = R/r. Default integration limit is sufficient for up to p = 10
      """
      if (p1+p2+p3)%2: return 0 # checks if the product of three polynomials is an odd function
      import math
      from scipy.special import hermite
      from scipy.integrate import quad
      def th0(m): return (R_ratio)**0.75/math.sqrt(m)
      def AngFunc(x,m,p): return hermite(p)(x/th0(m))*math.exp(-0.5*(x/th0(m))**2)/math.sqrt(2**(p+1)*math.factorial(p)*th0(m)*math.pi**1.5)
      def Core(x): return AngFunc(x,m1,p1)*AngFunc(x,m2,p2)*AngFunc(x,m3,p3)
      return 4*math.pi*quad(lambda x: Core(x), 0, limit*((max(p1,p2,p3)+1)*math.pi/min(m1,m2,m3)/2)**0.5)[0]

def AngOvlp4_new(R_ratio,m1,p1,m2,p2,m3,p3,m4,p4,limit=3):
      """
      Calculates angular overlap of four normalized spherical functions  using Breuing et al. Opt. Expr. 21 2013.
      R_ratio = R/r. Default integration limit is sufficient for up to p = 10
      """
      if (p1+p2+p3+p4)%2: return 0 # checks if the product of four polynomials is an odd function
      import math
      from scipy.special import hermite
      from scipy.integrate import quad
      def th0(m): return (R_ratio)**0.75/math.sqrt(m)
      def AngFunc(x,m,p): return hermite(p)(x/th0(m))*math.exp(-0.5*(x/th0(m))**2)/math.sqrt(2**(p+1)*math.factorial(p)*th0(m)*math.pi**1.5)
      def Core(x): return AngFunc(x,m1,p1)*AngFunc(x,m2,p2)*AngFunc(x,m3,p3)*AngFunc(x,m4,p4)
      return 4*math.pi*quad(lambda x: Core(x), 0, limit*((max(p1,p2,p3)+1)*math.pi/min(m1,m2,m3)/2)**0.5)[0]

def AngOvlp4(L1,p1,L2,p2,L3,p3,L4,p4,limit=3):
      """
      Calculates angular overlap of three spherical functions. Default integration limit is sufficient for up to p = 10
      """
      import math
      from scipy.special import hermite
      from scipy.integrate import quad
      def HG(p,L,x): return hermite(p)(x*math.sqrt(L))*math.exp(-L*x**2/2.0) # removed a large factor that would drop out in normalization anyway
      def HGnorm(p,L,x): return HG(p,L,x)*(4*math.pi*quad(lambda x: HG(p,L,x)**2, 0, limit*((p+1)*math.pi/L/2)**0.5)[0])**(-0.5)
      upper_limit=limit*((max(p1,p2,p3,p4)+1)*math.pi/min(L1,L2,L3,L4)/2)**0.5
      return 4*math.pi*quad(lambda x: HGnorm(p1,L1,x)*HGnorm(p2,L2,x)*HGnorm(p3,L3,x)*HGnorm(p4,L4,x), 0, upper_limit)[0]
	  
def FSR_simple(R,lda,T,pol,L,q = 1,MgO=0,MgO_th=5.0, n = n_LNO1,E = 0):
      import math
      from scipy.special import ai_zeros
      q -= 1
      AiRoots = -ai_zeros(40)[0]
      a = 29.9792/(2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E))
      b = 1+0.5*AiRoots[q]*((L/2.0+0.5)**(1.0/3.0)-(L/2.0-0.5)**(1.0/3.0))
      return a*b
	
def FSR_bulk(R,lda,T,pol,MgO=0,MgO_th=5.0,E = 0, n = n_LNO1):
      import math
      FSR = 29.9792/(2*math.pi*R*n(lda,T,pol,MgO,MgO_th,E))
      return FSR
	
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
				
def Get_frac_L_DS(R,r,lda,T,pol,q=1,p=0,n=n_LNO1,*args):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target wavelength lda
      and returs the fractional orbital number L
      MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in mm
      TM 	Er(o) pol_term numerator 1/n    pol=0 
      TE	Ez(e) pol_term numerator n       pol=1
      """
      import math
      A = 90
      if pol: A=0
      from scipy.special import ai_zeros
      q -= 1
      r = math.sqrt(r*R) # redefine the rim radius to spheroid semi-axis
      AiRoots = -ai_zeros(40)[0]
      freq_term = 2*math.pi*R*n(lda,A,T,*args)*1E7/lda
      geom_term = (2*p*(R-r)+R)/2.0/r
      pol_term = n(lda,A,T,*args)**(2*pol-1)/math.sqrt(n(lda,A,T,*args)**2-1)
      term5 = AiRoots[q]/12*(2*p*(R**3-r**3)+R**3)/r/r/r
      term5pol=AiRoots[q]/12*2*n(lda,A,T,*args)**(2*pol+1)\
                      *(2*n(lda,A,T,*args)**(4*pol-2)-4)*(n(lda,A,T,*args)**2-1)**(-3.0/2.0)
      mu = (R-r)/4.0/r
      term6 = (10-AiRoots[q]**3)/1400.0+(1+3*mu)*(2*p+1)**2*R*R*(r*r-R*R)/32.0/r/r/r/r
      L = freq_term
      dL = 1
      while math.fabs(dL) >= 1E-9:
            L1 = freq_term - AiRoots[q]*(L/2.0)**(1.0/3.0) - geom_term + pol_term\
                 - 3*AiRoots[q]**2*(L/2.0)**(-1.0/3.0)/20-(term5+term5pol)*(L/2.0)**(-2.0/3.0)-term6*2.0/L
            dL = L1-L
            L = L1
      return L
	
def Get_frac_L(R,r,lda,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = n_LNO1):
      """
      This iteratively solves the WGM dispersion equation for a spheroid at the target wavelength lda
      and returs the fractional orbital number L
      pol=0 "ordinary", pol=1 "extraordinary" MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in mm
      """
      import math
      #from WGM_lib import n_LNO1
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

def n_eff(wl,R,r,T,pol,q,p,MgO=0,MgO_th=5.0,E = 0,n = n_LNO1):
      '''
      Finds effective refraction index (material+geometrical) for a given WGM.
      R,r in cm, wl in nm. 
      '''
      import math
      from WGM_lib import Get_frac_L
      L = Get_frac_L(R,r,wl,T,pol,q,p,MgO,MgO_th,E,n)
      neff = L*wl*1e-7/(2*math.pi*R)
      return neff,L
	
def n_eff_FSR(wl,R,T,FSR,pol,MgO=0,MgO_th=5.0,E = 0,n = n_LNO1,n_disp=n_LNO1_disp):
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
	
def Phasematch(wlP,p1,q1,p2,q2,p3,q3,R,r,T,MgO=0,MgO_th=5.0,E = 0, n = n_LNO1):
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

def RadOvlp_new(R,r,m1,q1,m2,q2,m3,q3):
      """
      Calculates radial overlap of three normalized Airy functions using Breuing et al. Opt. Expr. 21 2013
      """
      import math
      from scipy.special import airy, ai_zeros
      from scipy.integrate import quad
      AiRoots =  -ai_zeros(40)[0]
      q1, q2, q3 = q1-1, q2-1, q3-1
      def xi(m): return R*(2*m*m)**(-1.0/3.0)/r
      def RadFunc(x,m,q): return airy(x/xi(m)-AiRoots[q])[0]
      N1 = R*r*r*quad(lambda x:  RadFunc(x,m1,q1)**2, 0, 1)[0]
      N2 = R*r*r*quad(lambda x:  RadFunc(x,m2,q2)**2, 0, 1)[0]
      N3 = R*r*r*quad(lambda x:  RadFunc(x,m3,q3)**2, 0, 1)[0]
      def Core(x): return RadFunc(x,m1,q1)*RadFunc(x,m2,q2)*RadFunc(x,m3,q3)
      return R*r*r*quad(lambda x: Core(x), 0, 1)[0]/math.sqrt(N1*N2*N3)

def RadOvlp4_new(R,r,m1,q1,m2,q2,m3,q3,m4,q4):
      """
      Calculates radial overlap of four normalized Airy functions using Breunig et al. Opt. Expr. 21 2013
      """
      import math
      from scipy.special import airy, ai_zeros
      from scipy.integrate import quad
      AiRoots =  -ai_zeros(40)[0]
      q1, q2, q3, q4 = q1-1, q2-1, q3-1, q4-1
      def xi(m): return R*(2*m*m)**(-1.0/3.0)/r
      def RadFunc(x,m,q): return airy(x/xi(m)-AiRoots[q])[0] # xi = u/r, c.f. Breunig et al.
      N1 = R*r*r*quad(lambda x:  RadFunc(x,m1,q1)**2, 0, 1)[0]
      N2 = R*r*r*quad(lambda x:  RadFunc(x,m2,q2)**2, 0, 1)[0]
      N3 = R*r*r*quad(lambda x:  RadFunc(x,m3,q3)**2, 0, 1)[0]
      N4 = R*r*r*quad(lambda x:  RadFunc(x,m4,q4)**2, 0, 1)[0]
      def Core(x): return RadFunc(x,m1,q1)*RadFunc(x,m2,q2)*RadFunc(x,m3,q3)*RadFunc(x,m4,q4)
      return R*r*r*quad(lambda x: Core(x), 0, 1)[0]/math.sqrt(N1*N2*N3*N4)

def RadOvlp4(R,L1,q1,L2,q2,L3,q3,L4,q4):
      """
      Calculates radial overlap of four normalized Airy functions
      """
      import math
      from scipy.special import airy, ai_zeros
      from scipy.integrate import quad
      AiRoots =  -ai_zeros(40)[0]
      qs = [q1-1, q2-1, q3-1, q4-1]
      Ls = [L1,L2,L3,L4]
      def nk(L,q): return L*(1+AiRoots[q]*(2*L**2)**(-1.0/3.0))
      def AiArg(x,L,q): return -AiRoots[q]*(L-nk(L+0.5,q)*x)/(L-nk(L+0.5,q))
      def RadFunc(x,L,q): return airy(AiArg(x,L,q))[0]/math.sqrt(x)
      L0 = min(Ls)
      q0 = max(qs)+1
      Rmin = 1-5*math.sqrt(q0)*L0**(-2.0/3.0)
      Ns = R**12
      for i in range(4):
            Ns*= quad(lambda x:  x*x*RadFunc(x,Ls[i],qs[i])**2, Rmin, 1)[0]
      def Core(x):
            C = x*x
            for i in range(4):
                  C*=RadFunc(x,Ls[i],qs[i])
            return C
      return R**3*quad(lambda x: Core(x), Rmin, 1)[0]/math.sqrt(Ns)

def SH_phasematch_blk(angle, wl0=1000,step=1.01, n = n_BBO):
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
	
def T_phasematch(L1,p1,q1,L2,p2,q2,L3,p3,q3,R,r,T0=70,MgO=0,MgO_th=5.0,n = n_LNO1,E = 0,Tstep=0.1):
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

def T_phasematch_blk(wl1,wl2, T0=70,MgO=0,MgO_th=5.0,E = 0,Tstep=0.1, n = n_LNO1):
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
	
def T_phasematch_DS(L1,p1,q1,L2,p2,q2,L3,p3,q3,R,r,T0=20,Tstep=0.1,n=n_LNO1,alphaT=1.7e-5,*args):
      '''
      Finds 1/wl2 + 1/wl2 = 1/wl3, L1+L2 -> L3 phase matching temperature down to the accuracy deltaf (GHz)
      based on the initial guess T0 with initial search step Tstep. R is at the returned T0!
      Requires importing Sellmeyer equations. This program is unaware of the orbital selection rules!
      '''
      import math
      from WGM_lib1 import WGM_freq
      def Df(T): # freqency detuning in GHz
            RR = R*(1+alphaT*(T-T0))
            rr = r*(1+alphaT*(T-T0))
            lda_1, freq_1 = WGM_freq(RR,rr,L1,T,0,q1,p1,n,*args)
            lda_2, freq_2 = WGM_freq(RR,rr,L2,T,0,q2,p2,n,*args)
            lda_3, freq_3 = WGM_freq(RR,rr,L3,T,1,q3,p3,n,*args)
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
	
def WGM_FSR(R,r,L,T,pol,q,p ,MgO,MgO_th=5.0,E = 0, n = n_LNO1):
     wlmode_minus, frmode_minus = WGM_freq(R,r,L-1.,T,pol,q,p,MgO,MgO_th,E,n) 
#     wlmode, frmode = WGM_freq(R,r,L,T,pol,q,p,MgO,MgO_th,E, n) 
     wlmode_plus, frmode_plus   = WGM_freq(R,r,L+1.,T,pol,q,p,MgO,MgO_th,E,n) 
     FSR = (frmode_plus - frmode_minus ) / 2.
     return FSR
	
def WGM_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = n_LNO1):
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
      #from WGM_lib import n_LNO1
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
	
def WGMxcut_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = n_LNO1):
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
      #from WGM_lib import n_LNO1
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
      
def WGM_freq_m(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,E = 0, n = n_LNO1):
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
      #from WGM_lib import n_LNO1
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
	


### --------------------------------------
### -------- INDIVIDUAL FUNCTIONS --------
### --------------------------------------

### --------- DMITRY'S FUNCTIONS ---------
### --------------------------------------
def PDC_findR_eff(wlp,signal,T,R_guess,R_range,qp,polp,qs,pols,qi,poli,Rratio,n=n_LNO1,*args):
      '''
      Finds the resonator radius for PDC phase-matching using effective refractive indices,
      around a guessed value, within R_range = [min,max]. All equatorial modes are assumed (p=0)
      '''
      import math
      from WGM_lib1 import n_eff
      idler = 1.0/(1.0/wlp-1.0/signal)
      Ap, As, Ai = 90,90,90
      if polp: Ap=0
      if pols: As=0
      if poli: Ai=0
      np = n(wlp,Ap,T)
      ns = n(signal,As,T)
      ni = n(idler,Ai,T)
      def npf(*args): return np # This way Sellmeyer eqs do not have to be evaluated at each step
      def nsf(*args): return ns
      def nif(*args): return ni
      def Dn(R): # index mismatch
          return  n_eff(wlp,R,R/Rratio,T,polp,qp,0,npf)[0]-wlp*n_eff(signal,R,R/Rratio,T,pols,qs,0,nsf)[0]/signal-wlp*n_eff(idler,R,R/Rratio,T,poli,qi,0,nif)[0]/idler
      R=R_guess
      R_step = 0.1*R
      dn = Dn(R)
      while math.fabs(dn)>1E-8: # delta n accuracy
            #print "R = ",R," dn = ", dn
            dn1 = Dn(R+R_step)
            ddndwl = (dn1-dn)/R_step
            dR = dn/ddndwl
            R-=dR
            dn = Dn(R)
            R_step = min(0.1*R,math.fabs(dR)) 
            if R <R_range[0]: return -1
            if R >R_range[1]: return -1
      return R
    
def PDC_findT_eff(wlp,signal,R,r,T_guess,T_range,qp,polp,qs,pols,qi,poli,alphaT,n=n_LNO1,*args):
      '''
      Finds PDC phase-matching temperature using effective refractive indices,
      around a guessed value, within T_range = [min,max]. All equatorial modes (p=0)
      R is at 20 C! alphaT is thermal expansion.
      '''
      import math
      from WGM_lib1 import n_eff
      def RT(t,x): return x*(1+alphaT*(t-20))
      idler = 1.0/(1.0/wlp-1.0/signal)
      def Dn(T): # index mismatch
          R1 = RT(T,R)
          r1 = RT(T,r)
          return n_eff(wlp,R1,r1,T,polp,qp,0,n)[0]-wlp*n_eff(signal,R1,r1,T,pols,qs,0,n)[0]/signal-wlp*n_eff(idler,R1,r1,T,poli,qi,0,n)[0]/idler
      T=T_guess
      T_step = 20
      dn = Dn(T)
      while math.fabs(dn)>1E-8: # delta n accuracy
            #print "R = ",R," dn = ", dn
            dn1 = Dn(T+T_step)
            ddndwl = (dn1-dn)/T_step
            dT = dn/ddndwl
            if math.fabs(dT)<math.fabs(T_step): T_step = dT
            else: T_step = 1.5*T_step*dT/math.fabs(dT) # increasing step is useful to get out of a loop
            T-=T_step
            dn = Dn(T)
            T_step = min(20,math.fabs(dT)) 
            if T <T_range[0]: return ""
            if R >T_range[1]: return ""
      if T: return T
      else: return 1e-9 # never return 0: can use boolean to identify no-solution.

def PDCwl_eff(lambda_p,wls_guess,wl_step,qp,polp,qs,pols,qi,poli,R,r,T,n=n_LNO1,*args):
      '''
      Finds phase-matched signal and idler wavelengths for a given pump using effective refractive indices,
      around a guessed value.
      All equatorial mods are assumed (p=0) THIS CAN MISS MULTIPLE SOLUTIONS!!
      '''
      import math
      from WGM_lib1 import n_eff
      def Dn(wls): # index mismatch
          wli=1.0/(1.0/lambda_p-1.0/wls)
          return n_eff(lambda_p,R,r,T,polp,qp,0,n)[0]-lambda_p*n_eff(wls,R,r,T,pols,qs,0,n)[0]/wls-lambda_p*n_eff(wli,R,r,T,poli,qi,0,n)[0]/wli
      dn = Dn(wls_guess)
      while math.fabs(dn)>1E-8: # delta n accuracy
            #print "wavelength = ",wls_guess+wl_step," dn = ", dn," wl step = ", wl_step
            try: dn1 = Dn(wls+wl_step)
            except: return 0,0
            dwl = wl_step*dn/(dn1-dn)
            if math.fabs(dwl)<math.fabs(wl_step): wl_step = dwl
            else: wl_step = 1.5*wl_step*dwl/math.fabs(dwl) # increasing step is useful to get out of a loop
            wls -= wl_step
            try: dn = Dn(wls)
            except: return 0,0
            if wls<lambda_p: return 0,0
      return wls,1.0/(1.0/lambda_p-1.0/wls)
    
def PDCwl_eff1(wl_p,qp,polp,qs,pols,qi,poli,R,r,T,n=n_LNO1,wl_p_offset=300, wl_step=100,*args):
      '''
      Finds phase-matched signal and idler wavelengths for a given pump using effective refractive indices.
      All equatorial mods are assumed (p=0). Search range is limited to above wl_p+wl_p_offset for both the signal and idler.
      '''
      import math
      from WGM_lib1 import n_eff
      def idlr(wls): return 1.0/(1.0/wl_p-1.0/wls)
      def Dn(wls): return n_eff(wl_p,R,r,T,polp,qp,0,n)[0]-wl_p*n_eff(wls,R,r,T,pols,qs,0,n)[0]/wls-wl_p*n_eff(idlr(wls),R,r,T,poli,qi,0,n)[0]/idlr(wls)
      def DDn(wls,step): return (Dn(wls+step)-Dn(wls))/step
      def find_signal(wls,step): # finds the closest solution 
        dn = Dn(wls)
        while math.fabs(dn)>1E-8: # delta n accuracy equiv ~ 1 MHz
            try: ddndwl = DDn(wls,step)
            except: return 0
            dwl = dn/ddndwl
            if math.fabs(dwl)<math.fabs(step): step = dwl
            else: step = 1.5*step*dwl/math.fabs(dwl) # increasing step is useful to get out of a loop
            wls -= step
            try: dn = Dn(wls)
            except: return 0
            if wls<wl_p: return 0
        return wls
      signals, idlers = [],[]
      step = wl_step
      wls = wl_p+wl_p_offset
      while Dn(wls)*DDn(wls,step)>0:
          wls+=step
          if idlr(wls)-wl_p_offset<wl_p: return signals, idlers # delta n never started to approach zero
      while idlr(wls)-wl_p_offset>wl_p: # main loop
          while Dn(wls)*Dn(wls+step)>0: # and DDn(wls,step)*DDn(wls+step,step)>0 :
              wls+=step
              if idlr(wls)-wl_p_offset<wl_p: return signals, idlers # delta n did not reached zero
          s = find_signal(wls,step)
          try: i = idlr(s)
          except:
              i = 0
              print "Solution exists between signal = ", wls, " nm and ", wls+wl_step, " nm, but I cannot find it!"
          signals.append(s)
          idlers.append(i)
          step = wl_step # resetting step after each found solution
          while Dn(s-step)*Dn(s+step)>0: # next solution is within step.
               step=step/2.0
          wls = s + step 
      return signals, idlers
	
def SH_findR_eff(wlp,T,R_guess,R_range,qp,polp,qs,pols,Rratio,n=n_LNO1,*args):
      '''
      Finds the resonator radius for SHG phase-matching using effective refractive indices,
      around a guessed value, within R_range = [min,max]. All equatorial modes are assumed (p=0)
      wlp is the pump (fundamental) wavelength
      '''
      import math
      from WGM_lib1 import n_eff
      Ap, As = 90,90
      if polp: Ap=0
      if pols: As=0
      np = n(wlp,Ap,T)
      ns = n(wlp/2.0,As,T)
      def npf(*args): return np # This way Sellmeyer eqs do not have to be evaluated at each step
      def nsf(*args): return ns
      def Dn(R): # index mismatch
          return  n_eff(wlp,R,R/float(Rratio),T,polp,qp,0,npf)[0]-n_eff(wlp/2.0,R,R/float(Rratio),T,pols,qs,0,nsf)[0]
      R=R_guess
      R_step = 0.1*R
      dn = Dn(R)
      while math.fabs(dn)>1E-8: # delta n accuracy
            #print "R = ",R," dn = ", dn
            dn1 = Dn(R+R_step)
            ddndwl = (dn1-dn)/R_step
            dR = dn/ddndwl
            R-=dR
            dn = Dn(R)
            R_step = min(0.1*R,math.fabs(dR)) 
            if R <R_range[0]: return -1
            if R >R_range[1]: return -1
      return R
	
def SHwl_eff(wlp_guess,wlp_range,wl_step,qp,polp,qs,pols,R,r,T,n=n_LNO1,*args):
      '''
      Finds phase-matched fundamental wavelengths for SHG using effective refractive indices,
      around a guessed value, within wlp_range = [min,max]. All equatorial modes are assumed (p=0)
      '''
      import math
      from WGM_lib1 import n_eff
      def Dn(wls): # index mismatch
          return  n_eff(2*wls,R,r,T,polp,qp,0,n)[0]-n_eff(wls,R,r,T,pols,qs,0,n)[0]
      wls=wlp_guess/2.0
      dn = Dn(wls)
      while math.fabs(dn)>1E-8: # delta n accuracy
            #print "wavelength = ",wls," dn = ", dn
            dn1 = Dn(wls+wl_step)
            ddndwl = (dn1-dn)/wl_step
            dwl = dn/ddndwl
            wls-=dwl
            if math.fabs(dwl)<math.fabs(wl_step): wl_step = math.fabs(dwl)
            if wls <wlp_range[0]/2.0: return 0
            if wls >wlp_range[1]/2.0: return 0
            try: dn = Dn(wls)
            except: return 0
      return 2*wls
	
def THG_phasematch(L1,pol1,pol3,p1,q1,p3,q3,n,R0=0.1,RRatio=6.0, *args):
      '''
      Finds the phase matched pump and third-harmonic Radius 
      '''
      import math
      from WGM_lib1 import WGM_freq
      def Df(R): # freqency detuning in GHz
            lda_1, freq_1 = WGM_freq(R,R/RRatio,L1,T,pol1,q1,p1,n,*args)
            lda_3, freq_3 = WGM_freq(R,R/RRatio,L3,T,pol3,q3,p3,n,*args)
            return 3*freq_3-freq_1
      R = R0
      df = Df(R)
      Rstep = R0*1e-2
      while math.fabs(df)>1E-6: #1 kHz accuracy
            #print "T = ",T," df = ", df, 
            df1 = Df(R+Rstep)
            ddfdt = (df1-df)/Rstep
            dR = df/ddfdt
            R -= dR
            Rstep = math.fabs(dR/2.0)   
            df = Df(R)
            #print " dT = ", dT, " T -> ",T, " df -> ", df
            if R<1e-3 or R> 1: break
      return R
	

### -------- GERHARD'S FUNCTIONS ---------
### --------------------------------------
def autocorr(x):
    result = np.correlate(x, x, mode='full')
    #result = np.convolve(x, x, mode='same')
    #result = np.convolve(x, x, mode='full')
    return result[result.size/2.:]

def Find_modeinsweep(fsweep,Lstart,R,r,lda,T,pol,q,p,MgO,MgO_th=5.0,E = 0, n = n_LNO1):
      """
      This iteratively reduces L until mode is in sweep range
      """
      import math
      wl, f = WGM_freq(R,r,Lstart,T,pol,q,p,MgO,MgO_th,E,n)
      FSR = WGM_FSR(R,r,Lstart,T,pol,q,p,MgO,MgO_th,E,n)
      #print 'n, FSR',n, FSR
      fdiff = f- fsweep
      f1 = f
      L1 = Lstart
      #x = 1
      while math.fabs(fdiff) >= FSR:
      #while math.fabs(4) >= x:
            L1 -= 1
            wl1, f1 = WGM_freq(R,r,L1,T,pol,q,p,MgO,MgO_th,E,n)
            fdiff = f1- fsweep
            #print 'f diff', fdiff, ' GHz'
            #x += 1
      return f1, L1
	
def gauss(x, *p):
    A, mu, sigma, y0 = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + y0

def getdiskparam(disknum=111):
    """
        dir=0 "ordinary" - big radius, dir=1 "extraordinary" - small radius
    """
    if disknum==103:
        from WGM_lib import n_LNO1
        R0_mm = 1.6
        r0_mm = 1./6.
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==111:
        from WGM_lib import n_LNO1
        R0_mm = 1.0
        r0_mm = 1./6.
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==106:
        from WGM_lib import n_LNO1
        R0_mm = 1.27
        #R0_mm = 0.1626
        r0_mm = R0_mm
        #r0_mm = 0.97
        MgO = 5.4 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]

    if disknum==666:
        from WGM_lib import n_LNO1
        R0_mm = 1.0
        r0_mm = R0_mm/6.
        r0_mm
        MgO = 5.7 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==966:
        from WGM_lib import n_LNO1
		#for JmO long paper
        R0_mm = 1.5
        r0_mm = .6
        MgO = 5.0
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    if disknum==9191:
        from WGM_lib import n_LNO1
		#for JmO long paper
        R0_mm = 2.4978
        r0_mm = .581000
        r0_mm = .65
        MgO = 5.0 #Mg conc
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
        
    if disknum==919:
        from WGM_lib import n_LNO1
        #R0_mm = 2.4601  MgO = 5.3 #Mg conc
		#MgO = 5.4 #Mg conc  R0_mm = 2.500
        #R0_mm = 2.499
        #R0_mm = 2.489
        #R0_mm = 2.486

        R0_mm = 2.4978
        #R0_mm = 2.
        #r0_mm = .581000
        r0_mm = .65
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]

    if disknum==966:
        from WGM_lib import n_LNO1
        R0_mm = 2.5
        r0_mm = .5
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]

    if disknum==977:
        #xcut disk
        from WGM_lib import n_LNO1
        R0_mm = 2.5
        r0_mm = .5
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        #tilted extension coefficients
        lamdir0 = [stretchx* (13.31  *10**(-6)+3.828 *10**(-6))/2., stretchx2*(0.03059*10**(-6)+stretchx2*0.00814*10**(-6))] 
        lamdir1 = [stretchx*  13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        
        
    if disknum==9193:
        from WGM_lib import n_LNO1

        R0_mm = 2.4978

        r0_mm = .1
        r0_mm = .581000
        r0_mm = 1.5
        r0_mm = R0_mm
        r0_mm = 4.
        
        MgO = 5.0 #Mg conc
        MgO = 5.8
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
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
        from WGM_lib import n_LNO1
        R0_mm = 2.42
        r0_mm = R0_mm/5.
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
                
    if disknum==2604:
        from WGM_lib import n_LNO1
        R0_mm = 0.75
        r0_mm = R0_mm/6.1
        MgO = 5.0 #Mg conc
        MgO_th = 5.0 #Mg threshold
        material = n_LNO1
        stretchx =  1
        stretchx2 = 1
        lamdir0 = [stretchx* 13.31  *10**(-6), stretchx2*0.03059*10**(-6)]
        lamdir1 = [stretchx*  3.828 *10**(-6), stretchx2*0.00814*10**(-6)]
        
    return R0_mm, r0_mm, material, MgO, MgO_th, lamdir0, lamdir1

def Get_Rfromf(fmeas,Lmeas,R,r,lda,T,pol,q,p,MgO,MgO_th=5.0,E = 0, n = n_LNO1):
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
           #print 'R loop is', R1, 'delta R is', deltaR
           R1 = R1 - deltaR					
      return R1, fmism

def lambda2nu(x):
    #for example: lambda in nm to frequency in GHz
    cms = 299792458 # speed of light
    return cms/x

def lorentz(x=range(1,100),x0=0,gamma=1,amp=1,y0=1):
    return amp / (1.0 + 4*((x - x0) / gamma)**2 ) + y0

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
	
def nonorthprism(alpha,beta,d,angtype='deg'): 
	#calculates distances in a triangle with given baseline and its two adjacent angles
	import math
	if angtype=='deg':
		alpha = alpha *math.pi / 180
		beta = beta *math.pi / 180
	h = d / (1 / math.tan(alpha) + 1 / math.tan(beta))	
	d1 = h / math.tan(alpha)
	d2 = h / math.tan(beta)
	return h, d1, d2

def nu2lambda(x):
    #for example: lambda in nm to frequency in GHz
    cms = 299792458 # speed of light
    return cms/x

def path_leaf(path):
    head, tail = ntpath.split(path)
    tail =  ('.').join(tail.split('.')[:-1])
    return tail or ntpath.basename(head)

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

def poly_iter(A, x):
    p = 0
    xn = 1
    for a in A:
        p += xn * a
        xn *= x
    return p

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
                
	#return xfreq_GHz, yfreq, expfilename
    return  x1, y1, y2, y3, y4, expfilename              

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

def tiltcorrect(ydata,windowlen=100,degree = 1):
    #import poly_iter
    x_ground = np.arange(len(ydata))
    x_ground_fit = smooth(x_ground,window_len=int(len(x_ground)/windowlen),window='hanning')
    ydata_fit    = smooth(ydata   ,window_len=int(len(ydata)   /windowlen),window='hanning')
    coeff = np.polyfit(x_ground_fit,ydata_fit, degree)
    
    #if degree == 1:
    #    ytilt = tuple(x*coeff[0] + coeff[1] for x in x_ground)
    #if degree == 2:
    #    ytilt = tuple(x*x*coeff[0] + x*coeff[1] + coeff[2]  for x in x_ground)
    #if degree == 3:
    #    ytilt = tuple(x*x*x*coeff[0] + x*x*coeff[1] + x*coeff[2] + coeff[3]  for x in x_ground)
    ytilt = np.array([], dtype='double')
    for xx in x_ground:
        ytilt = np.append(ytilt,poly_horner(np.flipud(coeff),xx))
    return  ydata - ytilt

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
      #from WGM_lib import n_LNO1
      from scipy.special import ai_zeros
      q -= 1
      AiRoots =  -ai_zeros(40)[0] # >0
      geom_term = p*(math.sqrt(Roverr)-1) + math.sqrt(Roverr)/2.0
      freq_term = 1
      pol_term = n**(-1+2*pol)/math.sqrt(n**2-1)
      
      f = freq_term * (L + AiRoots[q]*(L/2.0)**(1.0/3.0) + geom_term - pol_term + 3*AiRoots[q]**2*(L/2.0)**(-1.0/3.0)/20)
      
      return f
	

### --------- MICHL'S FUNCTIONS ----------
### --------------------------------------
def alert(text):
	print '---------------------------'
	print text
	print '---------------------------\n'
	
def coupling_Q_factor(radius,lda,distance,T=60,pol=0,MgO=5.0,MgO_th=5.0,nprism=ndiamond,ndisk=n_LNO1):
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
	
def show_func():
      print 'ntest(lda,T,pol=0,MgO=0,MgO_th=5.0,n = n_LNO1,E = 0)'
      print 'n_LNO1(lda,T,pol,MgO=0,MgO_th=5.0,n = n_LNO1,E = 0)'
      print 'nBBO(lda,T,pol=0,MgO=0,MgO_th=5.0,n = n_LNO1,E=0)'
      print 'nBBO1(lda,T,pol=0,MgO=0,MgO_th=5.0,n = n_LNO1,E=0)'
      print 'FSR_simple(R,lda,T,pol,L,q = 1,MgO=0,MgO_th=5.0,n = n_LNO1,E = 0)'
      print 'WGM_freq(R,r,L,T,pol,q = 1,p = 0,MgO=0,MgO_th=5.0,n = n_LNO1,E = 0)'
      print 'ndiamond(lda,T=20,pol=0)'
      print 'coupling_Q_factor(refr_diamond,refr_disk,radius,wavelength)'


###
### --------------------------------------
