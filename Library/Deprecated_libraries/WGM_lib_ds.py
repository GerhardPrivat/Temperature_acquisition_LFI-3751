#!/usr/bin/env mpython_q
def n_const(lda,A,*args):
    import math
    if A==0: return math.sqrt(4.66893170463187) #e
    else: return math.sqrt(5.00135666372516) #o
    
def nLNO1(lda,A,T,MgO=0,MgO_th=5.0,*args):
      """
      This is Sellmeyer eq. for Lithium Niobate from "Influence of the defect structure on the refractive indices
      of undoped and Mg-doped lithium niobate", U.Schlarb and K.Betzler, Phys. Rev. B 50, 751 (1994).
      Temperature in Centigrades, wavelength in nm.
      T is the temperature in C, A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      MgO=0 MgO concentration in % MgO_th the threshold concentration. External electric field E in V/cm
      """
      import math
      w0 = 223.219,218.203 # ordinary, extraordinary
      mu = 1.1082E-6,6.4047E-6
      A0 = 4.5312E-5,3.9466E-5
      A_NbLi = -1.4464E-8,23.727E-8
      A_Mg = -7.3548E-8,7.6243E-8
      A_IR = 3.6340E-8,3.0998E-8
      A_UV = 2.6613,2.6613
      B1 = A0[0]+A_Mg[0]*MgO,A0[1]+A_Mg[1]*MgO
      B = B1
      if MgO < MgO_th: B=B1[0]+(MgO_th-MgO)*A_NbLi[0],B1[1]+(MgO_th-MgO)*A_NbLi[1]
      def f(T): return (T + 273)**2 + 4.0238E5*(1/math.tanh(261.6/(T + 273)) - 1)
      no = math.sqrt(B[0]/((w0[0]+(f(T)-f(24.5))*mu[0])**(-2)-lda**(-2))-A_IR[0]*lda**2+A_UV[0])
      ne = math.sqrt(B[1]/((w0[1]+(f(T)-f(24.5))*mu[1])**(-2)-lda**(-2))-A_IR[1]*lda**2+A_UV[1])
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)
    
def nLTaO(lda,A,T,*args):
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

def nBBO(lda,A,*args):
      """
      This is Sellmeyer eq. for BBO from http://refractiveindex.info/
      (Handbook of Optics, 3rd edition, Vol. 4. McGraw-Hill 2009);
      Appl.Phys.A 52, 359-368 (1991) 
      A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      """
      import math 
      a = 2.7405,2.3730 # ordinary, extraordinary
      b = 0.0184,0.0128
      c = 0.0179,0.0156
      d = 0.0155,0.0044
      lda = lda/1000.0 # nm -> microns
      no = math.sqrt(a[0]-d[0]*lda*lda+b[0]/(lda*lda-c[0]))
      ne = math.sqrt(a[1]-d[1]*lda*lda+b[1]/(lda*lda-c[1]))
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)

def nYVO4(lda,A,*args):
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

def nBBO1(lda,A,*args):
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
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
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

def n_diamond(lda,*args):
      """
      This is Sellmeyer eq. for diamond from from F.Peter, Z Phys 15, 358 (1923)
      """
      import math 
      a = 0.3306*lda*lda/(lda*lda-175*175)
      b = 4.3356*lda*lda/(lda*lda-106*106)
      return math.sqrt(1+a+b)

def n_CaF2(lda,*args):
      """
      This is Sellmeyer eq. for fluorite from http://refractiveindex.info/?group=CRYSTALS&material=CaF2
      """
      import math 
      a = 0.5675888*lda*lda/(lda*lda-50.263605*50.263605)
      b = 0.4710914*lda*lda/(lda*lda-100.3909*100.3909)
      c = 3.8484723*lda*lda/(lda*lda-34649.040 *34649.040 )
      return math.sqrt(1+a+b+c)

def n_MgF2(lda,A,*args):
      """
      This is Sellmeyer eq. for magnesium fluoride from
      http://refractiveindex.info/?group=CRYSTALS&material=MgF2 at T=19C
      A is the angle (in degrees) between the optical axis and polarization:
      A = 0 extraordinary, A = 90 ordinary.
      """
      import math 
      a = 0.48755108,0.41344023 # ordinary, extraordinary
      b = 0.39875031,0.50497499
      c = 2.3120353,2.4904862
      d = 43.38408,36.84262 
      x = 94.61442,90.76162
      y = 23793.604,23771.995 
      no = math.sqrt(1+a[0]*lda*lda/(lda*lda-d[0]*d[0])+b[0]*lda*lda/(lda*lda-x[0]*x[0])+c[0]*lda*lda/(lda*lda-y[0]*y[0]))
      ne = math.sqrt(1+a[1]*lda*lda/(lda*lda-d[1]*d[1])+b[1]*lda*lda/(lda*lda-x[1]*x[1])+c[1]*lda*lda/(lda*lda-y[1]*y[1]))
      o = math.sin(math.pi*A/180)/no
      e = math.cos(math.pi*A/180)/ne
      return (o*o+e*e)**(-0.5)

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

def alpha_CaF2(T): # Thermal expansion coefficient of fluorite from Corning data sheet
    import math
    from lagrange_int import *
    temperatures = range(-180,220,20)
    alphas=[6.7,9.1,11.1,12.8,14.1,15.2,16.2,17.0,17.7,18.3,18.7,19.1,19.4,19.7,20.0,20.4,20.8,21.3,21.7,22.2]
    alpha,xintderiv,xerror_codex = lagrange_int(190,20,temperatures,alphas,3,0)
    return 1e-6*alpha
    
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
    
def n_fused_silica(lda,*args):
      """
      This is Sellmeyer eq. for fused silica http://refractiveindex.info/?group=GLASSES&material=F_SILICA
      """
      import math 
      a = 0.6961663*lda*lda/(lda*lda-68.4043*68.4043)
      b = 0.4079426*lda*lda/(lda*lda-116.2414*116.2414)
      c = 0.8974794*lda*lda/(lda*lda-9896.161*9896.161)
      return math.sqrt(1+a+b+c)

def FSR_simple(R,lda,A,T,L,q=1,n=nLNO1,*args): # returns FSR in GHz
      import math
      from scipy.special import ai_zeros
      q -= 1
      AiRoots = -ai_zeros(40)[0]
      a = 29.9792/(2*math.pi*R*n(lda,A,T,*args))
      b = 1+0.5*AiRoots[q]*((L/2.0+0.5)**(1.0/3.0)-(L/2.0-0.5)**(1.0/3.0))
      return a*b

def Get_frac_L(R,r,lda,T,pol,q=1,p=0,n=nLNO1,*args):
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
def WGM_freq(R,r,L,T,pol,q = 1,p = 0,n=nLNO1,*args):
      """
      This iteratively solves the WGM dispersion equation from
      Gorodetsky, Demchenko Proc. of SPIE Vol. 8236 823623-1 (2012)
      for a spheroid at the target
      orbital number L and returs the wavelength (nm) and frequency (GHz)
      MgO =0 MgO concentration in %
      MgO_th the threshold concentration. External electric field E in V/cm
      R, r are big and small radia in mm
      TM 	Er(o) pol_term numerator n*\chi=1/n    pol=0 
      TE	Ez(e) pol_term numerator n*\chi=n       pol=1
      """
      import math
      A = 90
      if pol: A=0 #e,z
      from scipy.special import ai_zeros
      q -= 1
      r = math.sqrt(r*R) # redefine the rim radius to spheroid semi-axis
      AiRoots =  -ai_zeros(40)[0] # >0
      geom_term = (2*p*(R-r)+R)/2.0/r
      term4 = 3*AiRoots[q]**2*(L/2.0)**(-1.0/3.0)/20
      term5 = AiRoots[q]/12*(L/2.0)**(-2.0/3.0)*(2*p*(R**3-r**3)+R**3)/r/r/r
      mu = (R-r)/4.0/r
      term6 = 2.0/L*((10-AiRoots[q]**3)/1400.0+(1+3*mu)*(2*p+1)**2*R*R*(r*r-R*R)/32.0/r/r/r/r)
      def nm2GHz(lda): return 2.99792E8/lda
      def wl1(lda):
            freq_term = 2*math.pi*R*n(lda,A,T,*args)*1E7
            pol_term = n(lda,A,T,*args)**(2*pol-1)/math.sqrt(n(lda,A,T,*args)**2-1)
            term5pol=AiRoots[q]/12*(L/2.0)**(-2.0/3.0)*2*n(lda,A,T,*args)**(2*pol+1)\
                      *(2*n(lda,A,T,*args)**(4*pol-4)-3)*(n(lda,A,T,*args)**2-1)**(-3.0/2.0)
            rhs = L+AiRoots[q]*(L/2.0)**(1.0/3.0)+geom_term-pol_term+term4+term5+term5pol+term6
            return freq_term/rhs
      lda0 = 2*math.pi*R*n(1000,A,T,*args)*1E7/L # initial wavelength guess based on n(1000 nm)
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
      Calculates radial overlap of three normalized Airy functions
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

def AngOvlp(L1,p1,L2,p2,L3,p3,limit=3):
      """
      Calculates angular overlap of three normalized spherical functions.
      Default integration limit is sufficient for up to p = 10
      """
      import math
      from scipy.special import hermite
      from scipy.integrate import quad
      #def HG(p,L,x): return (2**L*math.factorial(p))**-0.5*hermite(p)(x*math.sqrt(L))*(L/math.pi)**0.25*math.exp(-L*x**2/2.0)
      def HG(p,L,x): return hermite(p)(x*math.sqrt(L))*math.exp(-L*x**2/2.0) # removed a large factor that would drop out in normalization anyway
      def HGnorm(p,L,x): return HG(p,L,x)*(4*math.pi*quad(lambda x: HG(p,L,x)**2, 0, limit*((p+1)*math.pi/L/2)**0.5)[0])**(-0.5)
      return 4*math.pi*quad(lambda x: HGnorm(p1,L1,x)*HGnorm(p2,L2,x)*HGnorm(p3,L3,x), 0, limit*((max(p1,p2,p3)+1)*math.pi/min(L1,L2,L3)/2)**0.5)[0]

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

##############  In progress
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
###############################
def T_phasematch(L1,p1,q1,L2,p2,q2,L3,p3,q3,R,r,T0=20,Tstep=0.1,n=nLNO1,alphaT=1.7e-5,*args):
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
def SHwl_eff(wlp_guess,wlp_range,wl_step,qp,polp,qs,pols,R,r,T,n=nLNO1,*args):
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

def SH_findR_eff(wlp,T,R_guess,R_range,qp,polp,qs,pols,Rratio,n=nLNO1,*args):
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

def PDC_findR_eff(wlp,signal,T,R_guess,R_range,qp,polp,qs,pols,qi,poli,Rratio,n=nLNO1,*args):
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
    
def PDC_findT_eff(wlp,signal,R,r,T_guess,T_range,qp,polp,qs,pols,qi,poli,alphaT,n=nLNO1,*args):
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


def PDCwl_eff(lambda_p,wls_guess,wl_step,qp,polp,qs,pols,qi,poli,R,r,T,n=nLNO1,*args):
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
    
def PDCwl_eff1(wl_p,qp,polp,qs,pols,qi,poli,R,r,T,n=nLNO1,wl_p_offset=300, wl_step=100,*args):
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
      
def T_phasematch_blk(wl1,wl2,T0=70,Tstep=0.1, n = nLNO1,*args):
      '''
      Finds 1/wl2 + 1/wl2 = 1/wl3 (wl in nanometers), collinear bulk phase matching temperature 
      based on the initial guess T0 with initial search step Tstep
      Requires importing Sellmeyer equations. 
      '''
      import math
      wl3 = 1/(1.0/wl1+1.0/wl2)
      def Dk(T):
            k1 = 2*math.pi*n(wl1,90,T,*args)*1e7/wl1 # cm^-1
            k2 = 2*math.pi*n(wl2,90,T,*args)*1e7/wl2
            k3 = 2*math.pi*n(wl3,0,T,*args)*1e7/wl3
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

def SH_phasematch_blk(angle, wl0=1000,step=1.01, n = nBBO,*args):
      '''
      Finds the fundamental wavelength (in nm) for frequency doubling o+o->e at a given angle between the optical axis and the e.
      Requires importing Sellmeyer equations. 
      '''
      import math
      def Dk(wl):
            k1 = 2*math.pi*n(wl,90,*args)*1e7/wl # ordinary cm^-1
            k2 = 4*math.pi*n(wl/2.0,angle,*args)*1e7/wl # extraordinary with angle
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

def n_eff(wl,R,r,T,pol,q,p,n=nLNO1,*args):
      '''
      Finds effective refraction index (material+geometrical) for a given WGM.
      R,r in cm, wl in nm. pol=0 "ordinary", pol=1 "extraordinary"
      '''
      import math
      from WGM_lib1 import Get_frac_L
      L = Get_frac_L(R,r,wl,T,pol,q,p,n,*args)
      neff = (L-p)*wl*1e-7/(2*math.pi*R)
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


      