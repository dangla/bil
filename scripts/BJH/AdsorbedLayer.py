# Usage:
# python3 AdsorbedLayer.py

import numpy as np
import sys


class AdsorbedLayer:

    def __init__(self,T=293.15,K1=0.385e-9,K2=-0.189e-9,shuttleworth=4.97e-2,gas="water"):

        # Temperature
        self.temperature = T

        # Badam's constants
        self.K1 = K1
        self.K2 = K2

        # Shuttleworth term
        self.shuttleworth = shuttleworth

        # The nature of gas can be changed here.
        # Case of water
        if(gas == "water"):
            # 1. Surface tension
            Tc = 647.096
            self.surfacetension = 0.2358 * (1 - T/Tc)**(1.256) * (1 - 0.625*(1 - T/Tc))
        
            # 2. Molar volume of the gas
            self.molarvolume = 18.e-6
        else:
            print('gas %(name)s not available'%{'name':gas})
            sys.exit(1)
            return


    def thickness(self,hr):

        #hr = [float(i)/n for i in range(1,n)]

        t = self.K1 + self.K2 * np.log(-np.log(hr))

        #cv = Curve(hr,t,['Humidity','Thickness'])

        return t


    def poreradius(self,hr):

        # Capillary pressure
        R = 8.314         # Perfect gas constant
        molarvol = self.molarvolume
        T = self.temperature
        pc = - R * T / molarvol * np.log(hr)

        # Kelvin radius
        gamma = self.surfacetension
        rk = 2 * gamma / pc

        # Radius of pores
        rp = rk + self.thickness(hr)

        return rp


    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout
        
        # write header
        f.write('# Relative humidity  ' + 
                '  Thickness          ' + 
                '  Pore radius        ' + '\n')

        n = 100

        hr = [float(i)/n for i in range(1,n)]
        t  = self.thickness(hr)
        rp = self.poreradius(hr)

        # write columns
        col_format = " {:< 20e}"*3 + "\n"   # left-justified columns
        for y in zip(hr,t,rp):
            f.write(col_format.format(*y))

        #f.close()



import sys
import getopt
sys.path.append('../Curve')
from Curve import Curve
#===============================================================================
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'hT:1:2:s:',
                       ['help','Temperature=','K1=','K2=','sample='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)


    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            print('Usage: python3 %s -h -T<Temperature> -s<sample> -1<K1> -2<K2>' % sys.argv[0])
            print('\n')
            print('Options:')
            print('-s, --sample:          Number of sampling points')
            print('-T, --Temperature:     Temperature in K')
            print('-1, --K1:              Badam\'s constant')
            print('-2, --K2:              Badam\'s constant')
            print('-h, --help:            Display this help message')
            print('\n')
            sys.exit()
        elif opt in ('-s', '--sample'):
            sample = int(optval)
        elif opt in ('-T', '--Temperature'):
            T = int(optval)
        elif opt in ('-1', '--K1'):
            K1 = float(optval)
        elif opt in ('-2', '--K2'):
            K2 = float(optval)

    
    if(not "sample" in locals()):
        sample = 100
        
    if(not "T" in locals()):
        T = 293.15
        
    if(not "K1" in locals()):
        K1 = 0.385e-9
        
    if(not "K2" in locals()):
        K2 = -0.189e-9

    adslayer = AdsorbedLayer(T,K1,K2)

    hr = [float(i)/sample for i in range(1,sample)]
    t = adslayer.thickness(hr)
    cv = Curve(hr,t,['Humidity','Thickness'])
    cv.write()
    cv.plot()


