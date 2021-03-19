import numpy as np


class AdsorbedLayer:

    def __init__(self,T=293.15,K1=0.385e-9,K2=-0.189e-9,shuttleworth=4.97e-2):

        # Temperature
        self.temperature = T

        # Badam's constants
        self.K1 = K1
        self.K2 = K2

        # Surface tension
        Tc = 647.096
        self.surfacetension = 0.2358 * (1 - T/Tc)**(1.256) * (1 - 0.625*(1 - T/Tc))

        # Shuttleworth term
        self.shuttleworth = shuttleworth


    def thickness(self,hr):

        t = self.K1 + self.K2 * np.log(-np.log(hr))

        return t


    def poreradius(self,hr):

        # Capillary pressure
        R = 8.314         # Perfect gas constant
        molarvol = 18.e-6 # Molar volume of water
        T = self.temperature
        pc = - R * T / molarvol * np.log(hr)

        # Kelvin radius
        gamma = self.surfacetension
        rk = 2 * gamma / pc

        # Radius of (entry) pores
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
#===============================================================================
if __name__=="__main__":

    adslayer = AdsorbedLayer()
    adslayer.write()


