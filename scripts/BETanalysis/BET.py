# Not yet finished!
# Usage:
# python BET.py

import numpy as np
import sys
sys.path.append('../tools')
from Curve import Curve

class BET:
    
    def __init__(self,C=100,n_a_m=1,N=100):
        f = 1.091     # cst
        M = 14.e-3    # molar mass of nitrogen
        rho_l = 690.8 # density of liquid nitrogen
        molarvol = M / rho_l
        N_A = 6.02214086e23 # Avogadro constant
        molecularvol = molarvol / N_A
        # sigm = f * (molecularvol)**(2./3) # surface area occupied by a molecule of adsorbate
        sigm = 0.162e-18 # surface area occupied by the di-nitrogen at 77.4 K
        self.sigm = sigm
        # BET equation (h = p/p0 and N = maximum nb of layers): 
        # C = exp((E1-El)/RT)
        # E1 = heat of adsorption of the first adsorbed layer
        # El = heat of liquefaction of the adsorbate gas
        # n_a_m = number of molecule needed to cover the surface area with one layer per unit mass of the adsorbent
        # n_a / n_a_m = C * h / (1 - h) 
        #             * (1 - (N+1) * h**N + N * h**(N+1)) / (1 + (C - 1) * h - C * h**(N+1))
        # In case where N -> infinity ie h**N -> 0
        # n_a / n_a_m = C * h / (1 - h) / (1 + (C - 1) * h)
        # We cand measure n_am and C by fitting a line
        # h/(1 - h)/n_a = n_a_m/C * (1 + (C - 1) * h)
        # In case where N=1 we obtain the Langmuir equation.
        self.ssa = N_A * sigm * n_a_m
        self.C = C
        self.n_a_m = n_a_m
        self.N = N


    def adsorbedcontent(self,h):

        C = self.C
        n_a_m = self.n_a_m
        N = self.N
        if(N < 100): hN = h**N 
        else: hN = 0

        #n_a  = n_a_m * C * h / ((1 - h) * (1 + (C-1)*h))
        #n_a *= (1 - (N+1)*hN + N*h*hN) / (1 - C*h*hN / (1 + (C-1)*h))

        n_a  = n_a_m*C*h/(1-h) * (1 - (N+1)*hN + N*h*hN)/(1 + (C-1)*h - C*h*hN)

        return n_a


    def adsorptionisotherm(self,n=10):

        C = self.C
        n_a_m = self.n_a_m
        #N = self.N
        #hN = h**N

        #h = np.linspace(0,1,n)
        hr = [float(i)/n for i in range(0,n)]
        nads = []

        for i in range(0,n):
            h = hr[i]
            n_a = self.adsorbedcontent(h)
            nads.append(n_a)

        cv = Curve(hr,nads,['# Humidity   Adsorbed content'])

        return cv




import getopt
#===============================================================================
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],'f:',['file='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    #file_name = ""

    #for opt, optval in optlist:
    #    if opt in ('-f', '--file'):
    #        file_name = optval

    #if(not file_name):
    #    print 'Error: Must specify a file with -f: e.g -f vbb.txt'
    #    sys.exit(1)

    bet = BET(C=80,N=10,n_a_m=0.1)
    ads = bet.adsorptionisotherm(100)
    ads.write()

