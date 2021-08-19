# Usage:
# python BET.py

import numpy as np
import sys
sys.path.append('../tools')
from Curve import Curve

class BET:
    
    def __init__(self,C=100,n_a_m=1,N=100):
        #f = 1.091     # cst
        #M = 14.e-3    # molar mass of nitrogen
        #rho_l = 690.8 # density of liquid nitrogen
        #molarvol = M / rho_l
        #N_A = 6.02214086e23 # Avogadro constant
        #molecularvol = molarvol / N_A
        # sigm = f * (molecularvol)**(2./3) # surface area occupied by a molecule of adsorbate
        #sigm = 0.162e-18 # surface area occupied by the di-nitrogen at 77.4 K
        #self.sigm = sigm
        #self.ssa = N_A * sigm * n_a_m
        self.C = C
        self.n_a_m = n_a_m
        self.N = N


    def adsorbedcontent(self,h):

        C = self.C
        n_a_m = self.n_a_m
        N = self.N
        hN = h**N 
        
        #if(N > 100):
            #hN = 0
            #print('N = %d',N)

        # BET equation (input: h = p/p0)
        # N = maximum nb of layers
        # C = exp((E1-El)/RT)
        # E1 = heat of adsorption of the first adsorbed layer
        # El = heat of liquefaction of the adsorbate gas
        # n_a_m = number of molecule needed to cover the surface area with one layer
        # The number of layers, theta = n_a / n_a_m, is given by
        # theta = C*h/(1 - h)/(1 + (C - 1)*h - C*h*h**N)
        #          *(1 - (N+1-N*h)*h**N)
        # In case where N -> infinity ie h**N -> 0
        # theta = C*h/(1 - h)/(1 + (C - 1)*h)
        # We cand measure n_am and C by fitting a line
        # h/(1 - h)/n_a = n_a_m/C * (1 + (C - 1)*h)
        # In case where N=1 we obtain the Langmuir equation.

        theta  = C*h/(1 - h)/(1 + (C - 1)*h - C*h*hN)*(1 - (N+1-N*h)*hN)
        
        n_a  = n_a_m*theta

        return n_a


    def adsorptionisotherm(self,n=10):

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
        optlist,args = getopt.getopt(sys.argv[1:],
                       'hN:C:n:s:',
                       ['help','Number=','Constant=','n_a_m=','sample='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)


    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            print('Usage: python %s -h -N<Number> -s<sample> -C<Constant> -n<n_a_m>' % sys.argv[0])
            print('\n')
            print('Options:')
            print('-s, --sample:     Number of sampling points')
            print('-N, --Number:     Max number of layers')
            print('-C, --Constant:   BET constant i.e. C = exp((E1-El)/RT)')
            print('-n, --n_a_m:      Quantity of adsorbate needed to cover the surface of one layer')
            print('-h, --help:       Display this help message')
            print('\n')
            sys.exit()
        elif opt in ('-s', '--sample'):
            sample = int(optval)
        elif opt in ('-N', '--Number'):
            Number = int(optval)
        elif opt in ('-C', '--Constant'):
            Constant = float(optval)
        elif opt in ('-n', '--n_a_m'):
            n_a_m = float(optval)

    
    if(not "sample" in locals()):
        sample = 100
    
    if(not "Number" in locals()):
        Number = 100
        
    if(not "Constant" in locals()):
        Constant = 80
        
    if(not "n_a_m" in locals()):
        n_a_m = 1


    bet = BET(Constant,n_a_m,Number)
    ads = bet.adsorptionisotherm(sample)
    ads.write()

