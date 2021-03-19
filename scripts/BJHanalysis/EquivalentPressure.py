# Usage:
# python EquivalentPressure.py -f filename

import numpy as np
import sys
sys.path.append('../tools')
from Curve import Curve

class EquivalentPressure:
    def __init__(self,bjh,pg=0):

        # Adsorbed layer
        adsorbedlayer = bjh.adsorbedlayer

        # Humidity
        hr = bjh.humidity

        # Temperature
        T  = adsorbedlayer.temperature

        # Gas and liquid saturation degree
        sg = bjh.gassaturation
        sl = 1 - sg

        # Capillary pressure
        R = 8.314         # Perfect gas constant
        RT = R * T
        molarvol = 18.e-6 # Molar volume of water
        pc = - RT / molarvol * np.log(hr)

        # Liquid pressure
        pl = pg - pc

        # Surface area of pores per unit volume of pore (m2/m3)
        surfaceareaofpores = bjh.surfaceareaofpores

        # Adsorbed (volume) content per unit volume of pore (m3/m3)
        adscontent = bjh.adsorbedcontent
        
        # Average pressure term
        averagepressure = sl * pl + sg * pg

        # Interface energy term
        saturationcurve = Curve(sl,pc,['# Saturation Capillary pressure'])
        interfacepressure = 2./3 * saturationcurve.integrate().y

        # Adsorption (Bangham) term
        mu = RT / molarvol * np.log(hr)
        adscurve = Curve(mu,adscontent,['# Chemical potential Adsorbed content'])
        adsorptionpressure = 2./3 * adscurve.integrate().y

        # Shuttleworth term
        shuttleworthpressure = - 2./3 * adsorbedlayer.shuttleworth * surfaceareaofpores

        # backup
        self.humidity = hr
        self.average = averagepressure
        self.interface = interfacepressure
        self.adsorption = adsorptionpressure
        self.shuttleworth = shuttleworthpressure
        self.total = averagepressure + interfacepressure + adsorptionpressure + shuttleworthpressure


    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout
        
        # write header
        f.write('# Relative humidity  ' +
                '  Average pressure   ' +
                '  Interface energy   ' +
                '  Adsorption energy  ' +
                '  Shuttleworth term  ' +
                '  Total pressure     ' + '\n')

        # write columns
        col_format = " {:< 20e}"*6 + "\n"   # left-justified columns
	for y in zip(self.humidity,self.average,
                     self.interface,self.adsorption,
                     self.shuttleworth,self.total):
	    f.write(col_format.format(*y))

        #f.close()



import getopt
from BJH import BJH
from AdsorbedLayer import AdsorbedLayer
#===============================================================================
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],'f:',['file='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    file_name = ""

    for opt, optval in optlist:
        if opt in ('-f', '--file'):
            file_name = optval

    if(not file_name):
        print 'Error: Must specify a file with -f: e.g -f vbb.txt'
        sys.exit(1)

    adsorp  = Curve().read(file_name)
    adsorps = adsorp.sample(100)
    desorp  = adsorps.reverse()
    #desorp.write()

    adsorbedlayer = AdsorbedLayer(293.15)
    #adsorbedlayer.write()

    bjh = BJH(desorp,adsorbedlayer)

    equivalentpressure = EquivalentPressure(bjh)
    equivalentpressure.write()
    #equivalentpressure.write('Pi-'+file_name)



