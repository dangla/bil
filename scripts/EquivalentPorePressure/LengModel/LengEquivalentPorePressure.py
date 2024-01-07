# Usage:
# python3 LengEquivalentPorePressure.py -f filename -s 100 -r

import numpy as np
import sys
sys.path.append('../../Curve')
from Curve import Curve
sys.path.append('../../BJH')
from BJH import BJH


class LengEquivalentPorePressure:
    def __init__(self,bjh,pg=0):

        # Adsorbed layer
        adsorbedlayer = bjh.adsorbedlayer

        # Humidity
        hr = bjh.humidity

        # Temperature
        T  = adsorbedlayer.temperature
        
        # Physical constant
        R = 8.314         # Perfect gas constant
        RT = R * T

        # Gas and liquid saturation degree
        sg = bjh.gassaturation
        sl = 1 - sg
        
        # Molar volume of the gas
        molarvol = adsorbedlayer.molarvolume

        # Capillary pressure
        pc = - RT / molarvol * np.log(hr)

        # Liquid pressure
        pl = pg - pc

        # Surface area of pores per unit volume of pore (m2/m3)
        surfaceareaofpores = bjh.surfaceareaofpores
        
        # Average pressure term
        averagepressure = sl * pl + sg * pg
        
        # Surface stress term
        mu = RT * np.log(hr)
        beta = 1
        adsdensity = 1./molarvol * adsorbedlayer.thickness(hr)
        betagamma = Curve(mu,beta*adsdensity,['# beta x adsorbed density'])
        surfacestresspressure = 2./3. * betagamma.integrate().y * surfaceareaofpores
        
        #
        gamma = Curve(mu,adsdensity,['# adsorbed density'])
        derivativeofsurfaceareaofpores = bjh.derivativeofsurfaceareaofpores
        pctheory = - gamma.integrate().y * derivativeofsurfaceareaofpores

        # backup
        self.humidity = hr
        self.average = averagepressure
        self.surfacestress = surfacestresspressure
        self.total = averagepressure + surfacestresspressure
        self.capillarypressure = pc
        self.pctheory = pctheory


    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout
        
        # write header
        f.write('# Relative humidity  ' +
                '  Average pressure   ' +
                '  Surface stress     ' +
                '  Total pressure     ' + 
                '  Capillary pressure ' +
                '  Theoretical P_c    ' +
                '\n')

        # write columns
        col_format = " {:< 20e}"*6 + "\n"   # left-justified columns
        for y in zip(self.humidity,self.average,
                     self.surfacestress,self.total,
                     self.capillarypressure,self.pctheory):
            f.write(col_format.format(*y))

        #f.close()



import getopt
from BJH import BJH
from AdsorbedLayer import AdsorbedLayer
#===============================================================================
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'hf:s:r',
                       ['help','file=','sample=','reverse'])


    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)


    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            print('Usage: python3 %s -h -f<filename> -s<sample> -r' % sys.argv[0])
            print('\n')
            print('Options:')
            print('-s, --sample:     Number of sampling points')
            print('-r, --reverse:    Reverse the curve')
            print('-f, --file:       Sorption file name')
            print('-h, --help:       Display this help message')
            print('\n')
            sys.exit()
        elif opt in ('-f', '--file'):
            file_name = optval
        elif opt in ('-s', '--sample'):
            sample = int(optval)
        elif opt in ('-r', '--reverse'):
            reverse = True



    if(not "file_name" in locals()):
        print('Error: Must specify a file with -f: e.g -f vbb.txt')
        sys.exit(1)


    sorpcurve  = Curve().read(file_name)
    
    
    if("sample" in locals()):
        print('# Sample the curve with %(nb)d points'%{'nb':sample})
        scurve = sorpcurve.sample(sample)
        sorpcurve = scurve
        
    if("reverse" in locals()):
        print('# Reverse the curve')
        rcurve = sorpcurve.reverse()
        sorpcurve = rcurve
        

    adsorbedlayer = AdsorbedLayer(293.15)
    # idem as:
    #adsorbedlayer = AdsorbedLayer(293.15,0.385e-9,-0.189e-9,4.97e-2,"water")

    bjh = BJH(sorpcurve,adsorbedlayer)

    ET = LengEquivalentPorePressure(bjh)
    ET.write()
    #ET.write('Pi-'+file_name)



