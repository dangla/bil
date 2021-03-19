# Usage:
# python BJH.py -f filename

import numpy as np
from AdsorbedLayer import AdsorbedLayer


class BJH:

    def __init__(self,desorptioncurve=None,adsorbedlayer=AdsorbedLayer(293.15)):

        if(not desorptioncurve):
            return

        hr = desorptioncurve.x

        # Must be a desorption curve
        if(hr[0] < hr[-1]):
            return

        # The water content, sw, should be given in volume per unit volume 
        # of the porous space, i.e the measured liquid saturation degree.
        sw = desorptioncurve.y

        # Thickness of the film
        t = adsorbedlayer.thickness(hr)

        # Radius of (entry) pores
        rp = adsorbedlayer.poreradius(hr)

        # Mean radius of pores: rm[i] = 1/2 (rp[i-1] + rp[i])
        rm = []
        rm.append(rp[0])
        for i in range(1,len(hr)):
	    rm.append(0.5 * (rp[i] + rp[i-1]))

        rm = np.asarray(rm)


        # Input: volume of liquid
        VL = sw

        # Variations of the volume of liquid between two steps
        DeltaVL =[]
        # Assuming an initial non saturated state the first DeltaVL is defined
        # by the difference between 1 and the first VL or by 0 if
        # the initial state is saturated. In other words DeltaVL[0] is the variation
        # of the volume of liquid that would have been obtained if we had started
        # from a saturated state.
        DeltaVL.append(1 - VL[0])
        for i in range(1,len(hr)):
            DeltaVL.append(VL[i-1] - VL[i])


        # Variations of the volume of gas pores
        DeltaVP = []
        # Assuming an initial non saturated state, we start with a volume
        # of pore, DeltaVP[0], associated with DeltaVL[0] and taking into
        # account the adsorbed film on the surface of pore, as follows
        DeltaVP.append(DeltaVL[0] * ((rm[0]/(rm[0] - t[0]))**2))
        # Volume of gas pores
        VP = []
        VP.append(DeltaVP[0])

        # Variations of the surface area of gas pores
        DeltaAP = []
        # We assume that this first pore volume is cylinder-shaped.
        # This is a strong hypothesis if this volume is not small.
        # However the mean radius is large so this surface is small.
        DeltaAP.append(2 * DeltaVP[0] / rm[0])
        # Surface of gas pores
        AP=[]
        AP.append(DeltaAP[0])

        # Loop on the desorption steps (BJH method)
        for i in range(1,len(hr)):

            # BJH formula for the variation of the volume of pores: dvp
            b = 0
            for j in range(0,i):
                b += DeltaAP[j] * (1 - t[i]/rm[j])

            a = (rm[i]/(rm[i] - t[i]))**2
            deltat = t[i-1] - t[i]
            dvp = a * (DeltaVL[i] - deltat * b)

            # If all pores are emptied, BJH cannot apply anymore.
            # In such case we modify dvp accordingly.
            # (VP[-1] is the last value of VP already computed)
            if(VP[-1] + dvp > 1):
                dvp = 1 - VP[-1]

            # The surface of the pore, assuming cylinder-shaped pores
            dap = 2 * dvp / rm[i]

            DeltaVP.append(dvp)
            DeltaAP.append(dap)
            VP.append(VP[-1] + dvp)
            AP.append(AP[-1] + dap)

        VP = np.asarray(VP)
        AP = np.asarray(AP)


        # Surface area of pore wall per unit volume of material
        #TotalAreaOfPoreWall = AP[-1]


        # Backup
        self.humidity = hr
        self.poreradius = rp
        # volumes and surface per unit volume of porous space.
        self.surfaceareaofpores = AP
        self.adsorbedcontent = AP * t
        self.gassaturation = VP
        self.adsorbedlayer = adsorbedlayer


    def poresizedistribution(self):

        rp = self.poreradius
        vp = self.gassaturation

        psd = []
        psd.append(0)
        for i in range(1,len(rp)):
            psd.append((vp[i] - vp[i-1]) * rp[i] / (rp[i-1] - rp[i]))

        return psd



    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout
        
        # write header
        f.write('# Relative humidity  ' +
                '  Pore radius        ' +
                '  Volume of gas pores' +
                '  Surface of gas pores' +
                ' Pore size distribution' + '\n')


        psd = self.poresizedistribution()

        # write columns
        col_format = " {:< 20e}"*5 + "\n"   # left-justified columns
	for y in zip(self.humidity,self.poreradius,
                     self.gassaturation,self.surfaceareaofpores,psd):
	    f.write(col_format.format(*y))

        #f.close()



import sys
import getopt
sys.path.append('../tools')
from Curve import Curve
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
    bjh.write()
    #bjh.write('BJH-'+file_name)

    #SBJH = bjh.surfaceareaofpores[-1]
    #print "BJH surface area per unit of pore volume =",  SBJH,"(m2/m3)"
