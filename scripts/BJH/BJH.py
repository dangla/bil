# Usage:
# python3 BJH.py -f filename -s 100 -r

import sys
import numpy as np
from AdsorbedLayer import AdsorbedLayer


class BJH:

    def __init__(self,sorptioncurve=None,adsorbedlayer=AdsorbedLayer(293.15)):

        if(not sorptioncurve):
            print("Empty name for the desorption curve")
            sys.exit(1)
            return

        hr = sorptioncurve.x

        # Must be a desorption curve
        if(hr[0] < hr[-1]):
            print("Must be a desorption curve")
            sys.exit(1)
            return

        # Input: Sw = the water content
        # The water content, Sw, should be given in volume per unit volume 
        # of the porous space.
        # It takes into account the capillary water and the adsorbed water
        # on the surface of empty pores: Sw = SL + VADS.
        
        Sw = sorptioncurve.y

        # Thickness of the film
        
        t = adsorbedlayer.thickness(hr)

        # Radius of pores
        
        rp = adsorbedlayer.poreradius(hr)

        # Mean radius of pores: rm[i] = 1/2 (rp[i-1] + rp[i])
        
        rm = []
        rm.append(rp[0])
        for i in range(1,len(hr)):
            rm.append(0.5 * (rp[i] + rp[i-1]))

        rm = np.asarray(rm)



        # Variations of the volume of liquid between two steps
        
        DeltaSw =[]
        
        # Assuming an initial non saturated state the first DeltaSw is defined
        # by the difference between 1 and the first Sw or by 0 if
        # the initial state is saturated. In other words DeltaSw[0] is the variation
        # of the volume of liquid that would have been obtained if we had started
        # from a saturated state.
        
        DeltaSw.append(1 - Sw[0])
        for i in range(1,len(hr)):
            DeltaSw.append(Sw[i-1] - Sw[i])


        # Volume of gas pores, SG, and its variation
        
        SG = []
        DeltaSG = []
        
        # Assuming an initial non saturated state, we start with a volume
        # of pore, DeltaSG[0], associated with DeltaSw[0] and taking into
        # account the adsorbed film on the surface of pore, as follows
        
        DeltaSG.append(DeltaSw[0] * ((rm[0]/(rm[0] - t[0]))**2))
        SG.append(DeltaSG[0])

        # Surface of gas pores, AG, and its variation
        
        AG=[]
        DeltaAG = []
        
        # We assume that this first pore volume is cylinder-shaped.
        # This is a strong hypothesis if this volume is not small.
        # However the mean radius is large so this surface is small.
        
        DeltaAG.append(2 * DeltaSG[0] / rm[0])
        AG.append(DeltaAG[0])

        # Loop on the desorption steps (BJH method)
        for n in range(1,len(hr)):

            # The surface of the adorbate subjected to evaporation at step n
            # is firstly assessed to the surface of the gas pores, AG[n-1].
            
            adsurf = AG[n-1]
            
            # Because of the thickness of the layer it is smaller than 
            # the surface of the gas pores assumed as a series of cylinders.
            # So we correct it as follows

            for j in range(0,n):
                adsurf -= DeltaAG[j] * t[n] / rm[j]

            # BJH formula for the variation of the volume of gas pores: dsg
            
            a = (rm[n]/(rm[n] - t[n]))**2
            deltat = t[n-1] - t[n]
            dsg = a * (DeltaSw[n] - deltat * adsurf)
            
            # If dsg < 0 it's probably because the surface adsurf is too large.
            # So we cannot allow a negative dsg.
            
            if(dsg < 0):
                dsg = 0

            # If all pores are emptied, BJH cannot apply anymore.
            # In such case we modify dsg accordingly.
            
            if(SG[n-1] + dsg > 1):
                dsg = 1 - SG[n-1]

            # The surface of the pore, assuming cylinder-shaped pores
            
            dag = 2 * dsg / rm[n]

            DeltaSG.append(dsg)
            DeltaAG.append(dag)
            SG.append(SG[n-1] + dsg)
            AG.append(AG[n-1] + dag)

        SG = np.asarray(SG)
        AG = np.asarray(AG)


        # Backup
        self.humidity = hr
        self.poreradius = rp
        # volumes and surface per unit volume of porous space.
        self.surfaceareaofpores = AG
        self.gassaturation = SG
        
        SL = []
        VADS = []
        for i in range(0,len(SG)):
            SL.append(1 - SG[i])
            VADS.append(Sw[i] - 1 + SG[i])
        
        self.liquidsaturation = SL
        self.adsorbedcontent = VADS
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
            
        try:
            self.poreradius
        except:
            return
        
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



import getopt
sys.path.append('../Curve')
from Curve import Curve
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

    bjh = BJH(sorpcurve,adsorbedlayer)
    bjh.write()
    #bjh.write('BJH-'+file_name)

    #SBJH = bjh.surfaceareaofpores[-1]
    #print("BJH surface area per unit of pore volume =",  SBJH,"(m2/m3)")
