# Usage:
# python Curve.py -f filename

import numpy as np
import matplotlib.pyplot as plt

class Curve:
    
    def __init__(self,x=[],y=[],comments=[]):

        self.comments = comments
        self.x = x
        self.y = y


    def reverse(self):
        
        x_new = self.x[::-1]
        y_new = self.y[::-1]
        newcurve = Curve(x_new,y_new)
        return newcurve


    def sample(self,n):

        x_new = np.linspace(self.x[0],self.x[-1],n)
        y_new = np.interp(x_new,self.x,self.y)
        newcurve = Curve(x_new,y_new)
        return newcurve


    def integrate(self):

        x = self.x
        y = self.y
        integral = []

        integral.append(0)
        for i in range(1,len(x)):
            dx = x[i] - x[i-1]
            ds = 0.5 * (y[i] + y[i-1]) * dx
            s  = integral[-1] + ds
            integral.append(s)

        integral = np.asarray(integral)

        newcurve = Curve(x,integral)

        return newcurve
    

    def read(self,fname=None):

        # clear
        self.comments[:] = []
        self.x[:] = []
        self.y[:] = []
    
        # open file
        if(fname):
            f = open(fname,'r')
            contents = f.read()
            f.close()
        else:
            return self

        lines = contents.split('\n')
        for line  in lines:

            if(line.startswith('#')):

                toks = line.split()
                
                for i in range(1,len(toks)):
                    self.comments.append(toks[i])
            
            elif(line):
                
                tmp = line.split()
                x = float(tmp[0])
                y = float(tmp[1])
                self.x.append(x)
                self.y.append(y)

            else:
                break

        return self


    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout
        
        # write header
        f.write('# ')
        if(self.comments):
            f.write(' '.join([self.comments[i] for i in range(0,len(self.comments))]))
        f.write('\n')

        # write columns
        col_format = " {:< 20e}"*2 + "\n"   # left-justified columns
        for y in zip(self.x,self.y):
            f.write(col_format.format(*y))

        #f.close()
        
    
    def plot(self):
      
        x = self.x
        y = self.y
        
        lines = plt.plot(x,y)
        plt.setp(lines,color='b',linewidth=2.0,linestyle='-')
        plt.xlabel(self.comments[0])
        plt.ylabel(self.comments[1])
        #plt.axis([x[0],x[-1],y[0],y[-1]])
        plt.grid(True)
        plt.show()



import sys
import getopt
#===============================================================================
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'f:s:ri',
                       ['file=','sample=','reverse','integrate'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)


    for opt, optval in optlist:
        if opt in ('-f', '--file'):
            file_name = optval
        elif opt in ('-s', '--sample'):
            sample = int(optval)
        elif opt in ('-r', '--reverse'):
            reverse = True
        elif opt in ('-i', '--integrate'):
            integrate = True

    if(not "file_name" in locals()):
        print('Error: Must specify a file with -f')
        sys.exit(1)


    print('# The curve is %(file)s'%{'file':file_name})
    cv = Curve().read(file_name)


    if("sample" in locals()):
        print('# Sample the curve with %(nb)d points'%{'nb':sample})
        scv = cv.sample(sample)
        cv = scv

    if("reverse" in locals()):
        print('# Reverse the curve')
        rcv = cv.reverse()
        cv = rcv

    if("integrate" in locals()):
        print('# Integrate the curve')
        icv = cv.integrate()
        cv = icv
        
    
    # Output    
    cv.write()
    
    # Plot
    cv.plot()


