import numpy as np
from Curve import Curve

class Curves:
    
    def __init__(self,fname=None):

        self.comments = []
        self.curve = []
        self.n = 0
    
        # open file
        if(fname):
            f = open(fname,'r')
            contents = f.read()
            f.close()
        else:
            return

        lines = contents.split('\n')
        for line  in lines:

            if(line.startswith('#')):

                self.comments.append(line)
            
            elif(line):
                
                tmp = line.split()

                x = float(tmp[0])

                for i in range(0,len(tmp)-1):
                    if(i > len(self.curve)-1):
                        cv = Curve()
                        self.curve.append(cv)

                    y = float(tmp[i+1])
                    self.curve[i].x.append(x)
                    self.curve[i].y.append(y)

                n = len(tmp) - 1
                if(n > self.n):
                    self.n = n

            else:
                break


    def reverse(self):
        
        newcurves = Curves()
        newcurves.n = self.n

        for i in range(0,len(self.comments)):
            newcurves.comments.append(self.comments)

        for i in range(0,self.n):
            newcurves.curve.append(self.curve[i].reverse())

        return newcurves


    def sample(self,n):

        newcurves = Curves()
        newcurves.n = self.n

        for i in range(0,len(self.comments)):
            newcurves.comments.append(self.comments)

        for i in range(0,self.n):
            newcurves.curve.append(self.curve[i].sample(n))
        return newcurves


    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout

        n = self.n

        # write columns
        col_format = " {:< 20e}"*2 + "\n"   # left-justified columns
        for i in range(0,n):
            f.write('# Curve {0:d} \n'.format(i))
            for c in self.comments:
                f.write(str(c) + '\n')
	    for y in zip(self.curve[i].x,self.curve[i].y):
	        f.write(col_format.format(*y))
            f.write("\n")

        #f.close()



import sys
import getopt
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
        print 'Error: Must specify a file with -f'
        sys.exit(1)

    cv = Curves(file_name)
    cv.write()

    scv = cv.sample(100)
    scv.write()

    rcv = cv.reverse()
    rcv.write()


