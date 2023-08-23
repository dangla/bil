# Usage:
# python String.py -f filename

import re

class String:

    def __init__(self,contents=[]):
    
        self.contents = contents


    def read(self,fname=None):
    
        # open file and store its contents
        if(fname):
            f = open(fname,'r')
            self.contents = f.read()
            f.close()


    def removecomments(self):
    
        s = self.contents
    
        for x in re.findall(r'("[^\n]*"(?!\\))|(//[^\n]*$|/(?!\\)\*[\s\S]*?\*(?!\\)/)',s,8):
            s = s.replace(x[1],'')
        
        r = String(s)
        
        return r
    
    
    def write(self):
    
        print(self.contents)



import sys
import getopt
#===============================================================================
if __name__=="__main__":
    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'f:r',
                       ['file=','remove'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    for opt, optval in optlist:
        if opt in ('-f', '--file'):
            file_name = optval
        elif opt in ('-r', '--remove'):
            remove = True

    if(not "file_name" in locals()):
        print('Error: Must specify a file with -f')
        sys.exit(1)

    s = String()
    
    s.read(file_name)
    
    if("remove" in locals()):
        print('# Remove the comments')
        rs = s.removecomments()
        s = rs

    s.write()

