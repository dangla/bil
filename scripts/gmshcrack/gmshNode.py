#! /usr/bin/python

# Node
#=====

import sys

class Node:
    def __init__(self,ID=-1,coords=None,elements=None):
        self.ID = ID
        if coords is None:
            self.coords = []
        else:
            self.coords = coords
        if elements is None:
            self.elements = []
        else:
            self.elements = elements

    def write(self,f=sys.stdout):
        line = str(self.ID) + ' ' + ' '.join([str(i) for i in self.coords])
        f.write(line)

        if(len(self.elements) > 0):
            line = ' ' + '(' + ','.join([str(i) for i in self.elements]) + ')'
            f.write(line)

        f.write('\n')

  
#===============================================================================
if __name__=="__main__":
    id = 0
    coords = [0,1,2]

    n = Node(id,coords)

    n.write()
