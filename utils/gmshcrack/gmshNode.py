#! /usr/bin/python

# Node
#=====

import sys

class Node:
    def __init__(self,ID=-1,coords=[]):
        self.ID = ID
        self.coords=coords

    def write(self,f=sys.stdout):
        line = str(self.ID) + ' ' + ' '.join([str(i) for i in self.coords]) + '\n'
        f.write(line)

  
#===============================================================================
if __name__=="__main__":
    id = 0
    coords = [0,1,2]

    n = Node(id,coords)

    n.write()
