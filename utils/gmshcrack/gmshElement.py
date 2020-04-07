#! /usr/bin/python

# Element
#========

import sys

class Element:
    def __init__(self,ID=-1,eltype=-1,tags=[],nodes=[],phys_ID=-1,elem_ID=-1):
        self.ID = ID
        self.eltype = eltype
        self.tags = tags
        self.nodes = nodes
        self.phys_ID = phys_ID
        self.elem_ID = elem_ID

    def write(self,f=sys.stdout):
        line = str(self.ID) + ' ' + str(self.eltype) + ' ' + str(len(self.tags)) + ' ' \
             + ' '.join([str(i) for i in self.tags]) + ' ' \
             + ' '.join([str(i) for i in self.nodes]) + '\n'
        f.write(line)

  
#===============================================================================
if __name__=="__main__":
    id = 0
    eltp = 2
    ntags = 2
    tags = ['1','2']
    physid = int(tags[0])
    elemid = int(tags[1])
    el_nodes = ['1','2','3']

    e = Element(id,eltp,tags,el_nodes,physid,elemid)

    e.write()
