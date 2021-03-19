#! /usr/bin/python

# Element
#========

import sys

gmshDim = {15:0,1:1,8:1,2:2,9:2,3:2,10:2,16:2,4:3,11:3,5:3,12:3,17:3,6:3,13:3,18:3,7:3,14:3,19:3}
gmshNbOfVerticeNodes = {15:1,1:2,8:2,2:3,9:3,3:4,10:4,16:4,4:4,11:4,5:8,12:8,17:8,6:6,13:6,18:6,7:5,14:5,19:5}
gmshNbOfEdgeNodes    = {8:1,9:3,10:4,16:4,11:6,12:12,17:12,13:9,18:9,14:8,19:8}
gmshNbOfFaceNodes    = {10:1,12:6,14:1}
gmshNbOfVolumeNodes  = {12:1}
gmshNbOfNodes = {15:1,1:2,8:2+1,2:3,9:3+3,3:4,10:4+4+1,16:4+4,4:4,11:4+6,5:8,12:8+12+6+1,17:8+12,6:6,13:6+9+3,18:6+9,7:5,14:5+8+1,19:5+8}
gmshType = {(0,1):15,(1,2):1,(1,3):8,(2,3):2,(2,6):9,(2,4):3,(2,9):10,(2,8):16,(3,4):4,(3,10):11,(3,8):5,(3,27):12,(3,20):17,(3,6):6,(3,18):13,(3,15):18,(3,5):7,(3,14):14,(3,13):19}

class Element:
    def __init__(self,ID=-1,eltype=-1,tags=None,nodes=None,phys_ID=-1,elem_ID=-1):
        self.ID = ID
        self.eltype = eltype
        if tags is None:
            self.tags = []
        else:
            self.tags = tags
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes
        self.phys_ID = phys_ID
        self.elem_ID = elem_ID


    def dim(self):
        return gmshDim[self.eltype]


    def surfacenodeordering(self,node_mapping):
        # Return the nodes of a surface element. The latter should oriented correctly
        # and so the nodes should be given in the right order.
        # These nodes are the keys of node_mapping and their ordering depend on the nodes 
        # of the element they overlap. The orientation of the surface element must be
        # opposite to that of the element.

        n_surfel_nodes = len(node_mapping.keys())
        dim = self.dim()
        if not gmshType.has_key((dim-1,n_surfel_nodes)):
            return []

        surfel_nodes = []

        # first set the nodes of the surface element in the order so that they match 
        # those of the element in the inverse order of its numbering
        for i in range(0,len(self.nodes)):
            j = len(self.nodes) - i - 1
            ni = self.nodes[j]
            if node_mapping.has_key(ni):
                surfel_nodes.append(ni)


        #     type)   description
        if (self.eltype == 0):
            print('Not available type = %(type)d'%{'type':self.eltype})
            sys.exit(1)

        #     1)  2-node line
        #elif (self.eltype == 1):

        #     2)  3-node triangle
        elif (self.eltype == 2):
            if(node_mapping.has_key(self.nodes[0]) and node_mapping.has_key(self.nodes[2])):
                surfel_nodes = [self.nodes[0],self.nodes[2]]

        #     3)  4-node quadrangle
        elif (self.eltype == 3):
            if(node_mapping.has_key(self.nodes[0]) and node_mapping.has_key(self.nodes[3])):
                surfel_nodes = [self.nodes[0],self.nodes[3]]

        #     4)  4-node tetrahedron
        #elif (elf.eltype == 4):

        #     5)  8-node hexahedron
        #elif (elf.eltype == 5):

        #     6)  6-node prism
        #elif (elf.eltype == 6):

        #     7)  5-node pyramid
        #elif (self.eltype == 7):

        #     8)  3-node second order line
        #elif (elf.eltype == 8):

        #     9)  6-node second order triangle
        #elif (elf.eltype == 9):

        #     10) 9-node second order quadrangle
        #elif (elf.eltype == 10):

        #     15) 1-node point
        #elif (self.eltype == 15):
        else:
            print('Not available type = %(type)d'%{'type':self.eltype})
            sys.exit(1)


        return surfel_nodes


            
    def turnsurfaceintovolume(self,node_mapping):

        # Append nodes
        for n in self.nodes:
            if node_mapping.has_key(n):
                self.nodes.append(node_mapping[n])

        # Turn 2-node line (type = 1) into 
        if (self.eltype == 1):
            #     i)   4-node quadrangle  (type = 3)
            if (len(self.nodes) == 4):
                self.nodes = [self.nodes[0],self.nodes[1],self.nodes[3],self.nodes[2]]
                self.eltype = 3
            #     ii)  3-node triangle    (type = 2)
            elif (len(self.nodes) == 3):
                self.eltype = 2
            else:
                print('Not available nb of nodes = %(nb)d'%{'nb':len(self.nodes)})
                sys.exit(1)

        # Turn 3-node triangle (type = 2) into 
        elif (self.eltype == 2):
            #     i)   6-node prism       (type = 6)
            if (len(self.nodes) == 6):
                self.eltype = 6
            #     ii)  5-node pyramid     (type = 7)
            elif (len(self.nodes) == 5):
                if not node_mapping.has_key(self.nodes[0]):
                    self.nodes = [self.nodes[1],self.nodes[2],self.nodes[4],self.nodes[3],self.nodes[0]]
                elif not node_mapping.has_key(self.nodes[1]):
                    self.nodes = [self.nodes[0],self.nodes[2],self.nodes[4],self.nodes[3],self.nodes[1]]
                elif not node_mapping.has_key(self.nodes[2]):
                    self.nodes = [self.nodes[0],self.nodes[1],self.nodes[4],self.nodes[3],self.nodes[2]]
                self.eltype = 7
            #     iii) 4-node tetrahedron (type = 2)
            elif (len(self.nodes) == 4):
                self.eltype = 2
            else:
                print('Not available nb of nodes = %(nb)d'%{'nb':len(self.nodes)})
                sys.exit(1)

        # Turn 4-node quadrangle (type = 4) into 
        elif (self.eltype == 4):
            #     i)   8-node hexahedron  (type = 5)
            if (len(self.nodes) == 8):
                self.eltype = 5
            else:
                print('Not available nb of nodes = %(nb)d'%{'nb':len(self.nodes)})
                sys.exit(1)

        # Turn 1-node point (type = 15) into
        elif (self.eltype == 15):
            #     i)   2-node line    (type = 1)
            if (len(self.nodes) == 2):
                self.eltype = 1
            else:
                print('Not available nb of nodes = %(nb)d'%{'nb':len(self.nodes)})
                sys.exit(1)
            return

        else:
            print('Not available type = %(typ)d'%{'typ':self.eltype})
            sys.exit(1)


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
