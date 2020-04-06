#! /usr/bin/python

# Mesh
#=====

import sys
from gmshNode    import Node
from gmshElement import Element

class Mesh:
    
    def __init__(self,nodes=[], elements=[]):
        self.nodes = nodes
        self.elements = elements

    def nNodes(self):
        return len(self.nodes)
        
    def nElements(self):
        return len(self.elements)
    
    def write(self,fname=''):

        if fname:
            f = open(fname,'w')
        else:
            f = sys.stdout
        
        # write header
        f.write('$MeshFormat\n')
        f.write('2.2 0 8\n')
        f.write('$EndMeshFormat\n')
        
        # write nodes
        f.write('$Nodes\n')
        f.write(str(self.nNodes())+'\n')
        for n in self.nodes:
            n.write(f)
        f.write('$EndNodes\n')

        #write elements
        f.write('$Elements\n')
        f.write(str(self.nElements())+'\n')
        for e in self.elements:
            e.write(f)
        f.write('$EndElements')

        f.close()


    def read(self,fname):

        if not fname:
            print 'Error: Must specify mesh file'
            sys.exit(1)

        nodes = self.nodes
        elements = self.elements
    
        nnodes = -1
        nelements = -1

        # python 2.7 doesn't have enums...
        PARSING_NODES = 1
        PARSING_ELEMENTS = 2
        ACTION_UNSET = 3
    
        # set current action in file
        ACTION = ACTION_UNSET
    
        #open file
        f = open(fname, 'r')
        contents = f.read()
        f.close()

        lines = contents.split('\n')
        for line  in lines:
            if(ACTION == ACTION_UNSET):
                if(line == '$Nodes'):
                    ACTION = PARSING_NODES

                elif(line == '$Elements'):
                    ACTION = PARSING_ELEMENTS

                continue
            
            elif(ACTION == PARSING_NODES):
                if(line == '$EndNodes'):
                    ACTION = ACTION_UNSET
                    continue
                if(nnodes == -1):
                    nnodes = int(line)
                    continue
                
                tmp = line.split()
                i = int(tmp[0])
                coords = [float(t) for t in tmp[1:]]

                nodes.append(Node(i,coords))

            elif(ACTION == PARSING_ELEMENTS):
                if(line == '$EndElements'):
                    ACTION = ACTION_UNSET
                    continue

                if(nelements == -1):
                    nelements = int(line)
                    continue
                
                tmp = line.split()
                i = int(tmp[0])
                eltp = int(tmp[1])
                ntags = int(tmp[2])
                tags = [int(t) for t in tmp[3:3+ntags]]
                physid = int(tags[0])
                elemid = int(tags[1])
                el_nodes = [int(t) for t in tmp[3+ntags:]]

                elements.append(Element(i,eltp,tags,el_nodes,physid,elemid))


    def makecrackedmesh(self, crack_id, oneside_id, cracktip_id):
    
        # make set of node IDs for nodes on crack lip
        crack_nodes = []
        for e in self.elements:
            if (e.elem_ID == crack_id):
                for n in e.nodes:
                    crack_nodes.append(n)
    
        # remove node IDs for nodes on crack tip
        for e in self.elements:
            if (e.elem_ID == cracktip_id):
                for n in e.nodes:
                    if n in crack_nodes:
                        crack_nodes.remove(n)
    
        crack_nodes = [n for n in set(crack_nodes)]
        
        new_nodes = []
        for i in range(0,len(crack_nodes)):
            new_nodes.append(self.nNodes() + i + 1)
        
        node_mapping = dict([(crack_nodes[i],new_nodes[i]) for i in range(0,len(crack_nodes))])
    
        # set nodes in elementary_ID=oneside_ID to the newly created nodes
        for e in self.elements:
            if (e.elem_ID == oneside_id):
                for i in range(0,len(e.nodes)):
                    n = e.nodes[i]
                    if node_mapping.has_key(n):
                        e.nodes[i] = node_mapping[n]
            
        
        # Turn the surface crack elements into volume elements
        for e in self.elements:
            if (e.elem_ID == crack_id):
                nb_nodes = len(e.nodes)

                # Append nodes
                for i in range(0,nb_nodes):
                    n = e.nodes[i]
                    if node_mapping.has_key(n):
                        e.nodes.append(node_mapping[n])

                # Turn 2-node line (type = 1) into 
                if (e.eltype == 1):
                #     i)   4-node quadrangle  (type = 3)
                    if (len(e.nodes) == 4):
                        e.nodes = [e.nodes[0],e.nodes[1],e.nodes[3],e.nodes[2]]
                        e.eltype = 3
                #     ii)  3-node triangle    (type = 2)
                    elif (len(e.nodes) == 3):
                        e.eltype = 2
                    else:
                        print('Not available nb of nodes = %(nb)d'%{'nb':len(e.nodes)})
                        sys.exit(1)

                # Turn 3-node triangle (type = 2) into 
                elif (e.eltype == 2):
                #     i)   6-node prism       (type = 6)
                    if (len(e.nodes) == 6):
                        e.eltype = 6
                #     ii)  5-node pyramid     (type = 7)
                    elif (len(e.nodes) == 5):
                        if not node_mapping.has_key(e.nodes[0]):
                            e.nodes = [e.nodes[1],e.nodes[2],e.nodes[4],e.nodes[3],e.nodes[0]]
                        elif not node_mapping.has_key(e.nodes[1]):
                            e.nodes = [e.nodes[0],e.nodes[2],e.nodes[4],e.nodes[3],e.nodes[1]]
                        elif not node_mapping.has_key(e.nodes[2]):
                            e.nodes = [e.nodes[0],e.nodes[1],e.nodes[4],e.nodes[3],e.nodes[2]]
                        e.eltype = 7
                #     iii) 4-node tetrahedron (type = 2)
                    elif (len(e.nodes) == 4):
                        e.eltype = 2
                    else:
                        print('Not available nb of nodes = %(nb)d'%{'nb':len(e.nodes)})
                        sys.exit(1)

                # Turn 4-node quadrangle (type = 4) into 
                elif (e.eltype == 4):
                #     i)   8-node hexahedron  (type = 5)
                    if (len(e.nodes) == 8):
                        e.eltype = 5
                    else:
                        print('Not available nb of nodes = %(nb)d'%{'nb':len(e.nodes)})
                        sys.exit(1)

                else:
                    print('Not available type = %(typ)d'%{'typ':e.eltype})
                    sys.exit(1)
        
        
        # make new node objects with coordinates overlapping
        new_mesh_nodes = self.nodes
        for n in self.nodes:
            if (n.ID in crack_nodes):
                nodeid = node_mapping[n.ID]
                new_mesh_nodes.append(Node(nodeid,n.coords))
    
        return Mesh(new_mesh_nodes,self.elements)

  
#===============================================================================
if __name__=="__main__":
    nodes = []
    elements = []
    M = Mesh(nodes,elements)
    M.write()
