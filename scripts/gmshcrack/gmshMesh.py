#! /usr/bin/python

# Mesh
#=====

import sys
import getopt
import copy
from gmshNode    import Node
from gmshElement import Element

gmshType = {(0,1):15,(1,2):1,(1,3):8,(2,3):2,(2,6):9,(2,4):3,(2,9):10,(2,8):16,(3,4):4,(3,10):11,(3,8):5,(3,27):12,(3,20):17,(3,6):6,(3,18):13,(3,15):18,(3,5):7,(3,14):14,(3,13):19}

class Mesh:
    
    def __init__(self,nodes=None, elements=None):
        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

        if elements is None:
            self.elements = []
        else:
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
        f.write('$EndElements\n')

        if fname:
            f.close()


    def read(self,fname):

        if not fname:
            print('Error: Must specify mesh file')
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


    def clearnodeconnectivities(self):

        for n in self.nodes:
            n.elements[:] = []


    def setnodeconnectivities(self):

        # we first clear the node connectivities if needed
        self.clearnodeconnectivities()

        for e in self.elements:
            for i in range(0,len(e.nodes)):
                ni = e.nodes[i]
                self.nodes[ni-1].elements.append(e.ID)


    def cracknodes(self, crack_id, cracktip_id=None):
        # Return the set of node IDs for nodes belonging to elements 
        # (curves in 2D, surfaces in 3D) of physical IDs in the list
        # "crack_id[]" but excluding those belonging to elements 
        # (points in 2D, curves in 3D) of physical IDs in the list
        # "cracktip_id[]".
    
        crack_nodes = []
        for e in self.elements:
            if (e.elem_ID in crack_id):
                for n in e.nodes:
                    crack_nodes.append(n)
    
        # remove node IDs for nodes on crack tip
        for e in self.elements:
            if (e.elem_ID in cracktip_id):
                for n in e.nodes:
                    if n in crack_nodes:
                        crack_nodes.remove(n)

        # remove duplicates
        crack_nodes = [n for n in set(crack_nodes)]
    
        return crack_nodes


    def makecrackedmesh(self, crack_id, oneside_id, cracktip_id=None):
    
        # make set of node IDs for nodes on crack
        crack_nodes = self.cracknodes(crack_id,cracktip_id)
        
        # add overlapping nodes
        new_nodes = []
        for i in range(0,len(crack_nodes)):
            new_nodes.append(self.nNodes() + i + 1)
        
        node_mapping = dict([(crack_nodes[i],new_nodes[i]) for i in range(0,len(crack_nodes))])
    
        # set nodes for elements in oneside_ID to the newly created nodes
        # Pay attention: elements in the oneside_ID region must touch 
        # the opposite side of the zero-thickness element in the orientation
        # given by the numbering of the surface element.
        for e in self.elements:
            if (e.elem_ID in oneside_id):
                for i in range(0,len(e.nodes)):
                    n = e.nodes[i]
                    if node_mapping.has_key(n):
                        e.nodes[i] = node_mapping[n]
            
        
        # Turn the surface crack elements into volume elements
        for e in self.elements:
            if (e.elem_ID in crack_id):
                e.turnsurfaceintovolume(node_mapping)
        
        
        # make new node objects with coordinates overlapping
        new_mesh_nodes = self.nodes
        for n in self.nodes:
            if (n.ID in crack_nodes):
                nodeid = node_mapping[n.ID]
                new_mesh_nodes.append(Node(nodeid,n.coords))
    
        return Mesh(new_mesh_nodes,self.elements)


    def allelem_ID(self):

        elem_ID = set()
        for e in self.elements:
            elem_ID.add(e.elem_ID)

        return elem_ID
        


    def smash(self,elem_ID=None):
        # Break apart the elements of the mesh and 
        # thus build unconnected set of elements with
        # new connectivity and new nodes created.
        # Make this for element IDs in the list "elem_ID[]" only.
        # If "elem_Id[]" is empty consider all elements.
        # Return a mapping for overlapping nodes.

        if elem_ID is None: elem_ID = self.allelem_ID()


        # Make new node objects for the smashed mesh
        # with coordinates overlapping.
        new_mesh_nodes = self.nodes
        # For each node 'n' of the initial mesh, map 'n' 
        # to the set of overlapping nodes (new nodes)
        # including n itself.
        overlapping_nodes_map = {}

        # First deal with elements not in the smashed region IDs
        for e in self.elements:
            if e.elem_ID in elem_ID: continue
            for ni in e.nodes:
                if not overlapping_nodes_map.has_key(ni):
                    overlapping_nodes_map[ni] = set((ni,))

        # Second deal with element in the smashed region IDs
        for e in self.elements:
            if not e.elem_ID in elem_ID: continue
            for ni in e.nodes:
                if overlapping_nodes_map.has_key(ni):
                    n = self.nodes[ni-1]
                    new_node_id = len(new_mesh_nodes) + 1
                    new_mesh_nodes.append(Node(new_node_id,n.coords))
                    val = overlapping_nodes_map.get(ni)
                    val.add(new_node_id)
                else:
                    overlapping_nodes_map[ni] = set((ni,))


        # Make a deep copy of the dictionary in order to return it unchanged
        overlapping_nodes_map_copy = copy.deepcopy(overlapping_nodes_map)


        # Update the connectivity of elements to the new mesh nodes
        for e in self.elements:
            if e.elem_ID in elem_ID: continue
            for ni in e.nodes:
                overlapping_nodes_map[ni].discard(ni)

        for e in self.elements:
            if not e.elem_ID in elem_ID: continue
            for i in range(0,len(e.nodes)):
                ni = e.nodes[i]
                e.nodes[i] = overlapping_nodes_map[ni].pop()


        return overlapping_nodes_map_copy


    def makebrokenmesh(self,phys_ID=-1,elem_ID=None):
        # Break apart the elements of a mesh and
        # put the pieces back together while inserting 
        # zero-thickness elements in between.

        if elem_ID is None: elem_ID = self.allelem_ID()

        # Set the node connectivities
        self.setnodeconnectivities()

        # Map the two-element set (ei,ej) to their common nodes ([n])
        common_nodes_map = {}
        for e_i in self.elements:
            ei = e_i.ID
            for ni in e_i.nodes:
                n_i = self.nodes[ni-1]
                for ej in n_i.elements:
                    if(ei >= ej): continue

                    e_j = self.elements[ej-1]

                    if (e_i.elem_ID in elem_ID) or (e_j.elem_ID in elem_ID):
                        key_ij = frozenset((ei,ej))
                        if common_nodes_map.has_key(key_ij):
                            val_ij = common_nodes_map.get(key_ij)
                            val_ij.add(n_i.ID)
                        else:
                            common_nodes_map[key_ij] = set((n_i.ID,))


        # Smash the mesh and map the nodes to their sets of overlapping nodes
        overlapping_nodes_map = self.smash(elem_ID)


        # Make new zero-thickness interface elements between elements (ei,ej)
        interface_elements = []
        for e_i in self.elements:
            ei = e_i.ID
            for e_j in self.elements:
                ej = e_j.ID
                if(ei >= ej): continue

                key_ij = frozenset((ei,ej))
                if not common_nodes_map.has_key(key_ij): continue
                val_ij = common_nodes_map.get(key_ij)

                # choose the master element e_k as a volume one
                if(e_i.dim() >= e_j.dim()):
                    e_k = e_i
                    e_l = e_j
                else:
                    e_k = e_j
                    e_l = e_i

                # map the node of e_k to that of e_l which overlaps it.
                n_mapping = {}
                for n_ij in val_ij:
                    sk = overlapping_nodes_map[n_ij].intersection(set(e_k.nodes))
                    sl = overlapping_nodes_map[n_ij].intersection(set(e_l.nodes))
                    if(len(sk) != 1 and len(sl) != 1):
                        print('nodes of e_k: ',e_k.nodes)
                        print('nodes of e_l: ',e_l.nodes)
                        print('overlapping nodes for node ',n_ij,': ',overlapping_nodes_map[n_ij])
                        sys.exit(1)
                    nk = sk.pop()
                    nl = sl.pop()
                    if(nk == nl):
                        print('n_mapping: ',n_mapping)
                        sys.exit(1)
                    n_mapping[nk] = nl


                # create a surface element (submanifold) with the right node ordering
                surfel_nodes = e_k.surfacenodeordering(n_mapping)

                if(len(surfel_nodes) == 0): continue

                # make a submanifold element between ei and ej
                dim = e_k.dim()
                new_eij = len(self.elements) + len(interface_elements) + 1
                eltp = gmshType[(dim-1,len(surfel_nodes))]
                # assign the elementary index of the surface element
                #elid = max([e_i.elem_ID,e_j.elem_ID]) + 1
                #elid = min([e_i.elem_ID,e_j.elem_ID])
                elid = (e_i.elem_ID + e_j.elem_ID) / 2
                tags = [phys_ID,elid]
                new_e = Element(new_eij,eltp,tags,surfel_nodes,phys_ID,elid)
                interface_elements.append(new_e)

                # turn the surface element into volume element
                new_e.turnsurfaceintovolume(n_mapping)


        new_mesh_elements = self.elements + interface_elements
        M = Mesh(self.nodes,new_mesh_elements)

        # Clear the node connectivities
        M.clearnodeconnectivities()

        return M

  
#===============================================================================
if __name__=="__main__":


    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'hs:c:t:f:r:p:',
                       ['help','side_crack=','crack=','tip_crack=',
                        'file=','region=','physical'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)



    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-c', '--crack'):
            crack_id = set([int(i) for i in optval.split(',')])
        elif opt in ('-t', '--tip_crack'):
            cracktip_id = set([int(i) for i in optval.split(',')])
        elif opt in ('-s','--side_crack'):
            oneside_id = set([int(i) for i in optval.split(',')])
        elif opt in ('-f', '--file'):
            file_name = optval
        elif opt in ('-r', '--region'):
            region_id = set([int(i) for i in optval.split(',')])
        elif opt in ('-p', '--physical'):
            phys_id = int(optval)


    if(not "file_name" in locals()):
        print('Error: Must specify mesh file with -f')
        sys.exit(1)


    M = Mesh()

    M.read(file_name)

    M.write()

    if("crack_id" in locals()):
        if(not "cracktip_id" in locals()): cracktip_id = set()
        if(not "oneside_id" in locals()): oneside_id = set()
        M1 = M.makecrackedmesh(crack_id, oneside_id, cracktip_id)
    elif("region_id" in locals()):
        if(not "phys_id" in locals()): phys_id = -1
    
        #print('\ninitial mesh:')
        #print('nb of initial nodes: ',M.nNodes())

        #print('\nsmashed mesh:')
        #nodes_map = M.smash(region_id)
        #print('nb of nodes of the smashed mesh: ',M.nNodes())
        #M.write()
        #for i in nodes_map.keys():
        #    print('overlapping nodes for node ',i,': ',nodes_map[i])

        print('\nbroken mesh:')
        M1 = M.makebrokenmesh(phys_id, region_id)
        M1.setnodeconnectivities()
        M1.write()
    else:
        print('Error: Must specify crack ID with -c or region ID with -r')
        sys.exit(1)
    
