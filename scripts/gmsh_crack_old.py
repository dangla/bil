#! /usr/bin/python

import sys
import getopt

#===============================================================================
class Node:
    def __init__(self,ID=-1,coords=[]):
        self.ID = ID
        self.coords=coords

#===============================================================================
class Element:
    def __init__(self,ID=-1,eltype=-1,tags=[],nodes=[],phys_ID=-1):
        self.ID = ID
        self.eltype = eltype
        self.tags = tags
        self.nodes = nodes
        self.phys_ID = phys_ID

#===============================================================================
class Mesh:
    
    def __init__(self,nodes=[], elements=[]):
        self.nodes = nodes
        self.elements = elements

    def nNodes(self):
        return len(self.nodes)
        
    def nElements(self):
        return len(self.elements)
    
    def write(self,fname):
        f = open(fname,'w')
        
        # write header
        f.write('$MeshFormat\n')
        f.write('2.2 0 8\n')
        f.write('$EndMeshFormat\n')
        
        # write nodes
        f.write('$Nodes\n')
        f.write(str(self.nNodes())+'\n')
        for n in self.nodes:
            line = str(n.ID) + ' ' + ' '.join([str(i) for i in n.coords])
            f.write(line+'\n')
        f.write('$EndNodes\n')

        #write elements
        f.write('$Elements\n')
        f.write(str(self.nElements())+'\n')
        for e in self.elements:
            e.tags[0] = e.phys_ID
            line = str(e.ID) + ' ' + str(e.eltype) + ' ' +str(len(e.tags)) + ' ' + ' '.join([str(i) for i in e.tags]) + ' ' + ' '.join([str(i) for i in e.nodes]) + '\n'
            f.write(line)
        f.write('$EndElements')

        f.close()

#===============================================================================
def main():
    try:
        optlist,args = getopt.getopt(sys.argv[1:],'ht:b:c:f:',
                                     ['help','top=','bottom=','crack=','file='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    top_id = -1
    bottom_id = -1
    crack_id = -1
    file_name = ""

    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-c', '--crack'):
            crack_id = int(optval)
        elif opt in ('-t','--top'):
            top_id = int(optval)
        elif opt in ('-b','--bottom'):
            bottom_id = int(optval)
        elif opt in ('-f', '--file'):
            file_name = optval

    if(top_id == -1):
        print 'Error: Must specify top_id with -t'
        sys.exit(1)
    if(bottom_id == -1):
        print 'Error: Must specify bottom_id with -b'
        sys.exit(1)
    if(crack_id == -1):
        print 'Error: Must specify crack_id with -c'
        sys.exit(1)
    if(not file_name):
        print 'Error: Must specify mesh file with -f'
        sys.exit(1)

    M = read_mesh(file_name)
    
    M = make_crack(M, crack_id, top_id, bottom_id)
    
    M.write('cracked_'+file_name)
#===============================================================================
def usage():
    print 'Usage: python %s -c<crack_id> -t<top_id> -b<bottom_id>' % sys.argv[0]
    print 'Options:'
    print '-c, --crack: Physical ID of crack in .msh file'
    print '-t, --top: Physical ID of surface touching one side of crack'
    print '-b, --bottom: Physical ID of surface touching other side of crack'
    print '-f, --file: Mesh file name'
    print '-h, --help: Display this help message'

#===============================================================================
def read_mesh(fname):
    
    nodes = []
    elements = []
    
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
            el_nodes = [int(t) for t in tmp[3+ntags:]]
            physid = int(tags[0])

            elements.append(Element(i,eltp,tags,el_nodes,physid))

    return Mesh(nodes,elements)

#===============================================================================
def make_crack(M, crack_id, top_id, bottom_id):
    
    crack_elems = []
    for e in M.elements:
        if (int(e.phys_ID) == int(crack_id)):
            crack_elems.append(e)

    # make set of node IDs for nodes on crack
    crack_nodes = []
    for e in crack_elems:
        for n in e.nodes:
            crack_nodes.append(n)

    crack_nodes = [n for n in set(crack_nodes)]
    
    new_nodes = []
    for i in range(1,len(crack_nodes)+1):
        new_nodes.append(M.nNodes()+i)
    
    node_mapping = dict([(crack_nodes[i],new_nodes[i]) for i in range(0,len(crack_nodes))])

    # set nodes in phyisical_ID=top_ID to the newly created nodes
    for ei in range(0,len(M.elements)):
        if(M.elements[ei].phys_ID == top_id):
            for i in range(0,len(M.elements[ei].nodes)):
                if(node_mapping.has_key(M.elements[ei].nodes[i])):
                    M.elements[ei].nodes[i] = node_mapping[M.elements[ei].nodes[i]]
        
    
    # create new elements in crack
    new_elems = []
    for i in range(0,len(crack_elems)):
        e = crack_elems[i]
        n = [e.nodes[0], e.nodes[1], node_mapping[e.nodes[1]], node_mapping[e.nodes[0]]]
        newid = M.nElements() + i + 1
        eltp = 3 # <-- 4-node quad element type
        physid = crack_id+1 # put the "filler" elements at phys_ID = crack_id + 1 (no checking is done if this collides with the input!)
        tags = e.tags
        tags[0] = physid
        new_elems.append(Element(newid,eltp,tags,n,physid))
    
    
    # make new node objects with coordinates overlapping
    new_mesh_nodes = M.nodes
    for n in M.nodes:
        if (n.ID in crack_nodes):
            new_mesh_nodes.append(Node(node_mapping[n.ID],n.coords))

    # make new mesh with crack, without old crack elements
    new_mesh_elems = M.elements
    for e in new_elems:
        new_mesh_elems.append(e)
    
    return Mesh(new_mesh_nodes,new_mesh_elems)
    
#===============================================================================
if __name__=="__main__":
    main()
