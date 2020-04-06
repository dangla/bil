#! /usr/bin/python

import sys
import getopt
from gmshMesh import Mesh

#===============================================================================
def main():
    try:
        optlist,args = getopt.getopt(sys.argv[1:],'hs:c:t:f:',
                                     ['help','side_crack=','crack=','tip_crack=','file='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    oneside_id = -1
    crack_id = -1
    cracktip_id = -1
    file_name = ""

    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-c', '--crack'):
            crack_id = int(optval)
        elif opt in ('-t', '--tip_crack'):
            cracktip_id = int(optval)
        elif opt in ('-s','--side_crack'):
            oneside_id = int(optval)
        elif opt in ('-f', '--file'):
            file_name = optval

    if(not file_name):
        print 'Error: Must specify mesh file with -f'
        sys.exit(1)

    nodes = []
    elements = []
    M = Mesh(nodes,elements)

    M.read(file_name)
    
    M1 = M.makecrackedmesh(crack_id, oneside_id, cracktip_id)
    
    M1.write('cracked_'+file_name)


#===============================================================================
def usage():
    print 'Usage: python %s -f<filename> -c<crack_id> -t<cracktip_id> -o<oneside_id>' % sys.argv[0]
    print 'Options:'
    print '-c, --crack:      Elementary ID of crack in .msh file'
    print '-t, --tip_crack:  Elementary ID of crack tip in .msh file'
    print '-s, --side_crack: Elementary ID of surface touching one side of crack'
    print '-f, --file:       Mesh file name'
    print '-h, --help:       Display this help message'

    
#===============================================================================
if __name__=="__main__":
    main()
