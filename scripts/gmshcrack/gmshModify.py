#! /usr/bin/python
# Usage:
# python gmshModify.py -f filename -r region_id -p physical_id

import os
import sys
import getopt
from gmshMesh import Mesh

#===============================================================================
def main():
    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'hs:f:r:p:',
                       ['help','file=',
                       'region=','physical'])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
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
    

    if("region_id" in locals()):
        if(not "phys_id" in locals()): phys_id = -1
        M.changephysicalid(phys_id, region_id)
    else:
        print('Error: Must specify region ID with -r')
        sys.exit(1)
    

    M.write('modified_'+file_name)


#===============================================================================
def usage():
    print('Usage (with python 2 only):')
    print('python2 %s -f<filename> -r<region_ids> -p<physical_id>' % os.path.basename(sys.argv[0]))
    
    print('\n')

    print('Options:')
    print('-r, --region:        Comma separated elementary IDs of elements to be modified.')
    print('-p, --physical:      Physical ID to be modified.')
    print('-f, --file:          Mesh file name.')
    print('-h, --help:          Display this help message.')

    print('\n')

    print('Example:')
    print('python %s -f filename -r region_id -p physical_id' % os.path.basename(sys.argv[0]))

    
#===============================================================================
if __name__=="__main__":
    main()
