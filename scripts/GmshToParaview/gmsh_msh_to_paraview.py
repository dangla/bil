#! /usr/bin/python
# Usage:
# python gmsh_msh_to_paraview.py -i input.msh -o output.vtu


import sys
import getopt
import meshio


#===============================================================================
def main():
    try:
        optlist,args = getopt.getopt(sys.argv[1:],
                       'hi:o:',
                       ['help','input=','output='])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(2)

    for opt, optval in optlist:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            gmsh_file = optval
        elif opt in ('-o', '--output'):
            paraview_file = optval


    if(not "gmsh_file" in locals()):
        print('Error: Must specify mesh file with -i')
        sys.exit(1)


    if(not "paraview_file" in locals()):
        print('Error: Must specify paraview file with -o')
        sys.exit(1)


    gmsh_msh_to_paraview(gmsh_file, paraview_file)


#===============================================================================
def usage():
    print('Usage: python3 %s -i<gmsh_file> -o<paraview_file>' % sys.argv[0])

    print('\n')

    print('Options:')
    print('-i, --input:      Gmsh mesh file name.')
    print('-o, --output:     Paraview mesh file name.')
    print('-h, --help:       Display this help message.')

    # Usage example
    #gmsh_file = "input.msh"
    #paraview_file = "output.vtu"
    #gmsh1_msh_to_paraview(gmsh_file, paraview_file)
    #gmsh2_msh_to_paraview(gmsh_file, paraview_file)



def gmsh_msh_to_paraview(gmsh_pos_file, paraview_vtu_file):
    # Read Gmsh post-processing file
    mesh = meshio.read(gmsh_pos_file)

    # Write Paraview file
    meshio.write(paraview_vtu_file, mesh, file_format="vtk", binary=False)

    print("Conversion complete.")



#===============================================================================
if __name__=="__main__":
    main()
