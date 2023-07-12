#! /usr/bin/python
# Usage:
# python3 gmsh_pos_to_paraview.py -i input.pos -o output.vtu


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


    gmsh_pos_to_paraview(gmsh_file, paraview_file)


#===============================================================================
def usage():
    print('Usage: python3 %s -i<gmsh_file> -o<paraview_file>' % sys.argv[0])

    print('\n')

    print('Options:')
    print('-i, --input:      Gmsh pos file name.')
    print('-o, --output:     Paraview vtu file name.')
    print('-h, --help:       Display this help message.')

    # Usage example
    #gmsh_post_file = "input.pos"
    #paraview_file = "output.vtk"
    #gmsh_post_to_paraview(gmsh_post_file, paraview_file)



def gmsh_pos_to_paraview(gmsh_post_file, paraview_file):
    # Read the Gmsh post-processing file
    mesh = meshio.read(gmsh_post_file)
    
    # Extract node data
    node_data = mesh.point_data
    
    # Write the Paraview file
    meshio.write_points_cells(
        paraview_file,
        points=mesh.points,
        cells={"unstructured": mesh.cells["triangle"]},  # Assuming triangle cells
        point_data=node_data,
        file_format="vtk"
    )



#===============================================================================
if __name__=="__main__":
    main()

