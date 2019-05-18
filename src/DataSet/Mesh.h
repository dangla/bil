#ifndef MESH_H
#define MESH_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Mesh_s         ; typedef struct Mesh_s         Mesh_t ;


#include "Materials.h"
#include "Geometry.h"
#include "BConds.h"
#include "DataFile.h"
#include "Solutions.h"
#include "Solver.h"

extern Mesh_t*  (Mesh_Create)(DataFile_t*,Materials_t*,Geometry_t*) ;
extern char*    (Mesh_Scan)(Mesh_t*,char*) ;
extern void     (Mesh_SetMatrixPermutationNumbering)(Mesh_t*,BConds_t*,DataFile_t*) ;
extern void     (Mesh_ResetMatrixNumbering)(Mesh_t*) ;
extern void     (Mesh_WriteGraph)(Mesh_t*,const char*,const char*) ;
extern void     (Mesh_WriteInversePermutation)(Mesh_t*,const char*,const char*) ;
extern void     (Mesh_InitializeSolutionPointers)(Mesh_t*,Solutions_t*) ;
extern int      (Mesh_LoadCurrentSolution)(Mesh_t*,DataFile_t*,double*) ;
extern int      (Mesh_StoreCurrentSolution)(Mesh_t*,DataFile_t*,double) ;
extern void     (Mesh_SetCurrentUnknownsWithBoundaryConditions)(Mesh_t*,BConds_t*,double) ;
extern void     (Mesh_UpdateCurrentUnknowns)(Mesh_t*,Solver_t*) ;
extern void     (Mesh_CreateEquationContinuity)(Mesh_t*) ;


#include "Mry.h"

#define Mesh_New() \
        Mry_New(Mesh_t)


/* Some constants */
#define Mesh_MaxLengthOfKeyWord        (30)
#define Mesh_MaxLengthOfFileName       (60)
#define Mesh_MaxLengthOfTextLine       (500)



/* Accessors */
#define Mesh_GetGeometry(MSH)               ((MSH)->geometry)
#define Mesh_GetNodes(MSH)                  ((MSH)->nodes)
#define Mesh_GetElements(MSH)               ((MSH)->elements)



/* Access to the geometry characteristics */
#define Mesh_GetDimension(MSH) \
        Geometry_GetDimension(Mesh_GetGeometry(MSH))

#define Mesh_GetSymmetry(MSH) \
        Geometry_GetSymmetry(Mesh_GetGeometry(MSH))

#define Mesh_GetCoordinateSystem(MSH) \
        Geometry_GetCoordinateSystem(Mesh_GetGeometry(MSH))



/* Access to Nodes */
#define Mesh_GetNbOfNodes(MSH) \
        Nodes_GetNbOfNodes(Mesh_GetNodes(MSH))

#define Mesh_GetNode(MSH) \
        Nodes_GetNode(Mesh_GetNodes(MSH))



/* Access to the elements */
#define Mesh_GetNbOfElements(MSH) \
        Elements_GetNbOfElements(Mesh_GetElements(MSH))

#define Mesh_GetElement(MSH) \
        Elements_GetElement(Mesh_GetElements(MSH))



/* Access to the nb of matrix rows/columns */
#define Mesh_GetNbOfMatrixColumns(MSH) \
        Nodes_GetNbOfMatrixColumns(Mesh_GetNodes(MSH))

#define Mesh_GetNbOfMatrixRows(MSH) \
        Nodes_GetNbOfMatrixRows(Mesh_GetNodes(MSH))



/* Compute the nb of matrix entries */
#define Mesh_ComputeNbOfMatrixEntries(MSH) \
        Elements_ComputeNbOfMatrixEntries(Mesh_GetElements(MSH))




/* Periodicities */
#define Mesh_GetPeriodicities(MSH) \
        Geometry_GetPeriodicities(Mesh_GetGeometry(MSH))

#define Mesh_IsPeriodic(MSH) \
        Mesh_GetPeriodicities(MSH)



#include "Elements.h"
#include "Nodes.h"

struct Mesh_s {
  Geometry_t* geometry ;
  Elements_t* elements ;
  Nodes_t*    nodes ;
} ;

#endif
