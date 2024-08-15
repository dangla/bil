#ifndef ELEMENT_H
#define ELEMENT_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct Element_s ;
typedef struct Element_s  Element_t ;
typedef Element_t*        Element_tt ;


#include "Mesh.h"
#include "Node.h"
#include "ShapeFcts.h"
#include "IntFcts.h"
#include "Buffers.h"

extern Element_t*  (Element_New)(void) ;
extern void        (Element_CreateMore)(Element_t*,Buffers_t*,ShapeFcts_t*,IntFcts_t*) ;
extern void        (Element_Delete)(void*) ;
extern void        (Element_AllocateMicrostructureSolutions)       (Element_t*,Mesh_t*,const int) ;
extern double**    (Element_ComputePointerToCurrentNodalUnknowns)  (Element_t*) ;
extern double**    (Element_ComputePointerToPreviousNodalUnknowns) (Element_t*) ;
extern double*     (Element_ComputeCurrentNodalUnknowns)           (Element_t*) ;
extern double*     (Element_ComputePreviousNodalUnknowns)          (Element_t*) ;
extern double*     (Element_ComputeIncrementalNodalUnknowns)       (Element_t*) ;
extern double*     (Element_ComputeDeepNodalUnknowns)              (Element_t*,unsigned int) ;
extern double**    (Element_ComputePointerToNodalCoordinates)      (Element_t*) ;
extern double*     (Element_ComputeNodalCoordinates)               (Element_t*) ;
extern int         (Element_FindUnknownPositionIndex)              (Element_t*,const char*) ;
extern int         (Element_FindEquationPositionIndex)             (Element_t*,const char*) ;
extern double*     (Element_ComputeIncrementalImplicitTerms)       (Element_t*) ;
extern double*     (Element_ComputeNormalVector)                   (Element_t*,double*,int,const int) ;
extern double*     (Element_ComputeCoordinateInReferenceFrame)     (Element_t*,double*) ;
extern int         (Element_ComputeNbOfSolutions)                  (Element_t*) ;
extern int*        (Element_ComputeMatrixRowAndColumnIndices)      (Element_t*) ;
extern int*        (Element_ComputeSelectedMatrixRowAndColumnIndices)(Element_t*,const int) ;
extern double      (Element_ComputeSize)                           (Element_t*) ;
extern double*     (Element_ComputeSizes)                          (Element_t*) ;
extern int         (Element_ComputeNbOfMatrixEntries)              (Element_t*) ;
extern int         (Element_ComputeNbOfSelectedMatrixEntries)      (Element_t*,const int) ;
extern double*     (Element_ComputeJacobianMatrix)                 (Element_t*,double*,int,const int) ;
extern double      (Element_ComputeJacobianDeterminant)            (Element_t*,double*,int,const int) ;
extern double*     (Element_ComputeInverseJacobianMatrix)          (Element_t*,double*,int,const int) ;
extern int         (Element_OverlappingNode)                       (Element_t*,const int) ;
extern int         (Element_HasZeroThickness)                      (Element_t*) ;
extern int         (Element_NbOfOverlappingNodes)                  (Element_t*) ;
extern void        (Element_MakeUnknownContinuousAcrossZeroThicknessElement)(Element_t*,const char*);
extern void        (Element_MakeEquationContinuousAcrossZeroThicknessElement)(Element_t*,const char*);
extern int         (Element_FindNodeIndex)                         (Element_t*,const Node_t*) ;
extern double*     (Element_ComputeCoordinateVector)               (Element_t*,double*) ;
extern void        (Element_CopyCurrentSolutionIntoPreviousSolution)(Element_t*) ;

/* Synonyms */
#define  Element_ComputePointerToNodalUnknowns \
         Element_ComputePointerToCurrentNodalUnknowns



/* Some constants */
#define Element_MaxNbOfNodes  (8)

#define Element_MaxNbOfDOF \
        (Element_MaxNbOfNodes*Model_MaxNbOfEquations)

#define Element_SizeOfBuffer \
        (IntFct_MaxNbOfIntPoints*Element_MaxNbOfNodes*100*sizeof(double))
        
//#define Element_MaxLengthOfRegionName  (50)



/* Accessors */
#define Element_GetElementIndex(ELT)           ((ELT)->ElementIndex)
#define Element_GetDimension(ELT)              ((ELT)->Dimension)
#define Element_GetNbOfNodes(ELT)              ((ELT)->NbOfNodes)
#define Element_GetPointerToNode(ELT)          ((ELT)->PointerToNode)
//#define Element_GetRegionTag(ELT)              ((ELT)->RegionTag)
#define Element_GetRegion(ELT)                 ((ELT)->Region)
#define Element_GetMaterial(ELT)               ((ELT)->Material)
#define Element_GetMaterialIndex(ELT)          ((ELT)->MaterialIndex)
#define Element_GetShapeFct(ELT)               ((ELT)->ShapeFct)
#define Element_GetIntFct(ELT)                 ((ELT)->IntFct)
#define Element_GetUnknownPosition(ELT)        ((ELT)->UnknownPosition)
#define Element_GetEquationPosition(ELT)       ((ELT)->EquationPosition)
#define Element_GetBuffers(ELT)                ((ELT)->Buffers)
#define Element_GetSolutions(ELT)              ((ELT)->Solutions)
#define Element_GetMatrix(ELT)                 ((ELT)->Matrix)
#define Element_GetResidu(ELT)                 ((ELT)->Residu)



/* Buffer */
#define Element_GetBuffer(ELT) \
        Buffers_GetBufferOfCurrentThread(Element_GetBuffers(ELT))



/* Nb of (im/ex)plicit and constant terms */
#define Element_GetNbOfImplicitTerms(ELT) \
        ElementSol_GetNbOfImplicitTerms(Element_GetElementSol(ELT))

#define Element_GetNbOfExplicitTerms(ELT) \
        ElementSol_GetNbOfExplicitTerms(Element_GetElementSol(ELT))

#define Element_GetNbOfConstantTerms(ELT) \
        ElementSol_GetNbOfConstantTerms(Element_GetElementSol(ELT))



/* Access to solution */
#define Element_GetSolution(ELT) \
        Solutions_GetSolution(Element_GetSolutions(ELT))
        
#define Element_GetPreviousSolution(ELT) \
        Solution_GetPreviousSolution(Element_GetSolution(ELT))



/* Access to elementsol */
#define Element_GetElementSol(ELT) \
        (Solution_GetElementSol(Element_GetSolution(ELT)) + Element_GetElementIndex(ELT))

#define Element_GetPreviousElementSol(ELT) \
        (Solution_GetElementSol(Element_GetPreviousSolution(ELT)) + Element_GetElementIndex(ELT))



/* Access to (im/ex)plicit and constant terms */
#define Element_GetCurrentImplicitTerm(ELT) \
        ((double*) ElementSol_GetImplicitTerm(Element_GetElementSol(ELT)))

#define Element_GetCurrentExplicitTerm(ELT) \
        ((double*) ElementSol_GetExplicitTerm(Element_GetElementSol(ELT)))

#define Element_GetPreviousImplicitTerm(ELT) \
        ((double*) ElementSol_GetImplicitTerm(Element_GetPreviousElementSol(ELT)))

#define Element_GetPreviousExplicitTerm(ELT) \
        ((double*) ElementSol_GetExplicitTerm(Element_GetPreviousElementSol(ELT)))

#define Element_GetConstantTerm(ELT) \
        ((double*) ElementSol_GetConstantTerm(Element_GetElementSol(ELT)))

#define Element_GetImplicitTerm(ELT) \
        Element_GetCurrentImplicitTerm(ELT)
        
#define Element_GetExplicitTerm(ELT) \
        Element_GetCurrentExplicitTerm(ELT)




/* Access to generic data */
#define Element_GetCurrentImplicitGenericData(ELT) \
        ElementSol_GetImplicitGenericData(Element_GetElementSol(ELT))
        
#define Element_GetCurrentExplicitGenericData(ELT) \
        ElementSol_GetExplicitGenericData(Element_GetElementSol(ELT))
        
#define Element_GetPreviousImplicitGenericData(ELT) \
        ElementSol_GetImplicitGenericData(Element_GetPreviousElementSol(ELT))
        
#define Element_GetPreviousExplicitGenericData(ELT) \
        ElementSol_GetExplicitGenericData(Element_GetPreviousElementSol(ELT))
        
#define Element_GetConstantGenericData(ELT) \
        ElementSol_GetConstantGenericData(Element_GetElementSol(ELT))



/* Access to data */
#define Element_FindCurrentImplicitData(ELT,...) \
        ElementSol_FindImplicitData(Element_GetElementSol(ELT),__VA_ARGS__)
        
#define Element_FindCurrentExplicitData(ELT,...) \
        ElementSol_FindExplicitData(Element_GetElementSol(ELT),__VA_ARGS__)
        
#define Element_FindPreviousImplicitData(ELT,...) \
        ElementSol_FindImplicitData(Element_GetPreviousElementSol(ELT),__VA_ARGS__)
        
#define Element_FindPreviousExplicitData(ELT,...) \
        ElementSol_FindExplicitData(Element_GetPreviousElementSol(ELT),__VA_ARGS__)
        
#define Element_FindConstantData(ELT,...) \
        ElementSol_FindConstantData(Element_GetElementSol(ELT),__VA_ARGS__)



/* Access to nodes and node coordinates */
#define Element_GetNode(ELT,i) \
        (Element_GetPointerToNode(ELT)[i])

#define Element_GetNodeCoordinate(ELT,i) \
        Node_GetCoordinate(Element_GetNode(ELT,i))
        
#define Element_GetNodeIndex(ELT,i) \
        Node_GetNodeIndex(Element_GetNode(ELT,i))



/* Access to the values of nodal unknowns */
#define Element_GetNodalUnknownPosition(ELT,n,i) \
        (Element_GetUnknownPosition(ELT)[(n)*Element_GetNbOfEquations(ELT) + (i)])

#define Element_GetCurrentNodalUnknown(ELT,i) \
        Node_GetCurrentUnknown(Element_GetNode(ELT,i))

#define Element_GetPreviousNodalUnknown(ELT,i) \
        Node_GetPreviousUnknown(Element_GetNode(ELT,i))
        
#define Element_GetValueOfCurrentNodalUnknown(ELT,n,i) \
        (Element_GetCurrentNodalUnknown(ELT,n)[Element_GetNodalUnknownPosition(ELT,n,i)])
       
#define Element_GetValueOfPreviousNodalUnknown(ELT,n,i) \
        (Element_GetPreviousNodalUnknown(ELT,n)[Element_GetNodalUnknownPosition(ELT,n,i)])

#define Element_GetValueOfNodalUnknown(ELT,u,n,i) \
        ((u)[n][Element_GetNodalUnknownPosition(ELT,n,i)])



/* Objective values */
#define Element_GetObjectiveValue(ELT) \
        Model_GetObjectiveValue(Element_GetModel(ELT))



/* Geometrical characteristics */
#define Element_GetGeometry(ELT) \
        Model_GetGeometry(Element_GetModel(ELT))

#define Element_GetSymmetry(ELT) \
        Geometry_GetSymmetry(Element_GetGeometry(ELT))

#define Element_GetCoordinateSystem(ELT) \
        Geometry_GetCoordinateSystem(Element_GetGeometry(ELT))

#define Element_GetDimensionOfSpace(ELT) \
        Geometry_GetDimension(Element_GetGeometry(ELT))



/* Material properties */
#define Element_GetProperty(ELT) \
        Material_GetProperty(Element_GetMaterial(ELT))
        
#define Element_GetPropertyValue(ELT,S) \
        Material_GetPropertyValue(Element_GetMaterial(ELT),S)
        
#define Element_GetGenericProperty(ELT) \
        Material_GetGenericData(Element_GetMaterial(ELT))

#define Element_GetCurves(ELT) \
        Material_GetCurves(Element_GetMaterial(ELT))

#define Element_GetCurve(ELT) \
        Material_GetCurve(Element_GetMaterial(ELT))
        
#define Element_GetModel(ELT) \
        Material_GetModel(Element_GetMaterial(ELT))

#define Element_FindCurve(ELT,S) \
        Material_FindCurve(Element_GetMaterial(ELT),S)
        
#define Element_FindMaterialData(ELT,T,N) \
        Material_FindData(Element_GetMaterial(ELT),T,N)
        
#define Element_FindMaterialGenericData(ELT,T,N) \
        Material_FindGenericData(Element_GetMaterial(ELT),T,N)




/* Equations, unknowns and degrees of freedom (dof) per element */
#define Element_GetNbOfEquations(ELT) \
        Material_GetNbOfEquations(Element_GetMaterial(ELT))

#define Element_GetNbOfDOF(ELT) \
        (Element_GetNbOfNodes(ELT)*Element_GetNbOfEquations(ELT))

#define Element_GetNameOfEquation(ELT) \
        Material_GetNameOfEquation(Element_GetMaterial(ELT))

#define Element_GetNameOfUnknown(ELT) \
        Material_GetNameOfUnknown(Element_GetMaterial(ELT))



/* Methods and procedures */
#define Element_DefineProperties(ELT,...) \
        (Model_GetDefineElementProperties(Element_GetModel(ELT)))(ELT,__VA_ARGS__)
        
#define Element_ComputeInitialState(ELT,...) \
        (Model_GetComputeInitialState(Element_GetModel(ELT)))(ELT,__VA_ARGS__)
        
#define Element_ComputeExplicitTerms(ELT,...) \
        (Model_GetComputeExplicitTerms(Element_GetModel(ELT)))(ELT,__VA_ARGS__)

#define Element_ComputeImplicitTerms(ELT,...) \
        (Model_GetComputeImplicitTerms(Element_GetModel(ELT)))(ELT,__VA_ARGS__)

#define Element_ComputeMatrix(ELT,...) \
        (Model_GetComputeMatrix(Element_GetModel(ELT)))(ELT,__VA_ARGS__)

#define Element_ComputeResidu(ELT,...) \
        (Model_GetComputeResidu(Element_GetModel(ELT)))(ELT,__VA_ARGS__)

#define Element_ComputeOutputs(ELT,...) \
        (Model_GetComputeOutputs(Element_GetModel(ELT)))(ELT,__VA_ARGS__)

#define Element_ComputeLoads(ELT,...) \
        (Model_GetComputeLoads(Element_GetModel(ELT)))(ELT,__VA_ARGS__)



/* Access to local variables */
#define Element_GetNbOfVariables(ELT) \
        Model_GetNbOfVariables(Element_GetModel(ELT))
        
#define Element_GetVariable(ELT,n) \
        Model_GetVariable(Element_GetModel(ELT),n)
        
#define Element_GetVariableDerivative(ELT,n) \
        Model_GetVariableDerivative(Element_GetModel(ELT),n)


/* Submanifold */
#define Element_IsSubmanifold(ELT) \
        (Element_GetDimension(ELT) < Element_GetDimensionOfSpace(ELT))



/* Operations on buffer */
#define Element_AllocateInBuffer(ELT,sz) \
        Buffer_Allocate(Element_GetBuffer(ELT),(sz))
        
#define Element_FreeBuffer(ELT)  \
        Buffer_Free(Element_GetBuffer(ELT))
        
#define Element_FreeBufferFrom(ELT,p) \
        Buffer_FreeFrom(Element_GetBuffer(ELT),(char*) (p))



/* Access to datafile */
#define Element_GetDataFile(ELT) \
        Model_GetDataFile(Element_GetModel(ELT))


/* Access to dataset */
#define Element_GetDataSet(ELT) \
        DataFile_GetDataSet(Element_GetDataFile(ELT))


/* Access to module */
#define Element_GetModule(ELT) \
        DataSet_GetModule(Element_GetDataSet(ELT))


/* Method */
#include "String_.h"

#define Element_MethodIs(ELT,MTH) \
        String_Is(Material_GetMethod(Element_GetMaterial(ELT)),MTH)



/* Sequential index */
#define Element_GetSequentialIndexOfUnknown(ELT) \
        Model_GetSequentialIndexOfUnknown(Element_GetModel(ELT))

#define Element_GetCurrentSequentialIndex(ELT) \
        DataFile_GetSequentialIndex(Element_GetDataFile(ELT))
        
#define Element_GetNbOfSequences(ELT) \
        DataFile_GetNbOfSequences(Element_GetDataFile(ELT))

#define Element_EquationIsActive(ELT,I) \
        (( \
        (Element_GetCurrentSequentialIndex(ELT) == Element_GetNbOfSequences(ELT) - 1) \
        && \
        (Element_GetCurrentSequentialIndex(ELT) <= Element_GetSequentialIndexOfUnknown(ELT)[I]) \
        ) || ( \
        (Element_GetCurrentSequentialIndex(ELT) <= Element_GetNbOfSequences(ELT) - 1) \
        && \
        (Element_GetCurrentSequentialIndex(ELT) == Element_GetSequentialIndexOfUnknown(ELT)[I]) \
        ))

#define Element_CurrentSequentialIndexIs \
        Element_EquationIsActive



/* The time, time step and step index */
#define Element_GetTimeOf(ELT,SOL) \
        Solution_GetSequentialTime(SOL)[Element_GetCurrentSequentialIndex(ELT)]
        
#define Element_GetTimeStepOf(ELT,SOL) \
        Solution_GetSequentialTimeStep(SOL)[Element_GetCurrentSequentialIndex(ELT)]
        
#define Element_GetStepIndexOf(ELT,SOL) \
        Solution_GetSequentialStepIndex(SOL)[Element_GetCurrentSequentialIndex(ELT)]


#define Element_GetTime(ELT) \
        Element_GetTimeOf(ELT,Element_GetSolution(ELT))
        
#define Element_GetTimeStep(ELT) \
        Element_GetTimeStepOf(ELT,Element_GetSolution(ELT))
        
#define Element_GetStepIndex(ELT) \
        Element_GetStepIndexOf(ELT,Element_GetSolution(ELT))
        
#define Element_GetPreviousTime(ELT) \
        Element_GetTimeOf(ELT,Element_GetPreviousSolution(ELT))
        
#define Element_GetPreviousTimeStep(ELT) \
        Element_GetTimeStepOf(ELT,Element_GetPreviousSolution(ELT))
        
#define Element_GetPreviousStepIndex(ELT) \
        Element_GetStepIndexOf(ELT,Element_GetPreviousSolution(ELT))




#include "IntFct.h"

/* Loop on integration points */
#define Element_LoopOnIntegrationPoints(ELT,p) \
        for(p = 0 ; p < IntFct_GetNbOfPoints(Element_GetIntFct(ELT)) ; p++)
        

#include "DistributedMS.h"

/* Supporting processor */
#define Element_RankOfSupportingProcessor(ELT) \
        (Element_GetElementIndex(ELT) % DistributedMS_NbOfProcessors)


/* Region name */
#define Element_GetRegionName(ELT) \
        Region_GetRegionName(Element_GetRegion(ELT))


/* Set a default region to an element  */
#define Element_SetDefaultRegion(ELT,REGS,I) \
        do { \
          char Element_name[Region_MaxLengthOfRegionName] ; \
          sprintf(Element_name,"%d",I) ; \
          Element_GetRegion(ELT) = Regions_FindRegion(REGS,Element_name) ; \
        } while(0)
        

#include "ShapeFct.h"
#include "Material.h"
#include "Buffers.h"
//#include "ElementSol.h"
#include "Solutions.h"
#include "Region.h"


/* Minmize the size of the struct by re-organizing the order of its members*/
struct Element_s {
  Node_tt* PointerToNode ;
  Material_t* Material ;
  ShapeFct_t* ShapeFct ;      /* Shape function */
  IntFct_t* IntFct ;          /* interpolation function */
  Buffers_t*  Buffers ;
  Solutions_t* Solutions ;    /* Pointer to the global solutions */
  Region_t*  Region ;
  double* Matrix ;            /* Elementary matrix */
  double* Residu ;            /* Elementary residu */
  short int* UnknownPosition ;  /* local position of unknowns at nodes */
  short int* EquationPosition ; /* local position of equations at nodes */
  //int    RegionTag ;
  int    MaterialIndex ;
  unsigned short int Dimension ;    /* dimension of the element (0,1,2,3) */
  unsigned int ElementIndex ;
  unsigned short int NbOfNodes ;
  /* n_vi and n_ve must be kept for old methods! */
  unsigned int n_vi ;         /* Nb of implicit terms */
  unsigned int n_ve ;         /* Nb of explicit terms */
//  unsigned int nbofconstterms ;   /* Nb of constant terms */
} ;


/* Old notations which I try to eliminate little by little */
#define elem_t         Element_t
#define MAX_NOEUDS     Element_MaxNbOfNodes


#ifdef __CPLUSPLUS
}
#endif
#endif
