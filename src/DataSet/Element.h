#ifndef ELEMENT_H
#define ELEMENT_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Element_s      ; typedef struct Element_s      Element_t ;



extern double** (Element_ComputePointerToCurrentNodalUnknowns)(Element_t*) ;
extern double** (Element_ComputePointerToPreviousNodalUnknowns)(Element_t*) ;
extern double*  (Element_ComputeCurrentNodalUnknowns)(Element_t*) ;
extern double*  (Element_ComputePreviousNodalUnknowns)(Element_t*) ;
extern double*  (Element_ComputeIncrementalNodalUnknowns)(Element_t*) ;
extern double*  (Element_ComputeDeepNodalUnknowns)(Element_t*,unsigned int) ;
extern double** (Element_ComputePointerToNodalCoordinates)(Element_t*) ;
extern double*  (Element_ComputeNodalCoordinates)(Element_t*) ;
extern int      (Element_FindUnknownPositionIndex)(Element_t*,char*) ;
extern int      (Element_FindEquationPositionIndex)(Element_t*,char*) ;
extern double*  (Element_ComputeIncrementalImplicitTerms)(Element_t*) ;
extern double*  (Element_ComputeNormalVector)(Element_t*,double*,int) ;
extern double*  (Element_ComputeCoordinateInReferenceFrame)(Element_t*,double*) ;
extern int      (Element_ComputeNbOfSolutions)(Element_t*) ;


/* Synonyms */
#define  Element_ComputePointerToNodalUnknowns \
         Element_ComputePointerToCurrentNodalUnknowns



/* Some constants */
#define Element_MaxNbOfNodes  (8)

#define Element_MaxNbOfDOF \
        (Element_MaxNbOfNodes*Model_MaxNbOfEquations)

#define Element_SizeOfBuffer \
        (IntFct_MaxNbOfIntPoints*Element_MaxNbOfNodes*100*sizeof(double))



/* Accessors */
#define Element_GetElementIndex(ELT)           ((ELT)->index)
#define Element_GetDimension(ELT)              ((ELT)->dim)
#define Element_GetNbOfNodes(ELT)              ((ELT)->nn)
#define Element_GetPointerToNode(ELT)          ((ELT)->node)
#define Element_GetRegionIndex(ELT)            ((ELT)->reg)
#define Element_GetMaterial(ELT)               ((ELT)->mat)
#define Element_GetMaterialIndex(ELT)          ((ELT)->imat)
#define Element_GetShapeFct(ELT)               ((ELT)->shapefct)
#define Element_GetIntFct(ELT)                 ((ELT)->fi)
#define Element_GetUnknownPosition(ELT)        ((ELT)->pin)
#define Element_GetEquationPosition(ELT)       ((ELT)->peq)
//#define Element_GetNbOfImplicitTerms(ELT)      ((ELT)->n_vi)
//#define Element_GetNbOfExplicitTerms(ELT)      ((ELT)->n_ve)
//#define Element_GetNbOfConstantTerms(ELT)      ((ELT)->nbofconstterms)
#define Element_GetBuffer(ELT)                 ((ELT)->buffer)
#define Element_GetElementSol(ELT)             ((ELT)->sol)




/* Nb of (im/ex)plicit and constant terms */
#define Element_GetNbOfImplicitTerms(ELT) \
        ElementSol_GetNbOfImplicitTerms(Element_GetElementSol(ELT))

#define Element_GetNbOfExplicitTerms(ELT) \
        ElementSol_GetNbOfExplicitTerms(Element_GetElementSol(ELT))

#define Element_GetNbOfConstantTerms(ELT) \
        ElementSol_GetNbOfConstantTerms(Element_GetElementSol(ELT))



/* Access to previous elementsol */
#define Element_GetPreviousElementSol(ELT) \
        ElementSol_GetPreviousElementSol(Element_GetElementSol(ELT))



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



/* Access to nodes and node coordinates */
#define Element_GetNode(ELT,i) \
        (Element_GetPointerToNode(ELT)[i])

#define Element_GetNodeCoordinate(ELT,i) \
        Node_GetCoordinate(Element_GetNode(ELT,i))



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

#define Element_GetCurves(ELT) \
        Material_GetCurves(Element_GetMaterial(ELT))

#define Element_GetCurve(ELT) \
        Material_GetCurve(Element_GetMaterial(ELT))
        
#define Element_GetModel(ELT) \
        Material_GetModel(Element_GetMaterial(ELT))
        
#define Element_GetPropertyValue(ELT,S) \
        (Element_GetProperty(ELT)[Model_GetComputePropertyIndex(Element_GetModel(ELT))(S)])

#define Element_FindCurve(ELT,S) \
        Material_FindCurve(Element_GetMaterial(ELT),S)




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




#include "IntFct.h"
#include "ShapeFct.h"
#include "Material.h"
#include "Node.h"
#include "Buffer.h"
#include "ElementSol.h"


struct Element_s {            /* element */
  unsigned int index ;        /* element index */
  unsigned short int nn ;     /* nb of nodes */
  Node_t** node ;             /* pointers to nodes */
  int    reg ;                /* region index */
  int    imat ;               /* material index */
  Material_t* mat ;           /* material */
  ShapeFct_t* shapefct ;      /* Shape function */
  IntFct_t* fi ;              /* interpolation function */
  unsigned short int dim ;    /* dimension of the element (0,1,2,3) */
  short int* pin ;            /* local position of unknowns at nodes */
  short int* peq ;            /* local position of equations at nodes */
  /* n_vi and n_ve must be kept for old methods! */
  unsigned int n_vi ;         /* Nb of implicit terms */
  unsigned int n_ve ;         /* Nb of explicit terms */
//  unsigned int nbofconstterms ;   /* Nb of constant terms */
  ElementSol_t* sol ;         /* Element Solution */
  Buffer_t*   buffer ;        /* Buffer */
} ;


/* Old notations which I try to eliminate little by little */
#define elem_t         Element_t
#define MAX_NOEUDS     Element_MaxNbOfNodes

#endif
