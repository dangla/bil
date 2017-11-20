#ifndef FEM_H
#define FEM_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct FEM_s     ; typedef struct FEM_s     FEM_t ;


#include "Mesh.h"
#include "Element.h"
#include "IntFct.h"
#include "Load.h"

/* Extern Functions */

extern FEM_t*   FEM_GetInstance(Element_t*) ;

/* Matrices */
extern double*  FEM_ComputeElasticMatrix(FEM_t*,IntFct_t*,double*,int) ;
extern double*  FEM_ComputeMassMatrix(FEM_t*,IntFct_t*,double*,int) ;
extern double*  FEM_ComputeConductionMatrix(FEM_t*,IntFct_t*,double*,int) ;
extern double*  FEM_ComputeBiotMatrix(FEM_t*,IntFct_t*,double*,int) ;
extern double*  FEM_ComputePoroelasticMatrix(FEM_t*,IntFct_t*,double*,int,int) ;
extern void     FEM_TransformMatrixFromDegree2IntoDegree1(FEM_t*,int,int,double*) ;

/* Residus */
extern double*  FEM_ComputeSurfaceLoadResidu(FEM_t*,IntFct_t*,Load_t*,double,double) ;
extern double*  FEM_ComputeBodyForceResidu(FEM_t*,IntFct_t*,double*,int) ;
extern double*  FEM_ComputeStrainWorkResidu(FEM_t*,IntFct_t*,double*,int) ;
extern double*  FEM_ComputeFluxResidu(FEM_t*,IntFct_t*,double*,int) ;
extern void     FEM_TransformResiduFromDegree2IntoDegree1(FEM_t*,int,double*) ;

/* Compute Shape functions */
extern double*  FEM_ComputeIsoShapeFctInActualSpace(FEM_t*,double*) ;


/* Compute unknowns and unknown gradients: method 1 */
extern double   FEM_ComputeCurrentUnknown(FEM_t*,double*,int,int) ;
extern double   FEM_ComputePreviousUnknown(FEM_t*,double*,int,int) ;

extern double*  FEM_ComputeCurrentUnknownGradient(FEM_t*,double*,int,int) ;
extern double*  FEM_ComputePreviousUnknownGradient(FEM_t*,double*,int,int) ;

extern double*  FEM_ComputeCurrentLinearStrainTensor(FEM_t*,double*,double*,int,int) ;
extern double*  FEM_ComputePreviousLinearStrainTensor(FEM_t*,double*,double*,int,int) ;
extern double*  FEM_ComputeIncrementalLinearStrainTensor(FEM_t*,double*,double*,int,int) ;


/* Compute unknowns and unknown gradients: method 2 */
extern double*  FEM_ComputeUnknowns(FEM_t*,IntFct_t*,double**,int) ;
extern double*  FEM_ComputeUnknownGradients(FEM_t*,IntFct_t*,double**,int) ;
extern double*  FEM_ComputeLinearStrainTensors(FEM_t*,IntFct_t*,double**,int) ;

extern double   FEM_ComputeUnknown(FEM_t*,double**,IntFct_t*,int,int) ;
extern double*  FEM_ComputeUnknownGradient(FEM_t*,double**,IntFct_t*,int,int) ;
extern double*  FEM_ComputeLinearStrainTensor(FEM_t*,double**,IntFct_t*,int,int) ;
extern double*  FEM_ComputeDisplacementVector(FEM_t*,double**,IntFct_t*,int,int) ;


/* Compute element integration */
extern double   FEM_IntegrateOverElement(FEM_t*,IntFct_t*,double*,int) ;

/* Averaging */
extern void     FEM_AverageStresses(Mesh_t*,double*) ;


#include "Models.h"

#define FEM_MaxShift                             (81 + 18*Model_MaxNbOfEquations + 9*Model_MaxNbOfEquations*Model_MaxNbOfEquations)
#define FEM_MaxNbOfMatrices                      (4)
#define FEM_MaxSizeOfMatrix                      (Element_MaxNbOfDOF*Element_MaxNbOfDOF*sizeof(double))
#define FEM_SizeOfBuffer                         (FEM_MaxNbOfMatrices*FEM_MaxSizeOfMatrix)


#define FEM_GetElement(fem)                      ((fem)->el)
#define FEM_GetPointerToIntFct(fem)              ((fem)->pintfct)
#define FEM_GetInput(fem)                        ((fem)->input)
#define FEM_GetOutput(fem)                       ((fem)->output)
#define FEM_GetShiftOfInput(fem)                 ((fem)->shift)
#define FEM_GetBuffer(fem)                       ((fem)->buffer)


#define FEM_AllocateInBuffer(fem,sz)         (Buffer_Allocate(FEM_GetBuffer(fem),(sz)))
#define FEM_FreeBuffer(fem)                  (Buffer_Free(FEM_GetBuffer(fem)))
#define FEM_FreeBufferFrom(fem,p)            (Buffer_FreeFrom(FEM_GetBuffer(fem),(char*) (p)))


#include "Buffer.h"

struct FEM_s {                /* Finite Element Method */
  Element_t* el ;             /* Element */
  IntFct_t** pintfct ;        /* Pointer to interpolation functions */
  void*      input ;
  void*      output ;
  int        shift ;
  Buffer_t*  buffer ;         /* Buffer */
} ;

#endif
