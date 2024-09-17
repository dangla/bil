#ifndef FVM_H
#define FVM_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct FVM_s     ; typedef struct FVM_s     FVM_t ;


#include "Element.h"
#include "Load.h"
#include "IntFct.h"

/*  Typedef names of Methods */
typedef void     FVM_ComputeFluxes_t(FVM_t*,double*,double*,int,int) ;


extern FVM_t*     (FVM_GetInstance)(Element_t*) ;
extern void       (FVM_Delete)(void*) ;

extern double*    (FVM_ComputeMassMatrix)(FVM_t*,double*,int) ;
extern double*    (FVM_ComputeIsotropicConductionMatrix)(FVM_t*,double*,int) ;
extern double*    (FVM_ComputeMassAndIsotropicConductionMatrix)(FVM_t*,double*,int) ;

extern double*    (FVM_ComputeSurfaceLoadResidu)(FVM_t*,Load_t*,double,double) ;
extern double*    (FVM_ComputeBodyForceResidu)(FVM_t*,double*) ;
extern double*    (FVM_ComputeFluxResidu)(FVM_t*,double*) ;
extern double*    (FVM_ComputeMassAndFluxResidu)(FVM_t*,double*) ;
extern double*    (FVM_ComputeMassBalanceEquationResidu)(FVM_t*,double const*,double const*,double const) ;
//extern double*    (FVM_ComputeMassBalanceEquationResidu)(FVM_t*,double const*,double const*,double const,int const) ;


extern double*    (FVM_ComputeCellVolumes)(FVM_t*) ;
extern double*    (FVM_ComputeCellSurfaceAreas)(FVM_t*) ;
extern double*    (FVM_ComputeCellVolumesAndSurfaceAreas)(FVM_t*) ;
extern short int  (FVM_FindLocalCellIndex)(FVM_t*,double*) ;
extern double*    (FVM_ComputeIntercellDistances)(FVM_t*) ;


extern double*    (FVM_ComputeVariableFluxes)(FVM_t*,FVM_ComputeFluxes_t*,int,int) ;

extern double*    (FVM_ComputeTheNodalFluxVector)(FVM_t*,double*) ;

extern double     (FVM_AverageCurrentImplicitTerm)(Mesh_t*,const char*,const int,const int) ;
extern double     (FVM_AveragePreviousImplicitTerm)(Mesh_t*,const char*,const int,const int) ;

extern double*    (FVM_ComputeGradient)(FVM_t*,double*,IntFct_t*,int,int);




#define FVM_MaxNbOfNodes                         (Element_MaxNbOfNodes)
#define FVM_MaxNbOfDOF                           (Element_MaxNbOfDOF)
#define FVM_MaxShift                             (Model_MaxNbOfEquations*Model_MaxNbOfEquations)
#define FVM_MaxNbOfMatrices                      (4)

#define FVM_MaxSizeOfMatrix                      (FVM_MaxNbOfDOF*FVM_MaxNbOfDOF*sizeof(double))
#define FVM_MaxSizeOfOutput                      (FVM_MaxNbOfDOF*FVM_MaxNbOfDOF*sizeof(double))
#define FVM_MaxSizeOfInput                       (FVM_MaxNbOfNodes*FVM_MaxNbOfNodes*FVM_MaxShift*sizeof(double))
#define FVM_SizeOfBuffer                         (FVM_MaxNbOfMatrices*FVM_MaxSizeOfMatrix)





#define FVM_GetElement(fvm)                      ((fvm)->el)
#define FVM_GetInput(fvm)                        ((fvm)->input)
#define FVM_GetOutput(fvm)                       ((fvm)->output)
#define FVM_GetShiftOfInput(fvm)                 ((fvm)->shift)
#define FVM_GetBuffers(fvm)                      ((fvm)->buffers)
#define FVM_GetCellVolumes(fvm)                  ((fvm)->cellvolumes)
#define FVM_GetCellSurfaceAreas(fvm)             ((fvm)->cellsurfaceareas)
#define FVM_GetIntercellDistances(fvm)           ((fvm)->celldistances)




/* Buffer */
#define FVM_GetBuffer(fvm) \
        Buffers_GetBufferOfCurrentThread(FVM_GetBuffers(fvm))



#define FVM_AllocateInBuffer(fvm,sz) \
        (Buffer_Allocate(FVM_GetBuffer(fvm),(sz)))
        
#define FVM_FreeBuffer(fvm) \
        (Buffer_Free(FVM_GetBuffer(fvm)))
        
#define FVM_FreeBufferFrom(fvm,p) \
        (Buffer_FreeFrom(FVM_GetBuffer(fvm),(char*) (p)))



#include "Buffers.h"
#include "GenericObject.h"

struct FVM_s {                /* Finite Volume Method */
  Element_t* el ;             /* Element */
  void*      input ;          /* Input */
  void*      output ;         /* Output*/
  int        shift ;          /* Shift of input */
  Buffers_t* buffers ;        /* Buffer */
  double*    cellvolumes ;
  double*    cellsurfaceareas ;
  double*    celldistances ;
  GenericObject_Delete_t* Delete ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
