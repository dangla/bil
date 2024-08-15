#ifndef FEM2_H
#define FEM2_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct FEM2_s     ; typedef struct FEM2_s     FEM2_t ;


#include "Mesh.h"
#include "Solvers.h"
#include "DataSet.h"
#include "Solutions.h"

extern FEM2_t*  (FEM2_GetInstance)(DataSet_t*,Solvers_t*,Solutions_t*,Solutions_t*) ;
extern void     (FEM2_Delete)(void*) ;
extern int      (FEM2_HomogenizeTangentStiffnessTensor)(FEM2_t*,double,double,double*) ;
extern int      (FEM2_ComputeHomogenizedStressTensor)(FEM2_t*,double,double,double*,double*) ;
extern void     (FEM2_InitializeMicrostructureDataSet)(FEM2_t*) ;
//extern int      (FEM2_HomogenizeTangentStiffnessTensor1)(Mesh_t*,Solver_t*,double,double,double*) ;
//extern int      (FEM2_ComputeHomogenizedStressTensor1)(DataSet_t*,Solver_t*,double,double,Solutions_t*,Solutions_t*,double*,double*) ;



#define FEM2_SizeOfBuffer                        (1)


#define FEM2_GetDataSet(F2)                      ((F2)->dataset)
#define FEM2_GetSolvers(F2)                      ((F2)->solvers)
#define FEM2_GetCurrentSolutions(F2)             ((F2)->sols)
#define FEM2_GetPreviousSolutions(F2)            ((F2)->sols_n)
#define FEM2_GetBuffers(F2)                      ((F2)->buffers)


#define FEM2_GetSolver(F2) \
        Solvers_GetSolver(FEM2_GetSolvers(F2))



/* Buffer */
#define FEM2_GetBuffer(F2) \
        Buffers_GetBufferOfCurrentThread(FEM2_GetBuffers(F2))
        
        

#define FEM2_AllocateInBuffer(F2,sz) \
        (Buffer_Allocate(FEM2_GetBuffer(F2),(sz)))
        
#define FEM2_FreeBuffer(F2) \
        (Buffer_Free(FEM2_GetBuffer(F2)))
        
#define FEM2_FreeBufferFrom(F2,p) \
        (Buffer_FreeFrom(FEM2_GetBuffer(F2),(char*) (p)))



#include "Buffers.h"
#include "GenericObject.h"

struct FEM2_s {               /* (FEM)^2 */
  DataSet_t* dataset ;
  Solvers_t* solvers ;
  Solutions_t* sols ;
  Solutions_t* sols_n ;
  Buffers_t*  buffers ;         /* Buffer */
  GenericObject_Delete_t* Delete ;
} ;

#ifdef __CPLUSPLUS
}
#endif

#endif
