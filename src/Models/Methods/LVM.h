#ifndef LVM_H
#define LVM_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct LVM_s     ; typedef struct LVM_s     LVM_t ;


#include "Elements.h"
#include "IntFcts.h"


/*  Typedef names of Methods */
//typedef double*  LVM_ComputeVariables_t(LVM_t*,double,int) ;
typedef void     LVM_ComputeSecondaryVariables_t(Element_t*,double,double*) ;
typedef double*  LVM_ComputeFluxes_t(Element_t*,double*,double*,int,int) ;


extern LVM_t*     LVM_Create(void) ;
extern double*    LVM_ComputeVariableDerivatives(LVM_t*,LVM_ComputeSecondaryVariables_t*,double,double,int,int) ;
extern double*    LVM_ComputeVariableFluxes(LVM_t*,LVM_ComputeFluxes_t*,int,int) ;


#define LVM_MaxNbOfVariableVectors              (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))
#define LVM_MaxNbOfVariables                    (100)
#define LVM_MaxNbOfVariableFluxes               (10)

#define LVM_MaxSizeOfVariables                  (LVM_MaxNbOfVariables*sizeof(double))

#define LVM_MaxSizeOfOutput                      (LVM_MaxNbOfVariableVectors*LVM_MaxNbOfVariables*sizeof(double))
#define LVM_MaxSizeOfInput                       (0)
#define LVM_SizeOfBuffer                         (0)

#define LVM_GetElement(lvm)                          ((lvm)->el)
#define LVM_GetPointerToNodalUnknowns(lvm)           ((lvm)->u)
#define LVM_GetPointerToPreviousNodalUnknowns(lvm)   ((lvm)->u_n)
#define LVM_GetPreviousImplicitTerms(lvm)            ((lvm)->vim_n)
#define LVM_GetInput(lvm)                            ((lvm)->input)
#define LVM_GetOutput(lvm)                           ((lvm)->output)
#define LVM_GetBuffer(lvm)                           ((lvm)->buffer)
#define LVM_GetPointerToVariables(lvm)               ((lvm)->variables)
#define LVM_GetPointerToVariableDerivatives(lvm)     ((lvm)->varderiv)
#define LVM_GetPointerToVariableFluxes(lvm)          ((lvm)->varfluxes)


#define LVM_GetVariables(lvm,i)              (LVM_GetPointerToVariables(lvm)[i])
#define LVM_GetVariableDerivatives(lvm,i)    (LVM_GetPointerToVariableDerivatives(lvm)[i])
#define LVM_GetVariableFluxes(lvm,i)         (LVM_GetPointerToVariableFluxes(lvm)[i])


#define LVM_AllocateInBuffer(lvm,sz)         (Buffer_Allocate(LVM_GetBuffer(lvm),(sz)))
#define LVM_FreeBuffer(lvm)                  (Buffer_Free(LVM_GetBuffer(lvm)))
#define LVM_FreeBufferFrom(lvm,p)            (Buffer_FreeFrom(LVM_GetBuffer(lvm),(char*) (p)))



#include "Buffer.h"

struct LVM_s {                /* Local Secondary Variables */
  Element_t* el ;             /* Element */
  double**   u ;              /* Nodal unknowns */
  double**   u_n ;            /* Previous Nodal unknowns */
  double*    vim_n ;          /* Previous implicit terms */
  double**   variables ;      /* Pointer to variables */
  double**   varderiv ;       /* Pointer to variables derivatives */
  double**   varfluxes ;      /* Pointer to variables fluxes */
  void*      input ;          /* Input */
  void*      output ;         /* Output*/
  Buffer_t*  buffer ;         /* Buffer */
} ;

#endif
