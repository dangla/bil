#ifndef PCM_H
#define PCM_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct PCM_s     ; typedef struct PCM_s     PCM_t ;

/*  Typedef names of Methods */
typedef double*  PCM_ComputeComponents_t(PCM_t*,double,int) ;
typedef void     PCM_ComputeSecondaryComponents_t(PCM_t*,double,double*) ;
typedef double*  PCM_ComputeFluxes_t(PCM_t*,double*,int,int) ;



#include "Elements.h"
#include "IntFcts.h"

extern PCM_t*     PCM_GetInstance(Element_t*,double**,double**,double*) ;
extern double*    PCM_ComputeComponentDerivatives(PCM_t*,PCM_ComputeSecondaryComponents_t*,double,double,int,int) ;
extern double*    PCM_ComputeComponentFluxes(PCM_t*,PCM_ComputeFluxes_t*,int,int) ;


#define PCM_MaxNbOfComponentVectors              (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))
#define PCM_MaxNbOfComponents                    (100)
#define PCM_MaxNbOfComponentFluxes               (10)

#define PCM_MaxSizeOfComponents                  (PCM_MaxNbOfComponents*sizeof(double))

#define PCM_MaxSizeOfOutput                      (PCM_MaxNbOfComponentVectors*PCM_MaxNbOfComponents*sizeof(double))
#define PCM_MaxSizeOfInput                       (0)
#define PCM_SizeOfBuffer                         (0)

#define PCM_GetElement(pcm)                          ((pcm)->el)
#define PCM_GetPointerToNodalUnknowns(pcm)           ((pcm)->u)
#define PCM_GetPointerToPreviousNodalUnknowns(pcm)   ((pcm)->u_n)
#define PCM_GetPreviousImplicitTerms(pcm)            ((pcm)->vim)
#define PCM_GetInput(pcm)                            ((pcm)->input)
#define PCM_GetOutput(pcm)                           ((pcm)->output)
#define PCM_GetBuffer(pcm)                           ((pcm)->buffer)
#define PCM_GetPointerToComponents(pcm)              ((pcm)->components)
#define PCM_GetPointerToComponentDerivatives(pcm)    ((pcm)->compoderiv)
#define PCM_GetPointerToComponentFluxes(pcm)         ((pcm)->compfluxes)


#define PCM_GetComponents(pcm,i)              (PCM_GetPointerToComponents(pcm)[i])
#define PCM_GetComponentDerivatives(pcm,i)    (PCM_GetPointerToComponentDerivatives(pcm)[i])


#define PCM_AllocateInBuffer(pcm,sz)         (Buffer_Allocate(PCM_GetBuffer(pcm),(sz)))
#define PCM_FreeBuffer(pcm)                  (Buffer_Free(PCM_GetBuffer(pcm)))
#define PCM_FreeBufferFrom(pcm,p)            (Buffer_FreeFrom(PCM_GetBuffer(pcm),(char*) (p)))



#include "Buffer.h"

struct PCM_s {                /* PhysicoChemistry Method */
  Element_t* el ;             /* Element */
  double**   u ;              /* Nodal unknowns */
  double**   u_n ;            /* Previous Nodal unknowns */
  double*    vim ;            /* Previous implicit terms */
  double**   components ;     /* Pointer to components */
  double**   compoderiv ;     /* Pointer to component derivatives */
  double**   compfluxes ;     /* Pointer to component fluxes */
  void*      input ;          /* Input */
  void*      output ;         /* Output*/
  Buffer_t*  buffer ;         /* Buffer */
} ;

#endif
