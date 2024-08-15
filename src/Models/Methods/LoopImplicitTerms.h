#ifndef LOOPIMPLICITTERMS_H
#define LOOPIMPLICITTERMS_H


#include "Element.h"
#include "IntFct.h"

#include "Tuple.h"
#include "Arg.h"

#define LoopOnGaussIntFct(EL,P) \
        for(int P = 0 ; P < IntFct_GetNbOfPoints(Element_GetIntFct(EL)) ; P++)
        


#define LoopImplicitTerms(EL,CV,CV_Args,VIM,NVI,...) \
        do { \
          LoopOnGaussIntFct(EL,p) {\
            double* VIM = Element_GetCurrentImplicitTerm(EL) + p*NVI ; \
            double* x = CV Tuple_RCONS(CV_Args,p) ; \
            FEM_Store(x,__VA_ARGS__) ; \
          } \
        } while(0)



/* Implementation */
#define FEM_Store0(X,A,N) \
        do { \
          { \
            int i; \
            for(i = 0 ; i < N ; i++) A[i] = X[Utils_CAT(I_,A) + i] ; \
          } \
        } while(0)


#define FEM_Store(X,TA,TN) \
        Algos_SEPWITH(Algos_MAP3(Tuple_TUPLE(Utils_DUP(Tuple_LEN(TA),X,)),TA,TN,FEM_Store0),;)


#endif
