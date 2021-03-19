#ifndef LOOPIMPLICITTERMS_H
#define LOOPIMPLICITTERMS_H



#include "Tuple.h"
#include "Arg.h"


#define LoopImplicitTerms(CV,CV_Args,VIM,NVI,...) \
        do { \
          double* vim0  = Element_GetCurrentImplicitTerm(EL) ; \
          IntFct_t* intfct = Element_GetIntFct(EL) ; \
          int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ; \
          int    p ; \
          for(p = 0 ; p < NbOfIntPoints ; p++) { \
            double* x = CV Tuple_RCONS(CV_Args,p) ; \
            double* VIM = vim0 + p*NVI ; \
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
