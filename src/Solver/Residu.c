#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "Residu.h"




/* Extern functions */

Residu_t*   (Residu_Create)(const int n_col,const int n)
{
  Residu_t* residu = (Residu_t*) Mry_New(Residu_t) ;

  Residu_GetLengthOfRHS(residu) = n_col ;
  Residu_GetNbOfRHS(residu) = n ;


  /* Allocation of space for the right hand sides */
  {
    double* rhs = (double*) Mry_New(double[n*n_col]) ;
    
    Residu_GetRHS(residu) = rhs ;
  }
  
  
  /* Allocation of space for the solutions */
  {
    double* sol = (double*) Mry_New(double[n*n_col]) ;
    
    Residu_GetSolution(residu) = sol ;
  }
  
  
  return(residu) ;
}



void Residu_Delete(void* self)
{
  Residu_t** prs = (Residu_t**) self ;
  Residu_t* rs   = *prs ;
  
  free(Residu_GetRHS(rs)) ;
  free(Residu_GetSolution(rs)) ;
  
  free(rs) ;
}





void (Residu_AssembleElementResidu)(Residu_t* residu,Element_t* el,double* re)
{
  int  nn  = Element_GetNbOfNodes(el) ;
  int  neq = Element_GetNbOfEquations(el) ;
  double* r = Residu_GetRHS(residu) ;
  int imatrix = Residu_GetResiduIndex(residu) ;
            
    if(Element_GetMaterial(el)) {
      int i ;
            
      for(i = 0 ; i < nn ; i++) {
        Node_t* node_i = Element_GetNode(el,i) ;
        int j ;
              
        for(j = 0 ; j < neq ; j++) {
          int ij = i*neq + j ;
          int ii = Element_GetUnknownPosition(el)[ij] ;
                
          if(ii >= 0) {
            int k = Node_GetSelectedMatrixColumnIndexOf(node_i,ii,imatrix) ;
            if(k >= 0) r[k] += re[ij] ;
          }
        }
      }
    }
}
