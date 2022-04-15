#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Mry.h"
#include "Message.h"
#include "Math_.h"
#include "IntFcts.h"



/* Extern functions */

IntFcts_t*  (IntFcts_Create)(void)
{
  IntFcts_t* intfcts = (IntFcts_t*) Mry_New(IntFcts_t) ;


  /* Memory space for interpolation functions */
  {
    IntFct_t* intfct = (IntFct_t*) Mry_New(IntFct_t[IntFcts_MaxNbOfIntFcts]) ;
    
    IntFcts_GetIntFct(intfcts) = intfct ;
  }

  /* Nb of interpolation functions set to 0 */
  IntFcts_GetNbOfIntFcts(intfcts) = 0 ;
  
  return(intfcts) ;
}



void (IntFcts_Delete)(void* self)
{
  IntFcts_t* intfcts = (IntFcts_t*) self ;
  
  {
    int n = IntFcts_GetNbOfIntFcts(intfcts) ;
    IntFct_t* intfct = IntFcts_GetIntFct(intfcts) ;
    
    if(intfct) {
      int i ;
    
      for(i = 0 ; i < n ; i++) {
        IntFct_t* intfct_i = intfct + i ;
      
        IntFct_Delete(intfct_i) ;
      }
    
      free(intfct) ;
    }
  }
}



int IntFcts_FindIntFct(IntFcts_t* intfcts,int nn,int dim,const char* type)
/** Find an Interpolation Function class defined by 
 *  nn = nb of functions 
 *  dim = local dimension
 *  type = characteristic of this interpolation functions
 *  Return the index of the Int. Fct. */
{
  int n_fi = IntFcts_GetNbOfIntFcts(intfcts) ;
  int    i ;

  /* Does the function already exist? */
  for(i = 0 ; i < n_fi ; i++) {
    IntFct_t* fi  = IntFcts_GetIntFct(intfcts) + i ;
    int    i_nn   = IntFct_GetNbOfFunctions(fi) ;
    int    i_dim  = IntFct_GetDimension(fi) ;
    char*  i_type = IntFct_GetType(fi) ;

    if((nn == i_nn) && (dim == i_dim) && (strcmp(type,i_type) == 0)) {
      return(i) ;
    }
  }

  i = IntFcts_AddIntFct(intfcts,nn,dim,type) ;

  return(i) ;
}


int IntFcts_AddIntFct(IntFcts_t* intfcts,int nn,int dim,const char* type)
/** Add an Interpolation Function class defined by 
 *  nn = nb of functions 
 *  dim = local dimension
 *  type = characteristic of this interpolation functions
 *  Return the index of the Int. Fct. */
{
  IntFct_t* intfct = IntFcts_GetIntFct(intfcts) ;
  int    i ;

  /* We add a new function */
  IntFcts_GetNbOfIntFcts(intfcts) += 1 ;
  
  if(IntFcts_GetNbOfIntFcts(intfcts) > IntFcts_MaxNbOfIntFcts) {
    arret("IntFcts_AddIntFct(1): too many functions") ;
  }

  /* Index of the new function */
  i = IntFcts_GetNbOfIntFcts(intfcts) - 1 ;
  
  {
    IntFct_t* intfcti = IntFct_Create(nn,dim,type) ;
    
    intfct[i] = intfcti[0] ;
    free(intfcti) ;
  }
  
  return(i) ;
}
