#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Message.h"
#include "Tools/Math.h"
#include "IntFcts.h"



/* Extern functions */

IntFcts_t*  IntFcts_Create(void)
{
  IntFcts_t* intfcts = (IntFcts_t*) malloc(sizeof(IntFcts_t)) ;
  
  if(!intfcts) arret("IntFcts_Create") ;


  /* Memory space for interpolation functions */
  {
    size_t sz = IntFcts_MaxNbOfIntFcts*sizeof(IntFct_t) ;
    IntFct_t* intfct = (IntFct_t*) malloc(sz) ;
    
    if(!intfct) arret("IntFcts_Create(1)") ;
    
    IntFcts_GetIntFct(intfcts) = intfct ;
  }

  /* Nb of interpolation functions set to 0 */
  IntFcts_GetNbOfIntFcts(intfcts) = 0 ;
  
  return(intfcts) ;
}



int IntFcts_FindIntFct(IntFcts_t* intfcts,int nn,int dim,const char* type)
/** Find an Interpolation Function class defined by 
 *  nn = Nb of nodes 
 *  dim = local dimension
 *  type = characteristic of this interpolation functions
 *  Return the index of the Int. Fct. */
{
  int n_fi = IntFcts_GetNbOfIntFcts(intfcts) ;
  int    i ;

  /* Does the function already exist? */
  for(i = 0 ; i < n_fi ; i++) {
    IntFct_t* fi = IntFcts_GetIntFct(intfcts) + i ;
    int    i_nn = IntFct_GetNbOfNodes(fi) ;
    int    i_dim = IntFct_GetDimension(fi) ;
    char   *i_type = IntFct_GetType(fi) ;

    if((nn == i_nn) && (dim == i_dim) && (strcmp(type,i_type) == 0)) {
      return(i) ;
    }
  }

  i = IntFcts_AddIntFct(intfcts,nn,dim,type) ;

  return(i) ;
}


int IntFcts_AddIntFct(IntFcts_t* intfcts,int nn,int dim,const char* type)
/** Add an Interpolation Function class defined by 
 *  nn = Nb of nodes 
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
  
  IntFct_Create(intfct + i,nn,dim,type) ;
  
  return(i) ;
}
