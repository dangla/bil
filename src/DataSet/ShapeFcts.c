#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Mry.h"
#include "Message.h"
#include "Math_.h"
#include "ShapeFcts.h"




/* Extern functions */

ShapeFcts_t*  ShapeFcts_Create(void)
{
  ShapeFcts_t* shapefcts = (ShapeFcts_t*) Mry_New(ShapeFcts_t) ;


  /* Memory space for shape functions */
  {
    int n = ShapeFcts_MaxNbOfShapeFcts ;
    ShapeFct_t* shapefct = (ShapeFct_t*) Mry_New(ShapeFct_t[n]) ;
    
    ShapeFcts_GetShapeFct(shapefcts) = shapefct ;
  }

  /* Nb of shape functions set to 0 */
  ShapeFcts_GetNbOfShapeFcts(shapefcts) = 0 ;
  
  return(shapefcts) ;
}



void  (ShapeFcts_Delete)(void* self)
{
  ShapeFcts_t* shapefcts = (ShapeFcts_t*) self ;
  
  {
    int n = ShapeFcts_GetNbOfShapeFcts(shapefcts) ;
    ShapeFct_t* shapefct = ShapeFcts_GetShapeFct(shapefcts) ;
    
    if(shapefct) {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        ShapeFct_t* shapefct_i = shapefct + i ;
      
        ShapeFct_Delete(shapefct_i) ;
      }
    
      free(shapefct) ;
    }
  }
}



int ShapeFcts_FindShapeFct(ShapeFcts_t* shapefcts,int nn,int dim)
/** Find a Shape Function class defined by 
 *  nn = Nb of nodes 
 *  dim = local dimension
 *  Return the index of the shape function. */
{
  int n_fi = ShapeFcts_GetNbOfShapeFcts(shapefcts) ;
  int    i ;

  /* Does the function already exist? */
  for(i = 0 ; i < n_fi ; i++) {
    ShapeFct_t* fi = ShapeFcts_GetShapeFct(shapefcts) + i ;
    int    i_nn = ShapeFct_GetNbOfNodes(fi) ;
    int    i_dim = ShapeFct_GetDimension(fi) ;

    if((nn == i_nn) && (dim == i_dim)) {
      return(i) ; /* Existing function */
    }
  }

  i = ShapeFcts_AddShapeFct(shapefcts,nn,dim) ;

  return(i) ;
}


int ShapeFcts_AddShapeFct(ShapeFcts_t* shapefcts,int nn,int dim)
/** Add a Shape Function class defined by 
 *  nn = Nb of nodes 
 *  dim = local dimension
 *  Return the index of the shape function. */
{
  ShapeFct_t* shapefct = ShapeFcts_GetShapeFct(shapefcts) ;
  int    i ;

  /* We add a new function */
  ShapeFcts_GetNbOfShapeFcts(shapefcts) += 1 ;
  
  if(ShapeFcts_GetNbOfShapeFcts(shapefcts) > ShapeFcts_MaxNbOfShapeFcts) {
    arret("ShapeFcts_AddShapeFct(1): too many functions") ;
  }

  /* Index of the new function */
  i = ShapeFcts_GetNbOfShapeFcts(shapefcts) - 1 ;
  
  {
    ShapeFct_t* shapefcti = ShapeFct_Create(nn,dim) ;
    
    shapefct[i] = shapefcti[0] ;
    free(shapefcti) ;
  }
  
  return(i) ;
}
