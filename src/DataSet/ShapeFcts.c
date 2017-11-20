#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Message.h"
#include "Tools/Math.h"
#include "ShapeFcts.h"

static void   ShapeFct_Create(ShapeFct_t*,int,int) ;


/* Extern functions */

ShapeFcts_t*  ShapeFcts_Create(void)
{
  ShapeFcts_t* shapefcts = (ShapeFcts_t*) malloc(sizeof(ShapeFcts_t)) ;
  
  if(!shapefcts) arret("ShapeFcts_Create") ;


  /* Memory space for shape functions */
  {
    size_t sz = ShapeFcts_MaxNbOfShapeFcts*sizeof(ShapeFct_t) ;
    ShapeFct_t* shapefct = (ShapeFct_t*) malloc(sz) ;
    
    if(!shapefct) arret("ShapeFcts_Create(1)") ;
    
    ShapeFcts_GetShapeFct(shapefcts) = shapefct ;
  }

  /* Nb of shape functions set to 0 */
  ShapeFcts_GetNbOfShapeFcts(shapefcts) = 0 ;
  
  return(shapefcts) ;
}



int ShapeFcts_FindShapeFct(ShapeFcts_t *shapefcts,int nn,int dim)
/** Find a Shape Function class defined by 
 *  nn = Nb of nodes 
 *  dim = local dimension
 *  Return the index of the shape function. */
{
  int n_fi = ShapeFcts_GetNbOfShapeFcts(shapefcts) ;
  int    i ;

  /* Does the function already exist? */
  for(i = 0 ; i < n_fi ; i++) {
    ShapeFct_t *fi = ShapeFcts_GetShapeFct(shapefcts) + i ;
    int    i_nn = ShapeFct_GetNbOfFunctions(fi) ;
    int    i_dim = ShapeFct_GetDimension(fi) ;

    if((nn == i_nn) && (dim == i_dim)) {
      return(i) ; /* Existing function */
    }
  }

  i = ShapeFcts_AddShapeFct(shapefcts,nn,dim) ;

  return(i) ;
}


int ShapeFcts_AddShapeFct(ShapeFcts_t *shapefcts,int nn,int dim)
/** Add a Shape Function class defined by 
 *  nn = Nb of nodes 
 *  dim = local dimension
 *  Return the index of the shape function. */
{
  ShapeFct_t *shapefct = ShapeFcts_GetShapeFct(shapefcts) ;
  int    i ;

  /* We add a new function */
  ShapeFcts_GetNbOfShapeFcts(shapefcts) += 1 ;
  
  if(ShapeFcts_GetNbOfShapeFcts(shapefcts) > ShapeFcts_MaxNbOfShapeFcts) {
    arret("ShapeFcts_AddShapeFct(1): too many functions") ;
  }

  /* Index of the new function */
  i = ShapeFcts_GetNbOfShapeFcts(shapefcts) - 1 ;
  
  ShapeFct_Create(shapefct + i,nn,dim) ;
  
  return(i) ;
}


/* Intern functions */

void ShapeFct_Create(ShapeFct_t* shapefct,int nn,int dim)
{
  ShapeFct_GetDimension(shapefct) = dim ;
  ShapeFct_GetNbOfFunctions(shapefct) = nn ;
  
  {
    int k = dim + nn*(1 + dim) ;
    double* b = (double*) malloc(k*sizeof(double)) ;
    
    if(!b) arret("ShapeFct_Create(1)") ;
    
    ShapeFct_GetCoordinate(shapefct) = b ;
    ShapeFct_GetFunction(shapefct) = ShapeFct_GetCoordinate(shapefct) + dim ;
    ShapeFct_GetFunctionGradient(shapefct) = ShapeFct_GetFunction(shapefct) + nn ;
  }
  
  return ;
}
