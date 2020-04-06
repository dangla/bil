#include <stdlib.h>
#include "Message.h"
#include "Mry.h"
#include "Model.h"



#if 0
Node_t*  Node_Create(const int dim)
{
  Node_t* node = (Node_t*) Mry_New(Node_t) ;
  
  {
    double* x = (double*) Mry_New(double[dim]) ;
    
    Node_GetCoordinate(node) = x ;
  }
  
  return(node) ;
}
#endif
