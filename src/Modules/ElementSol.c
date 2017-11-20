#include <stdio.h>
#include "ElementSol.h"
#include "Message.h"
#include "GenericData.h"



/* Extern functions */

ElementSol_t* (ElementSol_GetDeepElementSol)(ElementSol_t* elementsol,unsigned int depth)
{
  while(depth--) elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
  return(elementsol) ;
}
