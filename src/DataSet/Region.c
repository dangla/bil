#include <stdlib.h>
#include <ctype.h>
#include "String_.h"
#include "Math_.h"
#include "Message.h"
#include "Mry.h"
#include "Region.h"


Region_t*  (Region_New)(void)
{
  Region_t* region = (Region_t*) Mry_New(Region_t) ;
  
  /* Allocation of space for the region name */
  {
    char* name = (char*) Mry_New(char[Region_MaxLengthOfRegionName]) ;
    
    Region_GetRegionName(region) = name ;
  }
  
  return(region) ;
}



void (Region_Delete)(void* self)
{
  Region_t* region = (Region_t*) self ;

  {
    char* name = Region_GetRegionName(region) ;
      
    if(name) {
      free(name) ;
    }
  }
}
