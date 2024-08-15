#include <stdlib.h>
#include "Message.h"
#include "DataFile.h"
#include "String_.h"
#include "Mry.h"
#include "Regions.h"



Regions_t*  (Regions_New)(void)
{
  Regions_t* regions = (Regions_t*) Mry_New(Regions_t) ;
  
  Regions_GetNbOfRegions(regions) = 0 ;
  
  
  /* Allocation of space for each region */
  {
    int n = Regions_MaxNbOfRegions ;
    Region_t* region = Mry_Create(Region_t,n,Region_New()) ;
    
    Regions_GetRegion(regions) = region ;
  }
  
  return(regions) ;
}



void (Regions_Delete)(void* self)
{
  Regions_t* regions = (Regions_t*) self ;

  {
    int    n = Regions_MaxNbOfRegions ;
    Region_t* region = Regions_GetRegion(regions) ;
      
    Mry_Delete(region,n,Region_Delete) ;
    free(region) ;
  }
}





int (Regions_FindRegionIndex)(Regions_t* regions,const char* name)
{
  Region_t* region = Regions_GetRegion(regions) ;
  int n_regions = Regions_GetNbOfRegions(regions) ;
  int j = 0 ;
  
  while(j < n_regions && String_IsNot(Region_GetRegionName(region + j),name)) j++ ;
  
  if(j < n_regions) {
    return(j) ;
  } else if(j < Regions_MaxNbOfRegions) {
    j = n_regions ;
    Regions_GetNbOfRegions(regions) = n_regions + 1 ;
    strncpy(Region_GetRegionName(region + j),name,Region_MaxLengthOfRegionName) ;
    
    return(j) ;
  }

  return(-1) ;
}




Region_t* (Regions_FindRegion)(Regions_t* regions,const char* name)
{
  int j = Regions_FindRegionIndex(regions,name) ;
  
  if(j >= 0) {
    Region_t* region = Regions_GetRegion(regions) ;
    
    return(region + j) ;
  }

  return(NULL) ;
}
