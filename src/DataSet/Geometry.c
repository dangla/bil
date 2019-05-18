#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include "Geometry.h"
#include "DataFile.h"
#include "Message.h"
#include "Periodicities.h"
#include "String.h"


/* Global functions */

Geometry_t*  Geometry_Create(DataFile_t* datafile)
{
  Geometry_t* geom = Geometry_New() ; ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"DIME,GEOM,Geometry",",",1) ;
  
  Message_Direct("Enter in %s","Geometry") ;
  Message_Direct("\n") ;
  
  {
    short int dim = 3 ;
    char* line = DataFile_GetCurrentPositionInFileContent(datafile) ;
    //char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
  
    if(line) line = String_FindAnyChar(line,"0123") ;
  
    if(line) {
      //DataFile_ScanAdv(datafile,"%d",&dim) ;
      dim  = atoi(line) ;
    }
  
    Geometry_GetDimension(geom) = dim ;
    
    if(line) line += strspn(line,"0123456789 ") ;
  
    /* The symmetry */

    if(dim > 0 && dim < 3) {
      char* pline = line ;
    
      if(!strncasecmp(pline,"plan",4))  {
        Geometry_SetPlaneSymmetry(geom) ;
      } else if(!strncasecmp(pline,"axis",4))  {
        Geometry_SetCylindricalSymmetry(geom) ;
      } else if(!strncasecmp(pline,"sphe",4))  {
        Geometry_SetSphericalSymmetry(geom) ;
      } else {
        arret("Geometry_Create: geometry not available\n\
        Available geometries are:\n\
        PLAN, AXIS, SPHE") ;
      }
    }
  }

  DataFile_CloseFile(datafile) ;
  
  
  Geometry_GetPeriodicities(geom) = Periodicities_Create(datafile) ;
  
  return(geom) ;
}



Geometry_t*  Geometry_New(void)
{
  Geometry_t* geom = (Geometry_t*) malloc(sizeof(Geometry_t)) ;
  
  assert(geom) ;
  
  Geometry_GetDimension(geom) = 3 ;
  Geometry_SetNoSymmetry(geom) ;
  Geometry_GetPeriodicities(geom) = NULL ;
  
  return(geom) ;
}




void Geometry_Delete(void* self)
{
  Geometry_t** pgeom = (Geometry_t**) self ;
  Geometry_t*   geom = *pgeom ;
  
  Periodicities_Delete(&Geometry_GetPeriodicities(geom)) ;
  
  free(geom) ;
  
  *pgeom = NULL ;
}
