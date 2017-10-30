#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "Geometry.h"
#include "DataFile.h"
#include "Message.h"
#include "Periodicities.h"


/* Global functions */

Geometry_t*  Geometry_Create(DataFile_t* datafile)
{
  Geometry_t* geom = (Geometry_t*) malloc(sizeof(Geometry_t)) ;
  
  if(!geom) arret("Geometry_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"DIME,GEOM,Geometry",",",1) ;
  
  Message_Direct("Enter in %s","Geometry") ;
  Message_Direct("\n") ;
  
  {
    short int dim = 3 ;
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
  
    if(line) line = strpbrk(line,"0123") ;
  
    if(line) dim  = atoi(line) ;
  
    Geometry_GetDimension(geom) = dim ;
    
    if(line) line += strspn(line,"0123456789 ") ;
  
    /* The symmetry */
    Geometry_SetNoSymmetry(geom) ;

    if(dim < 3) {
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
