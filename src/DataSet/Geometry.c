#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include "Geometry.h"
#include "DataFile.h"
#include "Message.h"
#include "Periodicities.h"
#include "Mry.h"
#include "String_.h"


/* Global functions */

Geometry_t*  (Geometry_New)(void)
{
  Geometry_t* geom = (Geometry_t*) Mry_New(Geometry_t) ;
  
  Geometry_GetDimension(geom) = 3 ;
  Geometry_SetNoSymmetry(geom) ;
  Geometry_GetPeriodicities(geom) = NULL ;
  
  return(geom) ;
}




void (Geometry_Delete)(void* self)
{
  Geometry_t* geom = (Geometry_t*) self ;
  
  {
    Periodicities_t* periodicities = Geometry_GetPeriodicities(geom) ;
    
    if(periodicities) {
      Periodicities_Delete(periodicities) ;
      free(periodicities) ;
    }
  }
}



Geometry_t*  (Geometry_Create)(DataFile_t* datafile)
{
  Geometry_t* geom = Geometry_New() ; ;
  
  
  Message_Direct("Enter in %s","Geometry") ;
  Message_Direct("\n") ;
  
  {
    short int dim = 3 ;
    char* filecontent = DataFile_GetFileContent(datafile) ;
    char* c  = String_FindToken(filecontent,"DIME,GEOM,Geometry",",") ;
    char* line = String_SkipLine(c) ;
  
    if(line) line = String_FindAnyChar(line,"0123") ;
  
    if(line) {
      dim  = atoi(line) ;
    }
  
    Geometry_GetDimension(geom) = dim ;
    
    if(line) line = String_SkipAnyChars(line,"0123456789 ") ;
  
    /* The symmetry */

    if(dim > 0 && dim < 3) {
      char* pline = line ;
    
      if(String_CaseIgnoredIs(pline,"plane",4))  {
        
        Geometry_SetPlaneSymmetry(geom) ;
        
      } else if(String_CaseIgnoredIs(pline,"axis",4))  {
        
        Geometry_SetCylindricalSymmetry(geom) ;
        
      } else if(String_CaseIgnoredIs(pline,"sphe",4))  {
        
        Geometry_SetSphericalSymmetry(geom) ;
        
      } else {
        
        arret("Geometry_Create: geometry not available\n\
        Available geometries are:\n\
        PLANE, AXIS, SPHE") ;
        
      }
    }
  }
  
  Geometry_GetPeriodicities(geom) = Periodicities_Create(datafile) ;
  
  return(geom) ;
}
