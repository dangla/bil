#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Mry.h"
#include "Message.h"
#include "Periodicity.h"



Periodicity_t* (Periodicity_New)(void)
{
  Periodicity_t* periodicity = (Periodicity_t*) Mry_New(Periodicity_t) ;

  {
    double* vector = (double*) Mry_New(double[3]) ;
    
    Periodicity_GetPeriodVector(periodicity) = vector ;
  }
  
  
  /* Allocation of space for the region name */
  {
    char* name = (char*) Mry_New(char[2*Periodicity_MaxLengthOfRegionName]) ;
    
    Periodicity_GetMasterRegionName(periodicity) = name ;
    Periodicity_GetSlaveRegionName(periodicity) = name + Periodicity_MaxLengthOfRegionName ;
  }
  
  return(periodicity) ;
}



void (Periodicity_Delete)(void* self)
{
  Periodicity_t* periodicity = (Periodicity_t*) self ;
  
  {
    double* vector = Periodicity_GetPeriodVector(periodicity) ;
    
    if(vector) {
      free(vector) ;
    }
  }
  
  {
    char* name = Periodicity_GetMasterRegionName(periodicity) ;
    
    if(name) {
      free(name) ;
    }
  }
}



void (Periodicity_Scan)(Periodicity_t* periodicity,DataFile_t* datafile)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
  
    /* Master region */
    {
      char name[Periodicity_MaxLengthOfRegionName] ;
      int n = String_FindAndScanExp(line,"MasterRegion",","," = %s",name) ;
        
      if(n) {
        Periodicity_GetMasterRegion(periodicity) = atoi(name) ;
        strncpy(Periodicity_GetMasterRegionName(periodicity),name,Periodicity_MaxLengthOfRegionName)  ;
      } else {
        arret("Periodicity_Scan: no master region") ;
      }
    }

    
    /* Slave region */
    {
      char name[Periodicity_MaxLengthOfRegionName] ;
      int n = String_FindAndScanExp(line,"SlaveRegion",","," = %s",name) ;
        
      if(n) {
        Periodicity_GetSlaveRegion(periodicity) = atoi(name) ;
        strncpy(Periodicity_GetSlaveRegionName(periodicity),name,Periodicity_MaxLengthOfRegionName)  ;
      } else {
        arret("Periodicity_Scan: no slave region") ;
      }
    }
      
      
    /* Period vector */
    {
      double x[3] ;
      int n = String_FindAndScanExp(line,"PeriodVector",","," = %lf %lf %lf",x,x+1,x+2) ;
        
      if(n) {
        double* vector = Periodicity_GetPeriodVector(periodicity) ;
        int    j ;
        
        for(j = 0 ; j < 3 ; j++) {
          vector[j] = x[j] ;
        }
      } else {
        arret("Periodicity_Scan: no period vector") ;
      }
    }
}
