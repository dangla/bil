#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Dates.h"



Dates_t*  Dates_Create(DataFile_t *datafile)
{
  FILE *ficd ;
  Dates_t *dates = (Dates_t*) malloc(sizeof(Dates_t)) ;
  int n_dates ;
  int   i ;
  
  if(!dates) arret("Dates_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"TEMP,DATE,Dates",",",1) ;
  
  Message_Direct("Enter in %s","Dates") ;
  Message_Direct("\n") ;
  
  ficd = DataFile_GetFileStream(datafile) ;
  
  fscanf(ficd,"%d",&n_dates) ;
  Dates_GetNbOfDates(dates) = n_dates ;

  {
    double *date = (double *) malloc(n_dates*sizeof(double)) ;
    if(!date) arret("Dates_Create (1)") ;
    Dates_GetDate(dates) = date ;
  }

  for(i = 0 ; i < n_dates ; i++) {
    fscanf(ficd,"%lf",Dates_GetDate(dates) + i) ;
  }
  
  return(dates) ;
}
