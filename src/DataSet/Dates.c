#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Dates.h"



static Date_t*  Date_Create(int) ;



Dates_t*  Dates_Create(DataFile_t* datafile)
{
  FILE* ficd ;
  Dates_t* dates = (Dates_t*) malloc(sizeof(Dates_t)) ;
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
    Date_t* date = Date_Create(n_dates) ;
    
    Dates_GetDate(dates) = date ;
  }

  for(i = 0 ; i < n_dates ; i++) {
    Date_t* date = Dates_GetDate(dates) + i ;
    double t ;
    
    fscanf(ficd,"%lf",&t) ;
    
    Date_GetTime(date) = t ;
  }
  
  return(dates) ;
}




Date_t*  Date_Create(int n_dates)
{
  Date_t* date = (Date_t*) malloc(n_dates*sizeof(Date_t)) ;
  
  if(!date) arret("Date_Create") ;
  
  return(date) ;
}
