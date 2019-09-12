#ifndef DATES_H
#define DATES_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Dates_s        ; typedef struct Dates_s        Dates_t ;

#include "DataFile.h"


extern Dates_t*  (Dates_New)    (const int) ;
extern Dates_t*  (Dates_Create) (DataFile_t*) ;



#define Dates_GetNbOfDates(DATES)       ((DATES)->n_dates)
#define Dates_GetDate(DATES)            ((DATES)->date)



#include "Date.h"

struct Dates_s {
  unsigned int n_dates ;      /* nb of dates */
  Date_t* date ;              /* date */
} ;

#endif
