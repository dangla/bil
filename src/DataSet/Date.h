#ifndef DATE_H
#define DATE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Date_s        ; typedef struct Date_s        Date_t ;


#define Date_GetTime(DATE)            ((DATE)->time)


struct Date_s {
  double time ;
} ;

#endif
