#ifndef MESSAGE_H
#define MESSAGE_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct Message_s      ; typedef struct Message_s      Message_t ;

/* Fonctions */
extern void   (Message_Delete)(void*) ;
extern void   (Message_RuntimeError0)(const char*,...) ;
extern void   (Message_FatalError0)(const char*,...) ;
extern void   (Message_InputError)(const char*,const int) ;
extern void   (Message_Warning)(const char*,...) ;
extern void   (Message_Info)(const char*,...) ;
extern int    (Message_Direct)(const char*,...) ;
extern void   (Message_Initialize)(Message_t*) ;
extern char*  (Message_LaunchDate)(void) ;
extern double (Message_CPUTime)(void) ;
extern int    (Message_SetVerbosity)(const int) ;


#define Message_GetLaunchClock(MSG)           ((MSG)->launchclock)
#define Message_GetLaunchTime(MSG)            ((MSG)->launchtime)
#define Message_GetLaunchDate(MSG)            ((MSG)->launchdate)
#define Message_GetVerbosity(MSG)             ((MSG)->verbosity)
#define Message_GetDelete(MSG)                ((MSG)->Delete)


#include <stdlib.h>
#include <stdio.h>

#define Message_PrintSourceLocation \
        fprintf(stdout,"\nAt %s, line %d",__FILE__,__LINE__)

#define Message_RuntimeError(...) \
        do { Message_PrintSourceLocation ; \
        Message_RuntimeError0(__VA_ARGS__) ; } while(0)

#define Message_FatalError(...) \
        do { Message_PrintSourceLocation ; \
        Message_FatalError0(__VA_ARGS__) ; } while(0)
                                           
#define Message_Exit \
        (exit(EXIT_SUCCESS))


#include <time.h>
#include "GenericObject.h"

struct Message_s {            /* message */
  clock_t  launchclock ;      /* Start up processor clock time */
  time_t*  launchtime ;       /* Start up time */
  char*    launchdate ;       /* Start up date */
  /* Verbosity level
   *   0:  silent except fatal error, 
   *   1:  +errors , 
   *   2:  +warnings, 
   *   3:  +infos ,
   *   4:  normal, 
   *   99: debug) */
  int  verbosity ;
  GenericObject_Delete_t* Delete ;
} ;




/* Old notations which should be eliminated */
#define arret          Message_FatalError
             

#ifdef __CPLUSPLUS
}
#endif             

#endif
