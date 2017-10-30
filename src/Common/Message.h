#ifndef MESSAGE_H
#define MESSAGE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Message_s      ; typedef struct Message_s      Message_t ;

/* Fonctions */
extern void   (Message_RuntimeError0)(const char*,...) ;
extern void   (Message_FatalError0)(const char*,...) ;
extern void   (Message_Warning)(const char*,...) ;
extern void   (Message_Info)(const char*,...) ;
extern int    (Message_Direct)(const char*,...) ;
extern void   (Message_Initialize)(void) ;
extern char*  (Message_LaunchDate)(void) ;
extern double (Message_CPUTime)(void) ;
extern int    (Message_SetVerbosity)(const int) ;


#define Message_GetLaunchClock(message)           ((message)->launchclock)
#define Message_GetLaunchTime(message)            ((message)->launchtime)
#define Message_GetLaunchDate(message)            ((message)->launchdate)
#define Message_GetVerbosity(message)             ((message)->verbosity)


#include <stdlib.h>
//#include <stdio.h>

#define Message_PrintSourceLocation      Message_Direct("\nAt %s, line %d",__FILE__,__LINE__)

#define Message_RuntimeError(...) {   \
        Message_PrintSourceLocation ; \
        Message_RuntimeError0(__VA_ARGS__) ; }

#define Message_FatalError(...) {     \
        Message_PrintSourceLocation ; \
        Message_FatalError0(__VA_ARGS__) ; }
                                           
#define Message_Exit        exit(EXIT_SUCCESS)


#include <time.h>

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
} ;




/* Old notations which should be eliminated */
#define arret          Message_RuntimeError
                          

#endif
