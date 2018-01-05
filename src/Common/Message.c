#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <time.h>


#include "Message.h"

//static Message_t* instancemsg = NULL ;
static Message_t msg1, *msg = &msg1 ;
static char launchdate[26] ;
static time_t now1, *now = &now1 ;


/* Extern functions */
/*
Message_t*  (Message_GetInstance)(void)
{
  if(!instancemsg) {
    instancemsg = Message_Create() ;
  }
  
  return(instancemsg) ;
}



Message_t*  (Message_Create)(void)
{
  Message_t* message = (Message_t*) malloc(sizeof(Message_t)) ;
  
  if(!message) {
    arret("Message_Create") ;
  }
  
  return(message) ;
}
*/



void (Message_Initialize)(void)
{
  clock_t start = clock() ;
  
  time(now) ;
  
  Message_GetLaunchClock(msg) = start ;
  Message_GetVerbosity(msg)  = 4 ;
  Message_GetLaunchTime(msg) = now ;
  Message_GetLaunchDate(msg) = launchdate ;
  strcpy(launchdate,ctime(now)) ;
}



void (Message_FatalError0)(const char *fmt, ...)
{
  fflush(stdout) ;
  
  fprintf(stderr,"\nBil fatal error...\n") ;
  
  {
    va_list args ;
    va_start(args,fmt) ;
    vfprintf(stderr,fmt,args) ;
    va_end(args) ;
  }
  
  fprintf(stderr,"\n...stop\n") ;
  fflush(stderr) ;
  
  exit(EXIT_SUCCESS) ;
}



void (Message_RuntimeError0)(const char *fmt, ...)
{
  if(Message_GetVerbosity(msg) < 1) {
    exit(EXIT_SUCCESS) ;
    return ;
  }
  
  fflush(stdout) ;
  
  fprintf(stderr,"\nBil runtime error...\n") ;
  
  {
    va_list args ;
    va_start(args,fmt) ;
    vfprintf(stderr,fmt,args) ;
    va_end(args) ;
  }
  
  fprintf(stderr,"\n...stop\n") ;
  fflush(stderr) ;
  
  exit(EXIT_SUCCESS) ;
}



void (Message_Warning)(const char *fmt, ...)
{
  if(Message_GetVerbosity(msg) < 2) return ;
  
  fflush(stdout) ;
  
  fprintf(stderr,"\n") ;
  fprintf(stderr,"Warning: ") ;
  
  {
    va_list args ;
    va_start(args,fmt) ;
    vfprintf(stderr,fmt,args) ;
    va_end(args) ;
  }
  
  fprintf(stderr,"\n") ;
  fflush(stderr) ;
}



void (Message_Info)(const char *fmt, ...)
{
  if(Message_GetVerbosity(msg) < 3) return ;
  
  fflush(stdout) ;
  
  fprintf(stdout,"\n") ;
  fprintf(stdout,"Info: ") ;
  
  {
    va_list args ;
    va_start(args,fmt) ;
    vfprintf(stdout,fmt,args) ;
    va_end(args) ;
  }
  
  fprintf(stdout,"\n") ;
  fflush(stdout) ;
}



int (Message_Direct)(const char *fmt, ...)
{
  int n ;
  if(Message_GetVerbosity(msg) < 4) return(0) ;
  
  fflush(stdout) ;
  
  {
    va_list args ;
    va_start(args,fmt) ;
    n = vfprintf(stdout,fmt,args) ;
    va_end(args) ;
  }
  
  fflush(stdout) ;
  return(n) ;
}


char* (Message_LaunchDate)(void)
{
  return Message_GetLaunchDate(msg) ;
}


double (Message_CPUTime)(void)
{
  double start = (double) Message_GetLaunchClock(msg) ;
  double end   = (double) clock() ;
  double t_cpu = end - start ;
  double elapsed = (t_cpu) / CLOCKS_PER_SEC ;
  
  return(elapsed) ;
}


int (Message_SetVerbosity)(const int newverb)
{
  int oldverb = Message_GetVerbosity(msg) ;
  
  Message_GetVerbosity(msg)  = newverb ;
  
  return(oldverb) ;
}
