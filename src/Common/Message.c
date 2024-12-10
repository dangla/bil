#include "Message.h"
#include "Session.h"
#include "GenericData.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <time.h>
#include <Mry.h>
#include <assert.h>
#include "DistributedMS.h"


static Message_t*  (Message_GetInstance)(void) ;
static Message_t*  (Message_Create)(void) ;




/* Extern functions */
Message_t*  (Message_GetInstance)(void)
{
  GenericData_t* gdat = Session_FindGenericData(Message_t,"Message") ;
  
  if(!gdat) {
    Message_t* msg = Message_Create() ;
    
    gdat = GenericData_Create(1,msg,Message_t,"Message") ;
    
    Session_AddGenericData(gdat) ;
    
    assert(gdat == Session_FindGenericData(Message_t,"Message")) ;
  }
  
  return((Message_t*) GenericData_GetData(gdat)) ;
}




Message_t*  (Message_Create)(void)
{
  Message_t* msg = (Message_t*) Mry_New(Message_t) ;
  
  {
    char* date = (char*) Mry_New(char[26]) ;
    
    Message_GetLaunchDate(msg) = date ;
  }
  
  {
    time_t* now = (time_t*) Mry_New(time_t) ;
    
    Message_GetLaunchTime(msg) = now ;
  }
  
  Message_GetDelete(msg) = &Message_Delete ;
  
  Message_Initialize(msg) ;
  
  return(msg) ;
}



void  (Message_Delete)(void* self)
{
  Message_t* msg = (Message_t*) self ;
  
  free(Message_GetLaunchDate(msg)) ;
  free(Message_GetLaunchTime(msg)) ;
  Message_GetDelete(msg) = NULL ;
}




void (Message_Initialize)(Message_t* msg)
{
  {
    clock_t start = clock() ;
  
    Message_GetLaunchClock(msg) = start ;
  }
  
  {
    time_t* now = Message_GetLaunchTime(msg) ;
  
    time(now) ;
  }
  
  {
    char* date = Message_GetLaunchDate(msg) ;
    time_t* now = Message_GetLaunchTime(msg) ;
    
    strcpy(date,ctime(now)) ;
  }
  
  Message_GetVerbosity(msg)  = 4 ;
}



void (Message_FatalError0)(const char* fmt, ...)
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



void (Message_RuntimeError0)(const char* fmt, ...)
{
  Message_t* msg = Message_GetInstance() ;
  
  if(!msg || Message_GetVerbosity(msg) < 1) {
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



void (Message_InputError)(const char* name,const int i)
{
  Message_t* msg = Message_GetInstance() ;
  
  if(!msg || Message_GetVerbosity(msg) < 1) {
    exit(EXIT_SUCCESS) ;
    return ;
  }
  
  fflush(stdout) ;
  
  fprintf(stderr,"\nBil input error...\n") ;
  
  {
    fprintf(stderr,"** On entry to %s, parameter number %d had an illegal value",name,i);
  }
  
  fprintf(stderr,"\n...stop\n") ;
  fflush(stderr) ;
  
  exit(EXIT_SUCCESS) ;
}



void (Message_Warning)(const char* fmt, ...)
{
  Message_t* msg = Message_GetInstance() ;
  
  if(DistributedMS_RankOfCallingProcess) return ;
  
  if(!msg || Message_GetVerbosity(msg) < 2) return ;
  
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



void (Message_Info)(const char* fmt, ...)
{
  Message_t* msg = Message_GetInstance() ;
  
  if(DistributedMS_RankOfCallingProcess) return ;
  
  if(!msg || Message_GetVerbosity(msg) < 3) return ;
  
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



int (Message_Direct)(const char* fmt, ...)
{
  Message_t* msg = Message_GetInstance() ;
  int n ;
  
  if(DistributedMS_RankOfCallingProcess) return(0) ;
  
  if(!msg || Message_GetVerbosity(msg) < 4) return(0) ;
  
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
  Message_t* msg = Message_GetInstance() ;
  
  return(Message_GetLaunchDate(msg)) ;
}


double (Message_CPUTime)(void)
{
  Message_t* msg = Message_GetInstance() ;
  double start = (double) Message_GetLaunchClock(msg) ;
  double end   = (double) clock() ;
  double t_cpu = end - start ;
  double elapsed = (t_cpu) / CLOCKS_PER_SEC ;
  
  return(elapsed) ;
}


int (Message_SetVerbosity)(const int newverb)
{
  Message_t* msg = Message_GetInstance() ;
  int oldverb = Message_GetVerbosity(msg) ;
  
  Message_GetVerbosity(msg)  = newverb ;
  
  return(oldverb) ;
}
