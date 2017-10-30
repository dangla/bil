#include <signal.h>
#include "Message.h"
#include "Exception.h"


static Exception_t* instanceexception = NULL ;

static Exception_t*  (Exception_Create)(void) ;
static void          (Exception_Initialize)(void) ;

static void InterruptHandler(int) ;
static void FloatingPointErrorHandler(int) ;
static void SegmentationFaultHandler(int) ;
static void SignalHandler(int) ;



/* Global functions */
Exception_t*  (Exception_GetInstance)(void)
{
  if(!instanceexception) {
    instanceexception = Exception_Create() ;
  }
  
  return(instanceexception) ;
}



/* Intern functions */
void (Exception_Initialize)(void)
{
  if(signal(SIGINT,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  if(signal(SIGFPE,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  if(signal(SIGSEGV,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
}



Exception_t*  (Exception_Create)(void)
{
  Exception_t* exception = (Exception_t*) malloc(sizeof(Exception_t)) ;
  
  if(!exception) {
    arret("Exception_Create") ;
  }
  
  Exception_Initialize() ;
  
  return(exception) ;
}



void (InterruptHandler)(int val)
{
  /* NOTE some versions of UNIX will reset signal to default
		 after each call. So for portability reset signal each time */
  if(signal(SIGINT,InterruptHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  Message_Warning("Interactive attention signal caught.\n") ;
  Exception_RestoreEnvironment(val) ;
  Message_FatalError("InterruptHandler") ;
}



void (FloatingPointErrorHandler)(int val)
{
  /* NOTE some versions of UNIX will reset signal to default
		 after each call. So for portability reset signal each time */
  if(signal(SIGFPE,FloatingPointErrorHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  Message_Warning("Floating-point error.\n") ;
  Message_FatalError("FloatingPointErrorHandler") ;
}



void (SegmentationFaultHandler)(int val)
{
  if(signal(SIGSEGV,SegmentationFaultHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  Message_Warning("Segmentation fault.\n") ;
  Message_FatalError("SegmentationFaultHandler") ;
}



void (SignalHandler)(int sigid)
{
  if(signal(sigid,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  switch(sigid) {
    case SIGSEGV:
  
      Message_Direct("Segmentation fault. ") ;
      Message_Direct("Possible sources:\n") ;
      Message_Direct("- Bad \"iperm\" file (remove it!)\n") ;
      break ;
    
    case SIGFPE:
  
      Message_Direct("Floating-point error.\n") ;
      break ;
      
    case SIGINT:
  
      Message_Direct("Program interrupted.\n") ;
      Exception_RestoreEnvironment(sigid) ;
      return ;
      
    default:
    
      break ;
  }
  
  Message_FatalError("SignalHandler: unexpected error") ;
}
