#include "Message.h"
#include "Exception.h"
#include "Session.h"
#include "GenericData.h"
#include "Mry.h"
#include <signal.h>
#include <assert.h>


//static Exception_t* instanceexception = NULL ;

static Exception_t*  (Exception_Create)(void) ;
static void          (Exception_Initialize)(void) ;

static void SignalHandler(int) ;



/* Global functions */
#if 0
//Exception_t*  (Exception_GetInstance0)(void)
{
  if(!instanceexception) {
    instanceexception = Exception_Create() ;
  }
  
  return(instanceexception) ;
}
#endif



Exception_t*  (Exception_GetInstance)(void)
{
  GenericData_t* gdat = Session_FindGenericData(Exception_t,"Exception") ;
  
  if(!gdat) {
    Exception_t* exc = Exception_Create() ;
    
    gdat = GenericData_Create(1,exc,Exception_t,"Exception") ;
    
    Session_AddGenericData(gdat) ;
    
    assert(gdat == Session_FindGenericData(Exception_t,"Exception")) ;
  }
  
  return((Exception_t*) GenericData_GetData(gdat)) ;
}



/* Intern functions */
Exception_t*  (Exception_Create)(void)
{
  Exception_t* exception = (Exception_t*) Mry_New(Exception_t) ;
  
  Exception_GetDelete(exception) = Exception_Delete ;
  
  Exception_Initialize() ;
  
  return(exception) ;
}



void  (Exception_Delete)(void* self)
{
  Exception_t** pexception = (Exception_t**) self ;
  
  free(*pexception) ;
  *pexception = NULL ;
}



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
  
  if(signal(SIGABRT,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  if(signal(SIGILL,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
  
  if(signal(SIGTERM,SignalHandler) == SIG_ERR) {
    Message_FatalError("An error occured while setting a signal handler.\n") ;
  }
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
      Message_Direct("- Mesh with physical index out of range\n") ;
      Exception_BackupAndTerminate ;
      return ;
    
    case SIGFPE:
  
      Message_Direct("Floating-point error.\n") ;
      Exception_BackupAndTerminate ;
      return ;
      
    case SIGINT:
  
      Message_Direct("Program interrupted.\n") ;
      Exception_BackupAndTerminate ;
      return ;
      
    case SIGABRT:
  
      Message_Direct("Abnormal termination.\n") ;
      Exception_BackupAndTerminate ;
      return ;
      
    case SIGILL:
  
      Message_Direct("Illegal operation.\n") ;
      Exception_BackupAndTerminate ;
      return ;
      
    case SIGTERM:
  
      Message_Direct("Termination request.\n") ;
      Exception_BackupAndTerminate ;
      return ;
      
    default:
    
      break ;
  }
  
  Message_FatalError("SignalHandler: unexpected error") ;
}




#if 0
static void InterruptHandler(int) ;
static void FloatingPointErrorHandler(int) ;
static void SegmentationFaultHandler(int) ;


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
#endif
