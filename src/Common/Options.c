#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Mry.h"
#include "Context.h"
#include "Message.h"
#include "Options.h"
#include "Modules.h"
#include "Threads.h"


/* Global functions */
static void   Options_SetDefault(Options_t*) ;
static void   Options_Initialize(Options_t*) ;



Options_t*  (Options_Create)(Context_t* ctx)
{
  Options_t* options = (Options_t*) Mry_New(Options_t) ;

  {
    int   max_mot_debug = Options_MaxLengthOfKeyWord ;
    char* c = (char*) Mry_New(char[5*max_mot_debug]) ;

    Options_GetPrintedInfos(options)      = c ;
    Options_GetResolutionMethod(options)  = (c += max_mot_debug) ;
    Options_GetPrintLevel(options)        = (c += max_mot_debug) ;
    Options_GetModule(options)            = (c += max_mot_debug) ;
    Options_GetNbOfThreads(options)       = (c += max_mot_debug) ;
  }
  
  
  Options_SetDefault(options) ;
  
  if(ctx) {
    Options_GetContext(options) = ctx ;
    Options_Initialize(options) ;
  }
  
  return(options) ;
}



void (Options_Delete)(void* self)
{
  Options_t* options = (Options_t*) self ;
  
  {
    char* c = Options_GetPrintedInfos(options) ;
    
    if(c) {
      free(c) ;
      Options_GetPrintedInfos(options) = NULL ;
    }
  }
}



/* Local functions */

void (Options_SetDefault)(Options_t* options)
/* Set default options */
{
  const char*  modulenames[] = {Modules_ListOfNames} ;
  const char*  defaultmodule = modulenames[0] ;
  
  strcpy(Options_GetPrintedInfos(options),"\0") ;
  strcpy(Options_GetResolutionMethod(options),"crout") ;
  strcpy(Options_GetPrintLevel(options),"1") ;
  strcpy(Options_GetModule(options),defaultmodule) ;
  Options_GetContext(options) = NULL ;
  strcpy(Options_GetNbOfThreads(options),"1") ;
}



void Options_Initialize(Options_t* options)
/* Set options from the command line arguments */
{
  Context_t* ctx = Options_GetContext(options) ;
  
  if(Context_GetSolver(ctx)) {
    char* mth = ((char**) Context_GetSolver(ctx))[1] ;
    
    //Options_GetResolutionMethod(options) = ((char**) Context_GetSolver(ctx))[1] ;
    strcpy(Options_GetResolutionMethod(options),mth) ;
  }
  
  if(Context_GetDebug(ctx)) {
    char* dbg = ((char**) Context_GetDebug(ctx))[1] ;
    
    //Options_GetPrintedInfos(options) = ((char**) Context_GetDebug(ctx))[1] ;
    strcpy(Options_GetPrintedInfos(options),dbg) ;
  }
  
  if(Context_GetGraph(ctx)) {
    Options_GetGraphMethod(options) = ((char**) Context_GetGraph(ctx))[1] ;
  }
  
  if(Context_GetElementOrdering(ctx)) {
    Options_GetElementOrderingMethod(options) = ((char**) Context_GetElementOrdering(ctx))[1] ;
  }
  
  if(Context_GetNodalOrdering(ctx)) {
    Options_GetNodalOrderingMethod(options) = ((char**) Context_GetNodalOrdering(ctx))[1] ;
  }
  
  if(Context_GetPrintLevel(ctx)) {
    char* level = ((char**) Context_GetPrintLevel(ctx))[1] ;
    
    //Options_GetPrintLevel(options) = ((char**) Context_GetPrintLevel(ctx))[1] ;
    strcpy(Options_GetPrintLevel(options),level) ;
  }
  
  if(Context_GetNbOfThreads(ctx)) {
    char* nthreads = ((char**) Context_GetNbOfThreads(ctx))[1] ;
    
    strcpy(Options_GetNbOfThreads(options),nthreads) ;
  }
  
  if(Context_GetUseModule(ctx)) {
    char* module = ((char**) Context_GetUseModule(ctx))[1] ;
    
    //Options_GetModule(options) = ((char**) Context_GetUseModule(ctx))[1] ;
    strcpy(Options_GetModule(options),module) ;
  }
  
  if(Context_GetPostProcessing(ctx)) {
    Options_GetPostProcessingMethod(options) = ((char**) Context_GetPostProcessing(ctx))[1] ;
  }

}
