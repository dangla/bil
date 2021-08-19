#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Mry.h"
#include "Context.h"
#include "Message.h"
#include "Options.h"
#include "Modules.h"


/* Global functions */
static void   Options_SetDefault(Options_t*) ;
static void   Options_Initialize(Options_t*) ;



Options_t*  (Options_Create)(Context_t* ctx)
{
  Options_t* options = (Options_t*) Mry_New(Options_t) ;

  {
    int   max_mot_debug = Options_MaxLengthOfKeyWord ;
    char* c = (char *) Mry_New(char[4*max_mot_debug]) ;

    Options_GetPrintedInfos(options)      = c ;
    Options_GetResolutionMethod(options)  = (c += max_mot_debug) ;
    Options_GetPrintLevel(options)        = (c += max_mot_debug) ;
    Options_GetModule(options)            = (c += max_mot_debug) ;
  }
  
  
  Options_SetDefault(options) ;
  
  if(ctx) {
    Options_GetContext(options) = ctx ;
    Options_Initialize(options) ;
  }
  
  return(options) ;
}



void Options_Delete(void* self)
{
  Options_t** poptions = (Options_t**) self ;
  Options_t*  options = *poptions ;
  
  Context_Delete(&(Options_GetContext(options))) ;
  free(Options_GetPrintedInfos(options)) ;
  free(options) ;
}



/* Local functions */

void Options_SetDefault(Options_t* options)
/* Set default options */
{
  const char*  modulenames[] = {Modules_ListOfNames} ;
  const char*  defaultmodule = modulenames[0] ;
  
  strcpy(Options_GetPrintedInfos(options),"\0") ;
  strcpy(Options_GetResolutionMethod(options),"crout") ;
  strcpy(Options_GetPrintLevel(options),"1") ;
  strcpy(Options_GetModule(options),defaultmodule) ;
  Options_GetContext(options) = NULL ;
}



void Options_Initialize(Options_t* options)
/* Set options from the command line arguments */
{
  Context_t* ctx = Options_GetContext(options) ;
  
  if(Context_GetSolver(ctx)) {
    Options_GetResolutionMethod(options) = ((char**) Context_GetSolver(ctx))[1] ;
  }
  
  if(Context_GetDebug(ctx)) {
    Options_GetPrintedInfos(options) = ((char**) Context_GetDebug(ctx))[1] ;
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
    Options_GetPrintLevel(options) = ((char**) Context_GetPrintLevel(ctx))[1] ;
  }
  
  if(Context_GetUseModule(ctx)) {
    Options_GetModule(options) = ((char**) Context_GetUseModule(ctx))[1] ;
  }
  
  if(Context_GetPostProcessing(ctx)) {
    Options_GetPostProcessingMethod(options) = ((char**) Context_GetPostProcessing(ctx))[1] ;
  }

}
