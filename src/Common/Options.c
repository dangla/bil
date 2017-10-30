#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Context.h"
#include "Message.h"
#include "Options.h"
#include "ListOfModules.h"


/* Global functions */
static void   Options_SetDefault(Options_t*) ;
static void   Options_SetFromContext(Options_t*) ;


Options_t*  Options_Create(void)
{
  Options_t *options = (Options_t*) malloc(sizeof(Options_t)) ;
  
  if(!options) arret("Options_Create(1)") ;

  {
    int   max_mot_debug = Options_MaxLengthOfKeyWord ;
    char* c = (char *) malloc(4*max_mot_debug*sizeof(char)) ;
  
    if(!c) arret("Options_Create(2)") ;

    Options_GetPrintData(options)         = c ;
    Options_GetResolutionMethod(options)  = (c += max_mot_debug) ;
    Options_GetPrintLevel(options)        = (c += max_mot_debug) ;
    Options_GetModule(options)            = (c += max_mot_debug) ;
  }
  
  
  Options_SetDefault(options) ;
  
  
  Options_SetFromContext(options) ;
  
  return(options) ;
}


/* Local functions */

void Options_SetDefault(Options_t *options)
/* Set default options */
{
  char  *modulenames[NB_MODULES] = {MODULENAMES} ;
  char  *defaultmodule = modulenames[0] ;
  
  strcpy(Options_GetPrintData(options),"\0") ;
  strcpy(Options_GetResolutionMethod(options),"crout") ;
  strcpy(Options_GetPrintLevel(options),"1") ;
  strcpy(Options_GetModule(options),defaultmodule) ;
  Options_GetInputFileName(options) = NULL ;
}


void Options_SetFromContext(Options_t *options)
/* Set options from the command line arguments */
{
  Context_t *ctx = Context_GetInstance() ;
  
  if(Context_GetSolver(ctx)) {
    Options_GetResolutionMethod(options) = ((char**) Context_GetSolver(ctx))[1] ;
  }
  
  if(Context_GetDebug(ctx)) {
    Options_GetPrintData(options) = ((char**) Context_GetDebug(ctx))[1] ;
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
  
  if(Context_GetInputFileName(ctx)) {
    Options_GetInputFileName(options) = ((char**) Context_GetInputFileName(ctx))[0] ;
  }
  
  if(Context_GetPostProcessing(ctx)) {
    Options_GetPostProcessingMethod(options) = ((char**) Context_GetPostProcessing(ctx))[1] ;
  }
  
  if(!Options_GetInputFileName(options)) {
    Message_FatalError("Missing file name") ;
  }
  
}
