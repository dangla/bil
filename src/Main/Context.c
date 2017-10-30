#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Help.h"
#include "BilVersion.h"
#include "BilInfo.h"
#include "Bil.h"
#include "Models.h"
#include "Modules.h"
#include "Message.h"
#include "Context.h"


static Context_t* instancectx = NULL ;

static Context_t* (Context_Create)(void) ;



/* Global functions */
Context_t*  (Context_GetInstance)(void)
{
  if(!instancectx) {
    instancectx = Context_Create() ;
  }
  
  return(instancectx) ;
}



/* Intern functions */
Context_t* (Context_Create)(void)
{
  Context_t* ctx = (Context_t*) calloc(1,sizeof(Context_t)) ;
  
  if(!ctx) {
    arret("Context_Create") ;
  }
  
  /* Initialization */
  /*
  Context_GetHelpOnline(ctx) = NULL ;
  Context_GetPrintInfo(ctx) = NULL ;
  Context_GetInversePermutation(ctx) = NULL ;
  Context_GetPrintModel(ctx) = NULL ;
  Context_GetPrintModule(ctx) = NULL ;
  Context_GetPostProcessing(ctx) = NULL ;
  Context_GetSolver(ctx) = NULL ;
  Context_GetReadOnly(ctx) = NULL ;
  Context_GetDebug(ctx) = NULL ;
  Context_GetPrintLevel(ctx) = NULL ;
  Context_GetUseModule(ctx) = NULL ;
  Context_GetGraph(ctx) = NULL ;
  Context_GetInputFileName(ctx) = NULL ;
  Context_GetMiscellaneous(ctx) = NULL ;
  Context_GetElementOrdering(ctx) = NULL ;
  Context_GetNodalOrdering(ctx) = NULL ;
  */
  
  return(ctx) ;
}
