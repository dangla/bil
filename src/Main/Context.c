#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "String_.h"
#include "Context.h"
#include "Mry.h"


static void          (Context_Initialize)(Context_t*) ;


/* Global functions */
Context_t* (Context_Create)(int argc,char** argv)
{
  Context_t* ctx = (Context_t*) Mry_New(Context_t) ;
  
  
  {
    CommandLine_t* cmd = CommandLine_Create(argc,argv) ;
    
    Context_GetCommandLine(ctx) = cmd ;
  }
  
  /* Initialization */
  if(argc > 0 && argv) {
    Context_Initialize(ctx) ;
  }
  
  
  {
    Options_t* opt = Options_Create(ctx) ;
    
    Context_GetOptions(ctx) = opt ;
  }
  
  return(ctx) ;
}



void (Context_Delete)(void* self)
{
  Context_t* ctx = (Context_t*) self ;
  
  {
    CommandLine_t* cmd = Context_GetCommandLine(ctx) ;
    
    if(cmd) {
      CommandLine_Delete(cmd) ;
      free(cmd) ;
      Context_GetCommandLine(ctx) = NULL ;
    }
  }
  
  {
    Options_t* opt = Context_GetOptions(ctx) ;
    
    if(opt) {
      Options_Delete(opt) ;
      free(opt) ;
      Context_GetOptions(ctx) = NULL ;
    }
  }
}



void (Context_Initialize)(Context_t* ctx)
{
  CommandLine_t* cmd = Context_GetCommandLine(ctx) ;
  int    argc = CommandLine_GetNbOfArg(cmd) ;
  char** argv = CommandLine_GetArg(cmd) ;
  int i ;
  

  if(argc == 1) {
    Context_GetPrintUsage(ctx) = (char**) argv ;
    return ;
  }
  
  
  /* Get line options */
  for(i = 1 ; i < argc ; i++) {
    
    if(argv[i][0] != '-') { /* File name */
      Context_GetInputFileName(ctx) = (char**) argv + i ;
      
    } else if(String_Is(argv[i],"-info",5)) {
      Context_GetPrintInfo(ctx) = (char**) argv + i ;
  
    } else if(String_Is(argv[i],"-help",5)) {
      Context_GetHelpOnline(ctx) = (char**) argv + i ;
      
    } else if(String_Is(argv[i],"-solver",strlen(argv[i]))) {
      Context_GetSolver(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing solver") ;
      }
      
      /* 
       * Skip two more entries if the following entry is either:
       * "-ff <factor>": the fill factor for multi-frontal methods
       **/
      {
        if(String_Is(argv[i + 1],"-ff")) {
          if(i + 2 < argc) {
            i += 2 ;
          } else {
            Message_FatalError("Missing fill factor") ;
          }
        }
      }
      
      /* 
       * Skip some more entries for PetscKSP methods
       **/
      if(String_Is(argv[i],"petscksp",strlen(argv[i]))) {
      
        /* 
         * Skip two more entries if the following entry is either:
         * "-ksp_XXX" <YYY>: option for PetscKSP methods
         **/
        while(String_Is(argv[i + 1],"-ksp_")) {
          if(i + 2 < argc) {
            i += 2 ;
          } else {
            Message_FatalError("Missing value for ksp option") ;
          }
        }
      
      
        /* 
         * Skip two more entries if the following entry is either:
         * "-pc_XXX" <YYY>: option for PetscKSP methods
         **/
      
        while(String_Is(argv[i + 1],"-pc_")) {
          if(i + 2 < argc) {
            i += 2 ;
          } else {
            Message_FatalError("Missing value for pc option") ;
          }
        }
      
      
        /* 
         * Skip two more entries if the following entry is either:
         * "-mat_XXX" <YYY>: option for PetscKSP methods
         **/
      
        while(String_Is(argv[i + 1],"-mat_")) {
          if(i + 2 < argc) {
            i += 2 ;
          } else {
            Message_FatalError("Missing value for mat option") ;
          }
        }
      
      
        /* 
         * Skip two more entries if the following entry is either:
         * "-log_XXX" <YYY>: option for PetscKSP methods
         **/
      
        while(String_Is(argv[i + 1],"-log_")) {
          if(i + 2 < argc) {
            i += 2 ;
          } else {
            Message_FatalError("Missing value for log option") ;
          }
        }
      
      
        /* 
         * Skip one more entry if the following entry is:
         * "-mpi_linear_solver_server": option for PetscKSP methods
         **/
      
        while(String_Is(argv[i + 1],"-mpi_linear_solver_server",5)) {
          i += 1 ;
        }
      }
    
    } else if(String_Is(argv[i],"-debug",strlen(argv[i]))) {
      Context_GetDebug(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing name of data to be printed") ;
      }

    } else if(String_Is(argv[i],"-level",strlen(argv[i]))) {
      Context_GetPrintLevel(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing level") ;
      }

    } else if(String_Is(argv[i],"-nthreads",strlen(argv[i]))) {
      Context_GetNbOfThreads(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing nb of requested threads") ;
      }

    } else if(String_Is(argv[i],"-with",strlen(argv[i]))) {
      Context_GetUseModule(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing module") ;
      }
      
      /* If the module is "Monolithic" we skip the next entry
       * which should be the nb of solutions requested. */
      {
        if(String_Is(argv[i],"Monolithic")) {
          if(i + 1 < argc) {
            i += 1 ;
          } else {
            Message_FatalError("Missing nb of solutions") ;
          }
        }
      }
      
      /* If the module is "SNIA" we skip the next entry
       * which should be the nb of sequences requested. */
      {
        if(String_Is(argv[i],"SNIA")) {
          if(i + 1 < argc) {
            i += 1 ;
          } else {
            Message_FatalError("Missing nb of sequences") ;
          }
        }
      }

    } else if(String_Is(argv[i],"-models",strlen(argv[i]))) {
      Context_GetPrintModel(ctx) = (char**) argv + i ;

    } else if(String_Is(argv[i],"-modules",strlen(argv[i]))) {
      Context_GetPrintModule(ctx) = (char**) argv + i ;

    } else if(String_Is(argv[i],"-readonly",strlen(argv[i]))) {
      Context_GetReadOnly(ctx) = (char**) argv + i ;

    } else if(String_Is(argv[i],"-create",strlen(argv[i]))) {
      Context_GetCreateInputFile(ctx) = (char**) argv + i ;

    } else if(String_Is(argv[i],"-graph",strlen(argv[i]))) {
      Context_GetGraph(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing graph method") ;
      }

    } else if(String_Is(argv[i],"-iperm",strlen(argv[i]))) {
      Context_GetInversePermutation(ctx) = (char**) argv + i ;

    } else if(String_Is(argv[i],"-eordering",strlen(argv[i]))) {
      Context_GetElementOrdering(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing element ordering method") ;
      }

    } else if(String_Is(argv[i],"-nordering",strlen(argv[i]))) {
      Context_GetNodalOrdering(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing nodal ordering method") ;
      }

    } else if(String_Is(argv[i],"-postprocessing",strlen(argv[i]))) {
      Context_GetPostProcessing(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing post-processing method") ;
      }

    } else if(String_Is(argv[i],"-miscellaneous",strlen(argv[i]))) {
      Context_GetMiscellaneous(ctx) = (char**) argv + i ;

    } else if(String_Is(argv[i],"-test",strlen(argv[i]))) {
      Context_GetTest(ctx) = (char**) argv + i ;
      
    } else {
      Message_FatalError("Unknown option") ;
    }
    
  }
  
}
