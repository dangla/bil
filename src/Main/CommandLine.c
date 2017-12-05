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
#include "CommandLine.h"



static void PrintUsage(char*) ;
static void PrintInfo(void) ;
static CommandLine_t cmd1 ;
static CommandLine_t* cmd = &cmd1 ;



/* Global functions */
void (CommandLine_Initialize)(int argc,char** argv)
{
  Context_t* ctx = Context_GetInstance() ;
  int i ;
  
  CommandLine_GetNbOfArg(cmd) = argc ;
  CommandLine_GetArg(cmd) = argv ;


  if(argc == 1) {
    PrintUsage(argv[0]) ;
    exit(EXIT_SUCCESS) ;
  }
  
  
  /* Get line options */
  for(i = 1 ; i < argc ; i++) {
    
    if(argv[i][0] != '-') { /* File name */
      Context_GetInputFileName(ctx) = (char**) argv + i ;
      
    } else if(strncmp(argv[i],"-info",5) == 0) {
      PrintInfo() ;
      exit(EXIT_SUCCESS) ;
  
    } else if(strncmp(argv[i],"-help",5) == 0) {
      Help_HelpOnline() ;
      exit(EXIT_SUCCESS) ;
      
    } else if(!strncmp(argv[i],"-solver",strlen(argv[i]))) {
      Context_GetSolver(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing solver") ;
      }
    
    } else if(strncmp(argv[i],"-debug",strlen(argv[i])) == 0) {
      Context_GetDebug(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing name of data to be printed") ;
      }

    } else if(strncmp(argv[i],"-level",strlen(argv[i])) == 0) {
      Context_GetPrintLevel(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing level") ;
      }

    } else if(strncmp(argv[i],"-with",strlen(argv[i])) == 0) {
      Context_GetUseModule(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing module") ;
      }

    } else if(strncmp(argv[i],"-models",strlen(argv[i])) == 0) {
      if(i + 1 < argc) {
        if(!strncmp(argv[i + 1],"all",strlen(argv[i + 1]))) {
          Models_Print(NULL,stdout) ;
        } else {
          char* codename = argv[i + 1] ;
          Models_Print(codename,stdout) ;
        }
      } else Models_Print(NULL,NULL) ;
      exit(EXIT_SUCCESS) ;

    } else if(strncmp(argv[i],"-modules",strlen(argv[i])) == 0) {
      if(i + 1 < argc) {
        if(!strncmp(argv[i + 1],"all",strlen(argv[i + 1]))) {
          Modules_Print(NULL) ;
        } else {
          char* codename = argv[i + 1] ;
          Modules_Print(codename) ;
        }
      } else Modules_Print(NULL) ;
      exit(EXIT_SUCCESS) ;

    } else if(strncmp(argv[i],"-readonly",strlen(argv[i])) == 0) {
      Context_GetReadOnly(ctx) = (char**) argv + i ;

    } else if(strncmp(argv[i],"-graph",strlen(argv[i])) == 0) {
      Context_GetGraph(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing graph method") ;
      }

    } else if(strncmp(argv[i],"-iperm",strlen(argv[i])) == 0) {
      Context_GetInversePermutation(ctx) = (char**) argv + i ;

    } else if(strncmp(argv[i],"-eordering",strlen(argv[i])) == 0) {
      Context_GetElementOrdering(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing element ordering method") ;
      }

    } else if(strncmp(argv[i],"-nordering",strlen(argv[i])) == 0) {
      Context_GetNodalOrdering(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing nodal ordering method") ;
      }

    } else if(strncmp(argv[i],"-postprocessing",strlen(argv[i])) == 0) {
      Context_GetPostProcessing(ctx) = (char**) argv + i ;
      if(i + 1 < argc) {
        i++ ;
      } else {
        Message_FatalError("Missing post-processing method") ;
      }

    } else if(strncmp(argv[i],"-miscellaneous",strlen(argv[i])) == 0) {
      Context_GetMiscellaneous(ctx) = (char**) argv + i ;
      
    } else {
      Message_FatalError("Unknown option") ;
    }
    
  }
  
}


CommandLine_t*  (CommandLine_GetInstance)(void)
{
  return cmd ;
}


/* 
 * Intern functions 
 */


void PrintUsage(char* path)
{
  Message_Direct("Usage: %s [options] file\n",path) ;
  Message_Direct("\n") ;
  Message_Direct("Options:\n") ;
  
  Message_Direct("  -debug \"input\"       Display the data structure of \"input\".\n") ;
  #if 0
  Message_Direct("                       Available inputs are:\n") ;
  Message_Direct("                       - geom   : geometry\n") ;
  Message_Direct("                       - mesh   : mesh\n") ;
  Message_Direct("                       - mate   : materials\n") ;
  Message_Direct("                       - field  : fields\n") ;
  Message_Direct("                       - init   : initialization\n") ;
  Message_Direct("                       - func   : time functions\n") ;
  Message_Direct("                       - bcond  : boundary conditions\n") ;
  Message_Direct("                       - load   : loads\n") ;
  Message_Direct("                       - poin   : points\n") ;
  Message_Direct("                       - obval  : objective variations\n") ;
  Message_Direct("                       - iter   : iterative process\n") ;
  Message_Direct("                       - time   : time steps\n") ;
  Message_Direct("                       - all    : all the dataset above\n") ;
  Message_Direct("                       - inter  : interpolation functions\n") ;
  Message_Direct("                       - continuity: equations and unknowns at nodes\n") ;
  Message_Direct("                       - matrix : the matrix\n") ;
  Message_Direct("                       - residu : the residu\n") ;
  #endif
  
  Message_Direct("  -graph \"fmt\"         Print the graph in \"file.graph\"\n") ;
  Message_Direct("                       in the format \"fmt\".\n") ;
  #if 0
  Message_Direct("                       Available formats are:\n") ;
  Message_Direct("                       - metis: for METIS (if installed),\n") ;
  Message_Direct("                       - hsl_mc40: for HSL_MC40,\n") ;
  #endif
  
  Message_Direct("  -help                Online help.\n") ;
  Message_Direct("  -info                Display general informations.\n") ;
  
  Message_Direct("  -iperm               Print the inverse permutation of nodes \n") ;
  Message_Direct("                       in \"file.graph.iperm\" by \"HSL_MC40\".\n") ;
  
  Message_Direct("  -eordering \"meth\"    Print the inverse permutation of elements\n") ;
  Message_Direct("                       in \"file.graph.iperm\" by the method \"meth\".\n") ;
  #if 0
  Message_Direct("                       Available methods are:\n") ;
  Message_Direct("                       - hsl_mc43\n") ;
  #endif
  
  Message_Direct("  -nordering \"meth\"    Print the inverse permutation of nodes \n") ;
  Message_Direct("                       in \"file.graph.iperm\" by the method \"meth\".\n") ;
  #if 0
  Message_Direct("                       Available methods are:\n") ;
  Message_Direct("                       - hsl_mc40 (same as -iperm)\n") ;
  #endif
  
  /*
  Message_Direct("  -verbosity \"degre\"   Define the verbosity level:\n") ;
  Message_Direct("                       \"degre\" = 0  : nothing\n") ;
  Message_Direct("                       \"degre\" = 1  : normal (default)\n") ;
  Message_Direct("                       \"degre\" = 2  : extra\n") ;
  */
  Message_Direct("  -level \"I\"           Define the screen output level \"I\".\n") ;
  
  Message_Direct("  -model               Display the available models.\n") ;
  
  Message_Direct("  -model \"mymodel\"     Display an example of material properties\n") ;
  Message_Direct("                       for the model \"mymodel\".\n") ;
  
  Message_Direct("  -modules             Display the available modules.\n") ;
  
  Message_Direct("  -solver \"meth\"       Use a solver defined by the method \"meth\".\n") ;
  #if 0
  Message_Direct("                       Available methods are:\n") ;
  Message_Direct("                       - crout : CROUT method (default),\n") ;
  Message_Direct("                       - slu   : SuperLU (if installed).\n") ;
  #endif
  
  Message_Direct("  -post \"fmt\"          Generates the post-processing files \n") ;
  Message_Direct("                       \"file.posI\" in the format \"fmt\".\n") ;
  #if 0
  Message_Direct("                       Available formats are:\n") ;
  Message_Direct("                       - GmshParsed: for GMSH parsed file format.\n") ;
  Message_Direct("                       - GmshASCII: for GMSH ASCII file format.\n") ;
  #endif
  
  Message_Direct("  -readonly            Read \"file\" only.\n") ;
  Message_Direct("  -with \"mod\"          Use the module \"mod\".\n") ;
}



void PrintInfo(void)
{
  char bil_progname[]  = BIL_PROGNAME ;
  char bil_copyright[] = BIL_COPYRIGHT ;
  char bil_version[]   = "Version        : " BIL_VERSION ;
  char bil_license[]   = "License        : " BIL_SHORT_LICENSE ;
  char bil_os[]        = "Build OS       : " BIL_OS ;
  char bil_date[]      = "Build date     : " BIL_DATE ;
  char bil_host[]      = "Build host     : " BIL_HOST ;
  char bil_packager[]  = "Packager       : " BIL_PACKAGER ;
  char bil_url[]       = "Web site       : " BIL_URL ;
  char bil_email[]     = "Contact        : " BIL_EMAIL ;
  
  Message_Direct("%s\n", bil_progname) ;
  Message_Direct("%s\n", bil_copyright) ;
  Message_Direct("%s\n", bil_version) ;
  Message_Direct("%s\n", bil_license) ;
  Message_Direct("%s\n", bil_os) ;
  Message_Direct("%s\n", bil_date) ;
  Message_Direct("%s\n", bil_host) ;
  Message_Direct("%s\n", bil_packager) ;
  Message_Direct("%s\n", bil_url) ;
  Message_Direct("%s\n", bil_email) ;
}
