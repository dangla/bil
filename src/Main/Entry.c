#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "DataSet.h"
#include "Message.h"
#include "Context.h"
#include "CommandLine.h"
#include "Help.h"
#include "Modules.h"
#include "OutputFiles.h"
#include "PosFilesForGMSH.h"
#include "Models.h"
#include "Entry.h"
#include "Exception.h"
#include "BilVersion.h"
#include "BilInfo.h"
#include "Session.h"
#include "Mry.h"
#include "TypeId.h"
#include "SharedMS.h"
#include "DistributedMS.h"



static void   (Entry_PrintUsage)(char*) ;
static void   (Entry_PrintInfo)(void) ;
static void   (Entry_CLI)(Entry_t*) ;


int (Entry_Main)(int argc,char** argv)
{
  #if DistributedMS_APIis(MPI)
    MPI_Init(&argc,&argv) ;
  #endif
  {
    Entry_t* entry = Entry_Create(argc,argv) ;
  
    Entry_Execute(entry) ;
  
    Entry_Delete(entry) ;
  }
  #if DistributedMS_APIis(MPI)
    MPI_Finalize() ;
  #endif
  
  return(0) ;
}



Entry_t*    (Entry_Create)(int argc,char** argv)
{
  Entry_t* entry = (Entry_t*) Mry_New(Entry_t) ;
  
  Session_Open() ;
  
  {
    Context_t* ctx = Context_Create(argc,argv) ;
    
    Entry_GetContext(entry) = ctx ;
  }
    
  Session_Close() ;
  
  return(entry) ;
}



void (Entry_Delete)(void* self)
{
  Entry_t* entry = (Entry_t*) self ;
  
  {
    Context_t* ctx = Entry_GetContext(entry) ;
    
    if(ctx) {
      Context_Delete(ctx) ;
      free(ctx) ;
      Entry_GetContext(entry) = NULL ;
    }
  }
}



int (Entry_Execute)(Entry_t* entry)
{
  int val = 0 ;
  
  Session_Open() ;
  
  {
    val = Exception_SaveEnvironment ;
    
    if(!val) {
      /* Command-line interface */
      Entry_CLI(entry) ;
      /* Graphical user interface */
      //Entry_GUI(entry) ; /* Not done yet */
    } else {
      Message_FatalError("An exception occurs with value %d\n",val) ;
    }
  }
    
  Session_Close() ;
  
  return(val) ;
}




/* 
 * Intern functions 
 */
 
#include "CementSolutionChemistry.h"
#include "HardenedCementChemistry.h"
#include "LogActivityOfWaterInBrine.h"

void Entry_CLI(Entry_t* entry)
/** Command-line interface */
{
  Context_t* ctx = Entry_GetContext(entry) ;

  Message_Info("Started on %s",Message_LaunchDate()) ;
  
  if(Context_IsPrintUsage(ctx)) {
    char** argv = (char**) Context_GetPrintUsage(ctx) ;
    
    Entry_PrintUsage(argv[0]) ;
    return ;
  }
  
  if(Context_IsPrintInfo(ctx)) {
    Entry_PrintInfo() ;
    return ;
  }
  
  if(Context_IsHelpOnline(ctx)) {
    Help_HelpOnline() ;
    return ;
  }
  
  if(Context_IsPrintModel(ctx)) {
    char** argv = (char**) Context_GetPrintModel(ctx) ;
    
    /* The standard requires that argv[argc] be a null pointer */
    if(argv[1]) {
      if(!strncmp(argv[1],"all",strlen(argv[1]))) {
        Models_Print(NULL,stdout) ;
      } else {
        char* codename = argv[1] ;
        Models_Print(codename,stdout) ;
      }
    } else {
      Models_Print(NULL,NULL) ;
    }
    return ;
  }
  
  if(Context_IsPrintModule(ctx)) {
    char** argv = (char**) Context_GetPrintModule(ctx) ;
    
    /* The standard requires that argv[argc] be a null pointer */
    if(argv[1]) {
      if(!strncmp(argv[1],"all",strlen(argv[1]))) {
        Modules_Print(NULL) ;
      } else {
        char* codename = argv[1] ;
        Modules_Print(codename) ;
      }
    } else Modules_Print(NULL) ;
    return ;
  }
  
  if(Context_IsReadOnly(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset = DataSet_Create(filename,options) ;
    
    DataSet_Delete(dataset) ;
    free(dataset) ;
    return ;
  }
  
  if(Context_IsCreateInputFile(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    FILE* str = fopen(filename,"r") ;
    
    if(!str) {
      Help_WriteData(filename) ;
      Message_Info("To start the computation, type bil %s\n",filename) ;
    } else {
      fclose(str) ;
      Message_Info("File %s already exists.\n",filename) ;
      Message_Exit ;
    }
    
    return ;
  }
  
  if(Context_IsGraph(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    char* method = Options_GetGraphMethod(options) ;
      
    Message_Direct("Graph (method %s)\n",method) ;
    Mesh_WriteGraph(mesh,filename,method) ;
    DataSet_Delete(dataset) ;
    free(dataset) ;
    return ;
  }
  
  if(Context_IsInversePermutation(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    char method[] = "hsl" ;
      
    Message_Direct("Inverse permutations (method %s)\n",method) ;
    Mesh_WriteInversePermutation(mesh,filename,method) ;
    DataSet_Delete(dataset) ;
    free(dataset) ;
    return ;
  }
  
  if(Context_IsNodalOrdering(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    char* method = Options_GetNodalOrderingMethod(options) ;
      
    Message_Direct("Nodal ordering (method %s)\n",method) ;
    Mesh_WriteInversePermutation(mesh,filename,method) ;
    DataSet_Delete(dataset) ;
    free(dataset) ;
    return ;
  }
  
  if(Context_IsElementOrdering(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    char* method = Options_GetElementOrderingMethod(options) ;
      
    Message_Direct("Element ordering (method %s)\n",method) ;
    Mesh_WriteInversePermutation(mesh,filename,method) ;
    DataSet_Delete(dataset) ;
    free(dataset) ;
    return ;
  }
  
  if(Context_IsPostProcessing(ctx)) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset =  DataSet_Create(filename,options) ;
    char* method = Options_GetPostProcessingMethod(options) ;
      
    Message_Direct("Post-processing\n") ;
    
    /* GMSH file formats */
    if(strncmp(method,"Gmsh",4) == 0) {
      PosFilesForGMSH_t* pf4gmsh = PosFilesForGMSH_Create(dataset) ;
      char* fmt = method + 4 ;
      
      /* GMSH Parsed file format */
      if(strncmp(fmt,"ParsedFileFormat",strlen(fmt)) == 0) {
        PosFilesForGMSH_ParsedFileFormat(pf4gmsh) ;
      /* GMSH ASCII file format */
      } else if(strncmp(fmt,"ASCIIFileFormat",strlen(fmt)) == 0) {
        PosFilesForGMSH_ASCIIFileFormat(pf4gmsh) ;
      } else {
        Message_FatalError("Format not available") ;
      }
      
      PosFilesForGMSH_Delete(pf4gmsh) ;
      free(pf4gmsh) ;
      
    } else {
      Message_FatalError("Format not available") ;
    }

    DataSet_Delete(dataset) ;
    free(dataset) ;
    return ;
  }
  
  if(Context_IsMiscellaneous(ctx)) {
    double temp = atof(((char**) Context_GetMiscellaneous(ctx))[1]) ;
    
    {
      HardenedCementChemistry_t* hcc = HardenedCementChemistry_Create() ;
    
      HardenedCementChemistry_SetRoomTemperature(hcc,temp) ;
      HardenedCementChemistry_PrintChemicalConstants(hcc) ;
      HardenedCementChemistry_Delete(hcc) ;
    }
    
    {
      CementSolutionChemistry_t* csc = CementSolutionChemistry_Create() ;
      
      CementSolutionChemistry_GetRoomTemperature(csc) = temp ;
      CementSolutionChemistry_PrintChemicalConstants(csc) ;
      CementSolutionChemistry_Delete(csc) ;
    }
    
    return ;
  }
  
  if(Context_IsTest(ctx)) {
    #if 0
    {
      LogActivityOfWaterInBrine_t* law = LogActivityOfWaterInBrine_Create() ;
      LogActivityOfWaterInBrine_Delete(law) ;
      return ;
    }
    #endif
    #if 0
    {
      char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
      Options_t* options = Context_GetOptions(ctx) ;
      DataSet_t* dataset = DataSet_Create1(filename,options) ;
      
      DataSet_Delete(dataset) ;
      free(dataset) ;
      return ;
    }
    #endif
    
    return ;
  }
  
  if(1) {
    char* filename = ((char**) Context_GetInputFileName(ctx))[0] ;
    Options_t* options = Context_GetOptions(ctx) ;
    DataSet_t* dataset =  DataSet_Create(filename,options) ;
    Module_t* module = DataSet_GetModule(dataset) ;

    #if SharedMS_APIisNot(None)
    {
      int nthreads = Options_NbOfThreads(options) ;
      
      SharedMS_SetTheNbOfThreads(nthreads) ;
      Message_Direct("Nb of requested threads: %d\n",nthreads) ;
      Message_Direct("Multithreading API: "Utils_STR(SharedMS_API)"\n") ;
    }
    #endif

    #if DistributedMS_APIisNot(None)
    {
      int nprocs = DistributedMS_NbOfProcessors ;
      
      Message_Direct("Nb of processors: %d\n",nprocs) ;
      Message_Direct("Distributed Memory API: "Utils_STR(DistributedMS_API)"\n") ;
    }
    #endif
  
    Message_Direct("Calculation of %s\n",filename) ;
    Message_Direct("Enter in %s\n",Module_GetCodeNameOfModule(module)) ;
    Module_ComputeProblem(module,dataset) ;
    Message_Direct("End of calculation\n") ;
    
    DataSet_Delete(dataset) ;
    free(dataset) ;
    
    Message_Info("CPU time %g seconds\n",Message_CPUTime()) ;
    return ;
  }
  
  return ;
}




void Entry_PrintUsage(char* path)
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
  Message_Direct("                       - ma38  : HSL-MA38 (if installed).\n") ;
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



void Entry_PrintInfo(void)
{
  char bil_progname[]  = BIL_PROGNAME ;
  char bil_copyright[] = BIL_COPYRIGHT ;
  char bil_version[]   = "Version           : " BIL_VERSION ;
  char bil_license[]   = "License           : " BIL_SHORT_LICENSE ;
  char bil_os[]        = "Build OS          : " BIL_OS ;
  char bil_date[]      = "Build date        : " BIL_DATE ;
  char bil_host[]      = "Build host        : " BIL_HOST ;
  char bil_packager[]  = "Packager          : " BIL_PACKAGER ;
  char bil_url[]       = "Web site          : " BIL_URL ;
  char bil_email[]     = "Contact           : " BIL_EMAIL ;
  char bil_mtapi[]     = "Multithreading API: " Utils_STR(SharedMS_API) ;
  char bil_nthreads[]  = "Nb of CPU threads : " ;
  int nthreads = SharedMS_NbOfLogicalCores ;
  
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
  Message_Direct("%s\n", bil_mtapi) ;
  Message_Direct("%s%d\n",bil_nthreads,nthreads) ;
}
