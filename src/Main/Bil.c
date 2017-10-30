#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "DataSet.h"
#include "Message.h"
#include "Context.h"
#include "CommandLine.h"
#include "Help.h"
#include "Modules.h"
#include "OutputFiles.h"
#include "DataBases/CementSolutionChemistry.h"
#include "DataBases/HardenedCementChemistry.h"
#include "Models.h"
#include "Bil.h"
#include "Exception.h"



static void   (Bil_Initialize)(int,char**) ;
static void   (Bil_CLI)(void) ;



int Bil_Main(int argc,char** argv)
{
  
  Bil_Initialize(argc,argv) ;
  
  {
    int val = Exception_SaveEnvironment ;
    
    if(!val) {
      Bil_CLI() ;
    } else {
      Message_Direct("An exception occurs with value %d\n",val) ;
    }
  
    return(val) ;
  }
  
  return(0) ;
}



void Bil_Initialize(int argc,char** argv)
{
  Message_Initialize() ;
  CommandLine_Initialize(argc,argv) ;
}


void Bil_CLI(void)
/** Command-line interface */
{
  Options_t* options = Options_Create() ;
  char* filename = Context_InputFileName ;

  Message_Info("Started on %s",Message_LaunchDate()) ;
  
  if(Context_IsReadOnly) {
    DataSet_Create(filename,options) ;
    return ;
  }
  
  if(Context_IsGraph) {
    DataSet_t* jdd =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(jdd) ;
    char* method = Options_GetGraphMethod(options) ;
      
    Message_Direct("Graph (method %s)\n",method) ;
    Mesh_WriteGraph(mesh,filename,method) ;
    return ;
  }
  
  if(Context_IsInversePermutation) {
    DataSet_t* jdd =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(jdd) ;
    char method[] = "hsl" ;
      
    Message_Direct("Inverse permutations (method %s)\n",method) ;
    Mesh_WriteInversePermutation(mesh,filename,method) ;
    return ;
  }
  
  if(Context_IsNodalOrdering) {
    DataSet_t* jdd =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(jdd) ;
    char* method = Options_GetNodalOrderingMethod(options) ;
      
    Message_Direct("Nodal ordering (method %s)\n",method) ;
    Mesh_WriteInversePermutation(mesh,filename,method) ;
    return ;
  }
  
  if(Context_IsElementOrdering) {
    DataSet_t* jdd =  DataSet_Create(filename,options) ;
    Mesh_t* mesh = DataSet_GetMesh(jdd) ;
    char* method = Options_GetElementOrderingMethod(options) ;
      
    Message_Direct("Element ordering (method %s)\n",method) ;
    Mesh_WriteInversePermutation(mesh,filename,method) ;
    return ;
  }
  
  if(Context_IsPostProcessing) {
    DataSet_t* jdd =  DataSet_Create(filename,options) ;
    int n_dates = Dates_GetNbOfDates(DataSet_GetDates(jdd)) ;
    int n_points = Points_GetNbOfPoints(DataSet_GetPoints(jdd)) ;
    OutputFiles_t* outputfiles = OutputFiles_Create(filename,n_dates,n_points) ;
    char* method = Options_GetPostProcessingMethod(options) ;
      
    Message_Direct("Post-processing\n") ;
    
    /* GMSH Parsed file format */
    if(strncmp(method,"GmshParsedFileFormat",strlen(method)) == 0) {
      OutputFiles_PostProcessForGmshParsedFileFormat(outputfiles,jdd) ;
      
    /* GMSH ASCII file format */
    } else if(strncmp(method,"GmshASCIIFileFormat",strlen(method)) == 0) {
      OutputFiles_PostProcessForGmshASCIIFileFormat(outputfiles,jdd) ;
      
    } else {
      Message_FatalError("Format not available") ;
    }

    OutputFiles_Delete(&outputfiles) ;
    return ;
  }
  
  if(Context_IsMiscellaneous) {
    double temp = atof(((char**) Context_IsMiscellaneous)[1]) ;
    
    HardenedCementChemistry_SetTemperature(temp) ;
    
    {
      HardenedCementChemistry_t* hcc = HardenedCementChemistry_GetInstance() ;
      
      HardenedCementChemistry_PrintChemicalConstants(hcc) ;
    }
    
    {
      CementSolutionChemistry_t* csc = CementSolutionChemistry_Create(1) ;
      
      CementSolutionChemistry_GetTemperature(csc) = temp ;
      CementSolutionChemistry_PrintChemicalConstants(csc) ;
    }
    
    
    //OutputFiles_t* of = OutputFiles_Create(filename,4,2) ;
    //OutputFiles_Delete(&of) ;
    return ;
  }
  
  if(1) {
    DataSet_t* jdd =  DataSet_Create(filename,options) ;
    char* codename = Options_GetModule(options) ;
    Modules_t* modules = DataSet_GetModules(jdd) ;
    Module_t* module_i = Modules_FindModule(modules,codename) ;
  
    Message_Direct("Calculation\n") ;
    
    Module_ComputeProblem(module_i,jdd) ;
    Message_Info("CPU time %g seconds\n",Message_CPUTime()) ;
    return ;
  }
  
  return ;
}
