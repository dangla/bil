#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Help.h"
#include "Message.h"
#include "DataFile.h"
#include "Geometry.h"
#include "Dates.h"
#include "Points.h"
#include "ObVals.h"
#include "TimeStep.h"
#include "IterProcess.h"
#include "Options.h"
#include "Modules.h"
#include "DataSet.h"
#include "Units.h"
#include "String_.h"
//#include "Parser.h"



/* Extern functions */

DataSet_t*  (DataSet_Create)(char* filename,Options_t* opt)
{
  DataSet_t* dataset = DataSet_New() ;
  char*   debug  = Options_GetPrintedInfos(opt) ;
  
  
  DataSet_GetOptions(dataset) = opt ;
  
  
  /* DataFile */
  {
    DataFile_t* datafile = DataFile_Create(filename) ;
    
    DataFile_RemoveComments(datafile) ;
    DataFile_GetParent(datafile) = dataset ;
  
    DataSet_GetDataFile(dataset) = datafile ;
  
    if(DataFile_DoesNotExist(datafile)) {
      Message_Info("File %s not found\n",filename) ;
      Message_Exit ;
    }
  }
  if(!strcmp(debug,"data")) DataSet_PrintData(dataset,debug) ;

  Message_Direct("Reading %s\n",filename) ;


  /* Units */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Units_t* units = Units_Create(datafile) ;
    
    DataSet_GetUnits(dataset) = units ;
  }
  
  
  /* Geometry */
  {
    DataFile_t* datafile = DataSet_GetDataFile(dataset) ;
    Geometry_t* geometry = Geometry_Create(datafile) ;
  
    DataSet_GetGeometry(dataset) = geometry ;
  }
  if(!strcmp(debug,"geom")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Fields */
  {
    DataFile_t*  datafile = DataSet_GetDataFile(dataset) ;
    Fields_t*    fields = Fields_Create(datafile) ;
  
    DataSet_GetFields(dataset) = fields ;
  }
  if(!strcmp(debug,"field")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Functions */
  {
    DataFile_t*  datafile = DataSet_GetDataFile(dataset) ;
    Functions_t* functions = Functions_Create(datafile) ;
  
    DataSet_GetFunctions(dataset) = functions ;
  }
  if(!strcmp(debug,"func")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Models */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Geometry_t*    geometry = DataSet_GetGeometry(dataset) ;
    Models_t*      models   = Models_Create(datafile,geometry) ;
  
    DataSet_GetModels(dataset) = models ;
  }
  if(!strcmp(debug,"model")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Materials */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Geometry_t*    geometry = DataSet_GetGeometry(dataset) ;
    Fields_t*      fields = DataSet_GetFields(dataset) ;
    Functions_t*   functions = DataSet_GetFunctions(dataset) ;
    Models_t*      models = DataSet_GetModels(dataset) ;
    Materials_t*   materials = Materials_Create(datafile,geometry,fields,functions,models) ;
  
    DataSet_GetMaterials(dataset) = materials ;
  }
  if(!strcmp(debug,"mate")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Mesh */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Geometry_t*    geometry = DataSet_GetGeometry(dataset) ;
    Materials_t*   materials = DataSet_GetMaterials(dataset) ;
    Mesh_t*        mesh = Mesh_Create(datafile,materials,geometry) ;
  
    DataSet_GetMesh(dataset) = mesh ;
  }
  if(!strcmp(debug,"mesh")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Initial conditions */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Fields_t*      fields = DataSet_GetFields(dataset) ;
    Functions_t*   functions = DataSet_GetFunctions(dataset) ;
    IConds_t*      iconds = IConds_Create(datafile,fields,functions) ;
  
    DataSet_GetIConds(dataset) = iconds ;
  }
  if(!strcmp(debug,"init")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Loads */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Fields_t*      fields = DataSet_GetFields(dataset) ;
    Functions_t*   functions = DataSet_GetFunctions(dataset) ;
    Loads_t*       loads = Loads_Create(datafile,fields,functions) ;
  
    DataSet_GetLoads(dataset) = loads ;
  }
  if(!strcmp(debug,"load")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Boundary conditions */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Fields_t*      fields = DataSet_GetFields(dataset) ;
    Functions_t*   functions = DataSet_GetFunctions(dataset) ;
    BConds_t*      bconds = BConds_Create(datafile,fields,functions) ;
  
    DataSet_GetBConds(dataset) = bconds ;
  }
  if(!strcmp(debug,"bcond")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Points */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Mesh_t*        mesh = DataSet_GetMesh(dataset) ;
    Points_t*      points = Points_Create(datafile,mesh) ;
  
    DataSet_GetPoints(dataset) = points ;
  }
  if(!strcmp(debug,"points")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Dates */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Dates_t*       dates = Dates_Create(datafile) ;
    
    DataSet_GetDates(dataset) = dates ;
  }
  if(!strcmp(debug,"dates")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Objective variations */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    Materials_t*   materials = DataSet_GetMaterials(dataset) ;
    Mesh_t*        mesh = DataSet_GetMesh(dataset) ;
    ObVals_t*      obvals = ObVals_Create(datafile,mesh,materials) ;
  
    DataSet_GetObVals(dataset) = obvals ;
  }
  if(!strcmp(debug,"obval")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Time steps */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    ObVals_t*      obvals = DataSet_GetObVals(dataset) ;
    TimeStep_t*    timestep = TimeStep_Create(datafile,obvals) ;
    
    DataSet_GetTimeStep(dataset) = timestep ;
  }
  if(!strcmp(debug,"time")) DataSet_PrintData(dataset,debug) ;
  
  
  /* Iterative process */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(dataset) ;
    ObVals_t*      obvals = DataSet_GetObVals(dataset) ;
    IterProcess_t* iterprocess = IterProcess_Create(datafile,obvals) ;
    
    DataSet_GetIterProcess(dataset) = iterprocess ;
  }
  if(!strcmp(debug,"iter")) DataSet_PrintData(dataset,debug) ;


  Message_Direct("End of reading %s\n",filename) ;
  Message_Direct("\n") ;
  
  
  /* Module */
  {
    Module_t*  module  = Module_New() ;
    char* codename = Options_GetModule(opt) ;

    Module_Initialize(module,codename) ;
    Module_GetNbOfSequences(module) = Options_GetNbOfSequences(opt) ;
    
    DataSet_GetModule(dataset) = module ;
    
    //DataSet_GetModule(dataset) = Modules_FindModule(modules,Options_GetModule(opt)) ;
  }
  if(!strcmp(debug,"module")) DataSet_PrintData(dataset,debug) ;
  
  
  
  /* Set up the system of equations before
   * passing through Elements_DefineProperties */
  {
    BConds_t*      bconds = DataSet_GetBConds(dataset) ;
    Mesh_t*        mesh = DataSet_GetMesh(dataset) ;
    
    Mesh_GetNbOfMatrices(mesh) = Options_GetNbOfSequences(opt) ;
  
    Mesh_SetMatrixRowColumnIndexes(mesh,bconds) ;
  }
  
  
  /* Other printings in debug mode */
  if(!strcmp(debug,"continuity")) DataSet_PrintData(dataset,debug) ;
  //if(!strcmp(debug,"numbering")) DataSet_PrintData(dataset,debug) ;
  //if(!strcmp(debug,"inter")) DataSet_PrintData(dataset,debug) ;
  if(!strcmp(debug,"all")) DataSet_PrintData(dataset,debug) ;
  
  return(dataset) ;
}



void (DataSet_Delete)(void* self)
{
  DataSet_t* dataset = (DataSet_t*) self ;
  
  
  {
    Units_t* units = DataSet_GetUnits(dataset) ;
    
    if(units) {
      Units_Delete(units) ;
      free(units) ;
    }
    
    DataSet_GetUnits(dataset) = NULL ;
  }
  
  {
    DataFile_t* datafile = DataSet_GetDataFile(dataset) ;
    
    if(datafile) {
      DataFile_Delete(datafile) ;
      free(datafile) ;
    }
    
    DataSet_GetDataFile(dataset) = NULL ;
  }
  
  {
    Geometry_t* geometry = DataSet_GetGeometry(dataset) ;
    
    if(geometry) {
      Geometry_Delete(geometry) ;
      free(geometry) ;
    }
    
    DataSet_GetGeometry(dataset) = NULL ;
  }
  
  {
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    
    if(mesh) {
      Mesh_Delete(mesh) ;
      free(mesh) ;
    }
    
    DataSet_GetMesh(dataset) = NULL ;
  }
  
  {
    Materials_t* materials = DataSet_GetMaterials(dataset) ;
    
    if(materials) {
      Materials_Delete(materials) ;
      free(materials) ;
    }
    
    DataSet_GetMaterials(dataset) = NULL ;
  }
  
  {
    Dates_t* dates = DataSet_GetDates(dataset) ;
    
    if(dates) {
      Dates_Delete(dates) ;
      free(dates) ;
    }
    
    DataSet_GetDates(dataset) = NULL ;
  }
  
  {
    Points_t* points = DataSet_GetPoints(dataset) ;
    
    if(points) {
      Points_Delete(points) ;
      free(points) ;
    }
    
    DataSet_GetPoints(dataset) = NULL ;
  }
  
  {
    IConds_t* iconds = DataSet_GetIConds(dataset) ;
    
    if(iconds) {
      IConds_Delete(iconds) ;
      free(iconds) ;
    }
    
    DataSet_GetIConds(dataset) = NULL ;
  }
  
  {
    BConds_t* bconds = DataSet_GetBConds(dataset) ;
    
    if(bconds) {
      BConds_Delete(bconds) ;
      free(bconds) ;
    }
    
    DataSet_GetBConds(dataset) = NULL ;
  }
  
  {
    Loads_t* loads = DataSet_GetLoads(dataset) ;
    
    if(loads) {
      Loads_Delete(loads) ;
      free(loads) ;
    }
    
    DataSet_GetLoads(dataset) = NULL ;
  }
  
  {
    Functions_t* functions = DataSet_GetFunctions(dataset) ;
    
    if(functions) {
      Functions_Delete(functions) ;
      free(functions) ;
    }
    
    DataSet_GetFunctions(dataset) = NULL ;
  }
  
  {
    Fields_t* fields = DataSet_GetFields(dataset) ;
    
    if(fields) {
      Fields_Delete(fields) ;
      free(fields) ;
    }
    
    DataSet_GetFields(dataset) = NULL ;
  }
  
  //IntFcts_t*     intfcts ; */
  
  {
    ObVals_t* obvals = DataSet_GetObVals(dataset) ;
    
    if(obvals) {
      ObVals_Delete(obvals) ;
      free(obvals) ;
    }
    
    DataSet_GetObVals(dataset) = NULL ;
  }
  
  //Models_t*      models ;
  
  {
    TimeStep_t* timestep = DataSet_GetTimeStep(dataset) ;
    
    if(timestep) {
      TimeStep_Delete(timestep) ;
      free(timestep) ;
    }
    
    DataSet_GetTimeStep(dataset) = NULL ;
  }
  
  {
    IterProcess_t* iterprocess = DataSet_GetIterProcess(dataset) ;
    
    if(iterprocess) {
      IterProcess_Delete(iterprocess) ;
      free(iterprocess) ;
    }
    
    DataSet_GetIterProcess(dataset) = NULL ;
  }
  
  {
    DataSet_GetOptions(dataset) = NULL ;
  }
  
  //Modules_t*     modules ;
  
  {
    Module_t* module = DataSet_GetModule(dataset) ;
    
    if(module) {
      Module_Delete(module) ;
      free(module) ;
    }
    
    DataSet_GetModule(dataset) = NULL ;
  }
}



#if 0
DataSet_t*  (DataSet_Create1)(char* filename,Options_t* opt)
{
  DataSet_t* dataset = DataSet_New() ;
  
  DataSet_GetOptions(dataset) = opt ;
  
  {
    DataFile_t* datafile = DataFile_Create(filename) ;
  
    DataSet_GetDataFile(dataset) = datafile ;
  
    if(DataFile_DoesNotExist(datafile)) {
      Help_WriteData(filename) ;
      Message_Info("To start the computation, type bil %s\n",filename) ;
      Message_Exit ;
    }
  }
  
  Parser_ParseFile(dataset) ;
  
  {
    char*   debug  = Options_GetPrintedInfos(opt) ;
    
    DataSet_PrintData(dataset,debug) ;
  }
  
  return(dataset) ;
}
#endif



/* Local functions */


#define DATASET       (dataset)
#define DATAFILE      DataSet_GetDataFile(DATASET)
#define GEOMETRY      DataSet_GetGeometry(DATASET)
#define MESH          DataSet_GetMesh(DATASET)
#define MODELS        DataSet_GetModels(DATASET)
#define MATERIALS     DataSet_GetMaterials(DATASET)
#define FIELDS        DataSet_GetFields(DATASET)
#define ICONDS        DataSet_GetIConds(DATASET)
#define BCONDS        DataSet_GetBConds(DATASET)
#define FUNCTIONS     DataSet_GetFunctions(DATASET)
#define LOADS         DataSet_GetLoads(DATASET)
#define POINTS        DataSet_GetPoints(DATASET)
#define DATES         DataSet_GetDates(DATASET)
#define OBVALS        DataSet_GetObVals(DATASET)
#define ITERPROCESS   DataSet_GetIterProcess(DATASET)
#define TIMESTEP      DataSet_GetTimeStep(DATASET)

#define N_EL          Mesh_GetNbOfElements(MESH)
#define EL            Mesh_GetElement(MESH)
#define ELTS          Mesh_GetElements(MESH)

#define INTFCTS       Elements_GetIntFcts(ELTS)

#define NOM           DataFile_GetFileName(DATAFILE)

//#define DIM           Mesh_GetDimension(MESH)
//#define SYMMETRY      Mesh_GetSymmetry(MESH) 
//#define COORSYS       Mesh_GetCoordinateSystem(MESH)

#define DIM           Geometry_GetDimension(GEOMETRY)
#define SYMMETRY      Geometry_GetSymmetry(GEOMETRY) 
#define COORSYS       Geometry_GetCoordinateSystem(GEOMETRY)

#define N_NO          Mesh_GetNbOfNodes(MESH)
#define NO            Mesh_GetNode(MESH)

#define N_MODELS      Models_GetNbOfModels(MODELS)
#define MODEL         Models_GetModel(MODELS)

#define N_MAT         Materials_GetNbOfMaterials(MATERIALS)
#define MAT           Materials_GetMaterial(MATERIALS)

#define N_CH          Fields_GetNbOfFields(FIELDS)
#define CH            Fields_GetField(FIELDS)

#define N_IC          IConds_GetNbOfIConds(ICONDS)
#define IC            IConds_GetICond(ICONDS)

#define N_FN          Functions_GetNbOfFunctions(FUNCTIONS)
#define FN            Functions_GetFunction(FUNCTIONS)

#define N_CL          BConds_GetNbOfBConds(BCONDS)
#define CL            BConds_GetBCond(BCONDS)

#define N_CG          Loads_GetNbOfLoads(LOADS)
#define CG            Loads_GetLoad(LOADS)

#define N_FI          IntFcts_GetNbOfIntFcts(INTFCTS)
#define FI            IntFcts_GetIntFct(INTFCTS)

#define N_POINTS      Points_GetNbOfPoints(POINTS)

#define N_DATES       Dates_GetNbOfDates(DATES)

#define N_OBJ         ObVals_GetNbOfObVals(OBVALS)
#define OBJ           ObVals_GetObVal(OBVALS)


#define PRINT(...) \
        fprintf(stdout,__VA_ARGS__)
        //Message_Direct(__VA_ARGS__)


void DataSet_PrintData(DataSet_t* dataset,char* mot)
{
  static int i_debug=0 ;
  
  if(!strcmp(mot,"\0")) return ;

  PRINT("\n") ;
  PRINT("debug(%d)\n",i_debug++) ;
  PRINT("-----\n") ;
  
  /* File content
   * ------------ */
  if(DataFile_GetFileContent(DataSet_GetDataFile(dataset)) && (!strncmp(mot,"data file content",4) || !strncmp(mot,"all",3))) {
    PRINT("\n") ;
    PRINT("Data file content:\n") ;
    
    PRINT("%s",DataFile_GetFileContent(DataSet_GetDataFile(dataset))) ;
    PRINT("\n") ;
  }

  /* Geometry
   * -------- */
  if(DataSet_GetGeometry(dataset) && (!strncmp(mot,"geometry",4) || !strncmp(mot,"all",3))) {
    PRINT("\n") ;
    PRINT("Geometry:\n") ;
    
    PRINT("\t Dimension = %dD\n",DIM) ;
    PRINT("\t Symmetry = ") ;
    
    if(0) {
      
    } else if(Symmetry_IsCylindrical(SYMMETRY)) {
      PRINT("Axisymmetrical\n") ;
      
    } else if(Symmetry_IsSpherical(SYMMETRY)) {
      PRINT("Spherical\n") ;

    } else if(Symmetry_IsPlane(SYMMETRY)) {
      PRINT("Plane\n") ;

    } else {
      PRINT("No symmetry\n") ;
    }
  }

  /* Mesh
   * ---- */
  if(DataSet_GetMesh(dataset) && (!strncmp(mot,"mesh",4) || !strncmp(mot,"all",3))) {
    int i ;
    int c1 = 14 ;
    int c2 = 30 ;
    int c3 = 45 ;
    
    PRINT("\n") ;
    PRINT("Mesh:\n") ;
    
    PRINT("\t Nodes:\n") ;
    PRINT("\t Nb of nodes = %d\n",N_NO) ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_t* node_i = NO + i ;
      int ne = Node_GetNbOfElements(node_i) ;
      int n = PRINT("\t no(%d)",i) ;
      int j ;
      
      while(n < c1) n += PRINT(" ") ;
      
      n += PRINT(":") ;
      
      for(j = 0 ; j < DIM ; j++) {
        n += PRINT(" % e",Node_GetCoordinate(node_i)[j]) ;
      }
      
      while(n < c3) n += PRINT(" ") ;
      
      if(ne) n += PRINT("  el(") ;
      
      for(j = 0 ; j < ne ; j++) {
        n += PRINT("%d",Element_GetElementIndex(Node_GetElement(node_i,j))) ;
        n += PRINT(((j < ne - 1) ? "," : ")")) ;
      }
      
      PRINT("\n") ;
    }
    
    PRINT("\n") ;
    PRINT("\t Elements:\n") ;
    PRINT("\t Nb of elements = %d\n",N_EL) ;
    
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_t* elt_i = EL + i ;
      int nn = Element_GetNbOfNodes(elt_i) ;
      int n = PRINT("\t el(%d)",i) ;
      int j ;
      
      while(n < c1) n += PRINT(" ") ;
      
      n += PRINT(":") ;
      
      n += PRINT("  reg(%s)",Element_GetRegionName(elt_i)) ;
      
      while(n < c2) n += PRINT(" ") ;
      
      n += PRINT("  mat(%d)",Element_GetMaterialIndex(elt_i)) ;
      
      while(n < c3) n += PRINT(" ") ;
      
      n += PRINT("  no(") ;
      
      for(j = 0 ; j < nn ; j++) {
        n += PRINT("%d",Node_GetNodeIndex(Element_GetNode(elt_i,j))) ;
        n += PRINT(((j < nn - 1) ? "," : ")")) ;
      }
      
      
      
      PRINT("\n") ;
    }
  }

  /* Models
   * --------- */
  if(DataSet_GetModels(dataset) && (!strncmp(mot,"model",5) || !strncmp(mot,"all",3))) {
    int c2 = 40 ;
    int i ;
    
    PRINT("\n") ;
    
    for(i = 0 ; i < (int) N_MODELS ; i++) {
      int nb_eqn = Model_GetNbOfEquations(MODEL + i) ;
      char* codename = Model_GetCodeNameOfModel(MODEL + i) ;
      char** name_eqn = Model_GetNameOfEquation(MODEL + i) ;
      char** name_unk = Model_GetNameOfUnknown(MODEL + i) ;
      int j ;
      
      PRINT("Model(%d) = %s\n",i,codename) ;
      
      PRINT("\n") ;
      
      PRINT("\t Equations:\n") ;
      PRINT("\t Nb of equations = %d\n",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        int n = PRINT("\t equation(%d): (%s)",j + 1,name_eqn[j]) ;
      
        while(n < c2) n += PRINT(" ") ;
        
        n += PRINT("unknown(%d): (%s)",j + 1,name_unk[j]) ;
        
        PRINT("\n") ;
      }
    }
  }

  /* Materials
   * --------- */
  if(DataSet_GetMaterials(dataset) && (!strncmp(mot,"material",3) || !strncmp(mot,"all",3))) {
    int i ;
    int c2 = 40 ;
    
    PRINT("\n") ;
    
    for(i = 0 ; i < (int) N_MAT ; i++) {
      int nb_pr = Material_GetNbOfProperties(MAT + i) ;
      int nb_eqn = Material_GetNbOfEquations(MAT + i) ;
      char** name_eqn = Material_GetNameOfEquation(MAT + i) ;
      char** name_unk = Material_GetNameOfUnknown(MAT + i) ;
      int nb_cv = Material_GetNbOfCurves(MAT + i) ;
      Curve_t* cv = Material_GetCurve(MAT + i) ;
      int j ;
      
      PRINT("Material(%d):\n",i) ;
      
      PRINT("\t Model = %s\n",Material_GetCodeNameOfModel(MAT + i)) ;
      
      PRINT("\n") ;
      
      PRINT("\t Equations:\n") ;
      PRINT("\t Nb of equations = %d\n",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        int n = PRINT("\t equation(%d): (%s)",j + 1,name_eqn[j]) ;
      
        while(n < c2) n += PRINT(" ") ;
        
        n += PRINT("unknown(%d): (%s)",j + 1,name_unk[j]) ;
        
        PRINT("\n") ;
      }
      
      PRINT("\n") ;
      
      PRINT("\t Properties:\n") ;
      PRINT("\t Nb of properties = %d\n",nb_pr) ;
      
      for(j = 0 ; j < nb_pr ; j++) {
        PRINT("\t prop(%d) = %e\n",j,Material_GetProperty(MAT + i)[j]) ;
      }
      
      PRINT("\n") ;
      
      PRINT("\t Curves:\n") ;
      PRINT("\t Nb of curves = %d\n",nb_cv) ;
      
      for(j = 0 ; j < nb_cv ; j++) {
        PRINT("\t curve(%d): np = %d\n",j + 1,Curve_GetNbOfPoints(cv + j)) ;
      }
    }
    
      
    PRINT("\n") ;
    {
      Models_t* usedmodels = Materials_GetUsedModels(MATERIALS) ;
      int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
      
      PRINT("Nb of used models = %d\n",n_usedmodels) ;
    
      for(i = 0 ; i < n_usedmodels ; i++) {
        Model_t* usedmodel = Models_GetModel(usedmodels) + i ;
      
        PRINT("\t Used model(%d): %s\n",i,Model_GetCodeNameOfModel(usedmodel)) ;
      }
    }
  }

  /* Continuity
   * ---------- */
  if(DataSet_GetMesh(dataset) && (!strncmp(mot,"continuity",3))) {
    int i ;
    
    PRINT("\n") ;
    PRINT("Continuity:\n") ;
    
    PRINT("\t Positions of unknowns and equations at nodes of elements\n") ;
    
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_t* elt_i = EL + i ;
      int nn = Element_GetNbOfNodes(elt_i) ;
      int neq = Element_GetNbOfEquations(elt_i) ;
      char** name_unk = Element_GetNameOfUnknown(elt_i) ;
      char** name_eqn = Element_GetNameOfEquation(elt_i) ;
      int j ;
      
      PRINT("\t el(%d): %d nodes\n",i,nn) ;
      
      PRINT("\t    %d unknowns\n",neq) ;
      
      for(j = 0 ; j < nn ; j++) {
        int k ;
        
        PRINT("\t    no(%d):",j) ;
        
        for(k = 0 ; k < neq ; k++) {
          PRINT(" %s(%d)",name_unk[k],Element_GetUnknownPosition(elt_i)[j*neq + k]) ;
        }
        
        PRINT("\n") ;
      }
      
      PRINT("\t    %d equations\n",neq) ;
      
      for(j = 0 ; j < nn ; j++) {
        int k ;
        
        PRINT("\t    no(%d):",j) ;
        
        for(k = 0 ; k < neq ; k++) {
          PRINT(" %s(%d)",name_eqn[k],Element_GetEquationPosition(elt_i)[j*neq + k]) ;
        }
        
        PRINT("\n") ;
      }
    }
    
    PRINT("\n") ;
    PRINT("\t Equations and unknowns at nodes:\n") ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_t* node_i = NO + i ;
      int nb_unk = Node_GetNbOfUnknowns(node_i) ;
      int nb_eqn = Node_GetNbOfEquations(node_i) ;
      int j ;
      
      PRINT("\t no(%d):\n",i) ;
      PRINT("\t    %d unknowns:",nb_unk) ;
      
      for(j = 0 ; j < nb_unk ; j++) {
        char* name = Node_GetNameOfUnknown(node_i)[j] ;
        
        PRINT(" %s",name) ;
      }
      
      PRINT("\n") ;
      PRINT("\t    %d equations:",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        char* name = Node_GetNameOfEquation(node_i)[j] ;
        
        PRINT(" %s",name) ;
      }
      PRINT("\n") ;
    }
  }

  /* Matrix numbering
   * ---------------- */
  if(DataSet_GetMesh(dataset) && !strncmp(mot,"numbering",3)) {
    int i ;
    
    PRINT("\n") ;
    PRINT("Matrix numbering:\n") ;
    
    PRINT("\n") ;
    PRINT("\t Matrix indexes of equations and unknowns at nodes:\n") ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_t* node_i = NO + i ;
      int nb_unk = Node_GetNbOfUnknowns(node_i) ;
      int nb_eqn = Node_GetNbOfEquations(node_i) ;
      int j ;
      
      PRINT("\t node(%d):\n",i) ;
      PRINT("\t    %d unknowns(col):",nb_unk) ;
      
      for(j = 0 ; j < nb_unk ; j++) {
        char* name = Node_GetNameOfUnknown(node_i)[j] ;
        int icol = Node_GetMatrixColumnIndex(node_i)[j] ;
        
        PRINT(" %s(%d)",name,icol) ;
      }
      
      PRINT("\n") ;
      PRINT("\t    %d equations(row):",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        char* name = Node_GetNameOfEquation(node_i)[j] ;
        int irow = Node_GetMatrixRowIndex(node_i)[j] ;
        
        PRINT(" %s(%d)",name,irow) ;
      }
      PRINT("\n") ;
    }
  }

  /* Functions
   * --------- */
  if(DataSet_GetFunctions(dataset) && (!strncmp(mot,"function",4) || !strncmp(mot,"all",3))) {
    int i ;
    
    PRINT("\n") ;
    
    for(i = 0 ; i < (int) N_FN ; i++) {
      int nb_pts = Function_GetNbOfPoints(FN + i) ;
      double* t  = Function_GetXValue(FN + i) ;
      double* f  = Function_GetFValue(FN + i) ;
      int j ;
      
      PRINT("Time Function(%d):\n",i) ;
      PRINT("\t Nb of points = %d\n",nb_pts) ;
      
      for(j = 0 ; j < nb_pts ; j++) {
        PRINT("\t F(%e) = %e\n",t[j],f[j]) ;
      }
    }
  }

  /* Fields
   * ------ */
  if(DataSet_GetFields(dataset) && (!strncmp(mot,"field",4) || !strncmp(mot,"all",3))) {
    int i ;
    
    PRINT("\n") ;
    
    for(i = 0 ; i < (int) N_CH ; i++) {
      char* type = Field_GetType(CH + i) ;
      
      PRINT("Field(%d):\n",i) ;
      PRINT("\t Type: %s\n",type) ;

      if(!strcmp(type,"affine")) {
        FieldAffine_t* affine =  (FieldAffine_t*) Field_GetFieldFormat(CH + i) ;
        double g[3] = {0.,0.,0.} ;
        double x[3] = {0.,0.,0.} ;
        int j ;
        
        PRINT("\t Value    = %e\n",FieldAffine_GetValue(affine)) ;
        
        PRINT("\t Gradient = ") ;
        
        for(j = 0 ; j < DIM ; j++) {
          g[j] = FieldAffine_GetGradient(affine)[j] ;
        }
        
        PRINT("(%e,%e,%e)",g[0],g[1],g[2]) ;
        
        PRINT("\n") ;
        
        PRINT("\t Point    = ") ;
        
        for(j = 0 ; j < DIM ; j++) x[j] = FieldAffine_GetCoordinate(affine)[j] ;
        
        PRINT("(%e,%e,%e)",x[0],x[1],x[2]) ;
        
        PRINT("\n") ;
        
      } else if(!strncmp(type,"grid",3)) {
        FieldGrid_t* grille =  (FieldGrid_t*) Field_GetFieldFormat(CH + i) ;
        int    n_x = FieldGrid_GetNbOfPointsAlongX(grille) ;
        int    n_y = FieldGrid_GetNbOfPointsAlongY(grille) ;
        int    n_z = FieldGrid_GetNbOfPointsAlongZ(grille) ;
        double* x  = FieldGrid_GetCoordinateAlongX(grille) ;
        double* y  = FieldGrid_GetCoordinateAlongY(grille) ;
        double* z  = FieldGrid_GetCoordinateAlongZ(grille) ;
        int u ;
        
        for(u = 0 ; u < n_x ; u++) {
          int j ;
          
          for(j = 0 ; j < n_y ; j++) {
            int k ;
            
            for(k = 0 ; k < n_z ; k++) {
              double* v = FieldGrid_GetValue(grille) ;
              
              PRINT("\t v(%e,%e,%e) = %e\n",x[u],y[j],z[k],v[(u) + (j)*n_x + (k)*n_x*n_y]) ;
            }
          }
        }

      } else if(!strcmp(type,"constant")) {
        FieldConstant_t* cst =  (FieldConstant_t*) Field_GetFieldFormat(CH + i) ;
        
        PRINT("\t Value    = %e\n",FieldConstant_GetValue(cst)) ;
        PRINT("\t RandomRange = %e\n",FieldConstant_GetRandomRangeLength(cst)) ;
        
        PRINT("\n") ;
        
      } else {
        arret("DataSet_PrintData: type de champ non connu") ;
      }
    }
  }

  /* Initial conditions
   * ------------------ */
  if(DataSet_GetIConds(dataset) && (!strncmp(mot,"initialization",3) || !strncmp(mot,"all",3))) {
    int i ;
    int c1 = 14 ;
    
    PRINT("\n") ;
    
    if(N_IC < 0) {
      char* nom = IConds_GetFileNameOfNodalValues(ICONDS) ;
      FILE*  fic_ini = fopen(nom,"r") ;
      
      if(!fic_ini) {
        arret("DataSet_PrintData: can't open file") ;
      }
    
      PRINT("\n") ;
      PRINT("Initialization of nodal unknowns from %s:\n",nom) ;
    
      for(i = 0 ; i < (int) N_NO ; i++) {
        Node_t* node_i = NO + i ;
        int neq = Node_GetNbOfEquations(node_i) ;
        int n = PRINT("\t no(%d)",i) ;
        int j ;
      
        while(n < c1) n += PRINT(" ") ;
      
        n += PRINT(":") ;
      
        for(j = 0 ; j < neq ; j++) {
          double u ;
          
          fscanf(fic_ini,"%le",&u) ;
          n += PRINT(" %d(%e)",j,u) ;
        }
      
        PRINT("\n") ;
        
      }
      
      fclose(fic_ini) ;
    }
    
    
    for(i = 0 ; i < (int) N_IC ; i++) {
      //int reg = ICond_GetRegionTag(IC + i) ;
      char* reg = ICond_GetRegionName(IC + i) ;
      char* name_unk =ICond_GetNameOfUnknown(IC + i) ;
      Field_t* ch = ICond_GetField(IC + i) ;
      Function_t* fn = ICond_GetFunction(IC + i) ;
      
      PRINT("Initial Condition(%d):\n",i) ;
      
      //PRINT("\t Region  = %d\n",reg) ;
      PRINT("\t Region  = %s\n",reg) ;
      
      PRINT("\t Unknown = %s\n",name_unk) ;
      
      if(ch) {
        int n = ch - CH ;
        
        PRINT("\t Field = %d (type %s)\n",n,Field_GetType(ch)) ;
        
      } else {
        PRINT("\t Natural initial condition (null)\n") ;
        
      }
      
      if(fn) {
        int n = fn - FN ;
        
        PRINT("\t Function = %d\n",n) ;
        
      } else {
        PRINT("\t Function unity (f(t) = 1)\n") ;
        
      }
    }
  }

  /* Boundary conditions
   * ------------------- */
  if(DataSet_GetBConds(dataset) && (!strncmp(mot,"bcondition",4) || !strncmp(mot,"all",3))) {
    int i ;
    
    PRINT("\n") ;
    
    for(i = 0 ; i < (int) N_CL ; i++) {
      //int reg = BCond_GetRegionTag(CL + i) ;
      char* reg = BCond_GetRegionName(CL + i) ;
      char* name_unk =BCond_GetNameOfUnknown(CL + i) ;
      int ich = BCond_GetFieldIndex(CL + i) ;
      int ifn = BCond_GetFunctionIndex(CL + i) ;
      
      PRINT("Boundary Condition(%d):\n",i) ;
      
      //PRINT("\t Region  = %d\n",reg) ;
      PRINT("\t Region  = %s\n",reg) ;
      
      PRINT("\t Unknown = %s\n",name_unk) ;
      
      if(ich >= 0 && ich < N_CH) {
        Field_t* ch = BCond_GetField(CL + i) ;
        
        PRINT("\t Field = %d (type %s)\n",ich,Field_GetType(ch)) ;
        
      } else {
        PRINT("\t Natural boundary condition (null)\n") ;
        
      }
      
      if(ifn >= 0 && ifn < N_FN) {
        
        PRINT("\t Function = %d\n",ifn) ;
        
      } else {
        PRINT("\t Function unity (f(t) = 1)\n") ;
        
      }
    }
  }

  /* Loads
   * ----- */
  if(DataSet_GetLoads(dataset) && (!strncmp(mot,"load",4) || !strncmp(mot,"all",3))) {
    int i ;
    
    PRINT("\n") ;
    
    for(i = 0 ; i < (int) N_CG ; i++) {
      //int reg = Load_GetRegionTag(CG + i) ;
      char* reg = Load_GetRegionName(CG + i) ;
      char* name_eqn = Load_GetNameOfEquation(CG + i) ;
      char* type = Load_GetType(CG + i) ;
      Field_t* ch = Load_GetField(CG + i) ;
      Function_t* fn = Load_GetFunction(CG + i) ;
      
      PRINT("Load(%d):\n",i) ;
      
      //PRINT("\t Region   = %d\n",reg) ;
      PRINT("\t Region   = %s\n",reg) ;
      PRINT("\t Equation = %s\n",name_eqn) ;
      PRINT("\t Type     = %s\n",type) ;
      
      if(ch) {
        int n = ch - CH ;
        
        PRINT("\t Field = %d (type %s)\n",n,Field_GetType(ch)) ;
        
      } else {
        PRINT("\t Natural load (null)\n") ;
        
      }
      
      if(fn) {
        int n = fn - FN ;
        
        PRINT("\t Function = %d\n",n) ;
        
      } else {
        PRINT("\t Function unity (f(t) = 1)\n") ;
        
      }
    }
  }

  /* Points
   * ------ */
  if(DataSet_GetPoints(dataset) && (!strncmp(mot,"points",4) || !strncmp(mot,"all",3))) {
    int n_points = N_POINTS ;
    Point_t* point = Points_GetPoint(POINTS) ;
    int i ;
    
    PRINT("\n") ;
    PRINT("Points:\n") ;
    
    PRINT("\t Nb of points = %d\n",n_points) ;
    
    for(i = 0 ; i < n_points ; i++) {
      double* coor = Point_GetCoordinate(point + i) ;
      double x = (DIM > 0) ? coor[0] : 0. ;
      double y = (DIM > 1) ? coor[1] : 0. ;
      double z = (DIM > 2) ? coor[2] : 0. ;
      Element_t* elt = Point_GetEnclosingElement(point + i) ;
      char* reg_el = (elt) ? Element_GetRegionName(elt) : NULL ;
      int index  = (elt) ? Element_GetElementIndex(elt) : -1 ;
      
      PRINT("\t Point(%d): ",i) ;
      
      PRINT("(%e,%e,%e)",x,y,z) ;
      
      if(reg_el) {
        PRINT(" in region %s",reg_el) ;
      }
      
      if(index >= 0) {
        PRINT(" in element %d",index) ;
      }
      
      PRINT("\n") ;
    }
  }

  /* Dates
   * ----- */
  if(DataSet_GetDates(dataset) && (!strncmp(mot,"dates",4) || !strncmp(mot,"all",3))) {
    int n_dates = N_DATES ;
    Date_t* date = Dates_GetDate(DATES) ;
    int i ;
    
    PRINT("\n") ;
    PRINT("Dates:\n") ;
    
    PRINT("\t Nb of dates = %d\n",n_dates) ;
    
    for(i = 0 ; i < n_dates ; i++) {
      double t = Date_GetTime(date + i) ;
      
      PRINT("\t Date(%d): ",i) ;
      
      PRINT("%e\n",t) ;
    }
  }

  /* Time steps
   * ---------- */
  if(DataSet_GetTimeStep(dataset) && (!strncmp(mot,"time",4) || !strncmp(mot,"all",3))) {
    PRINT("\n") ;
    PRINT("Time Step:\n") ;
    PRINT("\t Dtini = %e\n",TimeStep_GetInitialTimeStep(TIMESTEP)) ;
    PRINT("\t Dtmax = %e\n",TimeStep_GetMaximumTimeStep(TIMESTEP)) ;
    PRINT("\t Dtmin = %e\n",TimeStep_GetMinimumTimeStep(TIMESTEP)) ;
    PRINT("\t Max common ratio = %e\n",TimeStep_GetMaximumCommonRatio(TIMESTEP)) ;
    PRINT("\t Reduction factor = %e\n",TimeStep_GetReductionFactor(TIMESTEP)) ;
  }



  /* Iterative process
   * ----------------- */
  if(DataSet_GetIterProcess(dataset) && (!strncmp(mot,"iterations",4) || !strncmp(mot,"all",3))) {
    PRINT("\n") ;
    PRINT("Iterative Process:\n") ;
    PRINT("\t Nb of iterations = %d\n",IterProcess_GetNbOfIterations(ITERPROCESS)) ;
    PRINT("\t Tolerance = %e\n",IterProcess_GetTolerance(ITERPROCESS)) ;
    PRINT("\t Nb of repetitions = %d\n",IterProcess_GetNbOfRepetitions(ITERPROCESS)) ;
  }

  /* Objective variations
   * -------------------- */
  if(DataSet_GetObVals(dataset) && (!strncmp(mot,"obvariations",4) || !strncmp(mot,"all",3))) {
    int i ;
    
    PRINT("\n") ;
    PRINT("Objective values:\n") ;
    
    PRINT("\t Nb of objective values = %d\n",N_OBJ) ;
    
    for(i = 0 ; i < (int) N_OBJ ; i++) {
      PRINT("\t %s = %e",ObVal_GetNameOfUnknown(OBJ + i),ObVal_GetValue(OBJ + i)) ;
      PRINT(" , type = %c",ObVal_GetType(OBJ + i)) ;
      PRINT(" , relaxation factor = %e",ObVal_GetRelaxationFactor(OBJ + i)) ;
      PRINT("\n") ;
    }
  }

  /* Interpolation functions
   * ----------------------- */
  if(DataSet_GetMesh(dataset) && (!strncmp(mot,"interpolation",4))) {
    int i ;
    
    PRINT("\n") ;
    PRINT("Interpolation:\n") ;
    
    PRINT("\t Nb of interpolation functions = %d\n",N_FI) ;
    
    for(i = 0 ; i < (int) N_FI ; i++) {
      int np = IntFct_GetNbOfPoints(FI + i) ;
      int nn = IntFct_GetNbOfFunctions(FI + i) ;
      int dim = IntFct_GetDimension(FI + i) ;
      
      PRINT("\n") ;
      
      PRINT("\t Interpolation function %d\n",i) ;
      
      PRINT("\t Nb of integration points = %d",np) ;
      
      PRINT(", Dimension = %d\n",dim) ;
      
      if(np <= 0) continue ;
      
      
      PRINT("\t Point Coordinates:\n") ;
      
      {
        char axis[3] = {'x','y','z'} ;
        int k ;
        
        for(k = 0 ; k < dim ; k++) {
          int p ;
      
          PRINT("\t %c = ",axis[k]) ;
        
          for(p = 0 ; p < np ; p++) {
            double* ap = IntFct_GetCoordinatesAtPoint(FI + i,p) ;
            
            PRINT("% e ",ap[k]) ;
          }
        
          PRINT("\n") ;
        }
      }
      
      
      PRINT("\t Weights = ") ;
      
      {
        double* w = IntFct_GetWeight(FI + i) ;
        int p ;
        
        for(p = 0 ; p < np ; p++) {
          PRINT("%e ",w[p]) ;
        }
        
        PRINT("\n") ;
      }
      
      
      PRINT("\t Nb of functions = %d\n",nn) ;
      
      
      PRINT("\t Functions:\n") ;
      
      {
        int k ;
        
        for(k = 0 ; k < nn ; k++) {
          int p ;
          
          if(k == 0) {
            
            PRINT("\t hi = ") ;
            
            for(p = 0 ; p < np ; p++) {
              PRINT(" hi(pt %d)     ",p) ;
            }
            
            PRINT("\n") ;
          }
        
          PRINT("\t h%d = ",k) ;
        
          for(p = 0 ; p < np ; p++) {
            double* hp = IntFct_GetFunctionAtPoint(FI + i,p) ;
            
            PRINT("% e ",hp[k]) ;
          }
        
          PRINT("\n") ;
        }
      }
      
      
      
      PRINT("\t Function derivatives:\n") ;
      
      
#define DHP(n,i)  (dhp[(n)*dim+(i)])
      {
        int l ;
        
        for(l = 0 ; l < nn ; l++) {
          int k ;
          char axis[3] = {'x','y','z'} ;
          
          if(l == 0) {
            int p ;
            
            PRINT("\t hi,j = ") ;
            
            for(p = 0 ; p < np ; p++) {
              PRINT(" hi,j(pt %d)   ",p) ;
            }
            
            PRINT("\n") ;
          }
        
          for(k = 0 ; k < dim ; k++) {
            int p ;
        
            PRINT("\t h%d,%c = ",l,axis[k]) ;
          
            for(p = 0 ; p < np ; p++) {
              double* dhp = IntFct_GetFunctionGradientAtPoint(FI + i,p) ;
              
              PRINT("% e ",DHP(l,k)) ;
            }
        
            PRINT("\n") ;
          }
        }
      }
#undef DHJ
    }
  }

  fflush(stdout) ;
}
