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
#include "DataSet.h"
#include "Units.h"
#include "Parser.h"


static void   (DataSet_PrintData)(DataSet_t*,char*) ;



/* Extern functions */

DataSet_t*  (DataSet_Create)(char* filename,Options_t* opt)
{
  DataSet_t* jdd = DataSet_New() ;
  char*   debug  = Options_GetPrintedInfos(opt) ;
  
  //assert(jdd) ;
  
  DataSet_GetOptions(jdd) = opt ;
  
  {
    DataFile_t* datafile = DataFile_Create(filename) ;
  
    DataSet_GetDataFile(jdd) = datafile ;
  
    if(DataFile_DoesNotExist(datafile)) {
      Help_WriteData(filename) ;
      Message_Info("To start the computation, type bil %s\n",filename) ;
      exit(EXIT_SUCCESS) ;
    }
  }

  Message_Direct("Reading %s\n",filename) ;


  /* Units */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    
    Units_Create(datafile) ;
  }
  
  
  /* Geometry */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
  
    DataSet_GetGeometry(jdd) = Geometry_Create(datafile) ;
  }
  
  
  /* Materials */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Geometry_t*    geometry = DataSet_GetGeometry(jdd) ;
  
    DataSet_GetMaterials(jdd) = Materials_Create(datafile,geometry) ;
  }
  if(!strcmp(debug,"mate")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Mesh */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Geometry_t*    geometry = DataSet_GetGeometry(jdd) ;
    Materials_t*   materials = DataSet_GetMaterials(jdd) ;
  
    DataSet_GetMesh(jdd) = Mesh_Create(datafile,materials,geometry) ;
  }
  if(!strcmp(debug,"geom")) DataSet_PrintData(jdd,debug) ;
  if(!strcmp(debug,"mesh")) DataSet_PrintData(jdd,debug) ;
  if(!strcmp(debug,"continuity")) DataSet_PrintData(jdd,debug) ;
  if(!strcmp(debug,"inter")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Fields */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Geometry_t*    geometry = DataSet_GetGeometry(jdd) ;
    Materials_t*   materials = DataSet_GetMaterials(jdd) ;
  
    DataSet_GetFields(jdd) = Fields_Create(datafile,materials,geometry) ;
  }
  if(!strcmp(debug,"field")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Functions */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Materials_t*   materials = DataSet_GetMaterials(jdd) ;
  
    DataSet_GetFunctions(jdd) = Functions_Create(datafile,materials) ;
  }
  if(!strcmp(debug,"func")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Initial conditions */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Fields_t*      fields = DataSet_GetFields(jdd) ;
    Functions_t*   functions = DataSet_GetFunctions(jdd) ;
  
    DataSet_GetIConds(jdd) = IConds_Create(datafile,fields,functions) ;
  }
  if(!strcmp(debug,"init")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Loads */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Fields_t*      fields = DataSet_GetFields(jdd) ;
    Functions_t*   functions = DataSet_GetFunctions(jdd) ;
  
    DataSet_GetLoads(jdd) = Loads_Create(datafile,fields,functions) ;
  }
  if(!strcmp(debug,"load")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Boundary conditions */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Fields_t*      fields = DataSet_GetFields(jdd) ;
    Functions_t*   functions = DataSet_GetFunctions(jdd) ;
  
    DataSet_GetBConds(jdd) = BConds_Create(datafile,fields,functions) ;
  }
  if(!strcmp(debug,"bcond")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Matrix permutation numbering */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    BConds_t*      bconds = DataSet_GetBConds(jdd) ;
    Mesh_t*        mesh = DataSet_GetMesh(jdd) ;
  
    Mesh_SetMatrixPermutationNumbering(mesh,bconds,datafile) ;
  }
  if(!strcmp(debug,"numbering")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Points */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Mesh_t*        mesh = DataSet_GetMesh(jdd) ;
  
    DataSet_GetPoints(jdd) = Points_Create(datafile,mesh) ;
  }
  if(!strcmp(debug,"points")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Dates */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    
    DataSet_GetDates(jdd) = Dates_Create(datafile) ;
  }
  if(!strcmp(debug,"dates")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Objective variations */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    Materials_t*   materials = DataSet_GetMaterials(jdd) ;
    Mesh_t*        mesh = DataSet_GetMesh(jdd) ;
  
    DataSet_GetObVals(jdd) = ObVals_Create(datafile,mesh,materials) ;
  }
  if(!strcmp(debug,"obval")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Time steps */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    ObVals_t*      obvals = DataSet_GetObVals(jdd) ;
    
    DataSet_GetTimeStep(jdd) = TimeStep_Create(datafile,obvals) ;
  }
  if(!strcmp(debug,"time")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Iterative process */
  {
    DataFile_t*    datafile = DataSet_GetDataFile(jdd) ;
    ObVals_t*      obvals = DataSet_GetObVals(jdd) ;
    
    DataSet_GetIterProcess(jdd) = IterProcess_Create(datafile,obvals) ;
  }
  if(!strcmp(debug,"iter")) DataSet_PrintData(jdd,debug) ;
  
  
  /* Modules */
  DataSet_GetModules(jdd) = Modules_Create() ;
  if(!strcmp(debug,"module")) DataSet_PrintData(jdd,debug) ;
  
  
  if(!strcmp(debug,"all")) DataSet_PrintData(jdd,debug) ;


  Message_Direct("End of reading %s\n",filename) ;
  Message_Direct("\n") ;
  
  return(jdd) ;
}



DataSet_t*  (DataSet_Create1)(char* filename,Options_t* opt)
{
  DataSet_t* jdd = DataSet_New() ;
  
  DataSet_GetOptions(jdd) = opt ;
  
  {
    DataFile_t* datafile = DataFile_Create(filename) ;
  
    DataSet_GetDataFile(jdd) = datafile ;
  
    if(DataFile_DoesNotExist(datafile)) {
      Help_WriteData(filename) ;
      Message_Info("To start the computation, type bil %s\n",filename) ;
      exit(EXIT_SUCCESS) ;
    }
  }
  
  Parser_ParseFile(jdd) ;
  
  {
    char*   debug  = Options_GetPrintedInfos(opt) ;
    
    DataSet_PrintData(jdd,debug) ;
  }
  
  return(jdd) ;
}



/* Local functions */


#define DATASET       (jdd)
#define DATAFILE      DataSet_GetDataFile(DATASET)
#define MESH          DataSet_GetMesh(DATASET)
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

#define DIM           Mesh_GetDimension(MESH)
#define SYMMETRY      Mesh_GetSymmetry(MESH) 
#define COORSYS       Mesh_GetCoordinateSystem(MESH)

#define N_NO          Mesh_GetNbOfNodes(MESH)
#define NO            Mesh_GetNode(MESH)

#define N_MAT         (Materials_GetNbOfMaterials(MATERIALS))
#define MAT           (Materials_GetMaterial(MATERIALS))

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


void DataSet_PrintData(DataSet_t* jdd,char* mot)
{
  static int i_debug=0 ;
  
  if(!strcmp(mot,"\0")) return ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"debug(%d)\n",i_debug++) ;
  fprintf(stdout,"-----\n") ;

  /* Geometry
   * -------- */
  if(DataSet_GetGeometry(jdd) && (!strcmp(mot,"geom") || !strcmp(mot,"all"))) {
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Geometry:\n") ;
    
    fprintf(stdout,"\t Dimension = %dD\n",DIM) ;
    fprintf(stdout,"\t Symmetry = ") ;
    
    if(0) {
      
    } else if(Symmetry_IsCylindrical(SYMMETRY)) {
      fprintf(stdout,"Axisymmetrical\n") ;
      
    } else if(Symmetry_IsSpherical(SYMMETRY)) {
      fprintf(stdout,"Spherical\n") ;

    } else if(Symmetry_IsPlane(SYMMETRY)) {
      fprintf(stdout,"Plane\n") ;

    } else {
      fprintf(stdout,"No symmetry\n") ;
    }
  }

  /* Mesh
   * ---- */
  if(DataSet_GetMesh(jdd) && (!strcmp(mot,"mesh") || !strcmp(mot,"all"))) {
    int i ;
    int c1 = 14 ;
    int c2 = 30 ;
    int c3 = 45 ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Mesh:\n") ;
    
    fprintf(stdout,"\t Nodes:\n") ;
    fprintf(stdout,"\t Nb of nodes = %d\n",N_NO) ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_t* node_i = NO + i ;
      int n = fprintf(stdout,"\t no(%d)",i) ;
      int j ;
      
      while(n < c1) n += fprintf(stdout," ") ;
      
      n += fprintf(stdout,":") ;
      
      for(j = 0 ; j < DIM ; j++) {
        n += fprintf(stdout," % e",Node_GetCoordinate(node_i)[j]) ;
      }
      
      fprintf(stdout,"\n") ;
    }
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"\t Elements:\n") ;
    fprintf(stdout,"\t Nb of elements = %d\n",N_EL) ;
    
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_t* elt_i = EL + i ;
      int nn = Element_GetNbOfNodes(elt_i) ;
      int n = fprintf(stdout,"\t el(%d)",i) ;
      int j ;
      
      while(n < c1) n += fprintf(stdout," ") ;
      
      n += fprintf(stdout,":") ;
      
      n += fprintf(stdout,"  reg(%d)",Element_GetRegionIndex(elt_i)) ;
      
      while(n < c2) n += fprintf(stdout," ") ;
      
      n += fprintf(stdout,"  mat(%d)",Element_GetMaterialIndex(elt_i)) ;
      
      while(n < c3) n += fprintf(stdout," ") ;
      
      n += fprintf(stdout,"  no(") ;
      
      for(j = 0 ; j < nn ; j++) {
        n += fprintf(stdout,"%d",Node_GetNodeIndex(Element_GetNode(elt_i,j))) ;
        n += fprintf(stdout,((j < nn - 1) ? "," : ")")) ;
      }
      
      
      
      fprintf(stdout,"\n") ;
    }
  }

  /* Materials
   * --------- */
  if(DataSet_GetMaterials(jdd) && (!strcmp(mot,"mate") || !strcmp(mot,"all"))) {
    int i ;
    int c2 = 40 ;
    
    fprintf(stdout,"\n") ;
    
    for(i = 0 ; i < (int) N_MAT ; i++) {
      int nb_pr = Material_GetNbOfProperties(MAT + i) ;
      int nb_eqn = Material_GetNbOfEquations(MAT + i) ;
      char** name_eqn = Material_GetNameOfEquation(MAT + i) ;
      char** name_unk = Material_GetNameOfUnknown(MAT + i) ;
      int nb_cv = Material_GetNbOfCurves(MAT + i) ;
      Curve_t* cv = Material_GetCurve(MAT + i) ;
      int j ;
      
      fprintf(stdout,"Material(%d):\n",i) ;
      
      fprintf(stdout,"\t Model = %s\n",Material_GetCodeNameOfModel(MAT + i)) ;
      
      fprintf(stdout,"\n") ;
      
      fprintf(stdout,"\t Equations:\n") ;
      fprintf(stdout,"\t Nb of equations = %d\n",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        int n = fprintf(stdout,"\t equation(%d): (%s)",j + 1,name_eqn[j]) ;
      
        while(n < c2) n += fprintf(stdout," ") ;
        
        n += fprintf(stdout,"unknown(%d): (%s)",j + 1,name_unk[j]) ;
        
        fprintf(stdout,"\n") ;
      }
      
      fprintf(stdout,"\n") ;
      
      fprintf(stdout,"\t Properties:\n") ;
      fprintf(stdout,"\t Nb of properties = %d\n",nb_pr) ;
      
      for(j = 0 ; j < nb_pr ; j++) {
        fprintf(stdout,"\t prop(%d) = %e\n",j,Material_GetProperty(MAT + i)[j]) ;
      }
      
      fprintf(stdout,"\n") ;
      
      fprintf(stdout,"\t Curves:\n") ;
      fprintf(stdout,"\t Nb of curves = %d\n",nb_cv) ;
      
      for(j = 0 ; j < nb_cv ; j++) {
        fprintf(stdout,"\t curve(%d): np = %d\n",j + 1,Curve_GetNbOfPoints(cv + j)) ;
      }
    }
    
      
    fprintf(stdout,"\n") ;
    {
      Models_t* usedmodels = Materials_GetUsedModels(MATERIALS) ;
      int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
      
      fprintf(stdout,"Nb of used models = %d\n",n_usedmodels) ;
    
      for(i = 0 ; i < n_usedmodels ; i++) {
        Model_t* usedmodel = Models_GetModel(usedmodels) + i ;
      
        fprintf(stdout,"\t Used model(%d): %s\n",i,Model_GetCodeNameOfModel(usedmodel)) ;
      }
    }
  }

  /* Continuity
   * ---------- */
  if(DataSet_GetMesh(jdd) && (!strcmp(mot,"continuity"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Continuity:\n") ;
    
    fprintf(stdout,"\t Positions of unknowns and equations at nodes of elements\n") ;
    
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_t* elt_i = EL + i ;
      int nn = Element_GetNbOfNodes(elt_i) ;
      int neq = Element_GetNbOfEquations(elt_i) ;
      char** name_unk = Element_GetNameOfUnknown(elt_i) ;
      char** name_eqn = Element_GetNameOfEquation(elt_i) ;
      int j ;
      
      fprintf(stdout,"\t el(%d): %d nodes\n",i,nn) ;
      
      fprintf(stdout,"\t    %d unknowns\n",neq) ;
      
      for(j = 0 ; j < nn ; j++) {
        int k ;
        
        fprintf(stdout,"\t    no(%d):",j) ;
        
        for(k = 0 ; k < neq ; k++) {
          fprintf(stdout," %s(%d)",name_unk[k],Element_GetUnknownPosition(elt_i)[j*neq + k]) ;
        }
        
        fprintf(stdout,"\n") ;
      }
      
      fprintf(stdout,"\t    %d equations\n",neq) ;
      
      for(j = 0 ; j < nn ; j++) {
        int k ;
        
        fprintf(stdout,"\t    no(%d):",j) ;
        
        for(k = 0 ; k < neq ; k++) {
          fprintf(stdout," %s(%d)",name_eqn[k],Element_GetEquationPosition(elt_i)[j*neq + k]) ;
        }
        
        fprintf(stdout,"\n") ;
      }
    }
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"\t Equations and unknowns at nodes:\n") ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_t* node_i = NO + i ;
      int nb_unk = Node_GetNbOfUnknowns(node_i) ;
      int nb_eqn = Node_GetNbOfEquations(node_i) ;
      int j ;
      
      fprintf(stdout,"\t no(%d):\n",i) ;
      fprintf(stdout,"\t    %d unknowns:",nb_unk) ;
      
      for(j = 0 ; j < nb_unk ; j++) {
        char* name = Node_GetNameOfUnknown(node_i)[j] ;
        
        fprintf(stdout," %s",name) ;
      }
      
      fprintf(stdout,"\n") ;
      fprintf(stdout,"\t    %d equations:",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        char* name = Node_GetNameOfEquation(node_i)[j] ;
        
        fprintf(stdout," %s",name) ;
      }
      fprintf(stdout,"\n") ;
    }
  }

  /* Matrix numbering
   * ---------------- */
  if(DataSet_GetMesh(jdd) && !strcmp(mot,"numbering")) {
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Matrix numbering:\n") ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"\t Matrix indexes of equations and unknowns at nodes:\n") ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_t* node_i = NO + i ;
      int nb_unk = Node_GetNbOfUnknowns(node_i) ;
      int nb_eqn = Node_GetNbOfEquations(node_i) ;
      int j ;
      
      fprintf(stdout,"\t no(%d):\n",i) ;
      fprintf(stdout,"\t    %d unknowns(col):",nb_unk) ;
      
      for(j = 0 ; j < nb_unk ; j++) {
        char* name = Node_GetNameOfUnknown(node_i)[j] ;
        int icol = Node_GetMatrixColumnIndex(node_i)[j] ;
        
        fprintf(stdout," %s(%d)",name,icol) ;
      }
      
      fprintf(stdout,"\n") ;
      fprintf(stdout,"\t    %d equations(row):",nb_eqn) ;
      
      for(j = 0 ; j < nb_eqn ; j++) {
        char* name = Node_GetNameOfEquation(node_i)[j] ;
        int irow = Node_GetMatrixRowIndex(node_i)[j] ;
        
        fprintf(stdout," %s(%d)",name,irow) ;
      }
      fprintf(stdout,"\n") ;
    }
  }

  /* Functions
   * --------- */
  if(DataSet_GetFunctions(jdd) && (!strcmp(mot,"func") || !strcmp(mot,"all"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    
    for(i = 0 ; i < (int) N_FN ; i++) {
      int nb_pts = Function_GetNbOfPoints(FN + i) ;
      double* t  = Function_GetXValue(FN + i) ;
      double* f  = Function_GetFValue(FN + i) ;
      int j ;
      
      fprintf(stdout,"Function(%d):\n",i) ;
      fprintf(stdout,"\t Nb of points = %d\n",nb_pts) ;
      
      for(j = 0 ; j < nb_pts ; j++) {
        fprintf(stdout,"\t F(%e) = %e\n",t[j],f[j]) ;
      }
    }
  }

  /* Fields
   * ------ */
  if(DataSet_GetFields(jdd) && (!strcmp(mot,"field") || !strcmp(mot,"all"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    
    for(i = 0 ; i < (int) N_CH ; i++) {
      char* type = Field_GetType(CH + i) ;
      
      fprintf(stdout,"Field(%d):\n",i) ;
      fprintf(stdout,"\t Type: %s\n",type) ;

      if(!strcmp(type,"affine")) {
        FieldAffine_t* affine =  (FieldAffine_t*) Field_GetFieldFormat(CH + i) ;
        double g[3] = {0.,0.,0.} ;
        double x[3] = {0.,0.,0.} ;
        int j ;
        
        fprintf(stdout,"\t Value    = %e\n",FieldAffine_GetValue(affine)) ;
        
        fprintf(stdout,"\t Gradient = ") ;
        
        for(j = 0 ; j < DIM ; j++) {
          g[j] = FieldAffine_GetGradient(affine)[j] ;
        }
        
        fprintf(stdout,"(%e,%e,%e)",g[0],g[1],g[2]) ;
        
        fprintf(stdout,"\n") ;
        
        fprintf(stdout,"\t Point    = ") ;
        
        for(j = 0 ; j < DIM ; j++) x[j] = FieldAffine_GetCoordinate(affine)[j] ;
        
        fprintf(stdout,"(%e,%e,%e)",x[0],x[1],x[2]) ;
        
        fprintf(stdout,"\n") ;
        
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
              
              fprintf(stdout,"\t v(%e,%e,%e) = %e\n",x[u],y[j],z[k],v[(u) + (j)*n_x + (k)*n_x*n_y]) ;
            }
          }
        }

      } else if(!strcmp(type,"constant")) {
        FieldConstant_t* cst =  (FieldConstant_t*) Field_GetFieldFormat(CH + i) ;
        
        fprintf(stdout,"\t Value    = %e\n",FieldConstant_GetValue(cst)) ;
        
        fprintf(stdout,"\n") ;
        
      } else {
        arret("DataSet_PrintData: type de champ non connu") ;
      }
    }
  }

  /* Initial conditions
   * ------------------ */
  if(DataSet_GetIConds(jdd) && (!strcmp(mot,"init") || !strcmp(mot,"all"))) {
    int i ;
    int c1 = 14 ;
    
    fprintf(stdout,"\n") ;
    
    if(N_IC < 0) {
      char* nom = IConds_GetFileNameOfNodalValues(ICONDS) ;
      FILE*  fic_ini = fopen(nom,"r") ;
      
      if(!fic_ini) {
        arret("DataSet_PrintData: can't open file") ;
      }
    
      fprintf(stdout,"\n") ;
      fprintf(stdout,"Initialization of nodal unknowns from %s:\n",nom) ;
    
      for(i = 0 ; i < (int) N_NO ; i++) {
        Node_t* node_i = NO + i ;
        int neq = Node_GetNbOfEquations(node_i) ;
        int n = fprintf(stdout,"\t no(%d)",i) ;
        int j ;
      
        while(n < c1) n += fprintf(stdout," ") ;
      
        n += fprintf(stdout,":") ;
      
        for(j = 0 ; j < neq ; j++) {
          double u ;
          
          fscanf(fic_ini,"%le",&u) ;
          n += fprintf(stdout," %d(%e)",j,u) ;
        }
      
        fprintf(stdout,"\n") ;
        
      }
      
      fclose(fic_ini) ;
    }
    
    
    for(i = 0 ; i < (int) N_IC ; i++) {
      int reg = ICond_GetRegionIndex(IC + i) ;
      char* name_unk =ICond_GetNameOfUnknown(IC + i) ;
      Field_t* ch = ICond_GetField(IC + i) ;
      Function_t* fn = ICond_GetFunction(IC + i) ;
      
      fprintf(stdout,"Initial Condition(%d):\n",i) ;
      
      fprintf(stdout,"\t Region  = %d\n",reg) ;
      
      fprintf(stdout,"\t Unknown = %s\n",name_unk) ;
      
      if(ch) {
        int n = ch - CH ;
        
        fprintf(stdout,"\t Field = %d (type %s)\n",n,Field_GetType(ch)) ;
        
      } else {
        fprintf(stdout,"\t Natural initial condition (null)\n") ;
        
      }
      
      if(fn) {
        int n = fn - FN ;
        
        fprintf(stdout,"\t Function = %d\n",n) ;
        
      } else {
        fprintf(stdout,"\t Function unity (f(t) = 1)\n") ;
        
      }
    }
  }

  /* Boundary conditions
   * ------------------- */
  if(DataSet_GetBConds(jdd) && (!strcmp(mot,"bcond") || !strcmp(mot,"all"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    
    for(i = 0 ; i < (int) N_CL ; i++) {
      int reg = BCond_GetRegionIndex(CL + i) ;
      char* name_unk =BCond_GetNameOfUnknown(CL + i) ;
      Field_t* ch = BCond_GetField(CL + i) ;
      Function_t* fn = BCond_GetFunction(CL + i) ;
      
      fprintf(stdout,"Boundary Condition(%d):\n",i) ;
      
      fprintf(stdout,"\t Region  = %d\n",reg) ;
      
      fprintf(stdout,"\t Unknown = %s\n",name_unk) ;
      
      if(ch) {
        int n = ch - CH ;
        
        fprintf(stdout,"\t Field = %d (type %s)\n",n,Field_GetType(ch)) ;
        
      } else {
        fprintf(stdout,"\t Natural boundary condition (null)\n") ;
        
      }
      
      if(fn) {
        int n = fn - FN ;
        
        fprintf(stdout,"\t Function = %d\n",n) ;
        
      } else {
        fprintf(stdout,"\t Function unity (f(t) = 1)\n") ;
        
      }
    }
  }

  /* Loads
   * ----- */
  if(DataSet_GetLoads(jdd) && (!strcmp(mot,"load") || !strcmp(mot,"all"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    
    for(i = 0 ; i < (int) N_CG ; i++) {
      int reg = Load_GetRegionIndex(CG + i) ;
      char* name_eqn = Load_GetNameOfEquation(CG + i) ;
      char* type = Load_GetType(CG + i) ;
      Field_t* ch = Load_GetField(CG + i) ;
      Function_t* fn = Load_GetFunction(CG + i) ;
      
      fprintf(stdout,"Load(%d):\n",i) ;
      
      fprintf(stdout,"\t Region   = %d\n",reg) ;
      fprintf(stdout,"\t Equation = %s\n",name_eqn) ;
      fprintf(stdout,"\t Type     = %s\n",type) ;
      
      if(ch) {
        int n = ch - CH ;
        
        fprintf(stdout,"\t Field = %d (type %s)\n",n,Field_GetType(ch)) ;
        
      } else {
        fprintf(stdout,"\t Natural load (null)\n") ;
        
      }
      
      if(fn) {
        int n = fn - FN ;
        
        fprintf(stdout,"\t Function = %d\n",n) ;
        
      } else {
        fprintf(stdout,"\t Function unity (f(t) = 1)\n") ;
        
      }
    }
  }

  /* Points
   * ------ */
  if(DataSet_GetPoints(jdd) && (!strcmp(mot,"points") || !strcmp(mot,"all"))) {
    int n_points = N_POINTS ;
    Point_t* point = Points_GetPoint(POINTS) ;
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Points:\n") ;
    
    fprintf(stdout,"\t Nb of points = %d\n",n_points) ;
    
    for(i = 0 ; i < n_points ; i++) {
      double* coor = Point_GetCoordinate(point + i) ;
      double x = (DIM > 0) ? coor[0] : 0. ;
      double y = (DIM > 1) ? coor[1] : 0. ;
      double z = (DIM > 2) ? coor[2] : 0. ;
      Element_t* elt = Point_GetEnclosingElement(point + i) ;
      int reg = Element_GetRegionIndex(elt) ;
      
      fprintf(stdout,"\t Point(%d): ",i) ;
      
      fprintf(stdout,"(%e,%e,%e) ",x,y,z) ;
      
      fprintf(stdout,"in region %d\n",reg) ;
    }
  }

  /* Dates
   * ----- */
  if(DataSet_GetDates(jdd) && (!strcmp(mot,"dates") || !strcmp(mot,"all"))) {
    int n_dates = N_DATES ;
    Date_t* date = Dates_GetDate(DATES) ;
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Dates:\n") ;
    
    fprintf(stdout,"\t Nb of dates = %d\n",n_dates) ;
    
    for(i = 0 ; i < n_dates ; i++) {
      double t = Date_GetTime(date + i) ;
      
      fprintf(stdout,"\t Date(%d): ",i) ;
      
      fprintf(stdout,"%e\n",t) ;
    }
  }

  /* Time steps
   * ---------- */
  if(DataSet_GetTimeStep(jdd) && (!strcmp(mot,"time") || !strcmp(mot,"all"))) {
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Time Step:\n") ;
    fprintf(stdout,"\t Dtini = %e\n",TimeStep_GetInitialTimeStep(TIMESTEP)) ;
    fprintf(stdout,"\t Dtmax = %e\n",TimeStep_GetMaximumTimeStep(TIMESTEP)) ;
    fprintf(stdout,"\t Dtmin = %e\n",TimeStep_GetMinimumTimeStep(TIMESTEP)) ;
    fprintf(stdout,"\t Max common ratio = %e\n",TimeStep_GetMaximumCommonRatio(TIMESTEP)) ;
    fprintf(stdout,"\t Reduction factor = %e\n",TimeStep_GetReductionFactor(TIMESTEP)) ;
  }



  /* Iterative process
   * ----------------- */
  if(DataSet_GetIterProcess(jdd) && (!strcmp(mot,"iter") || !strcmp(mot,"all"))) {
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Iterative Process:\n") ;
    fprintf(stdout,"\t Nb of iterations = %d\n",IterProcess_GetNbOfIterations(ITERPROCESS)) ;
    fprintf(stdout,"\t Tolerance = %e\n",IterProcess_GetTolerance(ITERPROCESS)) ;
    fprintf(stdout,"\t Nb of repetitions = %d\n",IterProcess_GetNbOfRepetitions(ITERPROCESS)) ;
  }

  /* Objective variations
   * -------------------- */
  if(DataSet_GetObVals(jdd) && (!strcmp(mot,"obval") || !strcmp(mot,"all"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Objective values:\n") ;
    
    fprintf(stdout,"\t Nb of objective values = %d\n",N_OBJ) ;
    
    for(i = 0 ; i < (int) N_OBJ ; i++) {
      fprintf(stdout,"\t %s = %e",ObVal_GetNameOfUnknown(OBJ + i),ObVal_GetValue(OBJ + i)) ;
      fprintf(stdout," , type = %c",ObVal_GetType(OBJ + i)) ;
      fprintf(stdout," , relaxation factor = %e",ObVal_GetRelaxationFactor(OBJ + i)) ;
      fprintf(stdout,"\n") ;
    }
  }

  /* Interpolation functions
   * ----------------------- */
  if(DataSet_GetMesh(jdd) && (!strcmp(mot,"inter"))) {
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"Interpolation:\n") ;
    
    fprintf(stdout,"\t Nb of interpolation functions = %d\n",N_FI) ;
    
    for(i = 0 ; i < (int) N_FI ; i++) {
      int np = IntFct_GetNbOfPoints(FI + i) ;
      int nn = IntFct_GetNbOfFunctions(FI + i) ;
      int dim = IntFct_GetDimension(FI + i) ;
      
      fprintf(stdout,"\n") ;
      
      fprintf(stdout,"\t Interpolation function %d\n",i) ;
      
      fprintf(stdout,"\t Nb of integration points = %d",np) ;
      
      fprintf(stdout,", Dimension = %d\n",dim) ;
      
      if(np <= 0) continue ;
      
      
      fprintf(stdout,"\t Point Coordinates:\n") ;
      
      {
        double* a = IntFct_GetPointCoordinates(FI + i) ;
        char axis[3] = {'x','y','z'} ;
        int k ;
        
        for(k = 0 ; k < dim ; k++) {
          int p ;
      
          fprintf(stdout,"\t %c = ",axis[k]) ;
        
          for(p = 0 ; p < np ; p++) {
            double* ap = a + p*dim ;
            
            fprintf(stdout,"% e ",ap[k]) ;
          }
        
          fprintf(stdout,"\n") ;
        }
      }
      
      
      fprintf(stdout,"\t Weights = ") ;
      
      {
        double* w = IntFct_GetWeight(FI + i) ;
        int p ;
        
        for(p = 0 ; p < np ; p++) {
          fprintf(stdout,"%e ",w[p]) ;
        }
        
        fprintf(stdout,"\n") ;
      }
      
      
      fprintf(stdout,"\t Nb of functions = %d\n",nn) ;
      
      
      fprintf(stdout,"\t Functions:\n") ;
      
      {
        double* h = IntFct_GetFunction(FI + i) ;
        int k ;
        
        for(k = 0 ; k < nn ; k++) {
          int p ;
          
          if(k == 0) {
            
            fprintf(stdout,"\t hi = ") ;
            
            for(p = 0 ; p < np ; p++) {
              fprintf(stdout," hi(pt %d)     ",p) ;
            }
            
            fprintf(stdout,"\n") ;
          }
        
          fprintf(stdout,"\t h%d = ",k) ;
        
          for(p = 0 ; p < np ; p++) {
            double* hp = h + p*nn ;
            
            fprintf(stdout,"% e ",hp[k]) ;
          }
        
          fprintf(stdout,"\n") ;
        }
      }
      
      
      
      fprintf(stdout,"\t Function derivatives:\n") ;
      
      
#define DHP(n,i)  (dhp[(n)*dim+(i)])
      {
        double* dh = IntFct_GetFunctionGradient(FI + i) ;
        int l ;
        
        for(l = 0 ; l < nn ; l++) {
          int k ;
          char axis[3] = {'x','y','z'} ;
          
          if(l == 0) {
            int p ;
            
            fprintf(stdout,"\t hi,j = ") ;
            
            for(p = 0 ; p < np ; p++) {
              fprintf(stdout," hi,j(pt %d)   ",p) ;
            }
            
            fprintf(stdout,"\n") ;
          }
        
          for(k = 0 ; k < dim ; k++) {
            int p ;
        
            fprintf(stdout,"\t h%d,%c = ",l,axis[k]) ;
          
            for(p = 0 ; p < np ; p++) {
              double* dhp = dh + p*nn*dim ;
              
              fprintf(stdout,"% e ",DHP(l,k)) ;
            }
        
            fprintf(stdout,"\n") ;
          }
        }
      }
#undef DHJ
    }
  }

  fflush(stdout) ;
}
