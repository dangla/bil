#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "DataSet.h"
#include "Modules.h"

#define TITLE "Mechanics from a microstructure (2017)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (dim)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (12)
#define NVE     (0)
#define NV0     (9)

/* Equation index */
#define E_mec   (0)

/* Unknown index */
#define U_u     (0)

/* We define some names for implicit terms */
#define SIG           (vim + 0)
#define F_MASS        (vim + 9)

#define SIG_n         (vim_n + 0)

/* We define some names for explicit terms */

/* We define some names for constant terms */
#define SIG0          (v0  + 0)


/* Functions */
static Model_ComputePropertyIndex_t  pm ;
static void   GetProperties(Element_t*) ;

static int    ComputeTangentCoefficients(Element_t*,double,double,double*) ;
static int    ComputeMicrostructureMatrix(Mesh_t*,Solver_t*,double,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,Solutions_t*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,Solutions_t*,Solutions_t*) ;
//static double* ComputeVariableDerivatives(Element_t*,double,double,double*,Solutions_t*,Solutions_t*,double,int) ;

static void  ComputeMicrostructure(DataSet_t*,Solver_t*,double,double,Solutions_t*,Solutions_t*,double*,double*) ;
static void  InitializeMicrostructureDataSet(DataSet_t*) ;


static double* MacroGradient(Element_t*,double) ;
static double* MacroStrain(Element_t*,double) ;

/* Parameters */
//static double  gravity ;
//static double  rho_s ;
static double  macrogradient[9] ;
static double  macrostrain[9] ;


#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

#define ItIsPeriodic  (Geometry_IsPeriodic(Element_GetGeometry(el)))


#define NbOfVariables     (42)
//static double  Variable0[NbOfVariables] ;
//static double dVariable[NbOfVariables] ;
//static double  Variable2[NbOfVariables] ;
//static double*  Variable ;

#define I_U            (0)

#define I_EPS          (3)

#define I_SIG          (12)
#define I_EPS_n        (21)
#define I_Fmass        (30)
#define I_SIG_n        (33)



int pm(const char* s)
{
         if(!strcmp(s,"gravity")) {
    return(0) ;
  } else if(!strcmp(s,"rho_s")) {
    return(1) ;
  } else if(!strcmp(s,"macro-gradient")) {
    return(2) ;
  } else if(!strncmp(s,"macro-gradient_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(2 + 3*i + j) ;
  } else if(!strcmp(s,"macro-fctindex")) {
    return(11) ;
  } else if(!strncmp(s,"macro-fctindex_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(11 + 3*i + j) ;
  } else if(!strcmp(s,"charlen")) {
    return(20) ;
  } else return(-1) ;
}



double* MacroGradient(Element_t* el,double t)
{
  double* gradient = macrogradient ;
  double  f[9] = {0,0,0,0,0,0,0,0,0} ;
  
  {
    Functions_t* fcts = Material_GetFunctions(Element_GetMaterial(el)) ;
    Function_t*  fct = Functions_GetFunction(fcts) ;
    int nf = Functions_GetNbOfFunctions(fcts) ;
    double* fctindex = &GetProperty("macro-fctindex") ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      int idx = floor(fctindex[i] + 0.5) ;
      
      if(0 < idx && idx < nf + 1) {
        Function_t* macrogradfct = fct + idx - 1 ;
      
        f[i] = Function_ComputeValue(macrogradfct,t) ;
      }
    }
  }
  
  {
    double* g = &GetProperty("macro-gradient") ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      gradient[i] = g[i] * f[i] ;
    }
  }
  
  return(gradient) ;
}



double* MacroStrain(Element_t* el,double t)
{
  double* strain = macrostrain ;
  double* grd = MacroGradient(el,t) ;
  int i ;
    
  for(i = 0 ; i < 3 ; i++) {
    int j ;
      
    for(j = 0 ; j < 3 ; j++) {
      strain[3*i + j] = 0.5*(grd[3*i + j] + grd[3*j + i]) ;
    }
  }
  
  return(strain) ;
}



void GetProperties(Element_t* el)
{
  //gravity = GetProperty("gravity") ;
  //rho_s   = GetProperty("rho_s") ;
}



int SetModelProp(Model_t* model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  for(i = 0 ; i < dim ; i++) {
    char name_eqn[7] ;
    sprintf(name_eqn,"meca_%d",i + 1) ;
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn) ;
  }
  
  /** Names of the main unknowns */
  for(i = 0 ; i < dim ; i++) {
    char name_unk[4] ;
    sprintf(name_unk,"u_%d",i + 1) ;
    Model_CopyNameOfUnknown(model,U_u + i,name_unk) ;
  }
  
  Model_GetComputePropertyIndex(model) = pm ;
  
  Model_GetNbOfVariables(model) = NbOfVariables ;
  //Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int NbOfProp = 21 ;

  /* Par defaut tout a 0 */
  {
    int i ;
    
    for(i = 0 ; i < NbOfProp ; i++) {
      Material_GetProperty(mat)[i] = 0. ;
    }
  }
  
  /* Pre-initialization */
  {
    //Material_GetProperty(mat)[pm("charlen")] = 1 ;
  }
  

  Material_ScanProperties(mat,datafile,pm) ;


  /* The dataset and solver for the microstructure */
  {
    char* method = Material_GetMethod(mat) ;
    
    if(!strncmp(method,"Microstructure",14)) {
      char* p = strstr(method," ") ;
      char* cellname = p + strspn(p," ") ;
      Options_t* options = Options_Create(NULL) ;
      DataSet_t* jdd = DataSet_Create(cellname,options) ;
            
      InitializeMicrostructureDataSet(jdd) ;
      
      /* Store jdd in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,jdd,DataSet_t,"DataSet") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }

      
      /* The solver */
      {
        Mesh_t* mesh = DataSet_GetMesh(jdd) ;
        Solver_t* solver = Solver_Create(mesh,options,6) ;
        GenericData_t* gdat = GenericData_Create(1,solver,Solver_t,"Solver") ;
      
        Material_AppendGenericData(mat,gdat) ;
        
      }
    }
  }

  return(NbOfProp) ;
}


int PrintModelProp(Model_t* model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;

  printf("\n") ;
  printf("The system consists in dim equations\n") ;
  printf("\t The equilibrium equations (meca_1,meca_2,meca_3)\n") ;

  printf("\n") ;
  printf("The primary unknowns are\n") ;
  printf("\t The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n") ;
  printf("Example of input data\n\n") ;
  

  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # masse volumique du squelette sec\n") ;
  
  return(0) ;
}


//#define ElementSol_GetSolutions(ES) \
        ((Solutions_t*) GenericData_GetData(GenericData_GetNextGenericData(ElementSol_GetImplicitGenericData(ES))))
#define ElementSol_GetSolutions(ES) \
        ((Solutions_t*) GenericData_FindData(ElementSol_GetImplicitGenericData(ES),Solutions_t,"Solutions"))


/* The current and the previous solutions */
#define Element_GetSolutions(EL) \
        Element_GetCurrentSolutions(EL)

#define Element_GetCurrentSolutions(EL) \
        ElementSol_GetSolutions(Element_GetElementSol(EL))
        
#define Element_GetPreviousSolutions(EL) \
        ElementSol_GetSolutions(Element_GetPreviousElementSol(EL))
        
        

int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  

  /** The solutions from the microstructures */
  {
    /* Store "sols" for the whole history */
    {
      ElementSol_t* elementsol = Element_GetElementSol(el) ;

      if(elementsol) {
        DataSet_t* jdd = (DataSet_t*) Element_FindMaterialData(el,DataSet_t,"DataSet") ;
        Mesh_t* mesh = DataSet_GetMesh(jdd) ;
        const int nsol_micro = 2 ;
        
        do {
          Solutions_t* sols = Solutions_Create(NbOfIntPoints,mesh,nsol_micro) ;
  
          if(!sols) assert(sols) ;
    
          /* Merging explicit terms (not essential!) */
          {
            int i ;
    
            for(i = 0 ; i < NbOfIntPoints ; i++) {
              Solutions_MergeExplicitTerms(sols + i) ;
            }
          }
          
          {
            GenericData_t* gdat  = GenericData_Create(NbOfIntPoints,sols,Solutions_t,"Solutions") ;
          
            ElementSol_AddImplicitGenericData(elementsol,gdat) ;
          }
        
          elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
        } while(elementsol != Element_GetElementSol(el)) ;
      }
    }
  }
  
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
  
    {
      int i ;
      
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}



int ComputeInitialState(Element_t* el,double t)
{
  double* vim0  = Element_GetImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  Solutions_t* sols = Element_GetSolutions(el) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;
    
  {
    DataSet_t* jdd = (DataSet_t*) Element_FindMaterialData(el,DataSet_t,"DataSet") ;
    DataFile_t* df = DataSet_GetDataFile(jdd) ;
    DataFile_t* datafile = Element_GetDataFile(el) ;
    
    if(DataFile_ContextIsPartialInitialization(datafile)) {
      DataFile_ContextSetToPartialInitialization(df) ;
    } else {
      DataFile_ContextSetToFullInitialization(df) ;
    }
  }
  
  /* If there are initial displacements */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,sols,t,0,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
    
    }
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{
  double* vim0   = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  Solutions_t* sols_n = Element_GetPreviousSolutions(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;
    
  {
    DataSet_t* jdd = (DataSet_t*) Element_FindMaterialData(el,DataSet_t,"DataSet") ;
    DataFile_t* df = DataSet_GetDataFile(jdd) ;
    
    DataFile_ContextSetToPartialInitialization(df) ;
  }
  
  
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim_n,sols_n,t,dt,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
    
    }
  }
  
  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;


  /* Initialization */
  {
    double zero = 0. ;
    int    i ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;
  }


  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /*
    Input data
  */
  GetProperties(el) ;


  /*
  ** Tangent matrix
  */
  if(ndof) {
    double c[IntFct_MaxNbOfIntPoints*81] ;
    int dec = ComputeTangentCoefficients(el,t,dt,c) ;
    FEM_t* fem = FEM_GetInstance(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    double* kp = FEM_ComputeElasticMatrix(fem,intfct,c,dec) ;
    
    {
      int i ;
      
      for(i = 0 ; i < ndof*ndof ; i++) {
        k[i] = kp[i] ;
      }
    }
  }
  
  return(0) ;
#undef K
}




int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* vim = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialization */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;


  /* 1. Mechanics */
  
  /* 1.1 Stresses */
  {
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_mec + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  {
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_mec + dim - 1) -= -rbf[i] ;
    }
    
  }
  
  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 4 ;
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t* fem = FEM_GetInstance(el) ;

  //if(Element_IsSubmanifold(el)) return(0) ;

  /* Initialization */
  {
    int    i ;
    
    for(i = 0 ; i < NbOfOutputs ; i++) {
      Result_SetValuesToZero(r + i) ;
    }
  }
  
  /*
    Input data
  */
  GetProperties(el) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Displacement */
    double pdis[3] = {0,0,0} ;
    double dis[3] = {0,0,0} ;
    /* strains */
    double *eps ;
    /* stresses */
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    
    /* Displacements */
    {
      int    i ;
    
      for(i = 0 ; i < dim ; i++) {
        pdis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
        dis[i] = pdis[i] ;
      }

      if(ItIsPeriodic) {
      
        for(i = 0 ; i < 3 ; i++) {
          int j ;
        
          for(j = 0 ; j < 3 ; j++) {
            dis[i] += MacroGradient(el,t)[3*i + j] * s[j] ;
          }
        }
      }
    }
    
    
    /* Strains */
    {
      int    i ;
      
      eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
      
      if(ItIsPeriodic) {
      
        for(i = 0 ; i < 9 ; i++) {
          eps[i] += MacroStrain(el,t)[i] ;
        }
      }
    }
    
    /* Averaging */
    for(p = 0 ; p < np ; p++) {
      double* vim  = vim0 + p*NVI ;
      int j ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
    }
      
    {
      int i = 0 ;
      
      Result_Store(r + i++,dis   ,"Displacements",3) ;
      Result_Store(r + i++,sig   ,"Stresses",9) ;
      Result_Store(r + i++,pdis  ,"Perturbated-displacements",3) ;
      Result_Store(r + i++,eps   ,"Strains",9) ;
      
      if(i != NbOfOutputs) arret("ComputeOutputs") ;
    }
  }
  
  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(Element_t* el,double t,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
  DataSet_t* jdd = (DataSet_t*) Element_FindMaterialData(el,DataSet_t,"DataSet") ;
  DataFile_t* df = DataSet_GetDataFile(jdd) ;
  Solver_t* solver = (Solver_t*) Element_FindMaterialData(el,Solver_t,"Solver") ;
  Solutions_t* sols   = Element_GetCurrentSolutions(el) ;
  Mesh_t* mesh = DataSet_GetMesh(jdd) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 81 ;

  {
    int    p ;
    
    for(p = 0 ; p < np ; p++) {
      double* c0 = c + p*dec ;
    
      Mesh_InitializeSolutionPointers(mesh,sols + p) ;
      
      {
        Session_Open() ;
        Message_SetVerbosity(0) ;
        
        Message_Direct("\n") ;
        Message_Direct("Start a calculation of microstructure matrix: %s",DataFile_GetFileName(df)) ;
        Message_Direct("\n") ;
        
        ComputeMicrostructureMatrix(mesh,solver,t,dt,c0) ;
        
        Session_Close() ;
      }
      
#if 0
    {
      int i ;
      
      printf("\n") ;
      printf("4th rank stiffness tensor:\n") ;
      
      for(i = 0 ; i < 9 ; i++) {
        int j = i - (i/3)*3 ;
        
        printf("C%d%d--:",i/3 + 1,j + 1) ;
        
        for (j = 0 ; j < 9 ; j++) {
          printf(" % e",c0[i*9 + j]) ;
        }
        
        printf("\n") ;
      }
    }
#endif
    }
  }
  
  return(dec) ;
}




int ComputeMicrostructureMatrix(Mesh_t* mesh,Solver_t* solver,double t,double dt,double* c)
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  int dim = Mesh_GetDimension(mesh) ;
  Matrix_t* a = Solver_GetMatrix(solver) ;
  double*   b = Solver_GetRHS(solver) ;
  double*   u = Solver_GetSolution(solver) ;
  int ncol = Solver_GetNbOfColumns(solver) ;
  double*  pb[9] = {b,b+ncol,b+2*ncol,b+ncol,b+3*ncol,b+4*ncol,b+2*ncol,b+4*ncol,b+5*ncol} ;
  double*  pu[9] = {u,u+ncol,u+2*ncol,u+ncol,u+3*ncol,u+4*ncol,u+2*ncol,u+4*ncol,u+5*ncol} ;
  double    E[9] = {1,0.5,0.5,0.5,1,0.5,0.5,0.5,1} ;


  /* Initializations */
  {
    int j ;
    
    for(j = 0 ; j < 81 ; j++) c[j] = 0 ;
  }
  
    
  Matrix_SetValuesToZero(a) ;
    
  {
    int i ;
      
    for(i = 0 ; i < 6*ncol ; i++) {
      b[i] = 0. ;
    }
  }
  
  
  /* The matrix and the r.h.s. */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
#define NE (Element_MaxNbOfNodes*Model_MaxNbOfEquations)
    double ke[NE*NE] ;
#undef NE
    int ie ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Material_t* mat = Element_GetMaterial(el + ie) ;
    
      if(mat) {
      
        Element_FreeBuffer(el + ie) ;
        {
          int i = Element_ComputeMatrix(el + ie,t,dt,ke) ;
        
          if(i != 0) return(i) ;
        }
      
        Matrix_AssembleElementMatrix(a,el+ie,ke) ;
      
        {
          int  nn = Element_GetNbOfNodes(el + ie) ;
          int neq = Element_GetNbOfEquations(el + ie) ;
          int ndof = nn*neq ;
          int n ;
                
          for(n = 0 ; n < nn ; n++) {
            Node_t* node_n = Element_GetNode(el + ie,n) ;
            double* x_n = Node_GetCoordinate(node_n) ;
            int m ;
                
            for(m = 0 ; m < nn ; m++) {
              Node_t* node_m = Element_GetNode(el + ie,m) ;
              double* x_m = Node_GetCoordinate(node_m) ;
              int i ;
        
              for(i = 0 ; i < dim ; i++) {
                int ni = n*neq + i ;
                int jj_row = Element_GetEquationPosition(el + ie)[ni] ;
                
                /*  The r.h.s. stored in pb */
                if(jj_row >= 0) {
                  int row_i = Node_GetMatrixRowIndex(node_n)[jj_row] ;
                  int j ;
          
                  for(j = 0 ; j < dim ; j++) {
                    int mj = m*neq + j ;
                    int ij = ni * ndof + mj ;
                    int k ;
                      
                    for(k = j ; k < dim ; k++) {
                      int      jk = 3 * j + k ;
                      double* bjk = pb[jk] ;

                      bjk[row_i] -= ke[ij] * E[jk] * x_m[k] ;
                    }
                  }
                }
                
                /*  First part of C */
                {
                  int j ;
          
                  for(j = 0 ; j < dim ; j++) {
                    int mj = m*neq + j ;
                    int ij = ni * ndof + mj ;
                    int k ;
                      
                    for(k = 0 ; k < dim ; k++) {
                      int l ;
          
                      for(l = 0 ; l < dim ; l++) {
                        C(k,i,l,j) += x_n[k] * ke[ij] * x_m[l] ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  /* The solutions */
  {
    int k ;

    for(k = 0 ; k < dim ; k++) {
      int l ;
          
      for(l = k ; l < dim ; l++) {
        Solver_GetRHS(solver) = pb[3 * k + l] ;
        Solver_GetSolution(solver) = pu[3 * k + l] ;
        Solver_Solve(solver) ;
      }
    }
    
    Solver_GetRHS(solver) = b ;
    Solver_GetSolution(solver) = u ;
  }


  /* The macroscopic tangent stiffness matrix */
  {
    int i ;
        
    for(i = 0 ; i < dim ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        int k ;
          
        for(k = 0 ; k < dim ; k++) {
          int l ;
          
          for(l = 0 ; l < dim ; l++) {
            double  ub = 0. ;
            
            {
              double* uik = pu[3 * i + k] ;
              double* blj = pb[3 * l + j] ;
              int p ;
              
              for(p = 0 ; p < ncol ; p++) {
                ub += uik[p] * blj[p] ;
              }
            }
            
            {
              double  Eik =  E[3 * i + k] ;
              double  Elj =  E[3 * l + j] ;

              C(k,i,j,l) -= ub / (Eik * Elj) ;
            }
          }
        }
      }
    }
  }
  
  {
    double vol = FEM_ComputeVolume(mesh) ;
    int j ;
    
    for(j = 0 ; j < 81 ; j++) c[j] /= vol ;
  }
  
  return(0) ;
#undef C
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,Solutions_t* sols_n,double t,double dt,int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double*   x      = Model_GetVariable(model,p) ;
  //double*   x      = Variable ;
  
    
  /* Primary Variables */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }
  }
    
  /* Strains */
  {
    double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
    int    i ;
      
    if(ItIsPeriodic) {
      
      for(i = 0 ; i < 9 ; i++) {
        eps[i] += MacroStrain(el,t)[i] ;
      }
    }
    
    for(i = 0 ; i < 9 ; i++) {
      x[I_EPS + i] = eps[i] ;
    }
      
    FEM_FreeBufferFrom(fem,eps) ;
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    int    i ;
    
    /* Stresses, strains at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_u) ;
      double* vim_n = f_n + p*NVI ;
      
      if(ItIsPeriodic) {
        
        for(i = 0 ; i < 9 ; i++) {
          eps_n[i] += MacroStrain(el,t-dt)[i] ;
        }
      }
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS_n   + i] = eps_n[i] ;
        x[I_SIG_n   + i] = SIG_n[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
  }
  
  {
    Solutions_t* sols = Element_GetCurrentSolutions(el) ;
    
    ComputeSecondaryVariables(el,t,dt,x,sols_n + p,sols + p) ;
  }
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x,Solutions_t* sols_n,Solutions_t* sols)
{
  /* Strains */
  double* eps   =  x + I_EPS ;
  double* eps_n =  x + I_EPS_n ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = x + I_SIG ;
    double* sig_n = x + I_SIG_n ;
    
    
    {
      double  deps[9] ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Elastic trial stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
      
      /* Compute here the microstructure */
      {
        DataSet_t* jdd = (DataSet_t*) Element_FindMaterialData(el,DataSet_t,"DataSet") ;
        Solver_t* solver = (Solver_t*) Element_FindMaterialData(el,Solver_t,"Solver") ;
        DataFile_t* df = DataSet_GetDataFile(jdd) ;

        
        Session_Open() ;
        Message_SetVerbosity(0) ;
        
        Message_Direct("\n") ;
        Message_Direct("Start a calculation of microstructure: %s",DataFile_GetFileName(df)) ;
        Message_Direct("\n") ;
        
        ComputeMicrostructure(jdd,solver,t,dt,sols_n,sols,deps,sig) ;
        
        Session_Close() ;
      }

    }
  }
  
  
  /* Backup body force */
  {
    {
      double gravity = GetProperty("gravity") ;
      double rho_s   = GetProperty("rho_s") ;
      int dim = Element_GetDimensionOfSpace(el) ;
      double* f_mass = x + I_Fmass ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s)*gravity ;
    }
  }
}



void ComputeMicrostructure(DataSet_t* jdd,Solver_t* solver,double t,double dt,Solutions_t* sols_n,Solutions_t* sols,double* macrograd,double* sig)
{
  /* Set input data of the microstructure */
  {
    /* Update the time function */
    {
      Functions_t* fcts = DataSet_GetFunctions(jdd) ;
      int nfcts = 9 ;
      int j ;
      
      for(j = 0 ; j < nfcts ; j++) {
        Function_t* func = Functions_GetFunction(fcts) + j ;
        double* x = Function_GetXValue(func) ;
        double* f = Function_GetFValue(func) ;

        x[0] = t - dt ;
        x[1] = t ;
        f[0] = 0 ;
        f[1] = macrograd[j] ;
      }
    }
  }
    
  /* Compute the microstructure */
  {
    Modules_t* modules = DataSet_GetModules(jdd) ;
    Module_t* module_i = Modules_FindModule(modules,"Module1") ;
    Dates_t*   dates  = DataSet_GetDates(jdd) ;
    Date_t*    date   = Dates_GetDate(dates) ;
    TimeStep_t*  timestep  = DataSet_GetTimeStep(jdd) ;
    double dtini = TimeStep_GetInitialTimeStep(timestep) ;
    double t_n = Solutions_GetTime(sols_n) ;
    
    {
      /* t should be equal to t_n + dt. */
      if(fabs(t - dt - t_n) > 1.e-4*dtini) {
        printf("coucou") ;
        Message_FatalError("ComputeMicrostructure: t_n = %e ; t - dt = %e",t_n,t-dt) ;
      }
      
      /* There are 2 dates */
      Date_GetTime(date) = t_n ;
      Date_GetTime(date + 1) = t ;
    
      {
        Solution_t* sol   = Solutions_GetSolution(sols) ;
        Solution_t* sol_n = Solutions_GetSolution(sols_n) ;
        
        Solution_Copy(sol,sol_n) ;
      }
      
      {
        int i ;
        
        i = Module_SolveProblem(module_i,jdd,sols,solver,NULL) ;
        
        if(i < 0) {
          Message_Warning("ComputeMicrostructure: something went wrong") ;
          Exception_Interrupt ;
        }
      }
    }
  }

  /* Backup stresses as averaged stresses */
  {
    Mesh_t* mesh = DataSet_GetMesh(jdd) ;
    
    Mesh_InitializeSolutionPointers(mesh,sols) ;
    FEM_AverageStresses(mesh,sig) ;
  }
}






void InitializeMicrostructureDataSet(DataSet_t* jdd)
{
  /* Set input data of the microstructure */
  {
    /* Update the macro-gradient and the macro-fctindex */
    {
      Materials_t* mats = DataSet_GetMaterials(jdd) ;
      int nmats = Materials_GetNbOfMaterials(mats) ;
      int j ;
    
      for(j = 0 ; j < nmats ; j++) {
        Material_t* mat = Materials_GetMaterial(mats) + j ;
        Model_t* model = Material_GetModel(mat) ;
        Model_ComputePropertyIndex_t* pidx = Model_GetComputePropertyIndex(model) ;
        
        if(!pidx) {
          arret("InitializeMicrostructureDataSet(1): Model_GetComputePropertyIndex(model) undefined") ;
        }
        
        {
          double* grd = Material_GetProperty(mat) + pidx("macro-gradient") ;
          double* fid = Material_GetProperty(mat) + pidx("macro-fctindex") ;
          int i ;
    
          for(i = 0 ; i < 9 ; i++) {
            grd[i] = 1 ;
            fid[i] = i + 1 ;
          }
        }
      }
    }
    
    
    /* Check and update the function of time */
    {
      Functions_t* fcts = DataSet_GetFunctions(jdd) ;
      int nfcts = Functions_GetNbOfFunctions(fcts) ;
      int j ;
      
      if(nfcts < 9) {
        arret("InitializeMicrostructureDataSet(1): 9 functions are needed") ;
      }
      
      nfcts = 9 ;
      for(j = 0 ; j < nfcts ; j++) {
        Function_t* func = Functions_GetFunction(fcts) + j ;
        int npts = Function_GetNbOfPoints(func) ;
          
        if(npts < 2) {
          arret("InitializeMicrostructureDataSet(2): 2 points are needed") ;
        }
          
        {
          double* t = Function_GetXValue(func) ;
          double* f = Function_GetFValue(func) ;
            
          Function_GetNbOfPoints(func) = 2 ;
          t[0] = 0 ;
          t[1] = 1 ;
          f[0] = 0 ;
          f[1] = 1 ;
        }
      }
    }
    
    /* The dates */
    {
      {
        Dates_t* dates = DataSet_GetDates(jdd) ;
        int     nbofdates  = Dates_GetNbOfDates(dates) ;
          
        if(nbofdates < 2) {
          arret("InitializeMicrostructureDataSet(3): 2 dates are needed") ;
        }
      
        Dates_GetNbOfDates(dates) = 2 ;
      }
    }
  }
}



#if 0
double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,Solutions_t* sols_n,Solutions_t* sols,double dxi,int i)
{
  double* dx = dVariable ;
  double* x2 = Variable2 ;
  
  /* Primary Variables */
  {
    int j ;
    
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] = x[j] ;
      x2[j] = x[j] ;
    }
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  x2[i] -= dxi ;
  
  
  if(i >= I_EPS && i < I_EPS + 9) {
    /* Find indexes k and l: kl = 3*k + l */
    int kl = i - I_EPS ;
    int k = kl%3 ;
    int l = (kl - k)/3 ;
    int lk = 3*l + k ;
    int j = lk + I_EPS ;
    
    dx[j] += dxi ;
    x2[j] -= dxi ;
  }
  
  ComputeSecondaryVariables(el,t,dt,dx,sols_n,sols) ;
  ComputeSecondaryVariables(el,t,dt,x2,sols_n,sols) ;
  
  /* The numerical derivative as (f(x + dx) - f(x - dx))/(2dx) */
  {
    int j ;
    
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] -= x2[j] ;
      dx[j] /= (2*dxi) ;
    }
  }
  
  if(i >= I_EPS && i < I_EPS + 9) {
    int j ;
    
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] *= 0.5 ;
    }
  }

  return(dx) ;
}
#endif



#if 0
int ComputeTangentCoefficients1(FEM_t* fem,double t,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define C2(i,j)        (c1[(i)*9 +(j)])
  Element_t* el  = FEM_GetElement(fem) ;
  double deps[9] ;
  int    dec = 81 ;
  
  GetProperties(el) ;
  
  {
    int dim = Element_GetDimensionOfSpace(el) ;
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    double dxi[Model_MaxNbOfEquations] ;
    double chardis = 0 ;
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dxi[i] = 1.e-2*ObVal_GetValue(obval + i) ;
    }
    
    for(i = 0 ; i < dim ; i++) {
      double d = 1.e-2*ObVal_GetValue(obval + i) ;
      
      if(chardis < d) chardis = d ;
    }
    
    for(i = 0 ; i < 9 ; i++) {
      deps[i] = chardis/charlen ;
    }
  }
  
  /*
  Message_Direct("\n") ;
  Message_Direct("Start a stiffness matrix calculation") ;
  Message_Direct("\n") ;
  */


  {
    double* vim_n  = Element_GetPreviousImplicitTerm(el) ;
    double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
    double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
    Solutions_t* sols   = Element_GetCurrentSolutions(el) ;
    Solutions_t* sols_n = Element_GetPreviousSolutions(el) ;
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
  
    for(p = 0 ; p < np ; p++) {
      double* c0 = c + p*dec ;
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vim_n,sols_n,t,dt,p) ;
    
      /* Initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0 ;
      }
    

      /* Mechanics */
      {
        double* c1 = c0 ;
        int i ;
        
        for(i = 0 ; i < 81 ; i++) {
          c1[i] = 0 ;
        }
      
        /* Tangent stiffness matrix */
        {
          int k ;
        
          for(k = 0 ; k < 3 ; k++) {
            int l ;
            
            for(l = k ; l < 3 ; l++) {
              int kl = 3*k + l ;
              int lk = 3*l + k ;
              int m = I_EPS + kl ;
              double dxm = deps[kl] ;
              double* dx = ComputeVariableDerivatives(el,t,dt,x,sols_n + p,sols + p,dxm,m) ;
              int j ;
              
              for(j = 0 ; j < 9 ; j++) {
                C2(j,kl) = dx[I_SIG + (j)] ;
                C2(j,lk) = dx[I_SIG + (j)] ;
              }
            }
          }
        }
      }
    }
  }
  
  return(dec) ;
#undef C2
#undef C1
#undef T4
}
#endif
