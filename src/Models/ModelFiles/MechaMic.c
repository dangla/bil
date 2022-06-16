#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM2.h"
#include "FEM.h"
#include "DataSet.h"
#include "Solvers.h"
#include "Modules.h"

#define TITLE "Mechanics from a microstructure (2017)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (dim)

/* Equation index */
#define E_Mech   (0)

/* Unknown index */
#define U_dis     E_Mech



/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (12)
#define NVE     (0)
#define NV0     (9)

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

static double* ComputeVariables(Element_t*,double**,double**,double*,Solutions_t*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,double*,Solutions_t*,Solutions_t*) ;
//static double* ComputeVariableDerivatives(Element_t*,double,double,double*,Solutions_t*,Solutions_t*,double,int) ;



static double* MacroGradient(Element_t*,double) ;
static double* MacroStrain(Element_t*,double) ;

/* Parameters */
//static double  gravity ;
//static double  rho_s ;
static double  macrogradient[9] ;
static double  macrostrain[9] ;


#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

#define ItIsPeriodic  (Geometry_IsPeriodic(Element_GetGeometry(el)))



/* We define some indices for the local variables */
enum {
I_DIS = 0,
I_DIS2 =  I_DIS + 2,

I_EPS,
I_EPS8   = I_EPS  + 8,

I_SIG,
I_SIG8   = I_SIG  + 8,

I_F_MASS,
I_F_MASS2 = I_F_MASS + 2,
I_Last
} ;

#define NbOfVariables     (I_Last)
static double  Variable[NbOfVariables] ;
static double  Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;



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
  char i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  for(i = 0 ; i < dim ; i++) {
    char name_eqn[7] ;
    sprintf(name_eqn,"meca_%u",i + 1) ;
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn) ;
  }
  
  /** Names of the main unknowns */
  for(i = 0 ; i < dim ; i++) {
    char name_unk[4] ;
    sprintf(name_unk,"u_%u",i + 1) ;
    Model_CopyNameOfUnknown(model,U_dis + i,name_unk) ;
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
      DataSet_t* dataset = DataSet_Create(cellname,options) ;
      Mesh_t* mesh = DataSet_GetMesh(dataset) ;
      Solvers_t* solvers = Solvers_Create(mesh,options,6) ;
      
      /* Store options in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,options,Options_t,"Options") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }
      
      /* Store dataset in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,dataset,DataSet_t,"DataSet") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }

      
      /* Store solvers in mat */
      {
        GenericData_t* gdat = GenericData_Create(1,solvers,Solvers_t,"Solvers") ;
      
        Material_AppendGenericData(mat,gdat) ;
      }

      /* Initialize the dataset */
      {
        FEM2_t* fem2 = FEM2_GetInstance(dataset,NULL,NULL,NULL) ;
        
        FEM2_InitializeMicrostructureDataSet(fem2) ;
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



/* The current and the previous solutions */
#define Element_GetCurrentSolutions(EL) \
        ((Solutions_t*) Element_FindCurrentImplicitData(EL,Solutions_t,"Solutions"))

#define Element_GetPreviousSolutions(EL) \
        ((Solutions_t*) Element_FindPreviousImplicitData(EL,Solutions_t,"Solutions"))


/* The dataset */
#define Element_GetDataSet(EL) \
        ((DataSet_t*) Element_FindMaterialData(EL,DataSet_t,"DataSet"))


/* The solvers */
#define Element_GetSolvers(EL) \
        ((Solvers_t*) Element_FindMaterialData(EL,Solvers_t,"Solvers"))




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
    DataSet_t* dataset = Element_GetDataSet(el) ;
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    const int nsol_micro = 2 ;
    
    /* Store "sols" for the whole history */
    
    Element_AllocateSolutions(el,mesh,nsol_micro) ;
  }
  

  /** The solutions from the microstructures */
  {
    /* Merging explicit terms (not essential!) */
    #if 0
    {
      ElementSol_t* elementsol = Element_GetElementSol(el) ;

      if(elementsol) {
        do {
          Solutions_t* sols = ElementSol_FindImplicitData(elementsol,Solutions_t,"Solutions") ;
          
          {
            int i ;
    
            for(i = 0 ; i < NbOfIntPoints ; i++) {
              Solutions_MergeExplicitTerms(sols + i) ;
            }
          }
        
          elementsol = ElementSol_GetPreviousElementSol(elementsol) ;
        } while(elementsol != Element_GetElementSol(el)) ;
      }
    }
    #endif
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
  Solutions_t* sols = Element_GetCurrentSolutions(el) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;
    
  {
    DataSet_t* dataset = Element_GetDataSet(el) ;
    DataFile_t* df = DataSet_GetDataFile(dataset) ;
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
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_F_MASS + i] ;
    
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
    DataSet_t* dataset = Element_GetDataSet(el) ;
    DataFile_t* df = DataSet_GetDataFile(dataset) ;
    
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
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_F_MASS + i] ;
    
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
      
      for(j = 0 ; j < dim ; j++) R(i,E_Mech + j) -= rw[i*dim + j] ;
    }
    
  }
  
  /* 1.2 Body forces */
  {
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Mech + dim - 1) -= -rbf[i] ;
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
        pdis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_dis + i) ;
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
      
      eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_dis) ;
      
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
  DataSet_t* dataset = Element_GetDataSet(el) ;
  DataFile_t* df = DataSet_GetDataFile(dataset) ;
  Solvers_t* solvers = Element_GetSolvers(el) ;
  Solutions_t* sols   = Element_GetCurrentSolutions(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 81 ;

  {
    int    p ;
    
    for(p = 0 ; p < np ; p++) {
      double* c0 = c + p*dec ;
      
      {
        Session_Open() ;
        Message_SetVerbosity(0) ;
        
        Message_Direct("\n") ;
        Message_Direct("Start a calculation of microstructure matrix: %s",DataFile_GetFileName(df)) ;
        Message_Direct("\n") ;
        
        {
          FEM2_t* fem2 = FEM2_GetInstance(dataset,solvers,NULL,sols+p) ;
          int i = FEM2_HomogenizeTangentStiffnessTensor(fem2,t,dt,c0) ;
        
          if(i < 0) {
            Message_FatalError("ComputeTangentCoefficients: something went wrong") ;
          }
        }
        
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



double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,Solutions_t* sols_n,double t,double dt,int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  //double*   x      = Model_GetVariable(model,p) ;
  double*   x      = Variable ;
  double*   x_n    = Variable_n ;
  
    
  /* Primary Variables */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_DIS + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_dis + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_DIS + i] = 0 ;
    }
  }
    
  /* Strains */
  {
    double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_dis) ;
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
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_dis) ;
      double* vim_n = f_n + p*NVI ;
      
      if(ItIsPeriodic) {
        
        for(i = 0 ; i < 9 ; i++) {
          eps_n[i] += MacroStrain(el,t-dt)[i] ;
        }
      }
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_EPS   + i] = eps_n[i] ;
        x_n[I_SIG   + i] = SIG_n[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
  }
  
  {
    Solutions_t* sols = Element_GetCurrentSolutions(el) ;
    
    ComputeSecondaryVariables(el,t,dt,x_n,x,sols_n + p,sols + p) ;
  }
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x,Solutions_t* sols_n,Solutions_t* sols)
{
  /* Strains */
  double* eps   =  x   + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = x   + I_SIG ;
    double* sig_n = x_n + I_SIG ;
    
    
    {
      double  deps[9] ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Elastic trial stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
      
      /* Compute here the microstructure */
      {
        DataSet_t* dataset = Element_GetDataSet(el) ;
        Solvers_t* solvers = Element_GetSolvers(el) ;
        DataFile_t* df = DataSet_GetDataFile(dataset) ;

        
        Session_Open() ;
        Message_SetVerbosity(4) ;
        
        Message_Direct("\n") ;
        Message_Direct("Start a calculation of microstructure: %s",DataFile_GetFileName(df)) ;
        Message_Direct("\n") ;
        
        {
          FEM2_t* fem2 = FEM2_GetInstance(dataset,solvers,sols_n,sols) ;
          int j = FEM2_ComputeHomogenizedStressTensor(fem2,t,dt,deps,sig) ;
        
          if(j < 0) {
            Message_FatalError("ComputeSecondaryVariables: something went wrong") ;
          }
        }
        
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
      double* f_mass = x + I_F_MASS ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s)*gravity ;
    }
  }
}



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
