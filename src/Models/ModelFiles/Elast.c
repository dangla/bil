#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Common.h"
#include "FEM.h"
#include "DataSet.h"
#include "Modules.h"

#define TITLE   "Elasticity"
#define AUTHORS " "

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
#define SIG     (vim + 0)
#define F_MASS  (vim + 9)

/* We define some names for constant terms */
#define SIG0    (v0  + 0)


/* Functions */
static Model_ComputePropertyIndex_t  pm ;
static void    GetProperties(Element_t*,double) ;
static double* IsotropicElasticTensor(double,double,double*) ;
static double* MicrostructureElasticTensor(char*,double*) ;

static double* ComputeVariables(Element_t*,double**,double*,double,double,int) ;
static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;

static void    StressAveraging(Mesh_t*,double*) ;
static void    ComputeMicrostructure(DataSet_t*,double*,double*) ;

static double* MacroGradient(Element_t*,double) ;
static double* MacroStrain(Element_t*,double) ;


/* Material parameters */
static double  gravity ;
static double  rho_s ;
static double* sig0 ;
static double  macrogradient[9] ;
static double  macrostrain[9] ;
static double* cijkl ;

//#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
#define GetProperty(a)      Element_GetPropertyValue(el,a)

#define ItIsPeriodic  (Geometry_IsPeriodic(Element_GetGeometry(el)))


#define NbOfVariables  (24)

#define I_EPS          (3)
#define I_SIG          (12)
#define I_Fmass        (21)


int pm(const char* s)
{
       if(!strcmp(s,"gravity")) return(0) ;
  else if(!strcmp(s,"rho_s"))   return(1) ;
  else if(!strcmp(s,"poisson")) return(2) ;
  else if(!strcmp(s,"young"))   return(3) ;
  else if(!strcmp(s,"sig0"))    return(4) ;
  else if(!strncmp(s,"sig0_",5)) {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(4 + 3*i + j) ;
  } else if(!strcmp(s,"macro-gradient")) return(13) ;
  else if(!strncmp(s,"macro-gradient_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(13 + 3*i + j) ;
  } else if(!strcmp(s,"Cijkl")) {
    return(22) ;
  } else if(!strcmp(s,"macro-fctindex")) {
    return(103) ;
  } else if(!strncmp(s,"macro-fctindex_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(103 + 3*i + j) ;
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


void GetProperties(Element_t* el,double t)
{
  gravity = GetProperty("gravity") ;
  rho_s   = GetProperty("rho_s") ;
  sig0    = &GetProperty("sig0") ;
  cijkl   = &GetProperty("Cijkl") ;
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
  Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd */
{
  int NbOfProp = 112 ;
  int i ;

  /* Par defaut tout a 0 */
  for(i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;

  Material_ScanProperties(mat,datafile,pm) ;
  
  /* The 4th rank elastic tensor */
  {
    double* c = Material_GetProperty(mat) + pm("Cijkl") ;
    char* method = Material_GetMethod(mat) ;
    
    /* obtained from a microstructure */
    if(!strncmp(method,"Microstructure",14)) {
      char* p = strstr(method," ") ;
      char* cellname = p + strspn(p," ") ;
      
      MicrostructureElasticTensor(cellname,c) ;
      
    /* isotropic Hooke's law */
    } else {
      double young = Material_GetProperty(mat)[pm("young")] ;
      double poisson = Material_GetProperty(mat)[pm("poisson")] ;
    
      IsotropicElasticTensor(young,poisson,c) ;

    }

#if 0
    {
      printf("4th rank elastic tensor:\n") ;
      
      for(i = 0 ; i < 9 ; i++) {
        int j = i - (i/3)*3 ;
        
        printf("C%d%d--:",i/3 + 1,j + 1) ;
        
        for (j = 0 ; j < 9 ; j++) {
          printf(" % e",c[i*9 + j]) ;
        }
        
        printf("\n") ;
      }
    }
#endif
  }

  return(NbOfProp) ;
}



int PrintModelChar(Model_t* model,FILE* ficd)
/** Print the model characteristics */
{
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- Mechanical Equilibrium     (meca_[1,2,3])\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Displacements              (u_[1,2,3]) \n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  printf("gravity = 0          # gravity\n") ;
  printf("rho_s = 0            # mass density of the dry material\n") ;
  printf("young = 2713e6       # Young's modulus\n") ;
  printf("poisson = 0.339      # Poisson's ratio\n") ;
  printf("sig0_11 = 0          # initial stress 11 (any ij are allowed)\n") ;
  printf("macro-strain_23 = 1  # For periodic BC (any ij)\n") ;
  
  return(0) ;
}



int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads */
{
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FEM_t* fem = FEM_GetInstance(el) ;
  
  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,intfct,cg,t,dt) ;
    
    {
      int    i ;
      
      for(i = 0 ; i < NEQ*nn ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}



int ComputeInitialState(Element_t* el,double t)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 */ 
{
  double* vim0 = Element_GetImplicitTerm(el) ;
  double* v00  = Element_GetConstantTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  DataFile_t* datafile = Element_GetDataFile(el) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input Data
  */
  GetProperties(el,t) ;
    
    
  /* Pre-initialization */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* storage in v0 */
    {
      double* v0   = v00 + p*NV0 ;
      double* vim  = vim0 + p*NVI ;
      int    i ;
      
      /* Initial stresses */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
        for(i = 0 ; i < 9 ; i++) SIG0[i] = SIG[i] ;
      } else {
        for(i = 0 ; i < 9 ; i++) SIG0[i] = sig0[i] ;
      }
    }
  }
  
    
  /* If there are initial displacements */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,vim0,t,0,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
      
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
/** Compute the implicit terms */
{
  double* vim0 = Element_GetCurrentImplicitTerm(el) ;
  double* vex0 = Element_GetCurrentExplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el,t) ;
    
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,vim0,t,dt,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
    }
  }

  return(0) ;
}




int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int ndof = nn*NEQ ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input Data
  */
  GetProperties(el,t) ;
  
  /*
  ** Elastic Matrix
  */
  {
    FEM_t* fem = FEM_GetInstance(el) ;
    //double c[IntFct_MaxNbOfIntPoints*100] ;
    //int dec = Cijkl(fem,c) ;

    //double* kp = FEM_ComputeElasticMatrix(fem,intfct,c,dec) ;
    double* kp = FEM_ComputeElasticMatrix(fem,intfct,cijkl,0) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
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
  FEM_t* fem = FEM_GetInstance(el) ;
  int ndof = nn*NEQ ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
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
    
    for(i = 0 ; i < nn ; i++) R(i,E_mec + dim - 1) -= -rbf[i] ;
  }
  
  return(0) ;

#undef R
}




int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 3 ;
  double* vim0 = Element_GetCurrentImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t* fem = FEM_GetInstance(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el,t) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;    
    /* Variables */
    //double* x = ComputeVariables(el,u,vim0,t,0,p) ;
    
    double pdis[3] = {0,0,0} ;
    double dis[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    int    i ;
    
    /* Displacement */
    for(i = 0 ; i < dim ; i++) {
      pdis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
      dis[i] = pdis[i] ;
    }

    if(ItIsPeriodic) {
      for(i = 0 ; i < dim ; i++) {
        int j ;
        
        for(j = 0 ; j < dim ; j++) {
          dis[i] += MacroGradient(el,t)[3*i + j] * s[j] ;
        }
      }
    }
      
    /* Averaging */
    for(p = 0 ; p < np ; p++) {
      double* vim  = vim0 + p*NVI ;
      
      for(i = 0 ; i < 9 ; i++) sig[i] += SIG[i]/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,dis ,"Displacements",3) ;
    Result_Store(r + i++,sig ,"Stresses",9) ;
    Result_Store(r + i++,pdis,"Perturbated-displacements",3) ;
      
    if(i != NbOfOutputs) arret("ComputeOutputs") ;
  }

  return(NbOfOutputs) ;
}





double* ComputeVariables(Element_t* el,double** u,double* f_n,double t,double dt,int p)
{
  Model_t*  model  = Element_GetModel(el) ;
  double*   x      = Model_GetVariable(model,p) ;
  
  {
    IntFct_t* intfct = Element_GetIntFct(el) ;
    FEM_t*    fem    = FEM_GetInstance(el) ;
    int dim = Element_GetDimensionOfSpace(el) ;
    
    /* Primary Variables */
    {
      int i ;
    
      for(i = 0 ; i < dim ; i++) {
        x[U_u + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
      }
    
      for(i = dim ; i < 3 ; i++) {
        x[U_u + i] = 0 ;
      }
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;
      int i ;
      
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
    
    /* Initial stresses */
    {
      double* v00 = Element_GetConstantTerm(el) ;
      double* v0  = v00 + p*NV0 ;
      int i ;
  
      for(i = 0 ; i < 9 ; i++) {
        x[I_SIG + i] = SIG0[i] ;
      }
    }
  }
  
  /* Needed variables to compute secondary components */
    
  ComputeSecondaryVariables(el,dt,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double dt,double* x)
{
  /* Strains */
  double* eps =  x + I_EPS ;
    

  /* Backup stresses */
  {
    double* sig  = x + I_SIG ;
    int    i ;
    
    #define C(i,j)  (cijkl[(i)*9+(j)])
    for(i = 0 ; i < 9 ; i++) {
      int  j ;
      
      for(j = 0 ; j < 9 ; j++) {
        sig[i] += C(i,j)*eps[j] ;
      }
    }
    #undef C
  }
  
  /* Backup body forces */
  {
    int dim = Element_GetDimensionOfSpace(el) ;
    double* fmass = x + I_Fmass ;
    int    i ;
      
    for(i = 0 ; i < 3 ; i++) fmass[i] = 0. ;
      
    fmass[dim-1] = (rho_s)*gravity ;
  }
}



double* IsotropicElasticTensor(double Young,double Poisson,double* c)
/** Compute the 4th rank isotropic elastic tensor in c.
 *  Return c  */
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  double twomu   = Young/(1 + Poisson) ;
  double mu      = twomu/2 ;
  double lame    = twomu*Poisson/(1 - 2*Poisson) ;
   
  {
    int    i ;

    for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
    for(i = 0 ; i < 3 ; i++) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) {
        C(i,i,j,j) += lame ;
        C(i,j,i,j) += mu ;
        C(i,j,j,i) += mu ;
      }
    }
  }
  
  return(c) ;
#undef C
}



  
double* MicrostructureElasticTensor(char* cellname,double* c)
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  int verb = Message_SetVerbosity(0) ;
  Options_t* options = Options_Create() ;
  DataSet_t* jdd = DataSet_Create(cellname,options) ;
  
  {
    double grd[9] = {0,0,0,0,0,0,0,0,0} ;
    int k ;
    
    for(k = 0 ; k < 3 ; k++) {
      int l ;
      
      for(l = 0 ; l < 3 ; l++) {
        double sig[9] ;
        int i ;
          
        /* Set the macro gradient (k,l) as 1 */
        grd[k*3 + l] = 1 ;
          
        /* For eack (k,l) compute the stresses Sij = Cijkl */
        ComputeMicrostructure(jdd,grd,sig) ;

        for(i = 0 ; i < 3 ; i++) {
          int j ;
                
          for(j = 0 ; j < 3 ; j++) {
            C(i,j,k,l) = sig[i*3 + j] ;
          }
        }
        
        grd[k*3 + l] = 0 ;
      }
    }
  }
  
  Message_SetVerbosity(verb) ;
  
  return(c) ;
#undef C
}



void ComputeMicrostructure(DataSet_t* jdd,double* macrograd,double* sig)
{
  /* Set input data */
  {
    Materials_t* mats = DataSet_GetMaterials(jdd) ;
    int nmats = Materials_GetNbOfMaterials(mats) ;
    int j ;
    
    for(j = 0 ; j < nmats ; j++) {
      Material_t* mat = Materials_GetMaterial(mats) + j ;
      double* grd = Material_GetProperty(mat) + pm("macro-gradient") ;
      int i ;
    
      for(i = 0 ; i < 9 ; i++) {
        grd[i] = macrograd[i] ;
      }
    }
  }
    
  /* Compute the microstructure */
  {
    Modules_t* modules = DataSet_GetModules(jdd) ;
    Module_t* module_i = Modules_FindModule(modules,"Module1") ;
    
    Module_ComputeProblem(module_i,jdd) ;
  }

  /* Backup stresses as averaged stresses */
  {
    Mesh_t* mesh = DataSet_GetMesh(jdd) ;
    
    StressAveraging(mesh,sig) ;
  }
}



void   StressAveraging(Mesh_t* mesh,double* x)
{
  unsigned int nel = Mesh_GetNbOfElements(mesh) ;
  Element_t* el0 = Mesh_GetElement(mesh) ;
  double area = 0 ;
  int i ;
  
  /* The surface area */
  {
    unsigned int ie ;
    
    for(ie = 0 ; ie < nel ; ie++) {
      Element_t* el = el0 + ie ;
      IntFct_t* intfct = Element_GetIntFct(el) ;
      FEM_t*    fem    = FEM_GetInstance(el) ;
      double one = 1 ;
    
      if(Element_IsSubmanifold(el)) continue ;
      
      area +=  FEM_IntegrateOverElement(fem,intfct,&one,0) ;
    }
  }
  
  /* Stress integration */
  for(i = 0 ; i < 9 ; i++) {
    double sig = 0 ;
    unsigned int ie ;
    
    for(ie = 0 ; ie < nel ; ie++) {
      Element_t* el = el0 + ie ;
      double* vim = Element_GetCurrentImplicitTerm(el) ;
      IntFct_t* intfct = Element_GetIntFct(el) ;
      FEM_t*    fem    = FEM_GetInstance(el) ;
    
      if(Element_IsSubmanifold(el)) continue ;
    
      sig +=  FEM_IntegrateOverElement(fem,intfct,SIG + i,NVI) ;
    }
    
    /* Stress average */
    x[i] = sig/area ;
  }
}


/* Not used from here */
#if 0
int Cijkl(FEM_t* fem,double* c)
/*
**  Elastic matrix (c), return the shift (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
  Element_t* el = FEM_GetElement(fem) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  double twomu,lame,mu ;
  int    dec = 81 ;
  int    p ;
  
  twomu   = young/(1 + poisson) ;
  mu      = twomu/2 ;
  lame    = twomu*poisson/(1 - 2*poisson) ;
   
  for(p = 0 ; p < np ; p++) {
    int    i,j ;
    double* c1 = c + p*dec ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
      
    /* Mechanics */
    /* derivative of sig with respect to eps */
    for(i = 0 ; i < 3 ; i++) for(j = 0 ; j < 3 ; j++) {
      C1(i,i,j,j) += lame ;
      C1(i,j,i,j) += mu ;
      C1(i,j,j,i) += mu ;
    }

  }
  
  return(dec) ;
  
#undef C1
}
#endif
