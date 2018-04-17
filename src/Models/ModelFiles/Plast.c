#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"

#define TITLE "Plasticity with hardening (2017)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (dim)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (23)
#define NVE     (0)
#define NV0     (9)

/* Equation index */
#define E_mec   (0)

/* Unknown index */
#define U_u     (0)

/* We define some names for implicit terms */
#define SIG           (vim + 0)
#define F_MASS        (vim + 9)
#define EPS_P         (vim + 12)
#define GAM_P         (vim[21])
#define CRIT          (vim[22])

#define SIG_n         (vim_n + 0)
#define EPS_P_n       (vim_n + 12)
#define GAM_P_n       (vim_n[21])

/* We define some names for explicit terms */

/* We define some names for constant terms */
#define SIG0          (v0  + 0)


/* Functions */
static Model_ComputePropertyIndex_t  pm ;
static void   GetProperties(Element_t*,double) ;
static double* IsotropicElasticTensor(double,double,double*) ;
//static double* MicrostructureElasticTensor(char*,double*) ;

//static void  StressAveraging(Mesh_t*,double*) ;
//static void  ComputeMicrostructure(DataSet_t*,double*,double*) ;

static int    ComputeTangentCoefficients(FEM_t*,double,double*) ;
static double UpdateElastoplasticTensor(double*,double*,double*,double*,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
//static double* ComputeVariablesDerivatives(Element_t*,double,double*,double,int) ;

typedef double ReturnMapping_t(double*,double*,double*) ;
typedef double Criterion_t(const double*,const double,double*,double*,double*) ;

static ReturnMapping_t ReturnMapping_DruckerPrager ;
static Criterion_t     Criterion_DruckerPrager ;

static double* MacroGradient(Element_t*,double) ;
static double* MacroStrain(Element_t*,double) ;



#define Criterion       Criterion_DruckerPrager
#define ReturnMapping   ReturnMapping_DruckerPrager

//static double (*ReturnMapping[])(double*,double*,double*) = {ReturnMapping_DruckerPrager,ReturnMapping_CamClay} ;

/* Material parameters */
static double  gravity ;
static double  rho_s ;
static double  young,poisson ;
static double* sig0 ;
static double  cohesion,af,ad ;
static double  alpha,gam_R ;
static double  macrogradient[9] ;
static double  macrostrain[9] ;
static double* cijkl ;


//#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
#define GetProperty(a)      Element_GetPropertyValue(el,a)

#define ItIsPeriodic  (Geometry_IsPeriodic(Element_GetGeometry(el)))


#define NbOfVariables     (63)
static double  Variable[NbOfVariables] ;
//static double dVariable[NbOfVariables] ;

#define I_U            (0)

#define I_EPS          (3)

#define I_SIG          (12)
#define I_EPS_P        (21)
#define I_EPS_n        (30)
#define I_EPS_P_n      (39)
#define I_Fmass        (48)
#define I_GAM_P        (51)
#define I_CRIT         (52)
#define I_SIG_n        (53)
#define I_GAM_P_n      (62)



int pm(const char* s)
{
         if(!strcmp(s,"gravity")) {
    return(0) ;
  } else if(!strcmp(s,"rho_s")) {
    return(1) ;
  } else if(!strcmp(s,"poisson")) {
    return(2) ;
  } else if(!strcmp(s,"young")) {
    return(3) ;
  } else if(!strcmp(s,"sig0"))  {
    return(4) ;
  } else if(!strncmp(s,"sig0_",5)) {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(4 + 3*i + j) ;
  } else if(!strcmp(s,"macro-gradient")) {
    return(13) ;
  } else if(!strncmp(s,"macro-gradient_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(13 + 3*i + j) ;
  } else if(!strcmp(s,"Cijkl")) {
    return(22) ;
  } else if(!strcmp(s,"criterion")) {
    return(103) ;
    
    /* Drucker-Prager */
  } else if(!strcmp(s,"cohesion")) {
    return(104) ;
  } else if(!strcmp(s,"frottement")) {
    return(105) ;
  } else if(!strcmp(s,"dilatance"))  {
    return(106) ;
  } else if(!strcmp(s,"alpha"))      {
    return(107) ;
  } else if(!strcmp(s,"gamma_R"))    {
    return(108) ;
    
    /* Cam-clay */
  } else if(!strcmp(s,"kappa")) {
    return(104) ;
  } else if(!strcmp(s,"lambda")) {
    return(105) ;
  } else if(!strcmp(s,"M"))  {
    return(106) ;
  } else if(!strcmp(s,"phi"))      {
    return(107) ;
    
  } else if(!strcmp(s,"macro-fctindex")) {
    return(109) ;
  } else if(!strncmp(s,"macro-fctindex_",15)) {
    int i = (strlen(s) > 15) ? s[15] - '1' : 0 ;
    int j = (strlen(s) > 16) ? s[16] - '1' : 0 ;
    
    return(109 + 3*i + j) ;
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
  young   = GetProperty("young") ;
  poisson = GetProperty("poisson") ;
  
  cohesion       = GetProperty("cohesion") ;
  af      = GetProperty("frottement")*M_PI/180. ;
  ad      = GetProperty("dilatance")*M_PI/180. ;
  alpha   = GetProperty("alpha") ;
  gam_R   = GetProperty("gamma_R") ;
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
/** Read the material properties in the stream file ficd 
 *  Return the nb of (scalar) properties of the model */
{
  int NbOfProp = 118 ;
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
      
      //MicrostructureElasticTensor(cellname,c) ;
      
    /* isotropic Hooke's law */
    } else {
      young = Material_GetProperty(mat)[pm("young")] ;
      poisson = Material_GetProperty(mat)[pm("poisson")] ;
    
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


int PrintModelProp(Model_t* model,FILE *ficd)
/** Print the model properties 
 *  Return the nb of equations */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;

  printf("\n\
The system consists in (dim) equations\n\
\t 1. The equilibrium equation (meca_1,meca_2,meca_3)\n") ;

  printf("\n\
The primary unknowns are\n\
\t 1. The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n\
Example of input data\n\n") ;
  

  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # masse volumique du squelette sec\n") ;
  fprintf(ficd,"young = 5.8e+09   # module d\'Young\n") ;
  fprintf(ficd,"poisson = 0.3     # coefficient de Poisson\n") ;
  fprintf(ficd,"cohesion = 1e+06  # cohesion\n") ;
  fprintf(ficd,"frottement = 25   # frottement\n") ;
  fprintf(ficd,"dilatance = 25    # dilatance \n") ;
  fprintf(ficd,"sig0_ij = -11.5e6 # contrainte initiale sig0_ij\n") ;
  
  return(0) ;
}


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
    Input data
  */
  GetProperties(el,t) ;


  /* Pre-initialization */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      double* v0   = v00 + p*NV0 ;
      int    i ;
      
      /* Initial stresses */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
        for(i = 0 ; i < 9 ; i++) SIG0[i] = SIG[i] ;
      } else {
        for(i = 0 ; i < 9 ; i++) SIG0[i] = sig0[i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]  = sig0[i] ;
      }
      
      /* Initial plastic strains */
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0 ;
    
      GAM_P = 0 ;
    }
  }
  

  /* If there are initial displacements */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      GAM_P = x[I_GAM_P] ;
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
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n  = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
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
    double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;
    
      CRIT = x[I_CRIT] ;
      GAM_P = x[I_GAM_P] ;
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
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;


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
  GetProperties(el,t) ;


  /*
  ** Elastoplastic matrix
  */
  {
    double c[IntFct_MaxNbOfIntPoints*81] ;
    int dec = ComputeTangentCoefficients(fem,dt,c) ;
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
  int NbOfOutputs = 6 ;
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
  GetProperties(el,t) ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Displacement */
    double pdis[3] = {0,0,0} ;
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double gam_p = 0 ;
    double crit = 0 ;
    int    i ;
    
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
      int j ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) eps_p[j] += EPS_P[j]/np ;
      
      gam_p += GAM_P/np ;
      
      crit += CRIT/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,dis   ,"Displacements",3) ;
    Result_Store(r + i++,sig   ,"Stresses",9) ;
    Result_Store(r + i++,pdis  ,"Perturbated-displacements",3) ;
    Result_Store(r + i++,eps_p ,"Plastic-strains",9) ;
    Result_Store(r + i++,&gam_p,"Cumulative-plastic-shear-strain",1) ;
    Result_Store(r + i++,&crit ,"Yield function",1) ;
  }
  
  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(FEM_t* fem,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
  Element_t* el  = FEM_GetElement(fem) ;
  double*  vim0   = Element_GetCurrentImplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  int    dec = 81 ;
  int    p ;
  double zero = 0. ;
  
  /*
  ObVal_t* obval = Element_GetObjectiveValue(el) ;
  double dxi[Model_MaxNbOfEquations] ;
  
  {
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }
  */

  
  for(p = 0 ; p < np ; p++) {
    double* vim  = vim0 + p*NVI ;
    double* c0 = c + p*dec ;
    /* Variables */
    //double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
    
    /* Initialization */
    {
      int i ;
      
      for(i = 0 ; i < dec ; i++) c0[i] = zero ;
    }
    

    /* Mechanics */
    {
      double* c1 = c0 ;
      int i ;
        
      for(i = 0 ; i < 81 ; i++) {
        c1[i] = cijkl[i] ;
      }
      
      /* Tangent stiffness matrix */
      {
        /* Criterion */
        double crit = CRIT ;
        
        if(crit >= 0.) {
          
          /* Yield and potential function gradients */
          double dfsds[9] ;
          double dgsds[9] ;
          double fc[9] ;
          double cg[9] ;
          double hm ;
          double crit1 = Criterion(SIG,GAM_P,dfsds,dgsds,&hm) ;
          double fcg = UpdateElastoplasticTensor(dfsds,dgsds,fc,cg,hm,c1) ;
          
          if(fcg < 0) return(-1) ;
        }
      }
    }
  }
  
  return(dec) ;
#undef C1
#undef T4
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
//  double*   x      = Model_GetVariable(model,p) ;
  double*   x      = Variable ;
  
    
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
        x[I_EPS_P_n + i] = EPS_P_n[i] ;
      }
      
      x[I_GAM_P_n] = GAM_P_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
  }
    
  ComputeSecondaryVariables(el,dt,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double dt,double* x)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  x + I_EPS ;
  double* eps_n =  x + I_EPS_n ;
  /* Plastic strains */
  double* eps_p  = x + I_EPS_P ;
  double* eps_pn = x + I_EPS_P_n ;
    


  /* Backup stresses, plastic strains */
  {
    double* sig   = x + I_SIG ;
    double* sig_n = x + I_SIG_n ;
    double  gam_p = x[I_GAM_P_n] ;
    
    
    {
      double  deps[9] ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Elastic trial stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
    
      #define C(i,j)  (cijkl[(i)*9+(j)])
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
      
        for(j = 0 ; j < 9 ; j++) {
          sig[i] += C(i,j)*deps[j] ;
        }
      }
      #undef C
    
      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;
    
      /* Projection */
      {
        double crit = ReturnMapping(sig,eps_p,&gam_p) ;
        
        x[I_CRIT]  = crit ;
        x[I_GAM_P] = gam_p ;
      }
    }
  }
  
  
  /* Backup body force */
  {
    {
      double* f_mass = x + I_Fmass ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s)*gravity ;
    }
  }
}





#if 0
double* ComputeVariablesDerivatives(Element_t* el,double dt,double* x,double dxi,int i)
{
  double* dx = dVariable ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}
#endif



double UpdateElastoplasticTensor(double* dfsds,double* dgsds,double* fc,double* cg,double fcg0,double* c)
/** Update the 4th rank elastoplastic tensor in c.
 *  Inputs are: 
 *  dfsds is the yield function gradient
 *  dgsds is the potential function gradient
 *  Outputs are:
 *  fc(k,l) = dfsds(j,i) * C(i,j,k,l)
 *  cg(i,j) = C(i,j,k,l) * dgsds(l,k)
 *  fcg is updated as fcg += dfsds(j,i) * C(i,j,k,l) * dgsds(l,k))
 *  Tensor c is then updated as C(i,j,k,l) += cg(i,j) * fc(k,l) / fcg
 *  Return the inverse of fcg: 1/fcg */
{
#define C(i,j)  ((c)[(i)*9+(j)])
  double fcg = fcg0 ;
  
  /* Tangent elastoplastic tensor */
  {
      
    /* Criterion */
    {
      int i ;
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        fc[i] = 0. ;
        cg[i] = 0. ;
          
        for(j = 0 ; j < 9 ; j++) {
              
          fc[i] += dfsds[j]*C(j,i) ;
          cg[i] += C(i,j)*dgsds[j] ;
          fcg   += dfsds[i]*C(i,j)*dgsds[j] ;
        }
      }
        
      if(fcg > 0.) {
        fcg = 1./fcg ;
      } else {
            
        printf("\n") ;
        printf("dfsds = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dfsds[i]) ;
        }
        printf("\n") ;
        printf("dgsds = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dgsds[i]) ;
        }
        printf("\n") ;
        printf("fcg = %e\n",fcg) ;
        printf("\n") ;
          
        return(-1) ;
      }
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        for(j = 0 ; j < 9 ; j++) {
          C(i,j) -= cg[i]*fc[j]*fcg ;
        }
      }
    }
  }
  
  return(fcg) ;
#undef C
}



double Criterion_DruckerPrager(const double* sig,const double gam_p,double* dfsds,double* dgsds,double* hm)
/** Drucker-Prager criterion. 
 *  Inputs are: 
 *  the stresses (sig), the cumulative plastic shear strain (gam_p). 
 *  Parameters are:
 *  the friction angle (af), the dilatancy angle (ad) and the cohesion.
 *  On outputs the following values are modified:
 *  dfsds = derivative of the yield function wrt stresses
 *  dgsds = derivative of the potential function wrt stresses
 *  hm    = hardening modulus
 *  Return the value of the yield function. */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9] ;
  double crit ;
  double ff,dd,cc,cc0 ;
  int    i ;
  
  /*
    Input data
  */
  /*
  ff      = 3.*sin(af)/sqrt(3.+sin(af)*sin(af)) ;
  dd      = 3.*sin(ad)/sqrt(3.+sin(ad)*sin(ad)) ;
  cc      = 3.*cos(af)/sqrt(3.+sin(af)*sin(af))*c ;
  */
  ff      = 6.*sin(af)/(3. - sin(af)) ;
  dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  
  /*
    Cohesion
  */
  {
    //double c1 = (gam_p < gam_R) ? 1 - (1 - alpha)*gam_p/gam_R : alpha ;
    double c1 = 1 ;
    
    cc = cc0*c1*c1 ;
  }
  
  /*
    Yield function
  */ 
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q + ff*p - cc ;
  
  /*
    Deviatoric stresses
  */
  for(i = 0 ; i < 9 ; i++) {
    dev[i] = sig[i] - p*id[i] ;
  }
  
  /*
    Yield function gradient
  */
  if(q > 0.) {
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = 1.5*dev[i]/q + id[i]*ff/3. ;
    }
    
  } else {
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = id[i]*ff/3. ;
    }
  }
  
  /*
    Potential function gradient
  */
  
  /* Elastic case */
  if(crit <= 0.) {
    if(q > 0.) {
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    } else {
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = id[i]*dd/3. ;
      }
    }
  }
  
  /* Plastic case */
  if(crit > 0.) {
    double k   = young/(1. - 2.*poisson)/3. ;
    double dmu = young/(1.+poisson) ;
    double mu  = 0.5*dmu ;
    
    /* Smooth flow regime */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    /* Flow regime at the notch apex */
    } else {
      double dl = (ff*p - cc)/(k*ff*dd) ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = dev[i]/(dmu*dl) + id[i]*dd/3. ;
      }
    }
  }
  
  /* Hardening modulus */
  *hm = 0. ;
  /*
  if(gam_p < gam_R) {
    *hm = -2.*(1.-alpha)/gam_R*(1.-(1.-alpha)*gam_p/gam_R)*cc0 ;
    *hm *= sqrt(2*j2(dgsds)) ;
  }
  */
  
  return(crit) ;
}



double ReturnMapping_DruckerPrager(double* sig,double* eps_p,double* gam_p)
/** Drucker-Prager return mapping.
 *  Parameters are:
 *  the Young modulus (young),
 *  the Poisson's ratio (poisson),
 *  the friction angle (af), 
 *  the dilatancy angle (ad),
 *  the cohesion (cohesion).
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the cumulative plastic shear strain (gam_p).
 *  Return the value of the yield function. */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t,q_t,sdev_t[9] ;
  double crit ;
  double ff,dd,k,dmu,mu,cc0 ;
  int    i ;
  
  /*
    Data
  */
  dmu     = young/(1.+poisson) ;
  mu      = dmu/2. ;
  k       = young/(1. - 2.*poisson)/3. ;
  /*
  ff      = 3.*sin(af)/sqrt(3.+sin(af)*sin(af)) ;
  dd      = 3.*sin(ad)/sqrt(3.+sin(ad)*sin(ad)) ;
  cc0     = 3.*cos(af)/sqrt(3.+sin(af)*sin(af))*cohesion ;
  */
  ff      = 6.*sin(af)/(3. - sin(af)) ;
  dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  
  /*
    Trial stresses
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*j2(sig)) ;
  for(i = 0 ; i < 9 ; i++) {
    sdev_t[i] = sig[i] - p_t*id[i] ;
  }
  
  /*
    Criterion
  */ 
  {
    //double c1 = ((*gam_p) < gam_R) ? 1 - (1 - alpha)*(*gam_p)/gam_R : alpha ;
    double c1 = 1 ;
    double cc = cc0*c1*c1 ;
    
    crit = q_t + ff*p_t - cc ;
  }
  
  /*
    Return mapping: update plastic strains and stresses
  */
  if(crit > 0.) {
    double deps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double p = p_t ;
    double q = q_t ;
    double gam_pn = *gam_p ;
    double gam_p1 = *gam_p ;
    
    /* Smooth flow regime: assuming that q > 0 */
    if(q > 0) {
      double fcrit = crit ;
      int    nf    = 0 ;
      double dl    = 0 ;
      double tol   = 1.e-10 ;
      
      //while(fabs(fcrit) > tol*cc0) {
      while(fabs(fcrit/q) > tol) {
        double dqsdl ;
        double dpsdl ;
        double dccsdl ;
        double c1 ;
        double cc ;
        
	      /* Plastic strain increments */
	      for(i = 0 ; i < 9 ; i++) {
	        deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
	      }
        
	      /* Cumulative plastic shear strain */
	      gam_p1 = gam_pn + sqrt(2*j2(deps_p)) ;
        
	      /* p, q, cc */
	      //c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
	      c1 = 1 ;
	      cc = cc0*c1*c1 ;
	      q  = q_t - dl*3*mu ;
	      p  = p_t - dl*k*dd ;
        
	      /* dqsdl, dpsdl, dccsdl */
	      dqsdl = -3*mu ;
	      dpsdl = -k*dd ;
	      dccsdl = 0. ;
        /*
	      if(gam_p1 < gam_R) {
	        dccsdl = -2*(1 - alpha)/gam_R*c1*cc0 ;
	        dccsdl *= 1.5*sqrt(2*j2(sdev_t))/q_t ;
	      }
        */
        
        /* Criterion */
        fcrit = q + ff*p - cc ;
        
	      /* dl */
        {
          double df = dqsdl + ff*dpsdl - dccsdl ;
          
	        dl   -= fcrit/df ;
        }
        
	      if(nf++ > 20) {
	        printf("No convergence (ReturnMapping_DruckerPrager)") ;
	        exit(0) ;
	      }
      }
      
      /* Stresses */
      for(i = 0 ; i < 9 ; i++) {
	      sig[i] = sdev_t[i]*q/q_t + p*id[i] ;
      }
    }
    
    /* Flow regime at the notch apex */
    if(q <= 0.) {
      double c1 ;
      double cc ;
      double dl ;
      
      /* Deviatoric plastic strain increments */
      for(i = 0 ; i < 9 ; i++) deps_p[i]  = sdev_t[i]/dmu ;
        
      /* Cumulative plastic shear strain */
      gam_p1 = gam_pn + sqrt(2*j2(deps_p)) ;
      
      /* p, q, cc */
      //c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
      c1 = 1 ;
      cc = cc0*c1*c1 ;
      p  = cc/ff ;
      q  = 0. ;
      
      /* dl */
      dl   = (ff*p_t - cc)/(k*ff*dd) ;
      
      /* Plastic strain increments and stresses */
      for(i = 0 ; i < 9 ; i++) {
	      deps_p[i] = sdev_t[i]/dmu + dl*id[i]*dd/3. ;
	      sig[i]    = p*id[i] ;
      }
    }
    
      
    /* Total plastic strains */
    for(i = 0 ; i < 9 ; i++) eps_p[i] += deps_p[i] ;
    (*gam_p) = gam_p1 ;
  }
  
  return(crit) ;
}




#if 0
double Criterion_CamClay(double* sig,double pc,double* dfsds,double* dgsds,double* hm,Element_t *el)
/** Cam-Clay criterion */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit,m2 ;
  int    i ;
  /*
    Donnees
  */
  m       = GetProperty("M") ;
  m2      = m*m ;
  /* 
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Gradients
  */
  for(i = 0 ; i < 9 ; i++) {
    double dev = sig[i] - p*id[i] ;
    
    dfsds[i] = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    dgsds[i] = dfsds[i] ;
  }
  
  /* The hardening modulus */
  kappa   = GetProperty("kappa") ;
  lambda  = GetProperty("lambda") ;
  {
    double v = 1./(lambda - kappa) ;
    
    phi0    = GetProperty("phi") ;
    
    *hm = v/(1 - phi0)*p*(2*p + pc)*pc ;
  }
  return(crit) ;
}



double ReturnMapping_CamClay(double* sig,double* eps_p,double* p_co,Element_t *el)
/** Cam-Clay return mapping. Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (p_co),
 *  the porosity or the void ratio (phi0,e0).
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the hardening variable (p_co).
 *  Return the value of the yield function. 
 *  Algorithm from Borja & Lee 1990 modified by Dangla. */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t,q_t,p,q,pc,crit,m2,v,a ;
  double dl ;
  int    i ;
  /*
    Data
  */
  kappa   = GetProperty("kappa") ;
  mu      = GetProperty("mu") ;
  lambda  = GetProperty("lambda") ;
  m       = GetProperty("M") ;
  phi0    = GetProperty("phi") ;
  m2      = m*m ;
  v       = 1./(lambda - kappa) ;
  
  /* 
     The criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  pc   = *p_co ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
     Closest point projection algorithm.
   * Only one iterative loop is used to solve
                    q*q/m2 + p*(p + pc) = 0
     for p. The other variables (pc,q,dl) are expressed with p.
   */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    double tol = 1.e-8 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      double dfsdp  = 2*p + pc ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p ;
      double dpcsdp = -v*kappa*pc/p ;
      double dlsdp  = ((1 - phi0)*kappa/p - dl*(2+dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*dlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Variables (pc,dl,q) are explicit in p */
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      
      if(nf++ > 20) {
	      printf("pas de convergence (ReturnMapping_CamClay)") ;
	      exit(0) ;
      }
    }
  }
  
  /*
    Plastic stresses and strains
  */
  a = 1./(1 + 6*mu/m2*dl) ;
  
  for(i = 0 ; i < 9 ; i++) {
    double dev      = a*(sig[i] - p_t*id[i]) ;
    double dfsds    = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    
    sig[i]   = p*id[i] + dev ;
    eps_p[i] = dl*dfsds ;
  }
  
  /* Consolidation pressure */
  *p_co = pc ;
  return(crit) ;
}
#endif





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
