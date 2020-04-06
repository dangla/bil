#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"
#include "FEM.h"

#define TITLE   "Richards Equation (3D)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ     (1)
#define NVI     (4)
#define NVE     (1)
#define NV0     (0)

/* Equation index */
#define E_liq   (0)

/* Unknown index*/
#define U_P_L   (0)

/* We define some names for implicit terms */
#define M_L           (vim[0])
#define W_L           (vim + 1)

#define M_L_n         (vim_n[0])

/* We define some names for explicit terms */
#define K_L           (vex[0])


/* Material properties */
#define SATURATION_CURVE       (Element_GetCurve(el))
#define RELATIVEPERM_CURVE     (Element_GetCurve(el) + 1)

#define SATURATION(x)    (Curve_ComputeValue(SATURATION_CURVE,x))
#define DSATURATION(x)   (Curve_ComputeDerivative(SATURATION_CURVE,x))
#define RELATIVEPERM(x)  (Curve_ComputeValue(RELATIVEPERM_CURVE,x))



/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*) ;
static int    c10(FEM_t*,double*) ;
static int    k10(FEM_t*,double*) ;

static double* ComputeVariables(Element_t*,void*,void*,void*,const double,const double,const int) ;
static void    ComputeSecondaryVariables(Element_t*,const double,const double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,const double,const double,double*,const double,const int) ;


/* Parameters */
static double gravite,phi,rho_l,k_int,mu_l,p_g ;

  
enum {
I_P_L  = NEQ   ,
I_M_L          ,
I_S_L          ,
I_K_L          ,
I_GRD_P_L      ,
I_W_L = I_GRD_P_L + 3 ,
I_Last
} ;


#define NP               IntFct_MaxNbOfIntPoints
#define NbOfVariables    (I_Last)
static double Variables[NP][NbOfVariables] ;
static double dVariables[NbOfVariables] ;


int pm(const char *s) {
  if(strcmp(s,"gravite") == 0)    return (0) ;
  else if(strcmp(s,"phi") == 0)   return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0)  return (4) ;
  else if(strcmp(s,"p_g") == 0)   return (5) ;
  else return(-1) ;
}


void GetProperties(Element_t *el)
{
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  p_g     = GetProperty("p_g") ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_liq, "liq") ;

  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 7 ;

  Material_ScanProperties(mat,datafile,pm) ;
  
  return(NbOfProp) ;
}



int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("The equation to be solved is:\n") ;
  printf("\t- Mass balance of liquid (liq)\n") ;
  
  printf("\n") ;
  printf("The primary unknown is:\n") ;
  printf("\t- Liquid pressure        (p_l)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"gravite = -9.81  # La gravite\n") ;
  fprintf(ficd,"phi = 0.38       # La porosite\n") ;
  fprintf(ficd,"rho_l = 1000     # La masse volumique du fluide\n") ;
  fprintf(ficd,"k_int = 8.9e-12  # La permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001     # La viscosite du fluide\n") ;
  fprintf(ficd,"p_g = 1.e5       # La pression du gaz\n") ;
  fprintf(ficd,"Curves = billes  # Le nom du fichier p_c S_l k_rl\n") ;
  
  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  
  /* Continuity of pressure across zero-thickness element */
  {
    if(Element_HasZeroThickness(el)) {
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,U_P_L) ;
    }
  }
  
  {
    IntFct_t *intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  
    /** Define the length of tables */
    Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
    Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
    Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  }
  
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  IntFct_t *fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  
  {
    double *r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    int i ;
    
    for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}



int ComputeInitialState(Element_t *el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 */ 
{
  double* vim0 = Element_GetImplicitTerm(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  //if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Pre-initialization */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    double pl = FEM_ComputeUnknown(fem,u,intfct,p,U_P_L) ;
    double pc = p_g - pl ;
    
    /* transfert coefficient */
    double k_l = rho_l*k_int/mu_l*RELATIVEPERM(pc) ;
    
    /* storage in vex */
    double* vex = vex0 + p*NVE ;
    K_L = k_l ;
  }
  
  
    
  /* boucle sur les points d'integration */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    double* vim = vim0 + p*NVI ;
    double* vex = vex0 + p*NVE ;
    
    /* storage in vim */
    M_L = x[I_M_L] ;
    
    {
      int i ;
      double* w_l = x + I_W_L ;
      
      for(i = 0 ; i < 3 ; i++) W_L[i] = w_l[i] ;
    }
    
    /* storage in vex */
    K_L = x[I_K_L] ;
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t *el,double t)
/** Compute the explicit terms */
{
  double* vim0 = Element_GetPreviousImplicitTerm(el) ;
  double *vex0 = Element_GetExplicitTerm(el) ;
  double** u  = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  //if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    double* vex = vex0 + p*NVE ;
    
    /* pressures */
    double pl  = x[I_P_L] ;
    double pc  = p_g - pl ;
    
    /* transfert coefficient */
    double k_l = rho_l*k_int/mu_l*RELATIVEPERM(pc) ;
    
    /* storage in vex */
    K_L = k_l ;
  }
  return(0) ;
}






int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/** Compute the implicit terms */
{
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n= Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  //if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;
  
  
  
  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
    double* vim = vim0 + p*NVI ;
    
    /* storage in vim */
    M_L = x[I_M_L] ;
    
    {
      int i ;
      double* w_l = x + I_W_L ;
      
      for(i = 0 ; i < 3 ; i++) W_L[i] = w_l[i] ;
    }
  }

  return(0) ;
}






int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/** Compute the matrix (k) */
{
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double c[IntFct_MaxNbOfIntPoints*9] ;
  int i ;
  double zero = 0. ;

  
  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  //if(Element_IsSubmanifold(el)) return(0) ;


  /*
  ** Mass Matrix
  */
  {
    int dec = c10(fem,c) ;
    double *kp = FEM_ComputeMassMatrix(fem,intfct,c,dec) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  /*
  ** Conduction Matrix
  */
  {
    int dec = k10(fem,c) ;
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] += dt*kc[i] ;
    }
  }

  return(0) ;
}






int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  //if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /* 1 Accumulation Terms */
  {
    double *vim = vim_1 ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) {
      g1[i] = M_L - M_L_n ;
    }
    
    {
      double *ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
      for(i = 0 ; i < nn ; i++) r[i] -= ra[i] ;
    }
  }
  
  /* 2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_L,NVI) ;
  
    for(i = 0 ; i < nn ; i++) r[i] -= -dt*rf[i] ;
  }
  
  return(0) ;
}





int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 3 ;
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;

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
    //double* x = ComputeVariables(el,u,u,vim,t,0,p) ;
    /* pressures */
    double pl  = FEM_ComputeUnknown(fem,u,intfct,p,U_P_L) ;
    double pc  = p_g - pl ;
    
    /* saturation */
    double sl  = SATURATION(pc) ;
    
    double w_l[3] = {0,0,0} ;
    int    i ;
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) w_l[j] += W_L[j]/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&pl,"pressure",1) ;
    Result_Store(r + i++,w_l,"flow",3) ;
    Result_Store(r + i++,&sl,"saturation",1) ;
  }
  
  return(NbOfOutputs) ;
}




int c10(FEM_t *fem,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
  Element_t* el = FEM_GetElement(fem) ;
  //double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  //double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  //double** u_n= Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  
  int    dec = 1 ;
  int    p ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  for(p = 0 ; p < np ; p++) {
    double* c1 = c + p*dec ;
    /* pressures */
    double pl  = FEM_ComputeUnknown(fem,u,intfct,p,U_P_L) ;
    double pc  = p_g - pl ;
    
    /* saturation */
    double dslsdpc = DSATURATION(pc) ;
    
    c1[0] = -rho_l*phi*dslsdpc ;
  }
  
  return(dec) ;
}



int k10(FEM_t *fem,double *c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  double zero = 0. ;

  for(p = 0 ; p < np ; p++ , vex += NVE) {
    int i ;
    double* c1 = c + p*dec ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;
    
    /* Permeability tensor */
    c1[0] = K_L ;
    c1[4] = K_L ;
    c1[8] = K_L ;
  }
  
  return(dec) ;
}








double* ComputeVariables(Element_t* el,void* vu,void* vu_n,void* vf_n,const double t,const double dt,const int p)
/** This locally defined function compute the intern variables at
 *  the interpolation point p, from the nodal unknowns.
 *  Return a pointer on the locally defined array of the variables. */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  double** u   = (double**) vu ;
  //double** u_n = (double**) vu_n ;
  //double*  f_n = (double*)  vf_n ;
  double*  x   = Variables[p] ;
  
    
  /* Load the primary variables in x */
  {
    int    i ;
    
    /* Pressure */
    x[I_P_L] = FEM_ComputeUnknown(fem,u,intfct,p,U_P_L) ;
    
    /* Pressure gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_P_L) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_P_L + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    double* vex = vex0 + p*NVE ;
    
    /* transfert coefficient */
    x[I_K_L] = K_L ;
  }

    
  ComputeSecondaryVariables(el,t,dt,NULL,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,const double t,const double dt,double* x_n,double* x)
/** Compute the secondary variables from the primary ones. */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Pressures */
  double pl  = x[I_P_L] ;
  double pc  = p_g - pl ;
    
  /* saturation */
  double sl  = SATURATION(pc) ;
    
  /* liquid mass content */
  double m_l = rho_l*phi*sl ;
    
  /* transfert coefficient */
  double k_l = x[I_K_L] ;
    
  /* Pressure gradient */
  double* gpl = x + I_GRD_P_L ;
  double* w_l = x + I_W_L ;
  int i ;
  
  
  /* Backup */
    
  for(i = 0 ; i < 3 ; i++) w_l[i] = - k_l*gpl[i] ;
  w_l[dim-1] += k_l*rho_l*gravite ;
    
  x[I_M_L] = m_l ;
  
}


double* ComputeVariableDerivatives(Element_t* el,const double t,const double dt,double* x,const double dxi,const int i)
{
  double* dx = dVariables ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,t,dt,NULL,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}
