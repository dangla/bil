#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "FVM.h"

#define TITLE   "Two-phase flow (Liquid-Gas)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ     (2)
/* Nb of implicit terms */
#define NVI     (6)
/* Nb of explicit terms */
#define NVE     (2)
/* Nb of constant terms */
#define NV0     (0)


/* Equation index */
#define E_Liq (0)
#define E_Gas (1)

/* Unknown index */
#define U_P_L (0)
#define U_P_G (1)


#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])


/* We define some names for primary unknowns */
#define P_L(n)      (UNKNOWN(n,U_P_L))
#define P_G(n)      (UNKNOWN(n,U_P_G))

/* We define some names for implicit terms */
#define M_L(n)        (vim[(n)])
#define M_G(n)        (vim[(2+n)])
#define W_L           (vim[4])
#define W_G           (vim[5])

#define M_Ln(n)       (vim_n[(n)])
#define M_Gn(n)       (vim_n[(2+n)])

/* We define some names for explicit terms */
#define K_L           (vex[0])
#define K_G           (vex[1])


/* Material properties */
#define SATURATION_CURVE           (Element_GetCurve(el))
#define RELATIVEPERMLIQ_CURVE      (Element_GetCurve(el) + 1)
#define RELATIVEPERMGAS_CURVE      (Element_GetCurve(el) + 2)

#define SATURATION(x)    (saturation(x,p_c3,SATURATION_CURVE))
#define dSATURATION(x)   (dsaturation(x,p_c3,SATURATION_CURVE))
#define RelativePermeabilityToLiquid(x)  (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,x))
#define RelativePermeabilityToGas(x)  (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,x))



/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


/* Functions */
static int     pm(const char*) ;
static void    GetProperties(Element_t*) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static void    ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static double  saturation(double,double,Curve_t*) ;
static double  dsaturation(double,double,Curve_t*) ;


/* Parameters */
static double gravite,phi,rho_l,k_int,mu_l,mu_g,p_c3,MgsRT ;


#define NbOfComponents    (9)
static double Components[NbOfComponents] ;
static double dComponents[NbOfComponents] ;

#define I_P_L          (0)
#define I_P_G          (1)

#define I_M_L          (2)
#define I_M_G          (3)

#define I_RHO_G        (4)
#define I_S_L          (5)
#define I_H_L          (6)
#define I_H_G          (7)
#define I_COOR_X       (8)
  
  

#define NbOfComponentFluxes    (2)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_L           (0)
#define I_W_G           (1)


int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"phi") == 0) return (1) ;
  else if(strcmp(s,"rho_l") == 0) return (2) ;
  else if(strcmp(s,"k_int") == 0) return (3) ;
  else if(strcmp(s,"mu_l") == 0) return (4) ;
  else if(strcmp(s,"mu_g") == 0) return (5) ;
  else if(strcmp(s,"RT") == 0) return (6) ;
  else if(strcmp(s,"M_g") == 0) return (7) ;
  else if(strcmp(s,"p_c3") == 0) return (8) ;
  else return(-1) ;
}


void GetProperties(Element_t *el)
{
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  rho_l   = GetProperty("rho_l") ;
  k_int   = GetProperty("k_int") ;
  mu_l    = GetProperty("mu_l") ;
  mu_g    = GetProperty("mu_g") ;
  p_c3    = GetProperty("p_c3") ;
  
  {
    double M_g   = GetProperty("M_g") ;
    double RT    = GetProperty("RT") ;
    
    MgsRT   = M_g/RT ;
  }
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_Liq,"liq") ;
  Model_CopyNameOfEquation(model,E_Gas,"gas") ;

  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
  Model_CopyNameOfUnknown(model,U_P_G,"p_g") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 9 ;
  
  /* Self-initialization */
  {
    Material_GetProperty(mat)[pm("RT")] = 2436. ;
    Material_GetProperty(mat)[pm("M_g")] = 28.8e-3 ;
  }

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
  printf("\t- Mass balance of gas    (gas)\n") ;
  
  printf("\n") ;
  printf("The primary unknown is:\n") ;
  printf("\t- Liquid pressure        (p_l)\n") ;
  printf("\t- Gas pressure           (p_g)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"gravite = -9.81    # Gravity\n") ;
  fprintf(ficd,"phi = 0.3          # Porosity\n") ;
  fprintf(ficd,"rho_l = 1000       # mass density of liquid\n") ;
  fprintf(ficd,"k_int = 4.4e-13    # Intrinsic Permeability\n") ;
  fprintf(ficd,"mu_l = 1.e-3       # Liquid Viscosity\n") ;
  fprintf(ficd,"mu_g = 1.8e-5      # Gas Viscosity\n") ;
  fprintf(ficd,"M_g  = 28.8e-3     # Molar mass of gas\n") ;
  fprintf(ficd,"RT = 2436          # Product of R and T\n") ;
  fprintf(ficd,"p_c3 = 300         # Capillary pressure parameter\n") ;
  fprintf(ficd,"Curves = sol       # File name: p_c S_l k_rl k_rg\n") ;
  
  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  
  return(0) ;
}


int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  
  {
    double *r1 = FVM_ComputeSurfaceLoadResidu(fvm,cg,t,dt) ;
    int i ;
    
    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/* Initialise les variables du systeme (f,va) */ 
{
  double* vim = Element_GetImplicitTerm(el) ;
  double** u  = Element_ComputePointerToNodalUnknowns(el) ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Fluid mass contents */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,vim,0,i) ;
      
    /* Back up */
    M_L(i) = x[I_M_L] ;
    M_G(i) = x[I_M_G] ;
  }

  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,vim) ;

  /* Fluxes */
  ComputeFluxes(el,u) ;
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Explicites terms  */
{
  double*  f = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Implicit terms */
{
  double* vim = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Fluid mass contents */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,vim_n,dt,i) ;
      
    /* Back up */
    M_L(i) = x[I_M_L] ;
    M_G(i) = x[I_M_G] ;
    
    {
      double pg     = x[I_P_G] ;
    
      if(pg < 0.) return(-1) ;
    }
  }

  /* Fluxes */
  ComputeFluxes(el,u) ;
  
  return(0) ;
}


int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrix (k) */
{
#define K(i,j)    (k[(i)*2*NEQ + (j)])
  int    nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  double c[Element_MaxNbOfNodes*Element_MaxNbOfNodes*NEQ*NEQ] ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  TangentCoefficients(el,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }
  
  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *vim = Element_GetCurrentImplicitTerm(el) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double *volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Boundary Surface Area */
  {
    double *area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }
  
  /*
    Mass balance of liquid and gas
  */
  R(0,E_Liq)  -= volume[0]*(M_L(0)  - M_Ln(0))  + dt*surf*W_L ;
  R(1,E_Liq)  -= volume[1]*(M_L(1)  - M_Ln(1))  - dt*surf*W_L ;
  
  R(0,E_Gas)  -= volume[0]*(M_G(0)  - M_Gn(0))  + dt*surf*W_G ;
  R(1,E_Gas)  -= volume[1]*(M_G(1)  - M_Gn(1))  - dt*surf*W_G ;
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *vim = Element_GetCurrentImplicitTerm(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    nso = 5 ;
  int    i ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;

  /* initialisation */
  for(i = 0 ; i < nso ; i++) {
    Result_SetValuesToZero(r+i) ;
  }

  /* quantites exploitees */
  {
    int    j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Components */
    double *x = ComputeComponents(el,u,vim,0,j) ;
    double pl     = x[I_P_L] ;
    double pg     = x[I_P_G] ;
    double pc     = pg - pl ;
    double sl     = SATURATION(pc) ;
    /* flux */
    double wl  = W_L ;
    double wg  = W_G ;
    
    i = 0 ;
    Result_Store(r + i++,&pl,"p_l",1) ;
    Result_Store(r + i++,&pg,"p_g",1) ;
    Result_Store(r + i++,&wl,"w_l",1) ; 
    Result_Store(r + i++,&wg,"w_g",1) ; 
    Result_Store(r + i++,&sl,"saturation",1) ;
  }

  return(nso) ;
}


void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *vex = Element_GetExplicitTerm(el) ;
  int i ;
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) vex[i] = 0. ;
  
  
  /* Transfer coefficients */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;
    double pl    = x[I_P_L] ;
    double pg    = x[I_P_G] ;
    double pc    = pg - pl ;
    double k_rl  = RelativePermeabilityToLiquid(pc) ;
    double k_rg  = RelativePermeabilityToGas(pc) ;
    double rho_g = x[I_RHO_G] ;
    
    K_L += rho_l*k_int/mu_l*k_rl ;
    K_G += rho_g*k_int/mu_g*k_rg ;
  }
  
  /* Averaging */
  for(i = 0 ; i < NVE ; i++) vex[i] *= 0.5 ;
}




void ComputeFluxes(Element_t *el,double **u)
/* Fluxes */
{
  double *vim = Element_GetImplicitTerm(el) ;
  double *grd = dComponents ;


  {
    double *x1 = ComputeComponents(el,u,vim,0.,1) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] = x1[i] ;
  }
  {
    double *x0 = ComputeComponents(el,u,vim,0.,0) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] -= x0[i] ;
  }
  
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] /= dx ;
  }
  
  /* Fluxes */
  {
    double *w = Fluxes(el,grd) ;

    W_L     = w[I_W_L] ;
    W_G     = w[I_W_G] ;
  }
}


double* Fluxes(Element_t *el,double *grd)
{
  double* vex = Element_GetExplicitTerm(el) ;
  double* w   = ComponentFluxes ;

  /* Gradients */
  double grd_h_l       = grd[I_H_L] ;
  double grd_h_g       = grd[I_H_G] ;
    
    
  /* Flux */
  double w_l = - K_L*grd_h_l ;
  double w_g = - K_G*grd_h_g ;
   
  w[I_W_L ]  = w_l ;
  w[I_W_G ]  = w_g ;
   
  return(w) ;
}


double saturation(double pc,double pc3,Curve_t* curve)
/* Degre de saturation regularise autour de 1 */
{
  double* x = Curve_GetXRange(curve) ;
  double* y = Curve_GetYValue(curve) ;
  double pc1 = x[0] ;
  double sl ;
  
  if(pc >= pc3 || pc1 >= pc3) {
    sl = Curve_ComputeValue(curve,pc) ;
  } else {
    double sl1 = y[0] ;
    double sl3 = Curve_ComputeValue(curve,pc3) ;
    
    sl  = sl1 - (sl1 - sl3) * exp(-(pc - pc3)/(pc1 - pc3)) ;
  }
  
  return(sl) ;
}

double dsaturation(double pc,double pc3,Curve_t* curve)
{
  int   n_i = Curve_GetNbOfPoints(curve) - 1 ;
  double* x = Curve_GetXRange(curve) ;
  double dpc = (x[1] - x[0])/n_i ;
  double sl  = saturation(pc,pc3,curve) ;
  double sl1 = saturation(pc + dpc,pc3,curve) ;
  double dsl = sl1 - sl ;

  return(dsl/dpc) ;
}



int TangentCoefficients1(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *vex = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij   = dist[1] ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  

  /* Initialization */
  for(i = 0 ; i < nn*nn*NEQ*NEQ ; i++) c[i] = 0. ;
    
    
  /* Tangent coefficient */
  for(i = 0 ; i < 2 ; i++) {
    int j = 1 - i ;
    /* Content terms at node i */
    double *cii = c + (i*nn + i)*NEQ*NEQ ;
    /* Transfer terms from node i to node j */
    double *cij = c + (i*nn + j)*NEQ*NEQ ;
    double pl     = P_L(i) ;
    double pg     = P_G(i) ;
    double pc     = pg - pl ;
    double sl     = SATURATION(pc) ;
    double sg     = 1 - sl ;
    double rho_g  = pg*MgsRT ;
    double dslsdpc = dSATURATION(pc) ;
    double dsgsdpc = -dslsdpc ;
    
    /* Content term at node i */
    cii[E_Liq*NEQ + U_P_L] = rho_l*phi*(-dslsdpc) ;
    cii[E_Gas*NEQ + U_P_L] = rho_g*phi*(-dsgsdpc) ;
    
    cii[E_Liq*NEQ + U_P_G] = rho_l*phi*(dslsdpc) ;
    cii[E_Gas*NEQ + U_P_G] = rho_g*phi*(dsgsdpc) + MgsRT*phi*sg ;
    
    
    /* Transfer term from node i to node j */
    cij[E_Liq*NEQ + U_P_L] = dtdij*K_L ;
    cij[E_Gas*NEQ + U_P_L] = 0 ;
    
    cij[E_Liq*NEQ + U_P_G] = 0 ;
    cij[E_Gas*NEQ + U_P_G] = dtdij*K_G ;
  }

  return(dec) ;
}



int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij   = dist[1] ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < nn*nn*NEQ*NEQ ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < 2 ; i++) {
    int j = 1 - i ;
    /* Content terms at node i */
    double *cii = c + (i*nn + i)*NEQ*NEQ ;
    /* Transfer terms from node i to node j */
    double *cij = c + (i*nn + j)*NEQ*NEQ ;
    /* Components */
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    double dxi[NEQ] ;
    int k ;
    
    for(k = 0 ; k < NEQ ; k++) {
      /* dxi[k] =  1.e-1*ObVal_GetValue(obval + k) ; */
      dxi[k] =  1.e2 ;
    }
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      double *dw    = Fluxes(el,dx) ;
    
      cii[E_Liq*NEQ    + k] = dx[I_M_L] ;
      cii[E_Gas*NEQ    + k] = dx[I_M_G] ;
      
      cij[E_Liq*NEQ    + k] = - dtdij*dw[I_W_L] ;
      cij[E_Gas*NEQ    + k] = - dtdij*dw[I_W_G] ;
    }
  }

  return(dec) ;
}




double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *x = Components ;
  
  /* Primary Variables */
  x[I_P_L ] = P_L(n) ;
  x[I_P_G ] = P_G(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_COOR_X] = Element_ComputePointerToNodalCoordinates(el)[n][0] ;
    
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  dx[I_P_L ] = x[I_P_L ] ;
  dx[I_P_G ] = x[I_P_G ] ;
  
  /* Needed variables to compute secondary components */
  dx[I_COOR_X] = x[I_COOR_X] ;
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryComponents(el,dt,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfComponents ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}



void  ComputeSecondaryComponents(Element_t *el,double dt,double *x)
{
  double pl     = x[I_P_L] ;
  double pg     = x[I_P_G] ;
  double z      = x[I_COOR_X];
    
  /* Fluid components */
  double pc     = pg - pl ;
  double sl     = SATURATION(pc) ;
  double sg     = 1 - sl ;
  double rho_g  = pg*MgsRT ;
      
  double m_l    = rho_l*phi*sl ;
  double m_g    = rho_g*phi*sg ;


  /* Back up */
  
  /* Fluid components */
  x[I_S_L      ] = sl ;
  x[I_RHO_G    ] = rho_g ;
  x[I_H_L      ] = pl - rho_l*gravite*z ;
  x[I_H_G      ] = pg - rho_g*gravite*z ;
  
  /* Fluid mass contents */
  x[I_M_L      ] = m_l ;
  x[I_M_G      ] = m_g ;
}
