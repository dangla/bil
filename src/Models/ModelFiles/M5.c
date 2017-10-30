#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "FVM.h"


#define TITLE   "Drying-Wetting (1D isothermal case)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"

/* Nb of equations */
#define NEQ      (2)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI      (6)
#define NVE      (4)
#define NV0      (0)


/* Equation index */
#define E_Mass   (0)
#define E_Air    (1)

/* Unknown index */
#define U_P_L    (0)
#define U_P_G    (1)


#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])


/* We define some names for primary unknowns */
#define P_L(n)      (UNKNOWN(n,U_P_L))
#define P_G(n)      (UNKNOWN(n,U_P_G))

#define M_T(n)   (f[(n+0)])
#define M_A(n)   (f[(n+2)])
#define W_A      (f[(4)])
#define W_T      (f[(5)])

#define M_Tn(n)  (f_n[(n+0)])
#define M_An(n)  (f_n[(n+2)])

#define KD_L     (va[(0)])
#define KD_G     (va[(1)])
#define KC_V     (va[(2)])
#define KF_V     (va[(3)])



/* Water properties
 * ---------------- */
/* Molar mass */
#define M_h2o    (18.e-3)
/* Mass density of liquid water (kg/m3) */
#define rho_l    (1000.)
/* Molar volume of liquid water */
#define V_h2o    (18.e-6)
/* Mass density */
#define MassDensityOfWaterVapor(p_v)   (p_v*M_h2o/RT)
/* Vapor-Liquid Equilibrium */
#define RelativeHumidity(p_l)  (exp(V_h2o/RT*(p_l - p_l0)))
#define VaporPressure(p_l)     (p_v0*RelativeHumidity(p_l))
#define LiquidPressure(hr)     (p_l0 + RT/V_h2o*log(hr))


/* Dry air properties
 * ------------------ */
/* Molar mass */
#define M_air    (28.8e-3)
/* Mass density */
#define MassDensityOfDryAir(p_a)       (p_a*M_air/RT)


/* Material properties 
 * ------------------- */
#define SATURATION_CURVE           (Element_GetCurve(el) + 0)
#define RELATIVEPERMLIQ_CURVE      (Element_GetCurve(el) + 1)
#define RELATIVEPERMGAS_CURVE      (Element_GetCurve(el) + 2)

#define Saturation(pc)                   (saturation(pc,p_c3,SATURATION_CURVE))
#define RelativePermeabilityToLiquid(pc) (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,pc))
#define RelativePermeabilityToGas(pc)    (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,pc))
#define TortuosityToGas(f,sg)            (pow(f,aa)*pow(sg,bb))
#define aa                        (1.74)  /* 1/3 Millington, Thiery 1.74 */
#define bb                        (3.2)   /* 7/3 Millington, Thiery 3.2 */


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


/* Functions */
static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;
static int     TangentCoefficients(Element_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;

static double  saturation(double,double,Curve_t*) ;
static double  dsaturation(double,double,Curve_t*) ;


/* Parameters */
static double gravite ;
static double phi ;
static double kl_int ;
static double kg_int ;
static double p_c3 ;

static double p_l0 ;
static double p_g0 ;
static double p_v0 ;

static double RT ;
static double D_av0 ;
static double mu_l,mu_g ;



#include "DiffusionCoefficientOfMoleculeInAir.h"
#include "WaterVaporPressure.h"
#include "AtmosphericPressure.h"
#include "WaterViscosity.h"
#include "AirViscosity.h"
#include "PhysicalConstant.h"

void ComputePhysicoChemicalProperties(double TK)
{
  /* Diffusion Coefficient Of Molecules In Air (m2/s) */
  D_av0   = DiffusionCoefficientOfMoleculeInAir(H2O,TK) ;
  
  /* Water vapor pressure */
  p_v0    = WaterVaporPressure(TK) ;
  
  /* Reference pressures */
  p_l0    = AtmosphericPressure ;
  p_g0    = AtmosphericPressure ;
  
  /* Viscosities */
  mu_l    = WaterViscosity(TK) ;
  mu_g    = AirViscosity(TK) ;
  
  /* Physical constants */
  RT      = PhysicalConstant(PerfectGasConstant)*TK ;
}


#define NN                (2)
#define NbOfComponents    (19)
static double Components[NN][NbOfComponents] ;
static double dComponents[NbOfComponents] ;


#define I_M_L          (2)
#define I_M_V          (3)
#define I_M_A          (4)
#define I_M_G          (5)
#define I_M_T          (6)

#define I_RHO_V        (7)
#define I_RHO_A        (8)
#define I_RHO_G        (9)

#define I_S_L          (10)

#define I_H_L          (11)
#define I_H_G          (12)
#define I_COOR_X       (13)

#define I_P_L          (14)
#define I_P_G          (15)
#define I_P_V          (16)
#define I_P_A          (17)
#define I_C_V          (18)
  
  

#define NbOfComponentFluxes    (5)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_L           (0)
#define I_W_V           (1)
#define I_W_A           (2)
#define I_W_G           (3)
#define I_W_T           (4)



int pm(const char *s)
{
  if(!strcmp(s,"gravite"))     return (0) ;
  else if(!strcmp(s,"phi"))    return (1) ;
  else if(!strcmp(s,"k_int"))  return (2) ;
  else if(!strcmp(s,"kl_int")) return (2) ;
  else if(!strcmp(s,"kg_int")) return (3) ;
  else if(!strcmp(s,"p_c3"))   return (4) ;
  else if(!strcmp(s,"Temperature"))   return (5) ;
  else return(-1) ;
}


void GetProperties(Element_t *el)
{
  gravite = GetProperty("gravite") ;
  phi     = GetProperty("phi") ;
  kl_int  = GetProperty("kl_int") ;
  kg_int  = GetProperty("kg_int") ;
  p_c3    = GetProperty("p_c3") ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_Mass,"mass") ;
  Model_CopyNameOfEquation(model,E_Air,"air") ;

  Model_CopyNameOfUnknown(model,U_P_L,"p_l") ;
  Model_CopyNameOfUnknown(model,U_P_G,"p_g") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  NbOfProp = 6 ;
  
  
  {
    /* Self-initialization */
    Material_GetProperty(mat)[pm("kg_int")] = 0 ;
    Material_GetProperty(mat)[pm("Temperature")] = 293 ;
  }

  Material_ScanProperties(mat,datafile,pm) ;

  {
    kg_int = Material_GetProperty(mat)[pm("kg_int")] ;
    kl_int = Material_GetProperty(mat)[pm("kl_int")] ;

    if(kg_int == 0.) Material_GetProperty(mat)[pm("kg_int")] = kl_int ;
  }

  {
    double temperature = Material_GetProperty(mat)[pm("Temperature")] ;
    
    ComputePhysicoChemicalProperties(temperature) ;
  }
  
  return(NbOfProp) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n") ;
  printf("Equations to be solved are:\n") ;
  printf("\t- Total Mass balance     (mass)\n") ;
  printf("\t- Mass balance of air    (air)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Liquid pressure        (p_l)\n") ;
  printf("\t- Gas pressure           (p_g)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"gravite = -9.81    # Gravity\n") ;
  fprintf(ficd,"phi = 0.3          # Porosity\n") ;
  fprintf(ficd,"kl_int = 1.e-20    # Permeability to liquid water\n") ;
  fprintf(ficd,"kg_int = 1.e-18    # Permeability to air gas\n") ;
  fprintf(ficd,"p_c3 = 1.e6        # Capillary pressure parameter\n") ;
  fprintf(ficd,"Curves = my_file   # File name: p_c S_l k_rl k_rg\n") ;
  
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
  double*  f  = Element_GetImplicitTerm(el) ;
  double** u  = Element_ComputePointerToNodalUnknowns(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  
  
  /* Fluid mass contents */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;

    /* Back up */
    M_A(i) = x[I_M_A] ;
    M_T(i) = x[I_M_T] ;
  }

  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;

  /* Fluxes */
  {
    double* w = ComputeFluxes(el,u) ;

    W_A     = w[I_W_A] ;
    W_T     = w[I_W_T] ;
  }
  
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
  double* f  = Element_GetCurrentImplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* masses d'eau et d'air */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;

    /* Back up */
    M_A(i) = x[I_M_A] ;
    M_T(i) = x[I_M_T] ;
    
    {
      double p_a     = x[I_P_A] ;
    
      if(p_a < 0.) {
        printf("\n") ;
        printf("Air pressure = %e\n",p_a) ;
        return(1) ;
      }
    }
  }
  

  /* Fluxes */
  {
    double* w = ComputeFluxes(el,u) ;

    W_A     = w[I_W_A] ;
    W_T     = w[I_W_T] ;
  }
  
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
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
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
    Conservation of total mass
  */
  R(0,E_Mass)  -= volume[0]*(M_T(0)  - M_Tn(0))  + dt*surf*W_T ;
  R(1,E_Mass)  -= volume[1]*(M_T(1)  - M_Tn(1))  - dt*surf*W_T ;
  
  /*
    Conservation of air mass
  */
  R(0,E_Air)  -= volume[0]*(M_A(0)  - M_An(0))  + dt*surf*W_A ;
  R(1,E_Air)  -= volume[1]*(M_A(1)  - M_An(1))  - dt*surf*W_A ;
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int    nso = 8 ;
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
    double* x    = ComputeComponents(el,u,f,0,j) ;
    double  p_l  = x[I_P_L] ;
    double  p_v  = x[I_P_V] ;
    double  p_g  = x[I_P_G] ;
    double  p_c  = p_g - p_l ;
    double  s_l  = Saturation(p_c) ;
    /* Fluxes */
    double* w    = ComputeFluxes(el,u) ;
    double  w_l  = w[I_W_L] ;
    double  w_v  = w[I_W_V] ;
    double  w_g  = w[I_W_G] ;
    
    double h_r  = RelativeHumidity(p_l) ;
    
    i = 0 ;
    Result_Store(r + i++,&p_l,"p_l",1) ;
    Result_Store(r + i++,&p_g,"p_g",1) ;
    Result_Store(r + i++,&p_v,"p_v",1) ;
    Result_Store(r + i++,&w_l,"w_l",1) ;
    Result_Store(r + i++,&w_g,"w_g",1) ;
    Result_Store(r + i++,&w_v,"w_v",1) ;
    Result_Store(r + i++,&s_l,"saturation",1) ;
    Result_Store(r + i++,&h_r,"humidity",1) ;
  }

  return(nso) ;
}


void ComputeTransferCoefficients(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int i ;
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  
  /* Transfer coefficients */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;
    
    double p_l    = x[I_P_L] ;
    double p_g    = x[I_P_G] ;
    double p_c    = p_g - p_l ;

    double s_l    = x[I_S_L] ;
    double s_g    = 1 - s_l ;

    double rho_g  = x[I_RHO_G] ;
    
    double c_v    = x[I_C_V] ;
  
    double k_rl  = RelativePermeabilityToLiquid(p_c) ;
    double k_rg  = RelativePermeabilityToGas(p_c) ;
  
    double kh_l  = kl_int/mu_l*k_rl ;
    double kh_g  = kg_int/mu_g*k_rg ;
  
    double tau   = TortuosityToGas(phi,s_g) ;
  
    /* Darcy */
    double kd_l  = rho_l*kh_l ;           /* liquid */
    double kd_g  = rho_g*kh_g ;           /* gas */
    
    /* Fick */
    double D_av  = D_av0*p_g0/p_g ;
    double D_eff = phi*s_g*tau*D_av ;
    double kf_v  = D_eff ;                /* vapor */
    
    /* Back up */
    KD_L  += kd_l ;
    KD_G  += kd_g ;
    KF_V  += kf_v ;
    KC_V  += c_v ;
  }
  
  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;
}




double* ComputeFluxes(Element_t *el,double **u)
/* Fluxes */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *grd = dComponents ;


  {
    double *x1 = ComputeComponents(el,u,f,0.,1) ;
    int i ;
    
    for(i = 0 ; i < NbOfComponents ; i++)  grd[i] = x1[i] ;
  }
  {
    double *x0 = ComputeComponents(el,u,f,0.,0) ;
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
    
    return(w) ;
  }
}


double* Fluxes(Element_t *el,double *grd)
{
  double* va = Element_GetExplicitTerm(el) ;
  double* w   = ComponentFluxes ;

  /* Gradients */
  double grd_h_l    = grd[I_H_L] ;
  double grd_h_g    = grd[I_H_G] ;
  double grd_rho_v  = grd[I_RHO_V] ;
    
    
  /* Flux */
  double w_l = - KD_L*grd_h_l ;
  double w_g = - KD_G*grd_h_g ;
  double w_t =   w_l + w_g ;
  double j_v = - KF_V*grd_rho_v ;
  double c_v =   KC_V ;
  double w_v =   c_v*w_g + j_v ;
  double w_a =   w_g - w_v ;
   
   
  w[I_W_L ]  = w_l ;
  w[I_W_G ]  = w_g ;
  w[I_W_V ]  = w_v ;
  w[I_W_A ]  = w_a ;
  w[I_W_T ]  = w_t ;
   
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



int TangentCoefficients2(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij   = dist[1] ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < nn ; i++) {
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
      dxi[k] =  1.e-2*ObVal_GetValue(obval + k) ;
      /* dxi[k] =  1.e2 ; */
    }
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,k) ;
      double *dw    = Fluxes(el,dx) ;
    
      /* Content terms at node i */
      cii[E_Mass*NEQ   + k] = dx[I_M_T] ;
      cii[E_Air*NEQ    + k] = dx[I_M_A] ;
      
      /* Transfer terms from node i to node j */
      cij[E_Mass*NEQ   + k] = - dtdij*dw[I_W_T] ;
      cij[E_Air*NEQ    + k] = - dtdij*dw[I_W_A] ;
    }
  }

  return(dec) ;
}



int TangentCoefficients(Element_t *el,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij   = dist[1] ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  double dxi[NEQ] ;
  int    idi[NEQ] ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
    
  for(i = 0 ; i < NEQ ; i++) {
    dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    /* dxi[i] =  1.e2 ; */
  }
    
  idi[U_P_L  ] = U_P_L ; 
  idi[U_P_G  ] = U_P_G ; 
  
  
  for(i = 0 ; i < nn ; i++) {
    /* Components */
    double *x         = ComputeComponents(el,u,f_n,dt,i) ;
    int k ;
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      int    idk    = idi[k] ;
      double *dx    = ComputeComponentDerivatives(el,dt,x,dxk,idk) ;
      int j ;
      
      for(j = 0 ; j < nn ; j++) {
        if(j == i) {
          /* Content terms at node i */
          double *cii = c + (i*nn + i)*NEQ*NEQ ;
          
          cii[E_Mass*NEQ   + k] = dx[I_M_T] ;
          cii[E_Air*NEQ    + k] = dx[I_M_A] ;
        } else {
          /* Transfer terms from node i to node j */
          double *cij = c + (i*nn + j)*NEQ*NEQ ;
          double *dw  = Fluxes(el,dx) ;
          
          cij[E_Mass*NEQ   + k] = - dtdij*dw[I_W_T] ;
          cij[E_Air*NEQ    + k] = - dtdij*dw[I_W_A] ;
        }
      }
    }
  }

  return(dec) ;
}





double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *x = Components[n] ;
  
  /* Primary Variables */
  x[U_P_L ] = P_L(n) ;
  x[U_P_G ] = P_G(n) ;
  
  /* Needed variables to compute secondary components */
  x[I_COOR_X] = Element_GetNodeCoordinate(el,n)[0] ;
    
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  dx[U_P_L ] = x[U_P_L ] ;
  dx[U_P_G ] = x[U_P_G ] ;
  
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
  double p_l     = x[U_P_L] ;
  double p_g     = x[U_P_G] ;
  double z       = x[I_COOR_X];
    
  /* Fluid components */
  double p_v     = VaporPressure(p_l) ;
  double p_a     = p_g - p_v ;
  double p_c     = p_g - p_l ;

  double s_l     = Saturation(p_c) ;
  double s_g     = 1 - s_l ;

  double rho_v   = MassDensityOfWaterVapor(p_v) ;
  double rho_a   = MassDensityOfDryAir(p_a) ;
  double rho_g   = rho_v + rho_a ;
  
  double c_v     = rho_v/rho_g ;

  double m_l     = rho_l*phi*s_l ;
  double m_g     = rho_g*phi*s_g ;
  double m_v     = rho_v*phi*s_g ;
  double m_a     = rho_a*phi*s_g ;
  double m_t     = m_l + m_g ;
    
    
  /* Back up */
  
  /* Fluid components */
  x[I_P_L      ] = p_l ;
  x[I_P_G      ] = p_g ;
  x[I_P_V      ] = p_v ;
  x[I_P_A      ] = p_a ;
  x[I_S_L      ] = s_l ;
  x[I_RHO_V    ] = rho_v ;
  x[I_RHO_A    ] = rho_a ;
  x[I_RHO_G    ] = rho_g ;
  x[I_H_L      ] = p_l - rho_l*gravite*z ;
  x[I_H_G      ] = p_g - rho_g*gravite*z ;
  x[I_C_V      ] = c_v ;
  
  /* Fluid mass contents */
  x[I_M_L      ] = m_l ;
  x[I_M_V      ] = m_v ;
  x[I_M_A      ] = m_a ;
  x[I_M_G      ] = m_g ;
  x[I_M_T      ] = m_t ;
}
