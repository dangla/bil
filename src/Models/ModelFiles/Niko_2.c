#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"
#include "FEM.h"

#define TITLE   "Dual Porosity for Transversely Isotropic Coal Seam"
#define AUTHORS "Espinoza-Dangla"

#include "PredefinedMethods.h"

/* Macros */
#define NEQ     (1+dim)
#define NVI     (21)
#define NVE     (1)
#define NV0     (0)

#define E_co2   (dim)
#define E_mec   (0)

#define I_p_co2 (dim)
#define I_u     (0)

#define N_CO2   (vim[0])
#define W_CO2   (vim + 1)
#define SIG     (vim + 4)
#define F_MASS  (vim + 13)
#define PHI1    (vim[16])
#define N_CO2_M (vim[17])
#define N_CO2_m (vim[18])
#define TRE     (vim[19])
#define PHI_eul (vim[20])

#define N_CO2n  (vim_n[0])

#define K_CO2   (vex[0])

/* constantes physiques */
#define RT        (2436.)   /* produit de R=8.3143 et T=293 (Pa.m3/mole) */
/* masse molaire (kg/mole) */
#define M_CO2     (44.e-3)

/* Functions */
static int    pm(const char*) ;
static void   GetProperties(Element_t*) ;
static int    c69(FEM_t*,double*) ;
static int    k69(FEM_t*,double*) ;
static double kozeny_carman(double,double) ;
static double power10law(Element_t*,double*,double) ;

/* Curves */
#define CONCENTRATION_ADS(p)   Curve_ComputeValue(Element_GetCurve(el) + concentration_ads,p)
#define BULK_DENSITY(p)        Curve_ComputeValue(Element_GetCurve(el) + bulk_density,p)
#define PRESSURE_ADS(p)        Curve_ComputeValue(Element_GetCurve(el) + pressure_ads,p)
#define TANGENT_BIOT(p)        Curve_ComputeValue(Element_GetCurve(el) + tangent_biot,p)

#define DCONCENTRATION_ADS(p)  Curve_ComputeDerivative(Element_GetCurve(el) + concentration_ads,p)
#define DBULK_DENSITY(p)       Curve_ComputeDerivative(Element_GetCurve(el) + bulk_density,p)
#define DPRESSURE_ADS(p)       Curve_ComputeDerivative(Element_GetCurve(el) + pressure_ads,p)
#define DTANGENT_BIOT(p)       Curve_ComputeDerivative(Element_GetCurve(el) + tangent_biot,p)

#define KOZENY_CARMAN(phi)     kozeny_carman(phi0_1,phi)
#define POWER10LAW(sig,p)      power10law(el,sig,p)

#define AXIS(I)          (axis_##I)
#define STRESS(i,j)      (sig[3*(AXIS(i)) + (AXIS(j))])
#define STRAIN(i,j)      (eps[3*(AXIS(i)) + (AXIS(j))])

/* Parameters */
static double gravite ;
static double young,poisson,young_3,poisson_3,shear_3 ;
static short int axis_1,axis_2,axis_3 ;
static double K_m ;
static double k_int,mu_co2 ;
static double phi0_1,rho_s ;
static double p0_co2,sig0_11,sig0_22,sig0_33 ;
static int    concentration_ads = 0 ;
static int    bulk_density = 1 ;
static int    pressure_ads = 2 ;
static int    tangent_biot ;
static double alpha_h,alpha_v ;

#define GetProperty(a)   (Element_GetProperty(el)[pm(a)]) 

int pm(const char *s)
{
       if(!strcmp(s,"gravite"))   return (0) ;
  else if(!strcmp(s,"poisson"))   return (1) ;
  else if(!strcmp(s,"K_m"))       return (2) ;
  else if(!strcmp(s,"phi0_1"))    return (3) ;
  else if(!strcmp(s,"p0_co2"))    return (4) ;
  else if(!strcmp(s,"k_int"))     return (5) ;
  else if(!strcmp(s,"mu_co2"))    return (6) ;
  else if(!strcmp(s,"rho_s"))     return (7) ;
  else if(!strcmp(s,"sig0_11"))   return (8) ;
  else if(!strcmp(s,"sig0_22"))   return (9) ;
  else if(!strcmp(s,"sig0_33"))   return (10) ;
  else if(!strcmp(s,"young"))     return (11) ;
  else if(!strcmp(s,"young_3"))   return (12) ;
  else if(!strcmp(s,"poisson_3")) return (13) ;
  else if(!strcmp(s,"shear_3"))   return (14) ;
  else if(!strcmp(s,"axis_3"))    return (15) ;
  else if(!strcmp(s,"alpha_h"))   return (16) ;
  else if(!strcmp(s,"alpha_v"))   return (17) ;
  else return(-1) ;
}


void GetProperties(Element_t *el)
{
  gravite = GetProperty("gravite") ;
  young   = GetProperty("young") ;
  K_m     = GetProperty("K_m") ;
  poisson = GetProperty("poisson") ;
  young_3   = GetProperty("young_3") ;
  poisson_3 = GetProperty("poisson_3") ;
  shear_3   = GetProperty("shear_3") ;
  phi0_1  = GetProperty("phi0_1") ;
  k_int   = GetProperty("k_int") ;
  mu_co2  = GetProperty("mu_co2") ;
  rho_s   = GetProperty("rho_s") ;
  p0_co2  = GetProperty("p0_co2") ;
  sig0_11 = GetProperty("sig0_11") ;
  sig0_22 = GetProperty("sig0_22") ;
  sig0_33 = GetProperty("sig0_33") ;
  axis_3  = (short int) GetProperty("axis_3") - 1 ;
  axis_1  = (axis_3 + 1) % 3 ;
  axis_2  = (axis_3 + 2) % 3 ;
}


int SetModelProp(Model_t *model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  int i ;
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_co2,"co2") ;
  for(i = 0 ; i < dim ; i++) {
    char name_eqn[7] ;
    sprintf(name_eqn,"meca_%d",i + 1) ;
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn) ;
  }
  
  /** Names of the main unknowns */
  Model_CopyNameOfUnknown(model,I_p_co2,"p_co2") ;
  for(i = 0 ; i < dim ; i++) {
    char name_unk[4] ;
    sprintf(name_unk,"u_%d",i + 1) ;
    Model_CopyNameOfUnknown(model,I_u + i,name_unk) ;
  }
  
  return(0) ;
}



int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/** Read the material properties in the stream file ficd */
{
  int NbOfProp = 18 ;
  int i ;

  /* Par defaut tout a 0 */
  for(i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;
  { /* initialisation automatique (a completer) SAEID */
    Material_GetProperty(mat)[pm("K_m")]     = 2.e8 ;
    Material_GetProperty(mat)[pm("mu_co2")]  = 1.e-3 ;
    Material_GetProperty(mat)[pm("gravite")] = 0;
    Material_GetProperty(mat)[pm("rho_s")]   = 0;
    Material_GetProperty(mat)[pm("phi0_1")]  = 0.3;
    Material_GetProperty(mat)[pm("p0_co2")]  = 101325;
    Material_GetProperty(mat)[pm("k_int")]   = 1.e-16 ;
    Material_GetProperty(mat)[pm("sig0_11")] = 0 ;
    Material_GetProperty(mat)[pm("sig0_22")] = 0 ;
    Material_GetProperty(mat)[pm("sig0_33")] = 0 ;
    Material_GetProperty(mat)[pm("young")]   = 0.8e9 ;
    Material_GetProperty(mat)[pm("poisson")] = 0.25 ;
  }

  Material_ScanProperties(mat,datafile,pm) ;
  /* Create a new curve: the tangent Biot's coefficient 
   * by deriving adsorption induced pressure */
  {
    Curves_t *curves = Material_GetCurves(mat) ;
    Curve_t *curve = Material_GetCurve(mat) ;
    Curve_t *adpressure_curve = curve + pressure_ads ;
  
    tangent_biot = Curves_CreateDerivative(curves,adpressure_curve) ;
    
    {
      int NbOfCurves = Material_GetNbOfCurves(mat) ;
      
      if(NbOfCurves > Material_MaxNbOfCurves) {
        arret("ReadMatProp") ;
      }
      
      /* printf("NbOfCurves = %d\n",NbOfCurves) ; */
    }
  }

  return(NbOfProp) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/** Print the model characteristics */
{
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- Mass balance equation      (co2)\n") ;
  printf("\t- Mechanical Equilibrium     (meca_[1,2,3])\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- CO2 gas pressure           (p_co2)\n") ;
  printf("\t- Displacements              (u_[1,2,3]) \n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;
  
  fprintf(ficd,"gravite = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 0         # masse density of the dry material\n") ;
  fprintf(ficd,"K_m = 1.04e9      # Compression modulus of the coal matrix\n") ;
  fprintf(ficd,"phi0_1 = 0.032    # Macroporosity\n") ;
  fprintf(ficd,"p0_co2 = 1e5      # initial pressure of CO2\n") ;
  fprintf(ficd,"k_int = 2.6e-13   # intrinsic permeability\n") ;
  fprintf(ficd,"mu_co2 = 1.78e-5  # viscosity of CO2\n") ;
  fprintf(ficd,"sig0_11 = 0       # initial stress 11\n") ;
  fprintf(ficd,"sig0_22 = 0       # initial stress 22\n") ;
  fprintf(ficd,"sig0_33 = 0       # initial stress 33\n") ;
  fprintf(ficd,"young = 1.12e9    # Young's modulus of coal\n") ;
  fprintf(ficd,"poisson = 0.26    # Poisson's ratio of coal\n") ;
  fprintf(ficd,"young_3 = 2.71e9  # Young's modulus of coal in axis 3\n") ;
  fprintf(ficd,"poisson_3 = 0.26 # Poisson's ratio of coal in axis 3\n") ;
  fprintf(ficd,"shear_3 = 4.44e8  # Shear modulus in axis 3\n") ;
  fprintf(ficd,"axis_3 = 2        # Actual axis 3\n") ;
  fprintf(ficd,"alpha_h = 0.09e-6 # Coef of the horizontal stress dependent permeability\n") ;
  fprintf(ficd,"alpha_v = 0.02e-6 # Coef of the vertical stress dependent permeability\n") ;
  fprintf(ficd,"Curves = myfile   # Curves: p , n_ads , rho , ads_pressure\n") ;
  
  return(0) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
/** Define some properties attached to each element */
{
  IntFct_t *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  
  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ;
  
  return(0) ;
}


int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t *fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  
  {
    double *r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    
    /* hydraulic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_co2) {
      for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
    } else {
      for(i = 0 ; i < NEQ*nn ; i++) r[i] = r1[i] ;
    }
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
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double dmu1,dmu3,lam1,lam2,lam3 ;
  double b1,b3,N ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input Data
  */
  GetProperties(el) ;

  dmu1    = young/(1 + poisson) ;
  dmu3    = 2*shear_3 ;
  {
    double lamu1   = (1 - poisson)/young - 2*poisson_3*poisson_3/young_3 ;
    lam1 = 0.5/lamu1 - 0.5*dmu1 ;
    lam2 = poisson_3/lamu1 ;
    if(lam1 < 0) {
      printf("\n") ;
      printf("lambda_1 = %e should be positif\n",lam1) ;
      arret ("ComputeInitialState") ;
    }
  }
  lam3 = young_3 + 2*poisson_3*lam2 - dmu3 ;
  if(lam3 < 0) {
    printf("\n") ;
    printf("lambda_3 = %e should be positif\n",lam3) ;
    arret ("ComputeInitialState") ;
  }
  b1 = 1 - (2*lam1 + lam2 + dmu1)/(3*K_m) ;
  b3 = 1 - (2*lam2 + lam3 + dmu3)/(3*K_m) ;
  N  = (2*b1 + b3 - 3*phi0_1)/(3*K_m) ;
 
  /* Validation of initial parameters */
  if(3*phi0_1 > 2*b1 + b3) {
    printf("\n") ;
    printf("phi0_1 should be less than (2b1 + b3)/3\n") ; 
    printf("phi0_1 = %e\n",phi0_1) ;
    printf("b1 = %e , b3 = %e , (2b1 + b3)/3 = %f\n",b1,b3,(2*b1+b3)/3) ;  
    arret ("ComputeInitialState") ;
  }
    
  /* boucle sur les points d'integration */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vex += NVE) {
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* pressures */
    double p_co2  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_co2) ;
    double dp_co2 = p_co2 - p0_co2 ;
    
    /* adsorption induced pressure */
    double p_ads   = PRESSURE_ADS(p_co2) ;
    double p_ads0  = PRESSURE_ADS(p0_co2) ;
    double dp_ads  = p_ads - p_ads0 ;
    
    /* adsorbed concentration */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    
    /* bulk molar density */
    double rho_co2 = BULK_DENSITY(p_co2) ;
    
    /* apparent microporosity */
    double phi_m = c_co2/rho_co2 ;
    
    /* tangent Biot's coefficient */
    /* double b_m = phi_m*c_simu ; */
    double b_m = TANGENT_BIOT(p_co2) ;
    
    /* Strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
    
    /* Macroporosity */
    double dphi1 = b1*(STRAIN(1,1) + STRAIN(2,2)) + b3*STRAIN(3,3) + N*(dp_co2 - dp_ads) ;
    double phi1  = phi0_1 + dphi1 ;
    
    /* Stresses */
    double sig[9] ;
    
    STRESS(1,1) = sig0_11 + lam1*(STRAIN(1,1) + STRAIN(2,2)) + lam2*STRAIN(3,3) + dmu1*STRAIN(1,1) - b1*(dp_co2 - dp_ads)  - dp_ads ;
    STRESS(2,2) = sig0_22 + lam1*(STRAIN(1,1) + STRAIN(2,2)) + lam2*STRAIN(3,3) + dmu1*STRAIN(2,2) - b1*(dp_co2 - dp_ads)  - dp_ads ;
    STRESS(3,3) = sig0_33 + lam2*(STRAIN(1,1) + STRAIN(2,2)) + lam3*STRAIN(3,3) + dmu3*STRAIN(3,3) - b3*(dp_co2 - dp_ads)  - dp_ads ;
    STRESS(1,2) = dmu1*STRAIN(1,2) ; STRESS(2,1) = STRESS(1,2) ;
    STRESS(1,3) = dmu3*STRAIN(1,3) ; STRESS(3,1) = STRESS(1,3) ;
    STRESS(2,3) = dmu3*STRAIN(2,3) ; STRESS(3,2) = STRESS(2,3) ;
    {
      /* molar contents */
      /* in micropores */
      double n_co2_m = rho_co2*(phi_m*(1 - phi0_1) + b_m*(tre - dphi1)) ;
      
      /* in macropores */
      double n_co2_M = rho_co2*phi1 ;
      
      /* total */
      double n_co2   = n_co2_M + n_co2_m ;
      
      /* transport coefficient */
      /* double coeff_permeability = PermeabilityCoefficient(el,phi) ; */
      /* double coeff_permeability = KOZENY_CARMAN(dphi1) ; */
      double coeff_permeability = POWER10LAW(sig,p_co2) ;
      double k_co2 = rho_co2*k_int/mu_co2*coeff_permeability ;
      
      /* flux */
      double *grd_co2 = FEM_ComputeCurrentUnknownGradient(fem,dh,nn,I_p_co2) ;
      double w_co2[3] ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) w_co2[i] = - k_co2*grd_co2[i] ;
      w_co2[dim-1] += k_co2*rho_co2*gravite ;
      
      /* storage in vim */
      N_CO2 = n_co2 ;
      N_CO2_M = n_co2_M ;
      N_CO2_m = n_co2_m ;
      
      for(i = 0 ; i < 3 ; i++) W_CO2[i]  = w_co2[i] ;
      
      for(i = 0 ; i < 9 ; i++) SIG[i]    = sig[i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = 0. ;
      
      F_MASS[dim-1] = (rho_s + M_CO2*n_co2)*gravite ;
      PHI1 = phi1 ;
      TRE  = tre ;
      PHI_eul = phi1/(1 + tre) ; /* Eulerian porosity */
      
      /* storage in vex */
      K_CO2 = k_co2 ;
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t *el,double t)
/** Compute the explicit terms */
{
  double *vim = Element_GetPreviousImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double dmu1,dmu3,lam1,lam2,lam3 ;
  double b1,b3,N ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  dmu1    = young/(1 + poisson) ;
  dmu3    = 2*shear_3 ;
  {
    double lamu1   = (1 - poisson)/young - 2*poisson_3*poisson_3/young_3 ;
    lam1 = 0.5/lamu1 - 0.5*dmu1 ;
    lam2 = poisson_3/lamu1 ;
  }
  lam3 = young_3 + 2*poisson_3*lam2 - dmu3 ;
  b1 = 1 - (2*lam1 + lam2 + dmu1)/(3*K_m) ;
  b3 = 1 - (2*lam2 + lam3 + dmu3)/(3*K_m) ;
  N  = (2*b1 + b3 - 3*phi0_1)/(3*K_m) ;
   
  /* boucle sur les points d'integration */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vex += NVE) {
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* pressures */
    double p_co2  = FEM_ComputePreviousUnknown(fem,h,nn,I_p_co2) ;
    double dp_co2 = p_co2 - p0_co2 ;
    
    /* adsorption induced pressure */
    double p_ads   = PRESSURE_ADS(p_co2) ;
    double p_ads0  = PRESSURE_ADS(p0_co2) ;
    double dp_ads  = p_ads - p_ads0 ;
    
    /* molar density of the bulk */
    double rho_co2 = BULK_DENSITY(p_co2) ;
    
    /* Strains */
    double *eps = FEM_ComputePreviousLinearStrainTensor(fem,h,dh,nn,I_u) ;
    
    /* Macroporosity */
    double dphi1 = b1*(STRAIN(1,1) + STRAIN(2,2)) + b3*STRAIN(3,3) + N*(dp_co2 - dp_ads) ;
    
    /* transport coefficient */
    /* double coeff_permeability = PermeabilityCoefficient(el,phi) ; */
    /* double coeff_permeability = KOZENY_CARMAN(dphi1) ; */
    double coeff_permeability = POWER10LAW(SIG,p_co2) ;
    double k_co2 = rho_co2*k_int/mu_co2*coeff_permeability ;
    
    /* storage in vex */
    K_CO2 = k_co2 ;
  }
  return(0) ;
}



int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/** Compute the implicit terms */
{
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double dmu1,dmu3,lam1,lam2,lam3 ;
  double b1,b3,N ;
  int    p ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  dmu1    = young/(1 + poisson) ;
  dmu3    = 2*shear_3 ;
  {
    double lamu1   = (1 - poisson)/young - 2*poisson_3*poisson_3/young_3 ;
    lam1 = 0.5/lamu1 - 0.5*dmu1 ;
    lam2 = poisson_3/lamu1 ;
  }
  lam3 = young_3 + 2*poisson_3*lam2 - dmu3 ;
  b1 = 1 - (2*lam1 + lam2 + dmu1)/(3*K_m) ;
  b3 = 1 - (2*lam2 + lam3 + dmu3)/(3*K_m) ;
  N  = (2*b1 + b3 - 3*phi0_1)/(3*K_m) ;
    
  /* boucle sur les points d'integration */
  for(p = 0 ; p < NbOfIntPoints ; p++ , vim += NVI , vex += NVE) {
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* pressures */
    double p_co2  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_co2) ;
    double dp_co2 = p_co2 - p0_co2 ;
    
    /* adsorption induced pressure */
    double p_ads   = PRESSURE_ADS(p_co2) ;
    double p_ads0  = PRESSURE_ADS(p0_co2) ;
    double dp_ads  = p_ads - p_ads0 ;
    
    /* adsorbed concentration */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    
    /* bulk molar density */
    double rho_co2 = BULK_DENSITY(p_co2) ;
    
    /* apparent microporosity */
    double phi_m = c_co2/rho_co2 ;
    
    /* tangent Biot's coefficient */
    /* double b_m = phi_m*c_simu ; */
    double b_m = TANGENT_BIOT(p_co2) ;
    
    /* Strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
    
    /* Macroporosity */
    double dphi1 = b1*(STRAIN(1,1) + STRAIN(2,2)) + b3*STRAIN(3,3) + N*(dp_co2 - dp_ads) ;
    double phi1  = phi0_1 + dphi1 ;
    
    /* Stresses */
    double sig[9] ;
    
    STRESS(1,1) = sig0_11 ;
    STRESS(2,2) = sig0_22 ;
    STRESS(3,3) = sig0_33 ;
    
    STRESS(1,1) += lam1*(STRAIN(1,1) + STRAIN(2,2)) + lam2*STRAIN(3,3) + dmu1*STRAIN(1,1) ;
    STRESS(2,2) += lam1*(STRAIN(1,1) + STRAIN(2,2)) + lam2*STRAIN(3,3) + dmu1*STRAIN(2,2) ;
    STRESS(3,3) += lam2*(STRAIN(1,1) + STRAIN(2,2)) + lam3*STRAIN(3,3) + dmu3*STRAIN(3,3) ;
    STRESS(1,2) = dmu1*STRAIN(1,2) ; STRESS(2,1) = STRESS(1,2) ;
    STRESS(1,3) = dmu3*STRAIN(1,3) ; STRESS(3,1) = STRESS(1,3) ;
    STRESS(2,3) = dmu3*STRAIN(2,3) ; STRESS(3,2) = STRESS(2,3) ;
    
    STRESS(1,1) += - b1*(dp_co2 - dp_ads)  - dp_ads ;
    STRESS(2,2) += - b1*(dp_co2 - dp_ads)  - dp_ads ;
    STRESS(3,3) += - b3*(dp_co2 - dp_ads)  - dp_ads ;
    {
      /* molar content */
      /* in micropores */
      double n_co2_m = rho_co2*(phi_m*(1 - phi0_1) + b_m*(tre - dphi1)) ;
      
      /* in macropores */
      double n_co2_M = rho_co2*phi1 ;
      
      /* total */
      double n_co2   = n_co2_M + n_co2_m ;
      
      /* transport coefficient */
      double k_co2   = K_CO2 ;
      
      /* flux */
      double *grd_co2 = FEM_ComputeCurrentUnknownGradient(fem,dh,nn,I_p_co2) ;
      double w_co2[3] ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) w_co2[i] = - k_co2*grd_co2[i] ;
      w_co2[dim-1] += k_co2*rho_co2*gravite ;
      
      /* storage in vim */
      N_CO2 = n_co2 ;
      N_CO2_M = n_co2_M ;
      N_CO2_m = n_co2_m ;
      
      for(i = 0 ; i < 3 ; i++) W_CO2[i]  = w_co2[i] ;
      
      for(i = 0 ; i < 9 ; i++) SIG[i]    = sig[i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = 0. ;
      F_MASS[dim-1] = (rho_s + M_CO2*n_co2)*gravite ;
      
      PHI1 = phi1 ;
      TRE  = tre ;
      PHI_eul = phi1/(1 + tre) ; /* Eulerian porosity */
    }
    
    if(p_co2 < 0) {
      printf("p_co2 = %e\n",p_co2) ;
      return(-1) ;
    }
  }

  return(0) ;
}



int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  double c[IntFct_MaxNbOfIntPoints*100] ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
  ** Poroelastic Matrix
  */
  {
    int dec = c69(fem,c) ;
    double *kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  /*
  ** Conduction Matrix
  */
  {
    int dec = k69(fem,c) ;
    double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    int j ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_co2 + i*NEQ,I_p_co2 + j*NEQ) += dt*kc[i*nn + j] ;
    }
  }

  return(0) ;
#undef K
}



int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double *vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double *vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  double *vim = vim_1 ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < nn*NEQ ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /* 1. Mechanics */
  /* 1.1 Stresses */
  {
    double *rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    for(i = 0 ; i < nn ; i++) {
      int j ;
      for(j = 0 ; j < dim ; j++) R(i,E_mec + j) -= rw[i*dim + j] ;
    }
  }
  
  /* 1.2 Body forces */
  {
    double *rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_mec + dim - 1) -= -rbf[i] ;
  }

  /* 2. Hydraulic */
  /* 2.1 Accumulation terms */
  {
    double g1[IntFct_MaxNbOfIntPoints] ;
    double *ra ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) {
      g1[i] = N_CO2 - N_CO2n ;
    }
    
    ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_co2) -= ra[i] ;
  }
  
  /* 2.2 Transport terms */
  vim = vim_1 ;
  {
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_CO2,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_co2) -= -dt*rf[i] ;
  }
  
  return(0) ;

#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 15 ;
  double *vim = Element_GetImplicitTerm(el) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;

  {
    /* Interpolation functions at s */
    double *h_s = Element_ComputeIsoShapeFctInActualSpace(el,s) ;
    /* pressures */
    double p_co2  =  FEM_ComputeCurrentUnknown(fem,h_s,nn,I_p_co2) ;
    /* adsorbed concentration */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    /* molar density of the bulk */
    double rho_co2 = BULK_DENSITY(p_co2) ;
    /* adsorption induced pressure */
    double p_ads   = PRESSURE_ADS(p_co2) ;
    /* apparent microporosity */
    double phi_m = c_co2/rho_co2 ;
    /* tangent Biot's coefficient */
    double b_m = TANGENT_BIOT(p_co2) ;
    
    /* coefficient  */
    double c_simu = b_m/phi_m ;
    double u[3] = {0,0,0} ;
    double w_co2[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double kk0 = 0 ;
    double n_co2 = 0 ;
    double n_co2_M = 0 ;
    double n_co2_m = 0 ;
    double tre = 0 ;
    double phi = 0 ;
    double phi_eul = 0 ;
    int    i,p ;
    
    for(i = 0 ; i < dim ; i++) {
      u[i] = FEM_ComputeCurrentUnknown(fem,h_s,nn,I_u + i) ;
    }
      
    /* Averaging */
    for(p = 0 ; p < np ; p++ , vim += NVI , vex += NVE) {
      phi += PHI1/np ;
      for(i = 0 ; i < 9 ; i++) sig[i] += SIG[i]/np ;
      for(i = 0 ; i < 3 ; i++) w_co2[i] += W_CO2[i]/np ;
      /* kk0 += K_CO2/rho_co2/k_int*mu_co2/np ; */
      kk0 += POWER10LAW(SIG,p_co2)/np ;
      n_co2 += N_CO2/np ;
      n_co2_M += N_CO2_M/np ;
      n_co2_m += N_CO2_m/np ;
      tre += TRE/np ;
      phi_eul += PHI_eul/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&p_co2,"pressure_co2",1) ;
    Result_Store(r + i++,u,"displacements",3) ;
    Result_Store(r + i++,&c_co2,"adsorbed_concentration",1) ;
    Result_Store(r + i++,&p_ads,"adsorption-induced_pressure",1) ;
    Result_Store(r + i++,&rho_co2,"molar_density_co2",1) ;
    Result_Store(r + i++,&c_simu,"coeff",1) ;
    Result_Store(r + i++,&phi,"Macroporosity",1) ;
    Result_Store(r + i++,sig,"Stress_tensor",9) ;
    Result_Store(r + i++,w_co2,"flux",3) ;
    Result_Store(r + i++,&kk0,"permeability",1) ;
    Result_Store(r + i++,&n_co2,"total_molar_content",1) ;
    Result_Store(r + i++,&n_co2_M,"N_CO2_M",1) ;
    Result_Store(r + i++,&n_co2_m,"N_CO2_m",1) ;
    Result_Store(r + i++,&tre,"Volumetric_strain",1) ;
    Result_Store(r + i++,&phi_eul,"Eulerian_porosity",1) ;
      
    if(i != NbOfOutputs) arret("ComputeOutputs") ;
  }

  return(NbOfOutputs) ;
}

int c69(FEM_t *fem,double *c)
/*
**  Poroelastic matrix (c), return the shift (dec)
*/
{
#define C1(i,j,k,l)  (c1[((AXIS(i)*3+AXIS(j))*3+AXIS(k))*3+AXIS(l)])
#define B1(i,j)      (c1[AXIS(i)*3+AXIS(j)])
  Element_t *el = FEM_GetElement(fem) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int nn = IntFct_GetNbOfNodes(intfct) ;
  double dmu1,dmu3,lam1,lam2,lam3,mu1,mu3 ;
  double b1,b3,N ;
  int    dec = 100 ;
  int    p ;
  
  /*
    Donnees
  */
  young   = GetProperty("young") ;
  K_m     = GetProperty("K_m") ;
  poisson = GetProperty("poisson") ;
  young_3   = GetProperty("young_3") ;
  poisson_3 = GetProperty("poisson_3") ;
  shear_3   = GetProperty("shear_3") ;
  phi0_1  = GetProperty("phi0_1") ;
  p0_co2  = GetProperty("p0_co2") ;
  axis_3  = (short int) GetProperty("axis_3") - 1 ;
  axis_1  = (axis_3 + 1) % 3 ;
  axis_2  = (axis_3 + 2) % 3 ;

  dmu1    = young/(1 + poisson) ;
  dmu3    = 2*shear_3 ;
  mu1 = 0.5*dmu1 ;
  mu3 = 0.5*dmu3 ;
  {
    double lamu1   = (1 - poisson)/young - 2*poisson_3*poisson_3/young_3 ;
    lam1 = 0.5/lamu1 - 0.5*dmu1 ;
    lam2 = poisson_3/lamu1 ;
  }
  lam3 = young_3 + 2*poisson_3*lam2 - dmu3 ;
  b1 = 1 - (2*lam1 + lam2 + dmu1)/(3*K_m) ;
  b3 = 1 - (2*lam2 + lam3 + dmu3)/(3*K_m) ;
  N  = (2*b1 + b3 - 3*phi0_1)/(3*K_m) ;
   
  for(p = 0 ; p < np ; p++) {
    int    i ;
    double *c1 = c + p*dec ;
    
    /* interpolation functions */
    double *h  = IntFct_GetFunctionAtPoint(intfct,p) ;
    double *dh = IntFct_GetFunctionGradientAtPoint(intfct,p) ;
    
    /* pressures */
    double p_co2  = FEM_ComputeCurrentUnknown(fem,h,nn,I_p_co2) ;
    double dp_co2 = p_co2 - p0_co2 ;
    
    /* adsorption induced pressure */
    double p_ads   = PRESSURE_ADS(p_co2) ;
    double p_ads0  = PRESSURE_ADS(p0_co2) ;
    double dp_ads  = p_ads - p_ads0 ;
    
    /* molar density of the bulk */
    double rho_co2 = BULK_DENSITY(p_co2) ;
    
    /* adsorbed concentration */
    double c_co2   = CONCENTRATION_ADS(p_co2) ;
    
    /* apparent microporosity */
    double phi_m = c_co2/rho_co2 ;
    
    /* tangent Biot's coefficient */
    /* double b_m = phi_m*c_simu ; */
    double b_m = TANGENT_BIOT(p_co2) ;
    
    /* Strains */
    double *eps =  FEM_ComputeCurrentLinearStrainTensor(fem,h,dh,nn,I_u) ;
    double tre  = eps[0] + eps[4] + eps[8] ;
    
    /* Macroporosity */
    double dphi1 = b1*(STRAIN(1,1) + STRAIN(2,2)) + b3*STRAIN(3,3) + N*(dp_co2 - dp_ads) ;
    double phi1  = phi0_1 + dphi1 ;
     
    /* derivatives */
    /* ... with respect to p_co2 */
    double dp_adssdp_co2   = DPRESSURE_ADS(p_co2) ;
    double drho_co2sdp_co2 = DBULK_DENSITY(p_co2) ;
    double dphi1sdp_co2    = N*(1. - dp_adssdp_co2)  ;
    double dc_co2sdp_co2   = DCONCENTRATION_ADS(p_co2) ;
    double db_msdp_co2     = DTANGENT_BIOT(p_co2) ;
    double coef            = rho_co2*b_m ;
    /* double dcoefsdp_co2    = dc_co2sdp_co2*c_simu + c_co2*dc_simusdp_co2 ; */
    double dcoefsdp_co2    = drho_co2sdp_co2*b_m + rho_co2*db_msdp_co2 ;
    double dn_co2_Msdp_co2 = drho_co2sdp_co2*phi1 + rho_co2*dphi1sdp_co2 ;
    double dn_co2_msdp_co2 = dc_co2sdp_co2*(1 - phi0_1) + dcoefsdp_co2*(tre - dphi1) - coef*dphi1sdp_co2 ;
      
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
      
    /* Mechanics */
    /* derivative of sig with respect to eps */
    C1(1,1,1,1) = lam1 + dmu1 ; C1(1,1,2,2) = lam1 ; C1(1,1,3,3) = lam2 ;
    C1(2,2,1,1) = lam1 ; C1(2,2,2,2) = lam1 + dmu1 ; C1(2,2,3,3) = lam2 ;
    C1(3,3,1,1) = lam2 ; C1(3,3,2,2) = lam2 ; C1(3,3,3,3) = lam3 + dmu3 ;
    C1(1,2,1,2) = mu1 ; C1(1,2,2,1) = mu1 ; C1(2,1,1,2) = mu1 ; C1(2,1,2,1) = mu1 ;
    C1(1,3,1,3) = mu3 ; C1(1,3,3,1) = mu3 ; C1(3,1,1,3) = mu3 ; C1(3,1,3,1) = mu3 ;
    C1(2,3,2,3) = mu3 ; C1(2,3,3,2) = mu3 ; C1(3,2,2,3) = mu3 ; C1(3,2,3,2) = mu3 ;
    
    /* derivative of sig with respect to p_co2 */
    c1 += 81 ;
    B1(1,1) = - b1 -  (1 - b1)*dp_adssdp_co2 ;
    B1(2,2) = - b1 -  (1 - b1)*dp_adssdp_co2 ;
    B1(3,3) = - b3 -  (1 - b3)*dp_adssdp_co2 ;
      
    /* Hydraulic */
    /* derivative of n_co2 with respect to eps */
    c1 += 9 ;
    B1(1,1) = rho_co2*(b_m + b1*(1 - b_m)) ;
    B1(2,2) = rho_co2*(b_m + b1*(1 - b_m)) ;
    B1(3,3) = rho_co2*(b_m + b3*(1 - b_m)) ;
    
    /* derivative of n_co2 with respect to p_co2 */
    c1 += 9 ;
    c1[0] = dn_co2_Msdp_co2 + dn_co2_msdp_co2 ;
  }
  
  return(dec) ;
  
#undef C1
#undef B1
}


int k69(FEM_t *fem,double *c)
/*
**  Conduction matrix (c), return the shift (dec)
*/
{
  Element_t *el = FEM_GetElement(fem) ;
  double *vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  
  for(p = 0 ; p < np ; p++ , vex += NVE) {
    double *c1 = c + p*dec ;
    int    i ;
    
    /* initialisation */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
    
    /* Permeability tensor */
    c1[0] = K_CO2 ;
    c1[4] = K_CO2 ;
    c1[8] = K_CO2 ;
  }
  return(dec) ;
}


double kozeny_carman(double phi0,double dphi)
{
  /* linearization of (phi/phi0)^3*((1-phi0)/(1-phi))^2 */
  return(1 + (3/phi0 + 2/(1 - phi0))*dphi) ;
}

double power10law(Element_t *el,double *sig,double p_co2)
{
  double sig_h  = 0.5*(STRESS(1,1) + STRESS(2,2)) + p_co2 ;
  double sig0_h = 0.5*(sig0_11     + sig0_22    ) + p0_co2 ;
  double sig_v  = STRESS(3,3) + p_co2 ;
  double sig0_v = sig0_33     + p0_co2 ;
  double dsig_h = sig_h - sig0_h ;
  double dsig_v = sig_v - sig0_v ;
  
  alpha_h  = GetProperty("alpha_h") ;
  alpha_v  = GetProperty("alpha_v") ;
  
  return(pow(10,alpha_h*dsig_h + alpha_v*dsig_v)) ;
}
