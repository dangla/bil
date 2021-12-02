#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"

/* Choose the finite element method */
#include "FEM.h"
#include "Elasticity.h"

#define TITLE   "CO2 storage in dual porosity coal seam"
#define AUTHORS "Nikoosokhan-Espinoza-Dangla"

#include "PredefinedMethods.h"

/* Nb of equations of the model */
#define NEQ   (dim+1)


/* Indices of equations */
#define E_Mech   (0)
#define E_CO2    (dim)

/* Indices of unknowns (generic indices) */
#define U_Mech   E_Mech
#define U_CO2    E_CO2


/* Method chosen at compiling time.
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with specific modelings.
 * Uncomment/comment to let only one unknown per equation */
 
/* Mechanics */
#define U_dis      U_Mech

/* Mass balance of CO2 */
#define U_p_co2    U_CO2




/* Nb of terms per point */
#define NVI   (21)    /*  nb of implicit terms per point */
#define NVE   (1)     /*  nb of explicit terms per point */
#define NV0   (0)     /*  nb of constant terms per point */


/* We define some names for implicit terms (vim must be used as pointer below) */
#define N_CO2         (vim)[0]
#define N_CO2n        (vim_n)[0]
#define W_CO2         (vim + 1) /* this is a 3D vector */

#define SIG           (vim + 4) /* this a 3D tensor */

#define F_MASS        (vim + 13)
#define PHI_M         (vim[16])
#define N_CO2_M       (vim[17])
#define N_CO2_m       (vim[18])
#define TRE           (vim[19])
#define PHI_eul       (vim[20])



/* We define some names for explicit terms (vex must be used as pointer below) */
#define K_CO2         (vex[0])



/* We define some names for constant terms (v0 must be used as pointer below) */



/* Shorthands of some units */
#include "InternationalSystemOfUnits.h"

#define m        (InternationalSystemOfUnits_OneMeter)
#define m3       (m*m*m)
#define dm       (0.1*m)
#define cm       (0.01*m)
#define dm3      (dm*dm*dm)
#define cm3      (cm*cm*cm)
#define Pa       (InternationalSystemOfUnits_OnePascal)
#define MPa      (1.e6*Pa)
#define GPa      (1.e9*Pa)
#define Joule    (Pa*m3)
#define kg       (InternationalSystemOfUnits_OneKilogram)
#define mol      (InternationalSystemOfUnits_OneMole)
#define Kelvin   (InternationalSystemOfUnits_OneKelvin)
#define Watt     (InternationalSystemOfUnits_OneWatt)



/* Chosen temperature */
#define TEMPERATURE   293.



/* Physical constants
 * ------------------ */
#include "PhysicalConstant.h"
/* Perfect gas constant (J/mol/K) */
static double R_g ;
static double RT ;



/* Carbon dioxide properties
 * ------------------------- */
/* Molar mass (kg/mol) */
#define M_CO2     (44.e-3 * kg)



/* Material Properties
 * ------------------- */
//#define PermeabilityCoefficient(el,phi)   kozeny_carman(phi0_M,phi)
#define PermeabilityCoefficient(el,sig,p) power10law(el,sig,p)

#define AXIS(I)          (axis_##I)
#define STRESS(i,j)      (sig[3*(AXIS(i)) + (AXIS(j))])
#define STRESS0(i,j)     (sig0[3*(AXIS(i)) + (AXIS(j))])
#define STRAIN(i,j)      (eps[3*(AXIS(i)) + (AXIS(j))])



/* Adsorption properties 
 * --------------------- */
#define CONCENTRATION_ADS(p)   Curve_ComputeValue(adsorbedconcentrationcurve,p)
#define BULK_DENSITY(p)        Curve_ComputeValue(bulkdensitycurve,p)
#define PRESSURE_ADS(p)        Curve_ComputeValue(adsorptionstresscurve,p)
#define TANGENT_BIOT(p)        Curve_ComputeValue(tangentbiotcoefcurve,p)



/* The parameters below are read in the input data file */
static double gravity ;
static double young,poisson,young_3,poisson_3,shear_3 ;
static double b1,b3,N ;
static short int axis_1,axis_2,axis_3 ;
static double K_m ;
static double k_int,mu_co2 ;
static double phi0_M,rho_s ;
static double p0_co2 ;
static double alpha_h,alpha_v ;
static double* cijkl ;
static Elasticity_t* elasty ;
static double* sig0 ;
static Curve_t* adsorbedconcentrationcurve ;
static Curve_t* bulkdensitycurve ;
static Curve_t* adsorptionstresscurve ;
static Curve_t* tangentbiotcoefcurve ;



/* They are stored in the order specified in pm */
static int  pm(const char* s) ;
int pm(const char *s)
{
       if(!strcmp(s,"gravity"))   return (0) ;
  else if(!strcmp(s,"poisson"))   return (1) ;
  else if(!strcmp(s,"K_m"))       return (2) ;
  else if(!strcmp(s,"phi0_M"))    return (3) ;
  else if(!strcmp(s,"p0_co2"))    return (4) ;
  else if(!strcmp(s,"k_int"))     return (5) ;
  else if(!strcmp(s,"mu_co2"))    return (6) ;
  else if(!strcmp(s,"rho_s"))     return (7) ;
  else if(!strcmp(s,"young"))     return (8) ;
  else if(!strcmp(s,"young_3"))   return (9) ;
  else if(!strcmp(s,"poisson_3")) return (10) ;
  else if(!strcmp(s,"shear_3"))   return (11) ;
  else if(!strcmp(s,"axis_3"))    return (12) ;
  else if(!strcmp(s,"alpha_h"))   return (13) ;
  else if(!strcmp(s,"alpha_v"))   return (14) ;
  else if(!strcmp(s,"sig0"))       {
    return(15) ;
  } else if(!strncmp(s,"sig0_",5))   {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(15 + 3*i + j) ;
  } else return(-1) ;
}



/* They are retrieved automatically by calling the following function */
static void    GetProperties(Element_t*) ;
void GetProperties(Element_t* el)
{
/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

  gravity = GetProperty("gravity") ;
  young   = GetProperty("young") ;
  K_m     = GetProperty("K_m") ;
  poisson = GetProperty("poisson") ;
  young_3   = GetProperty("young_3") ;
  poisson_3 = GetProperty("poisson_3") ;
  shear_3   = GetProperty("shear_3") ;
  phi0_M  = GetProperty("phi0_M") ;
  k_int   = GetProperty("k_int") ;
  mu_co2  = GetProperty("mu_co2") ;
  rho_s   = GetProperty("rho_s") ;
  p0_co2  = GetProperty("p0_co2") ;
  sig0    = &GetProperty("sig0") ;
  axis_3  = (short int) GetProperty("axis_3") - 1 ;
  axis_1  = (axis_3 + 1) % 3 ;
  axis_2  = (axis_3 + 2) % 3 ;
  
  alpha_h  = GetProperty("alpha_h") ;
  alpha_v  = GetProperty("alpha_v") ;
  

  {
    double dmu1  = young/(1 + poisson) ;
    double dmu3  = 2*shear_3 ;
    double lamu1 = (1 - poisson)/young - 2*poisson_3*poisson_3/young_3 ;
    double lam1  = 0.5/lamu1 - 0.5*dmu1 ;
    double lam2  = poisson_3/lamu1 ;
    double lam3  = young_3 + 2*poisson_3*lam2 - dmu3 ;
  
    b1 = 1 - (2*lam1 + lam2 + dmu1)/(3*K_m) ;
    b3 = 1 - (2*lam2 + lam3 + dmu3)/(3*K_m) ;
    N  = (2*b1 + b3 - 3*phi0_M)/(3*K_m) ;
  }
  
  
  elasty  = Element_FindMaterialData(el,Elasticity_t,"Elasticity") ;
  cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
  
  adsorbedconcentrationcurve = Element_FindCurve(el,"n_co2") ;
  bulkdensitycurve = Element_FindCurve(el,"rho_co2") ;
  adsorptionstresscurve = Element_FindCurve(el,"stress_ads") ;
  tangentbiotcoefcurve = Element_FindCurve(el,"tangent_biot") ;
  
#undef GetProperty
}



/* Functions used below */
static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;


static int     ComputeTransferCoefficients(Element_t*,double,double*) ;

static int     ComputeTangentCoefficients(Element_t*,double,double,double*) ;

static void    ComputePhysicoChemicalProperties(void) ;

static double kozeny_carman(double,double) ;
static double power10law(Element_t*,double*,double) ;




void ComputePhysicoChemicalProperties(void)
{
/* Physical constants
 * ------------------ */
/* Perfect gas constant (J/mol/K) */
  R_g = PhysicalConstant(PerfectGasConstant) ;
  RT  = R_g * TEMPERATURE ;
}





/* We define some indices for the local variables */
enum {
I_DIS = 0,
I_DIS2 =  I_DIS + 2,

I_U_CO2,

I_N_CO2,
I_N_CO2_M,
I_N_CO2_m,

I_P_CO2,

I_EPS,
I_EPS8       = I_EPS        + 8,

I_SIG,
I_SIG8       = I_SIG        + 8,

I_W_CO2,
I_W_CO22  = I_W_CO2      + 2,

I_GRD_U_CO2,
I_GRD_U_CO22      = I_GRD_U_CO2  + 2,

I_K_CO2,
I_PHI_M,
I_TRE,
I_PHI_eul,
I_F_MASS,
I_Last
} ;



/* Locally defined intern variables  */
//#define NP               IntFct_MaxNbOfIntPoints
#define NbOfVariables    (I_Last)
static double Variable[NbOfVariables] ;
static double Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;




int SetModelProp(Model_t* model)
/** Set the model properties */
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */
  Model_CopyNameOfEquation(model,E_CO2,"co2") ;
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn[i]) ;
  }
  
  /** Names of the main unknowns */
#if defined (U_p_co2)
  Model_CopyNameOfUnknown(model,U_CO2,"p_co2") ;
#else
  #error "Undefined unknown"
#endif
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_dis + i,name_unk[i]) ;
  }
  
  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd */
{
  int NbOfProp = 24 ;
  int dim = Material_GetDimension(mat) ;

  /* By default */
  {
    int i ;
    
    for(i = 0 ; i < NbOfProp ; i++) {
      Material_GetProperty(mat)[i] = 0. ;
    }
  }
  
  ComputePhysicoChemicalProperties() ;
  
  /* Pre-initialization */
  {
    Material_GetProperty(mat)[pm("K_m")]     = 2.e8 ;
    Material_GetProperty(mat)[pm("mu_co2")]  = 1.e-3 ;
    Material_GetProperty(mat)[pm("phi0_M")]  = 0.3;
    Material_GetProperty(mat)[pm("p0_co2")]  = 101325;
    Material_GetProperty(mat)[pm("k_int")]   = 1.e-16 ;
    Material_GetProperty(mat)[pm("young")]   = 0.8e9 ;
    Material_GetProperty(mat)[pm("poisson")] = 0.25 ;
    
    Material_GetProperty(mat)[pm("axis_3")]  = dim ;
  }

  Material_ScanProperties(mat,datafile,pm) ;
  
  /* Create a new curve: the tangent Biot's coefficient 
   * by deriving adsorption induced pressure */
  {
    Curves_t* curves = Material_GetCurves(mat) ;
    Curve_t* curve = Material_GetCurve(mat) ;
    Curve_t* cv ;
    
    if((cv = Curves_FindCurve(curves,"stress_ads"))) {
      if(!Curves_FindCurve(curves,"tangent_biot")) {
        int i = Curves_CreateDerivative(curves,cv) ;
      
        if(i >= 0) {
          Curve_SetNameOfYAxis(curve+i,"tangent_biot") ;
        }
      }
    } else if((cv = Curves_FindCurve(curves,"tangent_biot"))) {
      if(!Curves_FindCurve(curves,"stress_ads")) {
        int i = Curves_CreateIntegral(curves,cv) ;
      
        if(i >= 0) {
          Curve_SetNameOfYAxis(curve+i,"stress_ads") ;
        }
      }
    } else {
      arret("ReadMatProp: curve missed") ;
    }
    
    {
      int NbOfCurves = Material_GetNbOfCurves(mat) ;
      
      if(NbOfCurves > Material_MaxNbOfCurves) {
        arret("ReadMatProp") ;
      }
      
      /* printf("NbOfCurves = %d\n",NbOfCurves) ; */
    }
  }
      
  /* Elasticity */
  {
    elasty = Elasticity_Create() ;
      
    Material_AppendData(mat,1,elasty,Elasticity_t,"Elasticity") ;
  }
  
  /* The 4th rank elastic tensor */
  {
    elasty = Material_FindData(mat,Elasticity_t,"Elasticity") ;

    young     = Material_GetProperty(mat)[pm("young")] ;
    poisson   = Material_GetProperty(mat)[pm("poisson")] ;
    young_3   = Material_GetProperty(mat)[pm("young_3")] ;
    poisson_3 = Material_GetProperty(mat)[pm("poisson_3")] ;
    shear_3   = Material_GetProperty(mat)[pm("shear_3")] ;
    axis_3    = Material_GetProperty(mat)[pm("axis_3")] - 1 ;
    
    
    /* Stiffness tensor */
    {
      /* Isotropic stiffness tensor */
      if(young_3 == 0) {
        Elasticity_SetToIsotropy(elasty) ;
        Elasticity_SetParameters(elasty,young,poisson) ;

        Material_GetProperty(mat)[pm("young_3")]   = young ;
        Material_GetProperty(mat)[pm("poisson_3")] = poisson ;
        Material_GetProperty(mat)[pm("shear_3")]   = young/(2 + 2*poisson) ;
        Material_GetProperty(mat)[pm("axis_3")]    = dim ;
      /* Transversely isotropic stiffness tensor */
      } else {
        Elasticity_SetToTransverselyIsotropy(elasty) ;
        Elasticity_SetParameters(elasty,young,poisson,young_3,poisson_3,shear_3,axis_3) ;
      }
      
      {
        double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
        Elasticity_ComputeStiffnessTensor(elasty,c) ;
      }
    }

#if 0
    {
      Elasticity_PrintStiffnessTensor(elasty) ;
    }
#endif
  }

  return(NbOfProp) ;
}


int PrintModelChar(Model_t* model,FILE *ficd)
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
  
  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 0         # masse density of the dry material\n") ;
  fprintf(ficd,"K_m = 1.04e9      # Compression modulus of the coal matrix\n") ;
  fprintf(ficd,"phi0_M = 0.032    # Macroporosity\n") ;
  fprintf(ficd,"p0_co2 = 1e5      # initial pressure of CO2\n") ;
  fprintf(ficd,"k_int = 2.6e-13   # intrinsic permeability\n") ;
  fprintf(ficd,"mu_co2 = 1.78e-5  # viscosity of CO2\n") ;
  fprintf(ficd,"sig0_ij = 0       # initial stress ij\n") ;
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


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;
  
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
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;
  
  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    
    /* hydraulic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_CO2) {
      for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r1[i] ;
    } else {
      for(i = 0 ; i < NEQ*nn ; i++) r[i] = r1[i] ;
    }
  }
  
  return(0) ;
}


int ComputeInitialState(Element_t* el)
/** Compute the initial state i.e. 
 *  the constant terms,
 *  the explicit terms,
 *  the implicit terms.
 */ 
{
  double* vim0 = Element_GetImplicitTerm(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input Data
  */
  GetProperties(el) ;

 
  /* Validation of initial parameters */
  if(3*phi0_M > 2*b1 + b3) {
    printf("\n") ;
    printf("phi0_M should be less than (2b1 + b3)/3\n") ; 
    printf("phi0_M = %e\n",phi0_M) ;
    printf("b1 = %e , b3 = %e , (2b1 + b3)/3 = %f\n",b1,b3,(2*b1+b3)/3) ;  
    arret ("ComputeInitialState") ;
  }
  
  
  /* Pre-initialization */
  {
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    
      double p_co2 = x[I_P_CO2] ;
      double rho_co2 = BULK_DENSITY(p_co2) ;
      double* sig = x + I_SIG ;
      
      /* transport coefficient */
      double coeff_permeability = PermeabilityCoefficient(el,sig,p_co2) ;
      double k_co2 = rho_co2*k_int/mu_co2*coeff_permeability ;
      
      /* storage in vex */
      {
        double* vex = vex0 + p*NVE ;
        
        K_CO2 = k_co2 ;
      }
    }
  }
  
  /* Loop on integration points */
  {
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
       
      /* storage in vim */
      {
        double* vim  = vim0 + p*NVI ;
        int i ;
        
        N_CO2   = x[I_N_CO2] ;
        N_CO2_M = x[I_N_CO2_M] ;
        N_CO2_m = x[I_N_CO2_m] ;
      
        for(i = 0 ; i < 3 ; i++) W_CO2[i]  = x[I_W_CO2 + i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
        for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_F_MASS + i] ;
      
        PHI_M = x[I_PHI_M] ;
        TRE   = x[I_TRE] ;
        PHI_eul = x[I_PHI_eul] ;
      }
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  double* vim0 = Element_GetPreviousImplicitTerm(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;

  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

   
  /* boucle sur les points d'integration */
  {
    IntFct_t* intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u_n,u_n,vim0,t,0,p) ;
    
      /* pressures */
      double p_co2  = x[I_P_CO2] ;
    
      /* molar density of the bulk */
      double rho_co2 = BULK_DENSITY(p_co2) ;
      
      /* stresses */
      double* sig = x + I_SIG ;
    
      /* transport coefficient */
      double coeff_permeability = PermeabilityCoefficient(el,sig,p_co2) ;
      double k_co2 = rho_co2*k_int/mu_co2*coeff_permeability ;
    
      /* storage in vex */
      {
        double* vex = vex0 + p*NVE ;
        
        K_CO2 = k_co2 ;
      }
    }
  }

  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the implicit terms */
{
  double* vim0  = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

    
  /* Loop on integration points */
  {
    IntFct_t* intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < NbOfIntPoints ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
       
      /* storage in vim */
      {
        double* vim = vim0 + p*NVI ;
        int i ;
        
        N_CO2   = x[I_N_CO2] ;
        N_CO2_M = x[I_N_CO2_M] ;
        N_CO2_m = x[I_N_CO2_m] ;
      
        for(i = 0 ; i < 3 ; i++) W_CO2[i]  = x[I_W_CO2 + i] ;
        for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
        for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_F_MASS + i] ;
      
        PHI_M = x[I_PHI_M] ;
        TRE  = x[I_TRE] ;
        PHI_eul = x[I_PHI_eul] ;
      }
    }
  }

  return(0) ;
}



int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
/** Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;
  double c[IntFct_MaxNbOfIntPoints*100] ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    We load some input data
  */
  GetProperties(el) ;
  
  /*
  ** Poroelastic Matrix
  */
  {
    int dec = ComputeTangentCoefficients(el,t,dt,c) ;
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,1) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {
      k[i] = kp[i] ;
    }
  }
  /*
  ** Conduction Matrix
  */
  {
    int dec = ComputeTransferCoefficients(el,dt,c) ;
    double* kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    int j ;
  
    for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
      K(E_CO2 + i*NEQ,U_CO2 + j*NEQ) += dt*kc[i*nn + j] ;
    }
  }

  return(0) ;
#undef K
}



int  ComputeResidu(Element_t* el,double t,double dt,double* r)
/** Comput the residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  double* vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < nn*NEQ ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;

  /* 1. Mechanics */
  /* 1.1 Stresses */
  {
    double* vim = vim_1 ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    for(i = 0 ; i < nn ; i++) {
      int j ;
      for(j = 0 ; j < dim ; j++) R(i,E_Mech + j) -= rw[i*dim + j] ;
    }
  }
  
  /* 1.2 Body forces */
  {
    double* vim = vim_1 ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_Mech + dim - 1) -= -rbf[i] ;
  }

  /* 2. Hydraulic */
  /* 2.1 Accumulation terms */
  {
    double g1[IntFct_MaxNbOfIntPoints] ;
    double* vim = vim_1 ;
    double* ra ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) {
      g1[i] = N_CO2 - N_CO2n ;
    }
    
    ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_CO2) -= ra[i] ;
  }
  
  /* 2.2 Transport terms */
  {
    double* vim = vim_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_CO2,NVI) ;
    for(i = 0 ; i < nn ; i++) R(i,E_CO2) -= -dt*rf[i] ;
  }
  
  return(0) ;

#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 15 ;
  double* vex  = Element_GetExplicitTerm(el) ;
  double* vim  = Element_GetCurrentImplicitTerm(el) ;
  double** u  = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  /* Initialization */
  {
    int i ;
    
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
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim,t,0,p) ;
    
    /* pressures */
    double p_co2  =  x[I_P_CO2] ;
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
    double w_co2[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double kk0 = 0 ;
    double n_co2 = 0 ;
    double n_co2_M = 0 ;
    double n_co2_m = 0 ;
    double tre = 0 ;
    double phi = 0 ;
    double phi_eul = 0 ;
    int    i ;
      
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      int    j ;
      
      phi += PHI_M/np ;
      
      for(j = 0 ; j < 9 ; j++) sig[j]   += SIG[j]/np ;
      for(j = 0 ; j < 3 ; j++) w_co2[j] += W_CO2[j]/np ;
      
      kk0 += PermeabilityCoefficient(el,SIG,p_co2)/np ;
      n_co2 += N_CO2/np ;
      n_co2_M += N_CO2_M/np ;
      n_co2_m += N_CO2_m/np ;
      tre += TRE/np ;
      phi_eul += PHI_eul/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,x + I_P_CO2,"pressure_co2",1) ;
    Result_Store(r + i++,x + I_DIS,"displacements",3) ;
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



int ComputeTangentCoefficients(Element_t* el,double t,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define C1(i,j,k,l)  (c1[((AXIS(i)*3+AXIS(j))*3+AXIS(k))*3+AXIS(l)])
#define B1(i,j)      (c1[AXIS(i)*3+AXIS(j)])
  double*  vim_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  double dui[NEQ] ;
  int    dec = 100 ;
  int    p ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  {
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dui[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }
  
  
  {
    for(p = 0 ; p < np ; p++) {
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;
      double* c0 = c + p*dec ;


      /* initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0. ;
      }
      
      
      /* The derivative of equations w.r.t unknowns */
      

      /* Derivatives w.r.t strains */
      {
        /* Tangent stiffness matrix */
        {
          double* c1 = c0 ;
          int i ;
          
          for(i = 0 ; i < 81 ; i++) {
            c1[i] = cijkl[i] ;
          }

          //Elasticity_ComputeStiffnessTensor(elasty,c1) ;
        }
        
        /* Coupling matrix */
        {
          double* c1 = c0 + 81 + 9 ;
          
          {
            double p_co2  = x[I_P_CO2] ;
            double rho_co2 = BULK_DENSITY(p_co2) ;
            double b_m = TANGENT_BIOT(p_co2) ;
            
            B1(1,1) = rho_co2*(b_m + b1*(1 - b_m)) ;
            B1(2,2) = rho_co2*(b_m + b1*(1 - b_m)) ;
            B1(3,3) = rho_co2*(b_m + b3*(1 - b_m)) ;
          }
        }
      }
      
      
      /* Derivatives w.r.t U_CO2 */
      {
        double  dp = dui[U_CO2] ;
        double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_U_CO2) ;
        
        /* Tangent Biot's coefficient */
        {
          double* dsigdp = dx + I_SIG ;
          double* c1 = c0 + 81 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* Storage coefficient */
        {
          double* c1 = c0 + 81 + 9 + 9 ;
        
          c1[0] = dx[I_N_CO2] ;
        }
      }
    }
  }
  
  return(dec) ;
  
#undef C1
#undef B1
}


int ComputeTransferCoefficients(Element_t* el,double dt,double* c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  double* vex = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 9 ;
  int    p ;
  
  for(p = 0 ; p < np ; p++ , vex += NVE) {
    double* c1 = c + p*dec ;
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





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
/** This locally defined function compute the intern variables at
 *  the interpolation point p, from the nodal unknowns.
 *  Return a pointer on the locally defined array of the variables. */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ; 
  double*  x   = Variable ;
  double*  x_n = Variable_n ;
  /* Variables is a locally defined array of array */
  //Model_t*  model  = Element_GetModel(el) ;
  //double*   x      = Model_GetVariable(model,p) ;
  
    
  /* Load the primary variables in x */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_DIS + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_Mech + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_DIS + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_Mech) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Pressure */
    x[I_U_CO2] = FEM_ComputeUnknown(fem,u,intfct,p,U_CO2) ;
    
    /* Pressure gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_CO2) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_U_CO2 + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    /* Transfer coefficient */
    {
      double* ve0 = Element_GetExplicitTerm(el) ;
      double* vex  = ve0 + p*NVE ;
      
      x[I_K_CO2]  = K_CO2 ;
    }
  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}






void  ComputeSecondaryVariables(Element_t* el,const double t,const double dt,double* x_n,double* x)
/** Compute the secondary variables from the primary ones. */
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Retrieve the primary variables from x */
  /* Strains */
  double* eps   = x + I_EPS ;
  
  /* pressures */
#if defined (U_p_co2)
  double p_co2  = x[I_U_CO2] ;
  double dp_co2 = p_co2 - p0_co2 ;
#endif
    
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
    
  /* Backup stresses */
  {
    double* sig = x + I_SIG ;
    
    {
      int    i ;
    
      /* Initialize stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig0[i] ;
      
      #define C(i,j)  (cijkl[(i)*9+(j)])
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
      
        for(j = 0 ; j < 9 ; j++) {
          sig[i] += C(i,j)*eps[j] ;
        }
      }
      #undef C
    
      STRESS(1,1) += - b1*(dp_co2 - dp_ads)  - dp_ads ;
      STRESS(2,2) += - b1*(dp_co2 - dp_ads)  - dp_ads ;
      STRESS(3,3) += - b3*(dp_co2 - dp_ads)  - dp_ads ;
    }
  }
  
  {
    double tre  = eps[0] + eps[4] + eps[8] ;
    
    /* Macroporosity */
    double dphi_M = b1*(STRAIN(1,1) + STRAIN(2,2)) + b3*STRAIN(3,3) + N*(dp_co2 - dp_ads) ;
    double phi_M  = phi0_M + dphi_M ;
    
    /* Molar content of CO2 */
    /* ... in micropores */
    double n_co2_m = rho_co2*(phi_m*(1 - phi0_M) + b_m*(tre - dphi_M)) ;
      
    /* ... in macropores */
    double n_co2_M = rho_co2*phi_M ;
      
    /* ... in total */
    double n_co2   = n_co2_M + n_co2_m ;
      
    /* transport coefficient */
    double k_co2   = x[I_K_CO2] ;
      
    /* flux */
    double* grd_co2 = x + I_GRD_U_CO2 ;
    double* w_co2   = x + I_W_CO2 ;
    double* f_mass  = x + I_F_MASS ;
    int i ;
      
    for(i = 0 ; i < 3 ; i++) w_co2[i] = - k_co2*grd_co2[i] ;
    w_co2[dim-1] += k_co2*rho_co2*gravity ;
      
    /* Backup */
    x[I_N_CO2  ] = n_co2 ;
    x[I_N_CO2_M] = n_co2_M ;
    x[I_N_CO2_m] = n_co2_m ;
      
    for(i = 0 ; i < 3 ; i++) f_mass[i] = 0. ;
    f_mass[dim-1] = (rho_s + M_CO2*n_co2)*gravity ;
      
    x[I_PHI_M] = phi_M ;
    x[I_TRE ]  = tre ;
    x[I_PHI_eul] = phi_M/(1 + tre) ; /* Eulerian porosity */
  }
    
  /* Backup */
  x[I_P_CO2   ] = p_co2 ;
  
}



double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{
  double*  x_n = Variable_n ;
  double* dx = dVariable ;
  int j ;
  
  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }
  
  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;
  
  ComputeSecondaryVariables(el,t,dt,x_n,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}




double kozeny_carman(double phi0,double phi)
{
  return(pow(phi/phi0,3)*pow((1-phi0)/(1-phi),2)) ;
  /* linearization of (phi/phi0)^3*((1-phi0)/(1-phi))^2 */
  //return(1 + (3/phi0 + 2/(1 - phi0))*dphi) ;
}

double power10law(Element_t* el,double* sig,double p_co2)
{
  double sig_h  = 0.5*(STRESS(1,1)  + STRESS(2,2) ) + p_co2 ;
  double sig0_h = 0.5*(STRESS0(1,1) + STRESS0(2,2)) + p0_co2 ;
  double sig_v  = STRESS(3,3)  + p_co2 ;
  double sig0_v = STRESS0(3,3) + p0_co2 ;
  double dsig_h = sig_h - sig0_h ;
  double dsig_v = sig_v - sig0_v ;
  
  return(pow(10,alpha_h*dsig_h + alpha_v*dsig_v)) ;
}
