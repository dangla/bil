#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
//#include <ctype.h>
#include "CommonModel.h"
#include "FEM.h"
#include "Elasticity.h"


#define TITLE   "Crystallization in porous media"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"



/* Indices of equations/unknowns */
/*
enum {
  E_Mech      ,
  E_Mass      ,
  E_Salt      ,
  E_Last
} ;
*/
#define  E_Mech      (0)
#define  E_Mass      (dim)
#define  E_Salt      (dim+1)
#define  E_Last(dim) ((dim)+2)


/* Nb of equations */
#define NbOfEquations      E_Last
#define NEQ                NbOfEquations(dim)





/* Method chosen at compiling time. 
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */


/* Mechanics: the displacements */
#define U_Dis       E_Mech

/* Mass: 
 * - the liquid pressure 
 * - The relative humidity */
#define U_P_L       E_Mass
//#define U_H_r       E_Mass


/* Salt: the molar concentration */
#define U_C_s       E_Salt




/* Nb of nodes (el must be used below) */
#define NbOfNodes               Element_GetNbOfNodes(el)


/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (18)
#define NVE     (4)
#define NV0     (9)



/* Names used for implicit terms */
/* Mechanics */
#define SIG           (vim + 0)
#define SIG_n         (vim_n + 0)
#define F_MASS        (vim + 9)

/* Total mass balance */
#define W_Tot         (vim + 10)
#define M_Tot         (vim + 13)[0]
#define M_Totn        (vim_n + 13)[0]

/* Salt mole balance */
#define W_Salt        (vim + 14)
#define N_Salt        (vim + 17)[0]
#define N_Saltn       (vim_n + 17)[0]



/* Names used for explicit terms */
#define KD_L          (vex + 0)[0]

#define KF_Vap        (vex + 1)[0]
#define KF_Salt       (vex + 2)[0]

#define KC_Salt       (vex + 3)[0]


/* We define some names for constant terms */
#define SIG0          (v0  + 0)






/* Water properties
 * ---------------- */
/* Molar mass */
#define M_H2O     (18.e-3)
/* Mass density of liquid water (kg/m3) */
#define rho_l0    (1000.)
/* Molar volume of liquid water */
#define V_H2O     (18.e-6)
/* Mass density */
#define MassDensityOfWaterVapor(p_v)   (p_v*M_H2O/RT)
/* Vapor-Liquid Equilibrium */
#define RelativeHumidity(p_l,lnaw)  (exp(V_H2O/RT*(p_l - p_l0)) + lnaw)
#define VaporPressure(p_l,lnaw)     (p_v0*RelativeHumidity(p_l,lnaw))
#define LiquidPressure(hr,lnaw)     (p_l0 + RT/V_H2O*(log(hr) - lnaw))
/* Logarithm of water activity in brine */
//#define LogActivityOfWater(c_s)     activite_w(c_s,TEMPERATURE)
#define LogActivityOfWater(c_s)     activite_w_ideal(c_s)



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
#define TORTUOSITYLIQ_CURVE        (Element_GetCurve(el) + 3)

#define SaturationDegree(pc)             (Curve_ComputeValue(SATURATION_CURVE,pc))
#define RelativePermeabilityToLiquid(sl) (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,sl))
#define RelativePermeabilityToGas(sg)    (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,sg))
#define TortuosityToGas(f,sg)            ((sg > 0) ? pow(f,aa)*pow(sg,bb) : 0)
#define aa                        (0.33)  /* 1/3 Millington, Thiery 1.74 */
#define bb                        (2.33)   /* 7/3 Millington, Thiery 3.2 */
//#define TortuosityToLiquid(f,sl)         (tortuosite_l(f)*pow(sl,4.5))
#define TortuosityToLiquid(f,sl)         ((sl > 0) ? (pow(f,0.33)*pow(sl,2.33)) : 0)
/*
#define TortuosityToLiquid(f,sl)         (0.00029*exp(9.95*f)/(1+625*((sl < 1) ? pow((1-sl),4) : 0)))
*/
    



#define GAMMA_LG        (0.07)


/* Types of salt */
#define NaCl    (0)
#define Na2SO4  (1)
#define KNO3    (2)

/* The chosen salt */
#define SALT     Na2SO4

/* Salt properties
 * --------------- */
#if SALT == NaCl
  /* Valencies */
  #define Z_A       (-1)
  #define Z_C       (+1)
  /* Stoichiometric coefficients */
  #define NU_A      1
  #define NU_C      1
  #define NU_H2O    0
  /* Partial molar volumes in liquid (m3/mole) */
  #define V_A       (2.52e-6)
  #define V_C       (1.87e-5)
  /* Solid molar volume of salt (m3/mole) */
  #define V_SALT    (24.5e-6)
  /* Molar mass */
  #define M_SALT    (58.45e-3)
  /* Molecular diffusion coefficients (m2/s) */
  #define D_A       (2.032e-9)
  #define D_C       (1.334e-9)
  /* Solubility (moles/m3) */
  #define K_SALT    (6.e3)
  /* Surface tensions (N/m) */
  #define GAMMA_CL  (0.1)
#elif SALT == Na2SO4
  #define Z_A       (-2)
  #define Z_C       (+1)
  #define NU_A      1
  #define NU_C      2
  #define NU_H2O    10
  #define V_A       (-8.334e-6)
  #define V_C       (1.87e-5)
  #define V_SALT    (220.e-6)
  #define M_SALT    (142.e-3)
  #define D_A       (1.065e-9)
  #define D_C       (1.334e-9)
  #define K_SALT    (1.24166e3)
  #define GAMMA_CL  (0.1)
#else
  #error "Salt not available"
#endif


#define NU_AC     (NU_A + NU_C)
#define V_AC      (NU_A*V_A + NU_C*V_C)
#define V_ACH     (NU_H2O*V_H2O + V_AC)
#define D_SALT    (Z_C - Z_A)*D_A*D_C/(Z_C*D_C - Z_A*D_A)
/* Logarithm of salt activity in brine  */
//#define LogActivityOfSalt(c_s)     activite_s(c_s,TEMPERATURE)
#define LogActivityOfSalt(c_s)     activite_s_ideal(c_s)

#define SolubilityProductOfSalt(p_l) \
        (K_SALT*exp(-(V_SALT-V_ACH)/RT*((p_l)-p0)))
       
#define SaturationIndexOfSalt(lnas,p_l) \
        (exp(lnas)/SolubilityProductOfSalt(p_l))
        
#define CrystallizationPressure(beta_s,p_l) \
        (p_l + RT*log(beta_s)/V_SALT)
        
#define ReactionGibbsEnergy(beta_s,p_l,kappa) \
        (-V_C*GAMMA_CL*kappa + (V_SALT-V_ACH)*((p_l)-p0) + RT*ln(beta_s))



#define TEMPERATURE         (293.)      /* Temperature (K) */


/* To retrieve the material properties */
#define GetProperty(a)      Element_GetPropertyValue(el,a)


/* Fonctions */
static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double,double*,double*) ;


static void    ComputePhysicoChemicalProperties(double) ;

static double  lna_i(double,double,double,double,double,double) ;
static double  lng_TQN(double,double,double,double,double,double,double,double) ;
static double  lng_LinLee(double,double,double,double,double,double) ;

static double  tortuosite_l(double) ;

static double  activite_w(double,double) ;
static double  activite_s(double,double) ;
static double  activite_w_ideal(double) ;
static double  activite_s_ideal(double) ;



/* Parameters */
static double* sig0 ;
static double  biot ;
static double  N ;

static double gravity ;
static double phi0 ;
static double kl_int ;

static double p0 = 0 ;

static double p_l0 ;
static double p_g0 = 0 ;
static double p_v0 ;

static double RT ;
static double D_av0 ;
static double mu_l ;
static double* cijkl ;
static Elasticity_t* elasty ;



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
  
  /* Physical constants */
  RT      = PhysicalConstant(PerfectGasConstant)*TK ;
}







enum {
I_U   = 5,
I_H_r = I_U + 3,
I_P_L,
I_P_C,
I_PP,
I_C_s,
I_EPS,
I_SIG = I_EPS + 9,
I_Fmass = I_SIG + 9,
I_W_Tot = I_Fmass + 3,
I_M_Tot = I_W_Tot + 3,
I_W_L,
I_M_L = I_W_L + 3,
I_W_Vap,
I_M_Vap = I_W_Vap + 3,
I_W_Salt,
I_N_Salt = I_W_Salt + 3,
I_GRD_Mass,
I_GRD_C_s = I_GRD_Mass + 3,
I_KD_L = I_GRD_C_s + 3,
I_KF_Vap,
I_KF_Salt,
I_KC_Salt,
I_S_L,
I_S_G,
I_S_C,
I_RHO_L,
I_Last
} ;


#define NbOfVariables    (I_Last)
static double Variable[2*NbOfVariables] ;
static double dVariable[NbOfVariables] ;
#define Variable_n(x)   ((x) + NbOfVariables)






int pm(const char *s)
{
         if(!strcmp(s,"gravity"))    { return (0) ;
  } else if(!strcmp(s,"young"))      { return (1) ;
  } else if(!strcmp(s,"poisson"))    { return (2) ;
  } else if(!strcmp(s,"porosity"))   { return (3) ;
  } else if(!strcmp(s,"rho_l"))      { return (4) ;
  } else if(!strcmp(s,"k_int"))      { return (5) ;
  } else if(!strcmp(s,"kl_int"))     { return (5) ;
  } else if(!strcmp(s,"mu_l"))       { return (6) ;
  } else if(!strcmp(s,"biot"))       { return (7) ;
  } else if(!strcmp(s,"N"))          { return (8) ;
  } else if(!strcmp(s,"sig0"))       {
    return(13) ;
  } else if(!strncmp(s,"sig0_",5))   {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    
    return(13 + 3*i + j) ;
  } else return (-1) ;
}


void GetProperties(Element_t* el)
{
  gravity = GetProperty("gravity") ;
  phi0    = GetProperty("porosity") ;
  kl_int  = GetProperty("kl_int") ;
  mu_l    = GetProperty("mu_l") ;
  //k_l     = GetProperty("k_l") ;
  //rho_s   = GetProperty("rho_s") ;
  p_l0    = GetProperty("p_l0") ;
  //b       = GetProperty("biot") ;
  N       = GetProperty("N") ;
  sig0    = &GetProperty("sig0") ;
  
  elasty  = (Elasticity_t*) Element_FindMaterialData(el,Elasticity_t,"Elasticity") ;
  cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
}


int SetModelProp(Model_t *model)
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"mech_1","mech_2","mech_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  

  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  
  
  /** Names of these equations */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn[i]) ;
  }
  Model_CopyNameOfEquation(model,E_Mass,"mass") ;
  Model_CopyNameOfEquation(model,E_Salt,"salt") ;



  /** Names of the main (nodal) unknowns */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,E_Mech + i,name_unk[i]) ;
  }
  
#if defined (U_H_r)
  Model_CopyNameOfUnknown(model,E_Mass,"h_r") ;
#elif defined (U_P_L)
  Model_CopyNameOfUnknown(model,E_Mass,"p_l") ;
#else
  #error "Ambiguous or undefined unknown"
#endif

#if defined (U_LogC_s)
  Model_CopyNameOfUnknown(model,E_Salt,"logc_s") ;
#elif defined (U_C_s)
  Model_CopyNameOfUnknown(model,E_Salt,"c_s") ;
#else
  #error "Ambiguous or undefined unknown"
#endif
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    NbOfProp = 6 ;


  Material_ScanProperties(mat,datafile,pm) ;

  {
    ComputePhysicoChemicalProperties(TEMPERATURE) ;
  }


#if 0
  { /* on stocke les courbes lna_w et lna_s */
    Curves_t* curves = Material_GetCurves(mat) ;
    int    n_points = 1000 ;
    double c_s1 = 1.e-10*K_SALT ;
    double c_s2 = 10*K_SALT ;
    Curve_t* curve_lnaw = Curve_Create(n_points) ;
    Curve_t* curve_lnas = Curve_Create(n_points) ;
    double* x_lnaw = Curve_GetXRange(curve_lnaw) ;
    double* x_lnas = Curve_GetXRange(curve_lnas) ;
      
    Curve_GetScaleType(curve_lnaw) = 'l' ;
    Curve_GetScaleType(curve_lnas) = 'l' ;
      
    x_lnaw[0] = c_s1 ;
    x_lnaw[1] = c_s2 ;
    x_lnas[0] = c_s1 ;
    x_lnas[1] = c_s2 ;
    
    Curves_Append(curves,curve_lnaw) ;
    Curves_Append(curves,curve_lnas) ;
      
    {
      double* x = Curve_CreateSamplingOfX(curve_lnaw) ;
      double* y_lnaw = Curve_GetYValue(curve_lnaw) ;
      double* y_lnas = Curve_GetYValue(curve_lnas) ;
      int i ;
        
      for(i = 0 ; i < n_points ; i++) {
        double c_s  = x[i] ;
          
        y_lnaw[i] = LogActivityOfWater(c_s) ;
        y_lnas[i] = LogActivityOfSalt(c_s) ;
      }
      free(x) ;
    }

    /* on met a jour le nb de proprietes */
    Material_GetNbOfCurves(mat) += 2 ;
  }
#endif

  return(NbOfProp) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  printf("\n") ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("Equations to be solved are:\n") ;
  printf("\t- Conservation of total mass content  (mass)\n") ;
  printf("\t- Conservation of air mass content    (air)\n") ;
  printf("\t- Conservation of salt mole content   (salt)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Relative humidity        (h_r)\n") ;
  printf("\t- Gas pressure             (p_g)\n") ;
  printf("\t- Salt concentration       (c_s)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

  fprintf(ficd,"Model = DWS1     # Model\n") ;
  fprintf(ficd,"porosite = 0.12  # Porosite\n") ;
  fprintf(ficd,"k_int = 1.e-20   # Permeabilite intrinseque (m2)\n") ;
  fprintf(ficd,"D_Cl = 6.25e-12  # Diffusion effective de Cl (m2/s)\n") ;
  fprintf(ficd,"r_d = 1.         # Rapport des tortuosites anions/cations\n") ;
  fprintf(ficd,"courbes = myfile # nom du fichier : p_c S_l k_rl tau_l\n") ;

  return(0) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element 
 *  Return 0 */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1 ;

  /** Define the length of tables */
  Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
  Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads 
 *  Return 0 if succeeded and -1 if failed */
{
  //int dim = Element_GetDimensionOfSpace(el) ;
  //IntFct_t* fi = Element_GetIntFct(el) ;
  //int nn = Element_GetNbOfNodes(el) ;
  //int ndof = nn*NEQ ;
  //FEM_t* fem = FEM_GetInstance(el) ;

  {
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
  double* vex0 = Element_GetExplicitTerm(el) ;
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
  GetProperties(el) ;
    
    
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
  
    
  /* If there are initial values */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;
    
    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) SIG[i]    = x[I_SIG + i] ;
      
      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;
      
      M_Tot = x[I_M_Tot] ;
      for(i = 0 ; i < 3 ; i++) W_Tot[i]  = x[I_W_Tot + i] ;
      
      N_Salt = x[I_N_Salt] ;
      for(i = 0 ; i < 3 ; i++) W_Salt[i] = x[I_W_Salt + i] ;
    }
    
    
    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
    
      /* saturation degrees */
      double s_l = x[I_S_L] ;
      double s_g = x[I_S_G] ;
    
      /* porosity */
      double phi = phi0 ;
    
      /* tortuosity */
      double tau_g  = TortuosityToGas(phi,s_g) ;
      double tau_l  = TortuosityToLiquid(phi,s_l) ;
      double rho_l  = x[I_RHO_L] ;
      double k_h    = rho_l*kl_int/mu_l*RelativePermeabilityToLiquid(s_l) ;
      double d_vap  = phi * s_g * tau_g * D_av0 ;
      double d_salt = phi * s_l * tau_l * D_SALT ;
    
      KD_L    = k_h ;
      KF_Vap  = d_vap ;
      KF_Salt = d_salt ;
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t* el,double t)
{
}


int  ComputeImplicitTerms(Element_t* el,double t,double dt)
{}


int  ComputeMatrix(Element_t* el,double t,double dt,double* k)
{}


int  ComputeResidu(Element_t* el,double t,double dt,double* r)
{}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 7 ;
  double* vex  = Element_GetExplicitTerm(el) ;
  double* vim  = Element_GetCurrentImplicitTerm(el) ;
  double** u   = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t* fem = FEM_GetInstance(el) ;

  if(Element_IsSubmanifold(el)) return(0) ;

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
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim,t,0,p) ;
    /* Liquid pressure and relative humidity */
#if 0
#if defined (U_P_L)
    double p_l    = FEM_ComputeUnknown(fem,u,intfct,p,E_Mass) ;
    double h_r    = RelativeHumidity(p_l,lna_w) ;
#elif defined (U_H_r)
    double h_r    = FEM_ComputeUnknown(fem,u,intfct,p,E_Mass) ;
    double p_l    = LiquidPressure(h_r,lna_w) ;
#endif
#endif
    double p_l    = x[I_P_L] ;
    double h_r    = x[I_H_r] ;
    /* Displacement */
    double dis[3] = {0,0,0} ;
    /* strains */
    double w_tot[3] = {0,0,0} ;
    double w_salt[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    int    i ;
    
    for(i = 0 ; i < dim ; i++) {
      //dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
      dis[i] = x[I_U + i] ;
    }
    
    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) w_tot[j] += W_Tot[j]/np ;
      
      for(j = 0 ; j < 3 ; j++) w_salt[j] += W_Salt[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
    }
      
    i = 0 ;
    Result_Store(r + i++,&p_l     ,"Pore pressure",1) ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,w_tot    ,"Fluid mass flow",3) ;
    Result_Store(r + i++,w_salt   ,"Salt molar flow",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
  }
  
  return(NbOfOutputs) ;
}



double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double*   x      = Variable ;
  double*   x_n    = Variable_n(x) ;
  
    
  /* Primary Variables */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,p,E_Mech + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,E_Mech) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Unknown of total mass balance equation */
    x[E_Mass] = FEM_ComputeUnknown(fem,u,intfct,p,E_Mass) ;
    
    /* Unknwon gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,E_Mass) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_Mass + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
    
    /* Salt concentration */
    x[I_C_s] = FEM_ComputeUnknown(fem,u,intfct,p,E_Salt) ;
    
    /* Salt concentration gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,E_Salt) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_C_s + i] = grd[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    int    i ;
    
    /* Stresses, strains at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,E_Mech) ;
      double* vim_n = f_n + p*NVI ;
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_EPS   + i] = eps_n[i] ;
        x_n[I_SIG   + i] = SIG_n[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    }
    
    /* Unknwon at previous time step */
    x_n[E_Mass] = FEM_ComputeUnknown(fem,u_n,intfct,p,E_Mass) ;
    
    /* Salt concentration at previous time step */
    x_n[I_C_s] = FEM_ComputeUnknown(fem,u_n,intfct,p,E_Salt) ;
    
    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + p*NVE ;
      
      x_n[I_KD_L]  = KD_L ;
    }
  }
    
  ComputeSecondaryVariables(el,t,dt,x_n,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x_n,double* x)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  /* Strains */
  double* eps   =  x   + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
  /* Salt concentration */
  double c_s    = x[I_C_s] ;
  /* Gas pressure */
  double p_g    = p_g0 ;

  /* Molar concentration of water */
  double c_w    = (1. - V_AC*c_s)/V_H2O ;
  
  /* Activities of water and salt */
  double lna_w  = LogActivityOfWater(c_s) ;
  double lna_s  = LogActivityOfSalt(c_s) ;
  
  /* Pressures and relative humidity */
#if defined (U_P_L)
  double p_l    = x[E_Mass] ;
  double h_r    = RelativeHumidity(p_l,lna_w) ;
#elif defined (U_H_r)
  double h_r    = x[E_Mass] ;
  double p_l    = LiquidPressure(h_r,lna_w) ;
#endif
  double p_v    = p_v0*h_r ;
  double p_a    = p_g - p_v ;
  
  /* Saturation index of salt */
  double beta_s = SaturationIndexOfSalt(lna_s,p_l) ;
  
  /* Crystallization pressure */
  double p_c    = CrystallizationPressure(beta_s,p_l) ;
  
  /* Mass densities */
  double rho_w  = M_H2O*c_w ;
  double rho_s  = M_SALT*c_s ;
  double rho_l  = rho_w + rho_s ;
  double rho_v  = MassDensityOfWaterVapor(p_v) ;
  double rho_a  = MassDensityOfDryAir(p_a) ;
  double rho_g  = rho_v + rho_a ;
  
  double c_v    = rho_v/rho_g ;
  
  /* Saturations */
  double s_l    = SaturationDegree(p_g - p_l) ;
  double s_c    = 0. ;
  double s_g    = 1 - s_l - s_c ;
  
  /* Porosity */
  double phi    = phi0 ;

  /* Mass/Molar contents */
  double n_s    = c_s*phi*s_l ;
  double m_l    = rho_l*phi*s_l ;
  double m_g    = rho_g*phi*s_g ;
  double m_v    = rho_v*phi*s_g ;
  double m_t    = m_l + m_g ;
    

  /* Backup 
   * ------ */

  /* Stresses */
  {
    double* sig   = x   + I_SIG ;
    double* sig_n = x_n + I_SIG ;
    
    {
      double  deps[9] ;
      int    i ;
      
      /* Incremental deformations */
      for(i = 0 ; i < 9 ; i++) deps[i] =  eps[i] - eps_n[i] ;
    
      /* Stresses */
      for(i = 0 ; i < 9 ; i++) sig[i] = sig_n[i] ;
      
      for(i = 0 ; i < 9 ; i++) sig[i] = 0 ;
      
      #define C(i,j)  (cijkl[(i)*9+(j)])
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
      
        for(j = 0 ; j < 9 ; j++) {
          //sig[i] += C(i,j)*deps[j] ;
          sig[i] += C(i,j)*eps[j] ;
        }
      }
      #undef C
    }
    
    {
      double  pp = p_l*s_l + p_g*s_g + p_c*s_c ;
      
      sig[0] += - biot*pp ;
      sig[4] += - biot*pp ;
      sig[8] += - biot*pp ;
    }
  }
    
  /* Body force */
  {
      double* f_mass = x + I_Fmass ;
      int i ;
      
      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_t)*gravity ;
  }
  
  
  /* Mass content */
  {
      x[I_M_L]    = m_l ;
      x[I_M_Vap]  = m_v ;
      x[I_M_Tot]  = m_l + m_v ;
      x[I_N_Salt] = n_s ;
  }
  
  
  
  
  /* Flows */
  {
    {
      /* Transfer coefficient */
      double k_h = x_n[I_KD_L] ;
      double kf_vap = x_n[I_KF_Vap] ;
      double kf_salt = x_n[I_KF_Salt] ;
      double kc_salt = x_n[I_KC_Salt] ;
    
      /* Gradients */
#if defined (U_P_L)
      double* gpl = x + I_GRD_Mass ;
#else
  #error "Not possible yet"
#endif
      double* grd = x + I_GRD_C_s ;
    
      /* Mass flows */
      double* w_liq  = x + I_W_L ;
      double* w_vap  = x + I_W_Vap ;
      double* w_tot  = x + I_W_Tot ;
      double* w_salt = x + I_W_Salt ;
      int i ;
    
    /* Liquid */
      for(i = 0 ; i < 3 ; i++) {
        w_liq[i] = - k_h * gpl[i] ;
      }
      w_liq[dim - 1] += k_h * rho_l * gravity ;
    
      /* Vapor */
      for(i = 0 ; i < 3 ; i++) {
        w_vap[i] = - kf_vap * gpl[i] ;
      }
    
      /* Total */
      for(i = 0 ; i < 3 ; i++) {
        w_tot[i] = w_liq[i] + w_vap[i] ;
      }
      
      /* Salt */
      for(i = 0 ; i < 3 ; i++) {
        w_salt[i] = kc_salt * w_liq[i] - kf_salt * grd[i] ;
      }
    }
  }
  
  /* Others */
  {
    x[I_RHO_L] = rho_l ;
    x[I_S_L  ] = s_l ;
    x[I_S_G  ] = s_g ;
    x[I_S_C  ] = s_c ;
  }
}





double activite_w(double c_s,double TK)
/* Water activity in brine */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  /* molalites (moles/mole) */
  double m_s   = (c_s > 0.) ? c_s/c_w : 0. ;
  double m_ani = NU_A*m_s ;
  double m_cat = NU_C*m_s ;
  /* force ionique */
  double I    = 0.5*(Z_A*Z_A*m_ani + Z_C*Z_C*m_cat) ;
  /* references */
  double T0   = 273.15 ;
  double b0   = sqrt(M_H2O),S0   = pow(M_H2O,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_cl_nacl = 1.827/b0,S_cl_nacl = 19.245/S0 ;
  double b_na_nacl = 4.352/b0,S_na_nacl = 26.448/S0 ;
  /* Na2SO4 (d'apres Lin et Lee) */
  double b_so4_na2so4 = 1.662/b0,S_so4_na2so4 = 0.022/S0 ;
  double b_na_na2so4 = 1.552/b0,S_na_na2so4 = 3.464/S0 ;

  double epsi = 0.0007*(TK - T0)*(TK - T0) - 0.3918*(TK - T0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*TK,1.5)/b0 ;

  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;

  switch(SALT) {
  case(NaCl) : {
    b_cat = b_na_nacl ;
    b_ani = b_cl_nacl ;
    S_cat = S_na_nacl ;
    S_ani = S_cl_nacl ;
    break ;
  } 
  case(Na2SO4) : {
    b_cat = b_na_na2so4 ;
    b_ani = b_so4_na2so4 ;
    S_cat = S_na_na2so4 ;
    S_ani = S_so4_na2so4 ;
    break ;
  }
  default : {
    arret("non prevu") ;
  }
  }

  if (I > 0.) {
    double lna_w ;
    lna_w = m_ani*lna_i(TK,I,Z_A,b_ani,S_ani,A)
          + m_cat*lna_i(TK,I,Z_C,b_cat,S_cat,A) ;
    return(lna_w) ;
  } else {
    return(0.) ;
  }
}


double activite_s(double c_s,double TK)
/* Salt activity in brine */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  /* molalites (moles/mole) */
  double m_s   = (c_s/c_w > DBL_MIN) ? c_s/c_w : DBL_MIN ; 
  double m_ani = NU_A*m_s ;
  double m_cat = NU_C*m_s ;
  /* force ionique */
  double I     = 0.5*(Z_A*Z_A*m_ani + Z_C*Z_C*m_cat) ; 
  /* references */
  double T0   = 273.15 ;
  double b0   = sqrt(M_H2O),S0   = pow(M_H2O,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_cl_nacl = 1.827/b0,S_cl_nacl = 19.245/S0 ;
  double b_na_nacl = 4.352/b0,S_na_nacl = 26.448/S0 ;
  /* Na2SO4 (d'apres Lin et Lee) */
  double b_so4_na2so4 = 1.662/b0,S_so4_na2so4 = 0.022/S0 ;
  double b_na_na2so4 = 1.552/b0,S_na_na2so4 = 3.464/S0 ;

  double epsi = 0.0007*(TK - T0)*(TK - T0) - 0.3918*(TK - T0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*TK,1.5)/b0 ;
  
  /* depend du sel */
  double b_cat ;
  double b_ani ;
  double S_cat ;
  double S_ani ;

  switch(SALT) {
  case(NaCl) : {
    b_cat = b_na_nacl ;
    b_ani = b_cl_nacl ;
    S_cat = S_na_nacl ;
    S_ani = S_cl_nacl ;
    break ;
  }
  case(Na2SO4) : {
    b_cat = b_na_na2so4 ;
    b_ani = b_so4_na2so4 ;
    S_cat = S_na_na2so4 ;
    S_ani = S_so4_na2so4 ;
    break ;
  }
  default : {
    arret("non prevu") ;
  }
  }

  {
    double m_t     = m_ani + m_cat ;
    double lna_w   = activite_w(c_s,TK) ;
    double lng_ani = lng_TQN(TK,I,Z_A,b_ani,S_ani,A,lna_w,m_t) ;
    double lng_cat = lng_TQN(TK,I,Z_C,b_cat,S_cat,A,lna_w,m_t) ;
    double lna_s   = (NU_A*lng_ani + NU_C*lng_cat) + NU_AC*log(m_s) ;
    return(lna_s) ;
  }
}


double activite_w_ideal(double c_s)
/* L'activite chimique de l'eau d'une solution ideale */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  double c_cat = NU_C*c_s ;
  double c_ani = NU_A*c_s ;
  double c_t   = c_w + c_cat + c_ani ;
  double lna_w = log(c_w/c_t) ;
  return(lna_w) ;
}


double activite_s_ideal(double c_s1)
/* L'activite chimique du sel d'une solution ideale */
{
  double c_s     = (c_s1 > DBL_MIN/V_H2O) ? c_s1 : DBL_MIN/V_H2O ; 
  double c_w     = (1. - V_AC*c_s)/V_H2O ;
  double c_cat   = NU_C*c_s ;
  double c_ani   = NU_A*c_s ;
  double c_t     = c_w + c_cat + c_ani ;
  double lna_cat = log(c_cat/c_t) ;
  double lna_ani = log(c_ani/c_t) ;
  double lna_s   = NU_C*lna_cat + NU_A*lna_ani ;
  return(lna_s) ;
}





double lng_LinLee(double T,double I,double z,double b,double S,double A)
/* Le log du coefficient d'activite d'un ion d'apres Lin & Lee */ 
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*(II/(1 + b*II) + 2*log(1 + b*II)/b) + S*pow(I,alpha)/T ;
  
  return(lng*z*z) ;
}

double lng_TQN(double T,double I,double z,double b,double S,double A,double lna_w,double m_t)
/* Le log du coefficient d'activite d'un ion (T.Q Nguyen) :
   lng_i = dGamma/dm_i = (dGamma/dm_i)_I - 0.5*z_i*z_i*(lna_w + m_t)/I 
   lna_w = - m_t - sum_i ( m_i*lng_i ) + Gamma */
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*2*log(1 + b*II)/b + S*pow(I,alpha)/(1+alpha)/T - 0.5*(lna_w + m_t)/I ;
  
  return(lng*z*z) ;
}

double lna_i(double T,double I,double z,double b,double S,double A)
/* Contribution de chaque ion au log de l'activite du solvant 
   lna_w = sum_i ( m_i*lna_i ) (T.Q Nguyen) */ 
{
  double alpha = 1.29,a1 = alpha/(1+alpha),II = sqrt(I) ;
  double lna ;
  
  lna = A*II/(1 + b*II) - a1*S*pow(I,alpha)/T ;
  
  return(-1 + lna*z*z) ;
}



double EquivalentePressure(Element_t* el,double p_l,double c_s)
{
  double p_g = p_g0 ;
  /* Saturations */
  double s_l    = SaturationDegree(p_g - p_l) ;
  double s_c    = 0. ;
  double s_g    = 1 - s_l - s_c ;
  double p_c    = 0 ;
  double pp     = s_l*p_l + s_g*p_g + s_c*p_c ;
  
  return(pp) ;
}
