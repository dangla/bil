#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*
#include <ctype.h>
*/
#include <float.h>
#include "CommonModel.h"
#include "FVM.h"


#define TITLE   "Drying-Wetting with Salt (only dissolved) (1D case)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"



/* Indices of equations/unknowns */
enum {
  E_Mass   ,
  /* Uncomment/Comment the next two lines to consider/suppress air */
  E_Air   ,
  #define E_Air E_Air
  E_Salt   ,
  E_Last
} ;


/* Nb of equations */
#define NEQ      (E_Last)



/* Value of the nodal unknown (u, u_n and el must be used below) */
#define UNKNOWN(n,i)     Element_GetValueOfNodalUnknown(el,u,n,i)
#define UNKNOWNn(n,i)    Element_GetValueOfNodalUnknown(el,u_n,n,i)


/* Generic names of nodal unknowns */
#define U_Mass(n)     (UNKNOWN(n,E_Mass))
#define Un_Mass(n)    (UNKNOWNn(n,E_Mass))

#define U_Air(n)      (UNKNOWN(n,E_Air))
#define Un_Air(n)     (UNKNOWNn(n,E_Air))

#define U_Salt(n)     (UNKNOWN(n,E_Salt))
#define Un_Salt(n)    (UNKNOWNn(n,E_Salt))



/* Method chosen at compiling time. 
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with a specific model.
 * Uncomment/comment to let only one unknown per equation */

/* Mass: the relative humidity */
#define U_H_r       U_Mass

/* Air: the gas pressure */
#ifdef E_Air
#define U_P_G       U_Air
#endif

/* Salt: the salt concentration */
//#define U_LogC_s    U_Salt
#define U_C_s       U_Salt



/* Names of nodal unknowns */
#if defined (U_LogC_s) && !defined (U_C_s)
  #define LogC_s(n)   U_Salt(n)
  #define LogC_sn(n)  Un_Salt(n)
  #define C_s(n)      (pow(10,LogC_s(n)))
  #define C_sn(n)     (pow(10,LogC_sn(n)))
#elif defined (U_C_s) && !defined (U_LogC_s)
  #define C_s(n)      U_Salt(n)
  #define C_sn(n)     Un_Salt(n)
  #define LogC_s(n)   (log10(C_s(n)))
  #define LogC_sn(n)  (log10(C_sn(n)))
#else
  #error "Ambiguous or undefined unknown"
#endif

#define H_r(n)      U_Mass(n)
#define P_G(n)      U_Air(n)





/* Nb of nodes (el must be used below) */
#define NbOfNodes               Element_GetNbOfNodes(el)



/* Nb of (im/ex)plicit terms and constant terms */
#define NbOfExplicitTerms       (6*NbOfNodes)
#define NbOfImplicitTerms       (3*NbOfNodes*NbOfNodes)
#define NbOfConstantTerms       (0)
#define NVI      (NbOfImplicitTerms)
#define NVE      (NbOfExplicitTerms)
#define NV0      (0)


/* Names used for implicit terms */
#define MassAndFlux(f,i,j)  ((f)[((i)*NbOfNodes + (j))])

#define MW_T          (f)
#define MW_Tn         (f_n)
#define M_T(i)        MassAndFlux(MW_T,i,i)
#define M_Tn(i)       MassAndFlux(MW_Tn,i,i)
#define W_T(i,j)      MassAndFlux(MW_T,i,j)


#define MW_A          (f   + NbOfNodes*NbOfNodes)
#define MW_An         (f_n + NbOfNodes*NbOfNodes)
#define M_A(i)        MassAndFlux(MW_A,i,i)
#define M_An(i)       MassAndFlux(MW_An,i,i)
#define W_A(i,j)      MassAndFlux(MW_A,i,j)


#define NW_S        (f   + 2*NbOfNodes*NbOfNodes)
#define NW_Sn       (f_n + 2*NbOfNodes*NbOfNodes)
#define N_S(i)        MassAndFlux(NW_S,i,i)
#define N_Sn(i)       MassAndFlux(NW_Sn,i,i)
#define W_S(i,j)      MassAndFlux(NW_S,i,j)


/* Names used for explicit terms */
#define TransferCoefficient(va,n)  ((va) + (n)*NbOfNodes)

#define KD_L          TransferCoefficient(va,0)
#define KD_G          TransferCoefficient(va,1)

#define KF_V          TransferCoefficient(va,2)
#define KF_S          TransferCoefficient(va,3)

#define KC_V          TransferCoefficient(va,4)
#define KC_S          TransferCoefficient(va,5)




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

#define SaturationDegree(pc)             (saturation(pc,p_c3,SATURATION_CURVE))
#define RelativePermeabilityToLiquid(pc) (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,pc))
#define RelativePermeabilityToGas(pc)    (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,pc))
#define TortuosityToGas(f,sg)            ((sg > 0) ? pow(f,aa)*pow(sg,bb) : 0)
#define aa                        (0.33)  /* 1/3 Millington, Thiery 1.74 */
#define bb                        (2.33)   /* 7/3 Millington, Thiery 3.2 */
//#define TortuosityToLiquid(f,sl)         (tortuosite_l(f)*pow(sl,4.5))
#define TortuosityToLiquid(f,sl)         ((sl > 0) ? (pow(f,0.33)*pow(sl,2.33)) : 0)
/*
#define TortuosityToLiquid(f,sl)         (0.00029*exp(9.95*f)/(1+625*((sl < 1) ? pow((1-sl),4) : 0)))
*/
    


/* Salt properties
 * --------------- */
/* Valencies */
#define z_cl      (-1.)
#define z_na      (1.)
#define z_so4     (-2.)

/* Stoichiometric coefficients */
#define nu_h2o_nacl   (0)
#define nu_na_nacl    (1)
#define nu_cl_nacl    (1)
#define nu_h2o_na2so4 (10)
#define nu_na_na2so4  (2)
#define nu_so4_na2so4 (1)

/* Partial molar volumes in liquid (m3/mole) */
#define v_na      (1.87e-5)
#define v_cl      (2.52e-6)
#define v_so4     (-8.334e-6)

/* Solid molar volumes volumes (m3/mole) */
#define v_nacl    (24.5e-6)
#define v_na2so4  (220.e-6)

/* Molar mass */
#define M_nacl    (58.45e-3)
#define M_na2so4  (142.e-3)

/* Molecular diffusion coefficients (m2/s) */
#define do_cl     (2.032e-9)
#define do_na     (1.334e-9)
#define do_so4    (1.065e-9)

/* Equilibrium constants */
#define K_nacl    (6.e3)      /* Solubilite de NaCl (moles/m3) */
#define K_na2so4  (1.24166e3) /* Solubilite de Na2SO4.10H2O (moles/m3) */

/* Surface tensions (N/m) */
#define Gamma_cl_nacl   (0.1)
#define Gamma_cl_na2so4 (0.1)

#define GAMMA_LG        (0.07)


/* Types of salt */
#define NaCl    (0)
#define Na2SO4  (1)
#define KNO3    (2)

/* The chosen salt */
#define SALT     Na2SO4

#if SALT == NaCl
#define Z_A       z_cl
#define Z_C       z_na
#define NU_A      nu_cl_nacl
#define NU_C      nu_na_nacl
#define NU_H2O    nu_h2o_nacl
#define V_A       v_cl
#define V_C       v_na
#define V_SALT    v_nacl
#define M_SALT    M_nacl
#define D_A       do_cl
#define D_C       do_na
#define K_SALT    K_nacl
#define GAMMA_CL  Gamma_cl_nacl
#elif SALT == Na2SO4
#define Z_A       z_so4
#define Z_C       z_na
#define NU_A      nu_so4_na2so4
#define NU_C      nu_na_na2so4
#define NU_H2O    nu_h2o_na2so4
#define V_A       v_so4
#define V_C       v_na
#define V_SALT    v_na2so4
#define M_SALT    M_na2so4
#define D_A       do_so4
#define D_C       do_na
#define K_SALT    K_na2so4
#define GAMMA_CL  Gamma_cl_na2so4
#else
#error "Salt not available"
#endif

#define NU_AC     (NU_A + NU_C)
#define V_AC      (NU_A*V_A + NU_C*V_C)
#define V_ACH     (NU_H2O*V_H2O + V_AC)
/* Logarithm of salt activity in brine  */
//#define LogActivityOfSalt(c_s)     activite_s(c_s,TEMPERATURE)
#define LogActivityOfSalt(c_s)     activite_s_ideal(c_s)



#define TEMPERATURE         (293.)      /* Temperature (K) */


/* Math constants */
#define Ln10      Math_Ln10


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])


/* Fonctions */
static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeVariables(Element_t*,double**,double*,double,double,int) ;
//static Model_ComputeSecondaryVariables_t    ComputeSecondaryVariables ;
static int     ComputeSecondaryVariables(Element_t*,double,double,double*) ;

static double* ComputeVariableDerivatives(Element_t*,double,double,double,int,int) ;

static void    ComputeTransferCoefficients(FVM_t*,double**,double*) ;
//static double* ComputeVariableFluxes(Element_t*,double**) ;
static FVM_ComputeFluxes_t    ComputeFluxes ;
static int     TangentCoefficients(FVM_t*,double,double*) ;

static void    ComputePhysicoChemicalProperties(double) ;

static double  lna_i(double,double,double,double,double,double) ;
static double  lng_TQN(double,double,double,double,double,double,double,double) ;
static double  lng_LinLee(double,double,double,double,double,double) ;

static double  tortuosite_l(double) ;

static double  saturation(double,double,Curve_t*) ;

static double  activite_w(double,double) ;
static double  activite_s(double,double) ;
static double  activite_w_ideal(double) ;
static double  activite_s_ideal(double) ;



/* Parameters */
static double phi0 ;
static double kl_int ;
static double kg_int ;
static double p_c3 ;

static double p_l0 ;
static double p_g0 ;
static double p_v0 ;

static double RT ;
static double D_av0 ;
static double mu_l ;
static double mu_g ;



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




enum {
I_M_T = NEQ    ,
I_M_A          ,
I_M_L          ,
I_M_V          ,
I_M_G          ,
I_N_S          ,


I_RHO_V        ,
I_RHO_A        ,
I_RHO_G        ,
I_RHO_L        ,
I_RHO_W        ,
I_RHO_S        ,

I_S_L          ,

I_H_L          ,
I_H_G          ,
I_COOR_X       ,

I_P_L          ,
I_P_G          ,
I_P_V          ,
I_P_A          ,

I_LNA_W        ,
I_LNA_S        ,

I_C_V          ,
I_C_W          ,
I_C_S          ,

I_H_R          ,
I_Last
} ;

//#define NN                (2)
#define NbOfVariables    (I_Last)
//static double Variables[NN][NbOfVariables] ;
//static double dVariables[NbOfVariables] ;
  
  


enum {
I_W_L           ,
I_W_V           ,
I_W_A           ,
I_W_G           ,
I_W_T           ,
I_W_S           ,
I_W_Ani         ,
I_W_Cat         ,
I_W_Last
} ;

#define NbOfVariableFluxes    (I_W_Last)
//static double VariableFluxes[NbOfVariableFluxes] ;



int pm(const char *s)
{
       if(!strcmp(s,"porosite")) return (0) ;
  else if(!strcmp(s,"k_int"))    return (1) ;
  else if(!strcmp(s,"kl_int"))   return (1) ;
  else if(!strcmp(s,"kg_int"))   return (2) ;
  else if(!strcmp(s,"p_c3"))     return (3) ;
  else return (-1) ;
}


void GetProperties(Element_t *el)
{
  phi0    = GetProperty("porosite") ;
  kl_int  = GetProperty("kl_int") ;
  kg_int  = GetProperty("kg_int") ;
  p_c3    = GetProperty("p_c3") ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_Mass,"mass") ;
#if defined (E_Air)
  Model_CopyNameOfEquation(model,E_Air ,"air") ;
#endif
  Model_CopyNameOfEquation(model,E_Salt,"salt") ;

  Model_CopyNameOfUnknown(model,E_Mass,"h_r") ;
#if defined (E_Air)
  Model_CopyNameOfUnknown(model,E_Air,"p_g") ;
#endif
#if defined (U_LogC_s) && !defined (U_C_s)
  Model_CopyNameOfUnknown(model,E_Salt,"logc_s") ;
#elif defined (U_C_s) && !defined (U_LogC_s)
  Model_CopyNameOfUnknown(model,E_Salt,"c_s") ;
#endif

  //Model_GetNbOfVariables(model) = NbOfVariables ;
  //Model_GetNbOfVariableFluxes(model) = NbOfVariableFluxes ;
  //Model_GetComputeSecondaryVariables(model) = ComputeSecondaryVariables ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    NbOfProp = 6 ;

  {
    int dim = Material_GetDimension(mat) ;
    
    if(dim > 1) arret("DWS1: dimension > 1 not available") ;
  }
  
  {
    /* Self-initialization */
    Material_GetProperty(mat)[pm("kg_int")] = 0 ;
  }

  Material_ScanProperties(mat,datafile,pm) ;

  {
    kg_int = Material_GetProperty(mat)[pm("kg_int")] ;
    kl_int = Material_GetProperty(mat)[pm("kl_int")] ;

    if(kg_int == 0.) Material_GetProperty(mat)[pm("kg_int")] = kl_int ;
  }

  {
    ComputePhysicoChemicalProperties(TEMPERATURE) ;
  }


#ifdef NOTDEFINED
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
  
  if(!ficd) return(NEQ) ;
  
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
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* Mass contents */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double *x = ComputeVariables(el,u,f,0,0,i) ;

    M_T(i)   = x[I_M_T] ;
    M_A(i)   = x[I_M_A] ;
    N_S(i)   = x[I_N_S] ;
  }

  /* Transfer coefficients */
  ComputeTransferCoefficients(fvm,u,f) ;

  /* Fluxes */
  {
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = FVM_ComputeVariableFluxes(fvm,ComputeFluxes,i,j) ;

        W_T(i,j)     = w[I_W_T] ;
        W_A(i,j)     = w[I_W_A] ;
        W_S(i,j)     = w[I_W_S] ;
        
        W_T(j,i)     = - w[I_W_T] ;
        W_A(j,i)     = - w[I_W_A] ;
        W_S(j,i)     = - w[I_W_S] ;
      }
    }
  }
  
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Explicites terms  */
{
  
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Transfer coefficients */
  {
    double*  f = Element_GetPreviousImplicitTerm(el) ;
    double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
    int nn = Element_GetNbOfNodes(el) ;
    FVM_t *fvm = FVM_GetInstance(el) ;
    int i ;
    
    for(i = 0 ; i < nn ; i++) {
      ComputeVariables(el,u,f,t,0,i) ;
    }
    
    ComputeTransferCoefficients(fvm,u,f) ;
  }

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Implicit terms */
{
  double* f   = Element_GetCurrentImplicitTerm(el) ;
  double* f_n = Element_GetPreviousImplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
   
  /* Mass contents */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    double *x = ComputeVariables(el,u,f_n,t,dt,i) ;

    M_T(i)   = x[I_M_T] ;
    M_A(i)   = x[I_M_A] ;
    N_S(i)   = x[I_N_S] ;
    
    {
      double h_r = x[I_H_R] ;
      double c_s = x[I_C_S] ;
      double c_w = x[I_C_W] ;

      //if(h_r <= 0 || c_s < 0. || c_w < 0.) {
      if(h_r <= 0 || c_w < 0.) {
        double xx = Element_ComputePointerToNodalCoordinates(el)[i][0] ;
        double s_l = x[I_S_L] ;
      
        printf("\n") ;
        printf("x       = %e\n",xx) ;
        printf("h_r     = %e\n",h_r) ;
        printf("c_s     = %e\n",c_s) ;
        printf("c_w     = %e\n",c_w) ;
        printf("s_l     = %e\n",s_l) ;
        return(-1) ;
      }
    }
  }
  

  /* Fluxes */
  {
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = i + 1 ; j < nn ; j++) {
        double* w = FVM_ComputeVariableFluxes(fvm,ComputeFluxes,i,j) ;

        W_T(i,j)     = w[I_W_T] ;
        W_A(i,j)     = w[I_W_A] ;
        W_S(i,j)     = w[I_W_S] ;
        
        W_T(j,i)     = - w[I_W_T] ;
        W_A(j,i)     = - w[I_W_A] ;
        W_S(j,i)     = - w[I_W_S] ;
      }
    }
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
  
  TangentCoefficients(fvm,dt,c) ;
  {
    double *km = FVM_ComputeMassAndIsotropicConductionMatrix(fvm,c,NEQ) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) k[i] = km[i] ;
  }
  

#if defined (U_LogC_s)
  {
    double** u = Element_ComputePointerToNodalUnknowns(el) ;
    
    for(i = 0 ; i < ndof ; i++){
      K(i,U_C_s)     *= Ln10*C_s(0) ;
      K(i,U_C_s+NEQ) *= Ln10*C_s(1) ;
    }
  }
#endif
  
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
  FVM_t* fvm = FVM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  

  /*
    Conservation of total mass
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,MW_T,MW_Tn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Mass) -= r1[i] ;
    }
  }
  
  /*
    Conservation of air mass
  */
#if defined (E_Air)
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,MW_A,MW_An,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Air) -= r1[i] ;
    }
  }
#endif

  /*
    Conservation of salt mass
  */
  {
    double* r1 = FVM_ComputeMassBalanceEquationResidu(fvm,NW_S,NW_Sn,dt) ;
      
    for(i = 0 ; i < nn ; i++) {
      R(i,E_Salt) -= r1[i] ;
    }
  }
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double*  f = Element_GetCurrentImplicitTerm(el) ;
  Model_t* model = Element_GetModel(el) ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    nso = 15 ;
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
  for(i = 0 ; i < nn ; i++) {
    ComputeVariables(el,u,f,t,0,i) ;
  }
  
  {
    int    j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Variables */
    //double* x    = ComputeVariables(el,u,f,t,0,j) ;
    double* x    = Model_GetVariable(model,j) ;
    /* Fluxes */
    double* w    = FVM_ComputeVariableFluxes(fvm,ComputeFluxes,0,1) ;
    
    /* quantites exploitees */
    i = 0 ;
    Result_Store(r + i++,(x + I_P_L  ),"Liquid pressure",1) ;
    Result_Store(r + i++,(x + I_P_G  ),"Gas pressure",1) ;
    Result_Store(r + i++,(x + I_P_V  ),"Vapor pressure",1) ;
    Result_Store(r + i++,(w + I_W_L  ),"Liquis mass flow",1) ;
    Result_Store(r + i++,(w + I_W_G  ),"Gas mass flow",1) ;
    Result_Store(r + i++,(w + I_W_V  ),"Vapor mass flow",1) ;
    Result_Store(r + i++,(x + I_S_L  ),"Liquid saturation",1) ;
    Result_Store(r + i++,(x + I_H_R  ),"Relative humidity",1) ;
    Result_Store(r + i++,(x + I_C_S  ),"Salt concentration",1) ;
    Result_Store(r + i++,(x + I_LNA_W),"Log of water activity",1) ;
    Result_Store(r + i++,(x + I_LNA_S),"Log of salt activity",1) ;
    {
      double c_s = x[I_C_S] ;
      double lna_w_ideal  = activite_w_ideal(c_s) ;
      
      Result_Store(r + i++,&lna_w_ideal,"Log of ideal water activity",1) ;
    }
    {
      double c_s = x[I_C_S] ;
      double lna_s_ideal  = activite_s_ideal(c_s) ;
      
      Result_Store(r + i++,&lna_s_ideal,"Log of ideal salt activity",1) ;
    }
    Result_Store(r + i++,(w + I_W_Ani),"Anion mass flow",1) ;
    Result_Store(r + i++,(w + I_W_Cat),"Cation mass flow",1) ;
    
    if(i != nso) {
      arret("DWS1") ;
    }
  }
  
  return(nso) ;
}




void ComputeTransferCoefficients(FVM_t* fvm,double **u,double *f)
/* Termes explicites (va)  */
{
  Element_t* el = FVM_GetElement(fvm) ;
  Model_t* model = Element_GetModel(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int i ;
  
  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  
  
  /* Transfer coefficients */
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    //double *x = ComputeVariables(el,u,f,0,0,i) ;
    double *x = Model_GetVariable(model,i) ;
    
    double p_l    = x[I_P_L] ;
    double p_g    = x[I_P_G] ;
    double p_c    = p_g - p_l ;

    double s_l    = x[I_S_L] ;
    double s_g    = 1 - s_l ;

    double rho_g  = x[I_RHO_G] ;
    double rho_l  = x[I_RHO_L] ;
    
    double c_v    = x[I_C_V] ;
    double c_s    = x[I_C_S] ;
  
    /* Relative permeabilities */
    double k_rl   = RelativePermeabilityToLiquid(p_c) ;
    double k_rg   = RelativePermeabilityToGas(p_c) ;
  
    /* Permeabilities */
    double kh_l   = kl_int/mu_l*k_rl ;
    double kh_g   = kg_int/mu_g*k_rg ;
  
    double phi    = phi0 ;
    
    /* tortuosites gaz et liquide*/
    double tau_g   = TortuosityToGas(phi,s_g) ;
    double tau_l   = TortuosityToLiquid(phi,s_l) ;
    double tau_ani = tau_l ;
    /*
    double tau_cat = tau_ani/r_d ;
    */
    double tau_cat = tau_ani ;
  
    /* Darcy */
    double kd_l  = rho_l*kh_l ;           /* liquid */
    double kd_g  = rho_g*kh_g ;           /* gas */
    
    /* Fick */
    double D_av  = D_av0*p_g0/p_g ;
    double D_eff = phi*s_g*tau_g*D_av ;
    double kf_v  = D_eff ;                /* vapor */
    double kf_a  = phi*s_l*tau_ani*D_A ;
    double kf_c  = phi*s_l*tau_cat*D_C ;
    double kf_s  = (Z_C - Z_A)*kf_c*kf_a/(Z_C*kf_c - Z_A*kf_a) ;
    
    
    
    /* Back up */
    KD_L[i]   = kd_l ;
    KD_G[i]   = kd_g ;
    
    KF_V[i]   = kf_v ;
    KF_S[i]   = kf_s ;
    KC_V[i]   = c_v ;
    KC_S[i]   = c_s/rho_l ;
  }
  
}



void ComputeFluxes(FVM_t* fvm,double* grd,double* w,int i,int j)
{
  Element_t* el = FVM_GetElement(fvm) ;
  double* va = Element_GetExplicitTerm(el) ;
  

  /* Gradients */
  /*
  double grd_h_l    = grd[I_H_L] ;
  double grd_h_g    = grd[I_H_G] ;
  */
  double grd_p_l    = grd[I_P_L] ;
  double grd_p_g    = grd[I_P_G] ;
  
  double grd_rho_v  = grd[I_RHO_V] ;
  
  double grd_c_s    = grd[I_C_S] ;
      
  /* Transfer terms */
  double kd_l   = 0.5 * (KD_L[i]   + KD_L[j]) ;
  double kd_g   = 0.5 * (KD_G[i]   + KD_G[j]) ;
  double kf_v   = 0.5 * (KF_V[i]   + KF_V[j]) ;
  double kc_v   = 0.5 * (KC_V[i]   + KC_V[j]) ;
  double kf_s   = 0.5 * (KF_S[i]   + KF_S[j]) ;
  double kc_s   = 0.5 * (KC_S[i]   + KC_S[j]) ;
    
    
  /* Fluxes */
  double w_l   = - kd_l*grd_p_l ;
  double w_g   = - kd_g*grd_p_g ;
  double w_t   =   w_l + w_g ;
  double j_v   = - kf_v*grd_rho_v ;
  double w_v   =   kc_v*w_g + j_v ;
  double w_a   =   w_g - w_v ;
 
  double j_s   = - kf_s*grd_c_s ;
  double w_s   =   kc_s*w_l + j_s ;
  double w_ani =   NU_A*w_s ;
  double w_cat =   NU_C*w_s ;
   
   
  w[I_W_L ]  = w_l ;
  w[I_W_G ]  = w_g ;
  w[I_W_V ]  = w_v ;
  w[I_W_A ]  = w_a ;
  w[I_W_T ]  = w_t ;
  
  w[I_W_S ]  = w_s ;
  w[I_W_Ani] = w_ani ;
  w[I_W_Cat] = w_cat ;
   
  //return(w) ;
}


double tortuosite_l(double phi)
/* After  
 * B.H. Oh and S.Y. Jang, Prediction of diffusivity of concrete based on simple analytic equations,
 * Cement and Concrete Research 34 (2004) 463–480.
 * */
{
  double phi_cap = phi/2  ;
  double phi_c   = 0.17 ; /*Percolation capilar porosity*/
  double n = 2.7 ;        /* OPC n = 2.7        , Silica fume n = 4.5 */
  double ds_norm = 5e-5 ; /* OPC ds_norm = 2e-4 , Silica fume ds_norm = 5e-5 */
  double dsn_norm = pow(ds_norm,1/n) ;
  double m_phi = 0.5*(dsn_norm + phi_cap/(1-phi_c)*(1 - dsn_norm) - phi_c/(1-phi_c)) ;
  double iff   =  pow(m_phi + sqrt(m_phi*m_phi +  dsn_norm*phi_c/(1-phi_c)),n) ;
    
  return(iff) ;
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



int TangentCoefficients(FVM_t* fvm,double dt,double *c)
/**  Tangent matrix coefficients (c) */
{
  Element_t* el = FVM_GetElement(fvm) ;
  Model_t* model = Element_GetModel(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  ObVal_t *obval = Element_GetObjectiveValue(el) ;
  double* dist = FVM_ComputeIntercellDistances(fvm) ;
  double dij   = dist[1] ;
  double dtdij = dt/dij ;
  int    dec = NEQ*NEQ ;
  double dxi[NEQ] ;
  int    i ;
  
  /* Initialization */
  for(i = 0 ; i < ndof*ndof ; i++) c[i] = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  {
    double*  f_n = Element_GetPreviousImplicitTerm(el) ;
    double** u   = Element_ComputePointerToNodalUnknowns(el) ;
    double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
    
    for(i = 0 ; i < nn ; i++) {
      ComputeVariables(el,u,f_n,0,dt,i) ;
    }
  }
  
    
  for(i = 0 ; i < NEQ ; i++) {
    dxi[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    /* dxi[i] =  1.e2 ; */
  }
  
  
  for(i = 0 ; i < nn ; i++) {
    /* Variables */
    //double *x         = ComputeVariables(el,u,f_n,t,dt,i) ;
    //double *x = Model_GetVariable(model,i) ;
    int k ;
  
    #if defined (U_LogC_s)
    dxi[U_LogC_s] =  1.e-2*ObVal_GetValue(obval + U_LogC_s)*C_sn(i) ;
    #endif
    
    for(k = 0 ; k < NEQ ; k++) {
      double dxk    = dxi[k] ;
      double *dx    = ComputeVariableDerivatives(el,0,dt,dxk,k,i) ;
      
      /* Content terms at node i */
      {
        double *cii = c + (i*nn + i)*NEQ*NEQ ;
          
        cii[E_Mass*NEQ   + k] = dx[I_M_T] ;
#if defined (E_Air)
        cii[E_Air*NEQ    + k] = dx[I_M_A] ;
#endif
        cii[E_Salt*NEQ   + k] = dx[I_N_S] ;
      }
      
      /* Transfer terms from node i to node j */
      {
        double* dw = Model_GetVariableFluxDerivative(model,i) ;
        int j ;
      
        for(j = 0 ; j < nn ; j++) {
          if(j != i) {
            double *cij = c + (i*nn + j)*NEQ*NEQ ;
        
            ComputeFluxes(fvm,dx,dw,i,j) ;
          
            cij[E_Mass*NEQ   + k] = - dtdij*dw[I_W_T] ;
#if defined (E_Air)
            cij[E_Air*NEQ    + k] = - dtdij*dw[I_W_A] ;
#endif
            cij[E_Salt*NEQ   + k] = - dtdij*dw[I_W_S] ;
          }
        }
      }
    }
  }

  return(dec) ;
}





double* ComputeVariables(Element_t* el,double **u,double *f_n,double t,double dt,int n)
{
  Model_t* model = Element_GetModel(el) ;
  double *x = Model_GetVariable(model,n) ;
  
  /* Primary Variables */
  x[E_Mass ] = H_r(n) ;
  x[E_Salt ] = C_s(n) ;
#if defined (E_Air)
  x[E_Air  ] = P_G(n) ;
#endif
  
  /* Needed variables to compute secondary components */
    
  ComputeSecondaryVariables(el,t,dt,x) ;
  return(x) ;
}




double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double dxi,int i,int n)
{
  Model_t* model = Element_GetModel(el) ;
  //int NbOfVariables = Model_GetNbOfVariables(model) ;
  
  if(NbOfVariables > Model_MaxNbOfVariables) {
    arret("ComputeVariableDerivatives") ;
  }
  
  {
    double* x  = Model_GetVariable(model,n) ;
    double* dx = Model_GetVariableDerivative(model,n) ;
    int j ;
  
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] = x[j] ;
    }
  
    dx[i] += dxi ;
  
    {
      //Model_ComputeSecondaryVariables_t* computesecondaryvariables = Model_GetComputeSecondaryVariables(model) ;
      
      ComputeSecondaryVariables(el,t,dt,dx) ;
    }
  
    for(j = 0 ; j < NbOfVariables ; j++) {
      dx[j] -= x[j] ;
      dx[j] /= dxi ;
    }

    return(dx) ;
  }
}



int  ComputeSecondaryVariables(Element_t *el,double t,double dt,double *x)
{
  double c_s    = x[E_Salt] ; 
  double h_r    = x[E_Mass] ;
#if defined (E_Air)
  double p_g    = x[E_Air] ;
#else
  double p_g    = p_g0 ;
#endif
    
  /* Water concentration */
  double c_w    = (1. - V_AC*c_s)/V_H2O ;
  
  /* Activities of water and salt */
  double lna_w  = LogActivityOfWater(c_s) ;
  double lna_s  = LogActivityOfSalt(c_s) ;
  
  /* Pressures */
  double p_l    = LiquidPressure(h_r,lna_w) ;
  double p_v    = p_v0*h_r ;
  double p_a    = p_g - p_v ;
  double p_c    = p_g - p_l ;
  
  /* Mass densities */
  double rho_w  = M_H2O*c_w ;
  double rho_s  = M_SALT*c_s ;
  double rho_l  = rho_w + rho_s ;
  double rho_v  = MassDensityOfWaterVapor(p_v) ;
  double rho_a  = MassDensityOfDryAir(p_a) ;
  double rho_g  = rho_v + rho_a ;
  
  double c_v    = rho_v/rho_g ;
  
  /* Saturations */
  double s_l    = SaturationDegree(p_c) ;
  double s_g    = 1 - s_l ;
  
  double phi    = phi0 ;

  /* Mass/Molar contents */
  double n_s    = c_s*phi*s_l ;
  double m_l    = rho_l*phi*s_l ;
  double m_g    = rho_g*phi*s_g ;
  double m_v    = rho_v*phi*s_g ;
  double m_a    = rho_a*phi*s_g ;
  double m_t    = m_l + m_g ;

  /* backup */
  x[I_P_L      ] = p_l ;
  x[I_P_G      ] = p_g ;
  x[I_P_V      ] = p_v ;
  x[I_P_A      ] = p_a ;
  x[I_S_L      ] = s_l ;
  x[I_LNA_W    ] = lna_w ;
  x[I_LNA_S    ] = lna_s ;
  x[I_RHO_V    ] = rho_v ;
  x[I_RHO_L    ] = rho_l ;
  x[I_RHO_G    ] = rho_g ;
  x[I_RHO_S    ] = rho_s ;
  x[I_C_V      ] = c_v ;
  x[I_C_W      ] = c_w ;
  x[I_C_S      ] = c_s ;
  x[I_H_R      ] = h_r ;
  
  /* Fluid contents */
  x[I_N_S      ] = n_s ;
  x[I_M_A      ] = m_a ;
  x[I_M_T      ] = m_t ;
  x[I_M_L      ] = m_l ;
  x[I_M_G      ] = m_g ;
  x[I_M_V      ] = m_v ;
  
  return(0) ;
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
