#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "Context.h"
#include "CommonModel.h"
#include "FEM.h"
#include "Plasticity.h"

#define TITLE "Transversely isotropic ThermoPoroplasticity with hardening and hydration process(2021)"
#define AUTHORS "Lecomte-Pierre-Braun-Dangla"

#include "PredefinedMethods.h"

//Ce fichier est un mouveau modèle développé par Lecomte et Pierre pour prendre en compte l'hydratation dans le modèle thermo-poro-plast
//Il s'appuie sur un fichier déja existant par Braun et Dangla

/* Nb of equations */
#define NEQ     (dim+2)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (34)
#define NVE     (6)
#define NV0     (12)

// la partie index sert à definir l'emplacement du stockage des variables

/* Equation index */
#define E_mec   (0)
#define E_liq   (dim)     // équation de fluide, probablement bilan de masse ?
#define E_the   (dim+1)   // équation thermique

/* Unknown index */
#define U_u     (0)       // inconnues en deplacement
#define U_p_l   (dim)     // inconnue sur la pression de pores
#define U_tem   (dim+1)   // inconnue thermique

/* We define some names for implicit terms */
/* Using vim[i] for scalars and (vim+i) for tables as they will be called with VAR[i] making them lvalue */
#define M_L           (vim[0])   // masse de liquide
#define W_L           (vim+1)   // flux de liquide imposé, vecteur
#define SIG           (vim+4)   // contraintes, tenseur
#define F_MASS        (vim+13)  // forces volumiques, vecteur
#define EPS_P         (vim+16)  // déformations plastiques, tenseur
#define HARDV         (vim[25])  // potentiel plastique
#define CRIT          (vim[26])  // critère de plasticité

#define S_TOT         (vim[27]) // entropy, scalar
#define W_THE         (vim+28)  // heat flux, vector
#define PHI           (vim[31]) // Porosity, scalar
#define ALPHA         (vim[32]) // Hydration Degree, scalar
#define S_F           (vim[33]) // Saturation degree, scalar


#define M_L_n         (vim_n[0])  // idem mais pour l'itération précédante
#define SIG_n         (vim_n+4)
#define EPS_P_n       (vim_n+16)
#define HARDV_n       (vim_n[25])
#define CRIT_n        (vim_n[26])

#define S_TOTn        (vim_n[27]) // entropy, scalar
#define PHI_n         (vim_n[31]) // Porosity, scalar
#define ALPHA_n       (vim_n[32]) // Hydration Degree, scalar
#define S_F_n         (vim_n[33]) // Saturation degree, scalar


/* We define some names for explicit terms */
#define K_L           (vex+0)     // permeabilité, vecteur
#define K_T           (vex+3)     // thermal conductivity, vecteur

/* We define some names for constant terms */
#define SIG0          (vc+0)    // initial stresses
#define P_L0          (vc[9])    // initial porepressure
#define TEMP0         (vc[10])   // initial temperature
#define T0            (vc[11])   // initial time relative to hydration


/* Functions */
static int    pm(const char *s) ;
static void   GetProperties(Element_t*) ;

static int    ComputeTangentCoefficients(FEM_t*,double,double,double*) ;
static int    ComputeTransferCoefficients(FEM_t*,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,double,int) ;
static double* ComputeVariableDerivatives(Element_t*,double,double,double*,double,int) ;

static void  ComputeSecondaryVariables(Element_t*,double,double,double*) ;
static double* ComputeVariablesDerivatives(Element_t*,double,double,double*,double,int) ;

/* Functions relative to Hydration */

double AlphaDot(double, double, double) ;
double SaturationDegree(double, double) ;
double Permeability(double) ;
double ModulusFromHydration(double) ;
double PoissonFromHydration(double) ;
double BiotCoefficientFromHydration(double) ;
double Chi(double) ;

#define ComputeFunctionGradients(...)  Plasticity_ComputeFunctionGradients(plasty,__VA_ARGS__)
#define ReturnMapping(...)             Plasticity_ReturnMapping(plasty,__VA_ARGS__)
#define CopyElasticTensor(...)         Plasticity_CopyElasticTensor(plasty,__VA_ARGS__)
#define UpdateElastoplasticTensor(...) Plasticity_UpdateElastoplasticTensor(plasty,__VA_ARGS__)
#define rotate(a,i,j)      ((a)[axis_v[i]*3+axis_v[j]]) //Note : on a besoin d'une rotation éventuelle pour aligner le repère d'anisotropie et le repère de travail
#define rotate1(a,i)      ((a)[axis_v[i]])


/* Parameters */
//Note : le fichier est isotrope transverse de base, on considère un cas isotrope pour débuter, mais on garde la formulation anisotrope qui servira à Maxime dans le futur
//général
static double  gravite ;      // intensité du champ de pesanteur
static double  rho_s ;        // masse volumique du squelette solide

//porosité
static double  phi0 ;         // porosité initiale
static double  b ;            // coefficient de biot
static double  b_3 ;          // coefficient de biot direction transverse
static double  N ;            // module de biot
static double  k_int_1;       // coefficient de perméabilité
static double  k_int_3 ;      // coefficient de perméabilité direction transverse
static double  beta ;         // variable d'ecrouissage ?
static double  beta_ksi_phi ;
static double  beta_ksi_w ;
static double  sig0 ;

//plasticité
static double  hardv0 ;       // potentiel plastique initial
static Plasticity_t* plasty ; // pointeur vers le modèle de plasticité ?

//thermique
static double lam_1 ;         // Thermal conductivity of bulk (W/m/K)
static double lam_3 ;         // Thermal conductivity of bulk (W/m/K)
static double alpha_s ;       // Thermal volumetric expansion coefficient of solid (1/K)
static double alpha_1 ;       // Thermal linear expansion coefficient of solid (1/K)
static double alpha_3 ;       // Thermal linear expansion coefficient of solid (1/K)
static double C_d ;           // Specific volumetric heat capacity (J/K/m3)

static double s_tot0 ;        // initial bulk entropy
static double K_d ;           // drained bulk modulus

//elasticite
static double* cijkl ;        // pointeur vers le tenseur de rigidité
static double young ;         // Module d'Young
static double young_3 ;       // Module d'Young direction transverse
static double poisson ;       // coefficient de poisson
static double poisson_3 ;     // coefiicient de poisson direction transverse
static double shear_3 ;       // module de cisaillement plan d'anisotropie

//anisotropie
static short int axis_1,axis_2,axis_3 ;  //axes d'anisotropie du materiau
static short int axis_v[3] = {0,0,0} ; // repère d'anisotropie, vecteur
static double  c_a_11,c_b_11,c_a_33,c_b_33,c_a_22,c_b_22  ; //what is this ? Peut-etre une matrice de rotation ?

//hydratation
static double t0 ; // initial time relative to hydration
static double beta_alpha ; //paramètre d'ecrouissage chimique
static double alpha_0 ; //alpha0 parameter for hydration kinetics
static double n_alpha ; //n parameter in alpha dot
static double p_alpha ; //p parameter in alpha dot
static double A_alpha ; // A parameter in alpha dot
static double Ea_alpha ; // E_a parameter in alpha dot
static double DeltaV_alpha ; // Delta V parameter in alpha dot
static double k_alpha_Sf ; // k_alpha parameter in S_f
static double n_Sf ; // n parameter in S_f
static double p_r0_Sf  ; // p_r^0 paramater in S_f


//water properties, fitted on IAWPS curves using 3rd order polynoms (in function of presssure and temperature)
static double rho_l0 ;            // densité du fluide
static double bulk_l ;            // bulk modulus of liquid phase
static double mu_l ;              // viscosité
static double s_l = 367;          //entropy at 25°C (J/kg/K) //! this is overwritten later in the code by the IAPWS properties
static double alpha_l = 207.e-6 ; //volumetric th expansion coeff of water at 20°C//! this is overwritten later in the code by the IAPWS properties

static double wf1[] = {1.00E+03,    4.87E-01,    -2.84E-02,    -3.95E-04,    -1.82E-03,    -0.004974,    -0.000001122,    0.00001615,    0.00001032};
static double wf2[] = {1.56E-03,    -1.91E-05,    1.53E-02,    -9.65E-07,    -9.90E-06,    -0.00002634,    1.972E-08,    2.318E-08,    4.362E-08};
static double wf3[] = {1.975,    5.34E-03,    1.37E-02,    6.79E-06,    6.98E-07,    -0.0001817,    -6.954E-08,    1.331E-07,    5.115E-07};
static double wf4[] = {-3.03E-05,    2.90E-06,    1.33E-05,    -5.70E-09,    -7.33E-08,    -8.702E-08,    1.238E-10,    2.771E-10,    3.188E-10};
static double wf5[] = {4.20E+00,    -3.84E-03,    -8.08E-04,    1.08E-05,    4.32E-05,    0.000007446,    -7.117E-08,    -2.739E-07,    2.531E-08};
static double wf6[] = {1.63E-03,    -1.11E-06,    -3.83E-05,    3.92E-09,    3.64E-08,    4.123E-07,    -4.878E-11,    -2.316E-10,    -1.652E-09};
static double wf7[] = {0.5589,    0.000485,    0.002186,    -3.391E-08,    -0.000001363,    -0.00001008,    -1.051E-09,    2.034E-08,    1.882E-09};

#define dens_w(T,p)     (wf1[0]+wf1[1]*p*1.e-6+wf1[2]*T+wf1[3]*p*1.e-6*p*1.e-6+wf1[4]*p*1.e-6*T+wf1[5]*T*T+wf1[6]*p*1.e-6*p*1.e-6*T+wf1[7]*p*1.e-6*T*T+wf1[8]*T*T*T)     //T in C, p in Pa, return density in kg/m3
#define enth_w(T,p)    ((wf2[0]+wf2[1]*p*1.e-6+wf2[2]*T+wf2[3]*p*1.e-6*p*1.e-6+wf2[4]*p*1.e-6*T+wf2[5]*T*T+wf2[6]*p*1.e-6*p*1.e-6*T+wf2[7]*p*1.e-6*T*T+wf2[8]*T*T*T)*1000)  //T in C, p in Pa, return entropy in Pa m3/kg/K
#define bulk_w(T,p)    ((wf3[0]+wf3[1]*p*1.e-6+wf3[2]*T+wf3[3]*p*1.e-6*p*1.e-6+wf3[4]*p*1.e-6*T+wf3[5]*T*T+wf3[6]*p*1.e-6*p*1.e-6*T+wf3[7]*p*1.e-6*T*T+wf3[8]*T*T*T)*1000*1.e6) //T in C, p in Pa, return bulk mod in Pa
#define alpha_w(T,p)   (wf4[0]+wf4[1]*p*1.e-6+wf4[2]*T+wf4[3]*p*1.e-6*p*1.e-6+wf4[4]*p*1.e-6*T+wf4[5]*T*T+wf4[6]*p*1.e-6*p*1.e-6*T+wf4[7]*p*1.e-6*T*T+wf4[8]*T*T*T)         // in C, p in Pa, return alpha in -/°C
#define visc_w(T,p)    ((wf6[0]+wf6[1]*p*1.e-6+wf6[2]*T+wf6[3]*p*1.e-6*p*1.e-6+wf6[4]*p*1.e-6*T+wf6[5]*T*T+wf6[6]*p*1.e-6*p*1.e-6*T+wf6[7]*p*1.e-6*T*T+wf6[8]*T*T*T)) //T in C, p in Pa, return visc in Pa s



#define GetProperty(a)      Element_GetPropertyValue(el,a)

#define NbOfVariables     (96)  //! always verify
static double Variable[NbOfVariables] ;
static double Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;

#define I_U            (0)   //deplacements, vecteur
#define I_P_L          (3)   //pression de fluide, scalaire
#define I_EPS          (4)   //déformations, tenseur
#define I_SIG          (13)  //contraintes, tenseur
#define I_EPS_P        (22)  //deformations plastiques, tenseur
#define I_EPS_n        (31)  //déformation itération precedante
#define I_EPS_P_n      (40)  //déformation plastique itération precedante
#define I_Fmass        (49)  //forces volumiques, vecteur
#define I_M_L          (52)  //masse de liquide, scalaire
#define I_W_L          (53)  //flux de fluide, vecteur
#define I_HARDV        (56)  //potentiel plastique, scalaire
#define I_CRIT         (57)  //critère de plasticité, scalaire
#define I_RHO_L        (58)  //masse volumique du fluide, scalaire
#define I_PHI          (59)  //porosité, scalaire
#define I_P_L_n        (60)  //pression de fluide itération precedante
#define I_K_H          (61)  //coeff de perméabilité, vecteur
#define I_GRD_P_L      (64)  //gradient de pression de liquide ?
#define I_SIG_n        (67)  //contraintes itération precedante
#define I_HARDV_n      (77)  //potentiel plast itération precedante
#define I_S_TOT        (78)  //total entropy
#define I_S_L          (79)  //liquid entropy
#define I_TEM           (80)  //Temperature
#define I_W_THE         (81)  //Heat flux
#define I_GRD_TEM       (84)  //Temperature gradient
#define I_KTH           (87)  //Thermal conductivity, vector
#define I_TEM_n         (90)  //Temperature at previous timestep
#define I_PHI_n         (91)  //Porosity at previous timestep
#define I_ALPHA         (92)  //degré d'hydratation
#define I_ALPHA_n       (93)  //degre hydra iteration precedante
#define I_S_F           (94)  //degre de saturation
#define I_S_F_n         (95)  //degre de saturation t-1




int pm(const char *s)
{
         if(!strcmp(s,"gravity"))    { return (0) ;
  } else if(!strcmp(s,"young"))      { return (1) ;
  } else if(!strcmp(s,"poisson"))    { return (2) ;
  } else if(!strcmp(s,"porosity"))   { return (3) ;
  } else if(!strcmp(s,"rho_l"))      { return (4) ; //TODO OBSOLETE
  } else if(!strcmp(s,"k_int_1"))    { return (5) ;
  } else if(!strcmp(s,"mu_l"))       { return (6) ;
  } else if(!strcmp(s,"b"))          { return (7) ;
  } else if(!strcmp(s,"N"))          { return (8) ;
  } else if(!strcmp(s,"rho_s"))      { return (9) ;
  } else if(!strcmp(s,"beta"))       { return (10) ;
  } else if(!strcmp(s,"p_l0"))       { return (11) ; //TODO OBSOLETE
  } else if(!strcmp(s,"bulk_l"))     { return (12) ; //TODO OBSOLETE
  } else if(!strcmp(s,"sig0"))       { return(13) ;
  } else if(!strncmp(s,"sig0_",5))   {
    int i = (strlen(s) > 5) ? s[5] - '1' : 0 ;
    int j = (strlen(s) > 6) ? s[6] - '1' : 0 ;
    return(13 + 3*i + j) ;

  } else if(!strcmp(s,"harv0")) {
    return(22) ;

    /* Drucker-Prager */
  } else if(!strcmp(s,"initial cumulative plastic shear strain")) {
    return(22) ;
  } else if(!strcmp(s,"cohesion"))   {
    return (23) ;
  } else if(!strcmp(s,"friction"))   {
    return (24) ;
  } else if(!strcmp(s,"dilatancy"))  {
    return (25) ;

    /* Cam-clay */
  } else if(!strcmp(s,"initial plastic void ratio")) {
    return(22) ;
  } else if(!strcmp(s,"slope of swelling line")) {
    return(23) ;
  } else if(!strcmp(s,"slope of virgin consolidation line")) {
    return(24) ;
  } else if(!strcmp(s,"slope of critical state line"))  {
    return(25) ;
  } else if(!strcmp(s,"initial pre-consolidation pressure")) {
    return(26) ;
  } else if(!strcmp(s,"initial void ratio")) {
    return(27) ;

  } else if(!strcmp(s,"lam_1")) {     /* ADDED Thermal conductivity of bulk (W/m/K) */
    return(28) ;
  } else if(!strcmp(s,"C_d")) {       /* ADDED Specific volumetric heat capacity (J/K/m3) */
    return(29) ;
  } else if(!strcmp(s,"alpha_s")) {   /* ADDED Thermal volumetric expansion coefficient of solid (1/K) */ //TODO OBSOLETE
    return(30) ;
  } else if(!strcmp(s,"TEMP0")) {       /* ADDED initial temperature (C)*/  //TODO OBSOLETE
    return(31) ;
  } else if(!strcmp(s,"s_tot0")) {        /* ADDED initial bulk entropy  */
    return(32) ;
  } else if(!strcmp(s,"c_a_11")) {        /* ADDED  grad*/
    return(33) ;
  } else if(!strcmp(s,"c_b_11")) {        /* ADDED   const*/
    return(34) ;
  } else if(!strcmp(s,"c_a_33")) {        /* ADDED  grad*/
    return(35) ;
  } else if(!strcmp(s,"c_b_33")) {        /* ADDED   const*/
    return(36) ;
  }
  else if(!strcmp(s,"young_3"))   return (37) ;// ADDED
  else if(!strcmp(s,"poisson_3")) return (38) ;// ADDED
  else if(!strcmp(s,"shear_3"))   return (39) ;// ADDED
  else if(!strcmp(s,"axis_3"))    return (40) ;// ADDED
  else if(!strcmp(s,"b_3"))       return (41) ;// ADDED
  else if(!strcmp(s,"alpha_1"))   return (42) ;// ADDED
  else if(!strcmp(s,"alpha_3"))   return (43) ;// ADDED
  else if(!strcmp(s,"k_int_3"))   return (44) ;// ADDED
  else if(!strcmp(s,"lam_3"))   return (45) ;// ADDED
  else if(!strcmp(s,"ACC_k"))   return (46) ;// ADDED -----------------ACC
  else if(!strcmp(s,"ACC_M"))   return (47) ;// ADDED
  else if(!strcmp(s,"ACC_N"))   return (48) ;// ADDED
  else if(!strcmp(s,"initial_isotropic_tensile_limit"))   return (49) ;// ADDED
  else if(!strcmp(s,"initial_preconsolidation_pressure"))   return (50) ;// ADDED
  else if(!strcmp(s,"volumetric_strain_hardening_parameter"))   return (51) ;
  else if(!strcmp(s,"c_a_22")) return(52) ;     /* ADDED  grad*/
  else if(!strcmp(s,"c_b_22")) return(53) ;         /* ADDED   const*/
  else if(!strcmp(s,"thermal_hardening_parameter"))   return (54) ;


  else if(!strcmp(s,"t0"))   return (55) ;
  else if(!strcmp(s,"hydration_hardening_parameter"))   return (56) ;
  else if(!strcmp(s,"beta_ksi_phi"))   return (57) ;
  else if(!strcmp(s,"beta_ksi_w"))   return (58) ;
  else if(!strcmp(s,"alpha_0"))   return (59) ;
  else if(!strcmp(s,"n_alpha"))   return (60) ;
  else if(!strcmp(s,"p_alpha"))   return (61) ;
  else if(!strcmp(s,"A_alpha"))   return (62) ;
  else if(!strcmp(s,"Ea_alpha"))   return (63) ;
  else if(!strcmp(s,"DeltaV_alpha"))   return (64) ;
  else if(!strcmp(s,"k_alpha_Sf"))   return (65) ;
  else if(!strcmp(s,"n_Sf"))   return (66) ;
  else if(!strcmp(s,"p_r0_Sf"))   return (67) ;
  else if(!strcmp(s,"sig0"))   return (68) ;

  //! verify always with NbOfProp in ReadMatProp below
  else return(-1) ;
}


/* Functions relative to hydration */

double AlphaDot(double alpha, double P, double T)
{
  return( A_alpha * exp( - (Ea_alpha  + DeltaV_alpha * P)/8.31/T ) * pow(1-alpha,n_alpha) * pow( -log(1-alpha),p_alpha) ) ;
}

double SaturationDegree(double p_l, double alpha)
{
  return( (p_l < 0.) ? pow( 1 + pow( -p_l / p_r0_Sf ,n_Sf ) * exp(- k_alpha_Sf * n_Sf * alpha), (1-n_Sf)/n_Sf ) : 1. ) ;
}

double Permeability(double s_f)
{
  return( (s_f < 1.) ? pow( s_f, 0.5 ) * pow( 1 - pow( 1 - pow( s_f, n_Sf / (n_Sf - 1) ), (n_Sf - 1)/n_Sf ), 2) : 1. ) ;
}

double ModulusFromHydration(double alpha)
{
  return( (alpha > 0.03) ? (alpha-0.03)/0.97 : 5.e-2 ) ;
}

double PoissonFromHydration(double alpha)
{
  return( alpha > 0.03 ? 0.5 - 0.3*(alpha-0.03)/0.97 : 0.49 ) ;
}

double BiotCoefficientFromHydration(double alpha)
{
  return( alpha > 0.03 ? 1.-0.4*(alpha-0.03)/0.97 : 1. ) ;
}

double Chi(double Sf)
{
  return( Sf ) ;
}


void GetProperties(Element_t* el)
{
  gravite = GetProperty("gravity") ;
  phi0    = GetProperty("porosity") ;
  k_int_1   = GetProperty("k_int_1") ;
  k_int_3   = GetProperty("k_int_3") ;
  mu_l    = GetProperty("mu_l") ;
  rho_l0  = GetProperty("rho_l") ;         //TODO OBSOLETE
  bulk_l     = GetProperty("bulk_l") ;
  rho_s   = GetProperty("rho_s") ;
  //p_l0    = GetProperty("p_l0") ;
  b       = GetProperty("b") ;
  b_3       = GetProperty("b_3") ;
  N       = GetProperty("N") ;
  beta    = GetProperty("beta") ;
  //sig0    = &GetProperty("sig0") ;
  hardv0  = GetProperty("initial_preconsolidation_pressure") ;

  young   = GetProperty("young") ;      //ADDED only to calculate bulk modulus
  poisson = GetProperty("poisson") ;    //ADDED only to calculate bulk modulus
  K_d = young / 3 / (1 - 2 * poisson) ;       //ADDED only to calculate bulk modulus
  young_3   = GetProperty("young_3") ;//ADDED
  poisson_3 = GetProperty("poisson_3") ;//ADDED
  shear_3   = GetProperty("shear_3") ;//ADDED
  axis_3  = (short int) GetProperty("axis_3") - 1 ;//ADDED
  axis_1  = (axis_3 + 1) % 3 ;//ADDED
  axis_2  = (axis_3 + 2) % 3 ;//ADDED

  axis_v[0] = axis_1;
  axis_v[1] = axis_2;
  axis_v[2] = axis_3;

  lam_1     = GetProperty("lam_1");       /* ADDED Thermal conductivity of bulk (W/m/K) */
  lam_3     = GetProperty("lam_3");       /* ADDED Thermal conductivity of bulk (W/m/K) */
  C_d       = GetProperty("C_d");         /* ADDED Specific volumetric heat capacity (J/K/m3) */
  alpha_1   = GetProperty("alpha_1") ;    /* ADDED Thermal linear expansion coefficient of solid (1/K) */
  alpha_3   = GetProperty("alpha_3") ;    /* ADDED Thermal linear expansion coefficient of solid (1/K) */
  alpha_s   = 2*alpha_1 + alpha_3 ;       /* ADDED Thermal volumetric expansion coefficient of solid (1/K) */
  //TEMP0         = GetProperty("TEMP0") ;            /* ADDED initial temperature (C)*/
  s_tot0    = GetProperty("s_tot0") ;       /* ADDED initial bulk entropy  */
  c_a_11  = GetProperty("c_a_11") ;
  c_b_11  = GetProperty("c_b_11") ;
  c_a_33  = GetProperty("c_a_33") ;
  c_b_33  = GetProperty("c_b_33") ;
  c_a_22  = GetProperty("c_a_22") ;
  c_b_22  = GetProperty("c_b_22") ;

  t0  = GetProperty("t0") ;  // initial time relative to hydration
  beta_alpha = GetProperty("hydration_hardening_parameter") ;
  beta_ksi_phi = GetProperty("beta_ksi_phi") ; // Porosity evolution with hydration
  beta_ksi_w = GetProperty("beta_ksi_w") ; //  Consumed water volume with hydration
  alpha_0 = GetProperty("alpha_0") ; // alpha0 parameter for hydration kinetics
  n_alpha = GetProperty("n_alpha") ; // n paramater in alpha dot
  p_alpha = GetProperty("p_alpha") ; // p parameter in alpha dot
  A_alpha = GetProperty("A_alpha") ; // A parameter in alpha dot
  Ea_alpha = GetProperty("Ea_alpha") ; // E_a parameter in alpha dot
  DeltaV_alpha = GetProperty("DeltaV_alpha") ; // Delta V parameter in alpha dot
  k_alpha_Sf = GetProperty("k_alpha_Sf") ; // k_alpha parameter in S_f
  n_Sf = GetProperty("n_Sf") ; // n parameter in S_f
  p_r0_Sf = GetProperty("p_r0_Sf") ; // p_r^0 paramater in S_f
  sig0 = GetProperty("sig0") ;

  plasty = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
  }
}



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
  Model_CopyNameOfEquation(model,E_liq,"liq") ;
  Model_CopyNameOfEquation(model,E_the,"the") ;     /* ADDED */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_mec + i,name_eqn[i]) ;
  }

  /** Names of the main (nodal) unknowns */
  Model_CopyNameOfUnknown(model,U_p_l,"p_l") ;
  Model_CopyNameOfUnknown(model,U_tem,"tem") ;  /* ADDED */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_u + i,name_unk[i]) ;
  }

  Model_GetComputePropertyIndex(model) = pm ;

  return(0) ;
}



int ReadMatProp(Material_t* mat,DataFile_t* datafile)
/** Read the material properties in the stream file ficd
 *  Return the nb of (scalar) properties of the model */
{
  int  NbOfProp = 69 ;                              //! verify always with above
  int dim = Material_GetDimension(mat) ;
  int i ;

  /* Par defaut tout a 0 */
  for(i = 0 ; i < NbOfProp ; i++) Material_GetProperty(mat)[i] = 0. ;

    /* Pre-initialization */
  {

    Material_GetProperty(mat)[pm("axis_3")]  = dim ; //TODO this is needed?
  }


  Material_ScanProperties(mat,datafile,pm) ;


  /* Plasticity */
  {
    plasty = Plasticity_Create() ;

    Material_AppendData(mat,1,plasty,Plasticity_t,"Plasticity") ;
  }

  /* Elastic and plastic properties */
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;

    {
      /* Elasticity */ //ADDED transverse isotropy from CO2coal
      {
        young     = Material_GetProperty(mat)[pm("young")] ;
        poisson   = Material_GetProperty(mat)[pm("poisson")] ;
        young_3   = Material_GetProperty(mat)[pm("young_3")] ;
        poisson_3 = Material_GetProperty(mat)[pm("poisson_3")] ;
        shear_3   = Material_GetProperty(mat)[pm("shear_3")] ;
        axis_3    = Material_GetProperty(mat)[pm("axis_3")] - 1 ;

        //Elasticity_SetToIsotropy(elasty) ; //CHANGED to transversely isotropic
        //Elasticity_SetParameters(elasty,young,poisson) ;//CHANGED to transversely isotropic

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
          Elasticity_SetToTransverseIsotropy(elasty) ;
          Elasticity_SetParameters(elasty,young,poisson,young_3,poisson_3,shear_3,axis_3) ;
        }

        {
          double* c = Elasticity_GetStiffnessTensor(elasty) ;

          Elasticity_ComputeStiffnessTensor(elasty,c) ;
        }
      }


      /* ACC */
      {
        double ACC_k  = Material_GetPropertyValue(mat,"ACC_k") ;
        double ACC_M = Material_GetPropertyValue(mat,"ACC_M") ;
        double ACC_N      = Material_GetPropertyValue(mat,"ACC_N") ;
        double ACC_pt0    = Material_GetPropertyValue(mat,"initial_isotropic_tensile_limit") ;
        double ACC_pc0    = Material_GetPropertyValue(mat,"initial_preconsolidation_pressure") ;
        double ACC_beta_eps   = Material_GetPropertyValue(mat,"volumetric_strain_hardening_parameter") ;
        double ACC_beta_TEMP0  = Material_GetPropertyValue(mat,"thermal_hardening_parameter") ;
        double ACC_beta_hyd = Material_GetPropertyValue(mat,"hydration_hardening_parameter") ;


        Plasticity_SetTo(plasty,ACCPierre) ;
        Plasticity_SetParameters(plasty,ACC_k,ACC_M,ACC_N,ACC_pt0,ACC_pc0,ACC_beta_eps,ACC_beta_TEMP0,ACC_beta_hyd) ;

      }

    }

  }

  return(NbOfProp) ;
}


int PrintModelProp(Model_t* model,FILE *ficd)
/** Print the model properties
 *  Return the nb of equations */
{

  int dim = Model_GetDimension(model) ;

  printf(TITLE) ;
  printf("\n") ;

  if(!ficd) return(0) ;

  printf("\n") ;

  printf("This model consists in %d equations:\n",NEQ) ;
  printf("\t 1. Conservation of fluid mass   (mass)\n") ;
  printf("\t 2. Entropy balance        (the)\n") ;          /* ADDED */
  printf("\t 3. Mechanical equilibrium   (meca1,meca_2,meca_3)\n") ;

  printf("\n") ;

  printf("The primary unknowns are:\n") ;
  printf("\t 1. The liquid pressure     (p_l)\n") ;
  printf("\t 2. The temperature         (tem)\n") ;         /* ADDED */
  printf("\t 3. The displacement vector (u_1,u_2,u_3)\n") ;

  printf("\n\
  Example of input data\n\n") ;


  fprintf(ficd,"gravity = 0       # gravity\n") ;
  fprintf(ficd,"rho_s = 2350      # mass density of solid skeleton\n") ;
  fprintf(ficd,"young = 5.8e+09   # Young's modulus\n") ;
  fprintf(ficd,"poisson = 0.3     # Poisson's ratio\n") ;
  fprintf(ficd,"young_3 = 2.71e9  # Young's modulus of coal in axis 3\n") ; //ADDED
  fprintf(ficd,"poisson_3 = 0.26  # Poisson's ratio of coal in axis 3\n") ;  //ADDED
  fprintf(ficd,"shear_3 = 4.44e8  # Shear modulus in axis 3\n") ;           //ADDED
  fprintf(ficd,"axis_3 = 2        # Actual axis 3\n") ;                     //ADDED
  fprintf(ficd,"porosity = 0.15   # porosity\n") ;
  fprintf(ficd,"rho_l = 1000      # mass density of fluid\n") ;
  fprintf(ficd,"p_l0 = 4.7e+06    # initial pressure of fluid\n") ;   //TODO OBSOLETE
  fprintf(ficd,"p_g = 0           # gas pressure\n") ;
  fprintf(ficd,"bulk_l = 2e+09    # bulk modulus of fluid\n") ;
  fprintf(ficd,"k_int = 1e-19     # intrinsic permeability\n") ;
  fprintf(ficd,"mu_l = 0.001      # viscosity of liquid\n") ;
  fprintf(ficd,"b = 0.8           # Biot's coefficient\n") ;
  fprintf(ficd,"b_3 = 0.8         # Biot's coefficient\n") ;
  fprintf(ficd,"N = 4.e-11        # compressibility of pores\n") ;
  fprintf(ficd,"cohesion = 1e+06  # cohesion\n") ;
  fprintf(ficd,"friction = 25     # friction angle\n") ;
  fprintf(ficd,"dilatancy = 25    # dilatancy angle \n") ;
  fprintf(ficd,"beta = 0.8        # plastic Biot's coefficient\n") ;
  //fprintf(ficd,"sig0_ij = -11.5e6 # initial stress sig0_ij\n") ;
  fprintf(ficd,"lam_d = 1.        # thermal conductivity of bulk (W/m/K) \n") ;           /* ADDED */
  fprintf(ficd,"C_d = 2e+06       # volumetric heat of bulk (J/m3/K) \n") ;               /* ADDED */
  fprintf(ficd,"alpha_s = 54.e-6  # thermal volumetric expansion of the solid matrix (1/K) \n") ;   /* ADDED */
  fprintf(ficd,"TEMP0 = 25        # initial temperature (C) \n") ;                        /* ADDED */  //TODO OBSOLETE
  fprintf(ficd,"t0 = 14400        # initial time relative to hydration \n") ;
  fprintf(ficd,"beta_alpha = 7.e+6 # Hydration hardening coefficient (Pa) \n") ;



  fprintf(ficd,"Curves = my_file  # file name: p_c S_l k_rl\n") ;

  return(0) ;
}


int DefineElementProp(Element_t* el,IntFcts_t* intfcts)
/** Define some properties attached to each element
 *  Return 0 */

{
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  /* Continuity of pressure across zero-thickness element */
  {
    if(Element_HasZeroThickness(el)) {
      //Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,U_p_l) ;
      Element_MakeUnknownContinuousAcrossZeroThicknessElement(el,"p_l") ;
      Element_MakeEquationContinuousAcrossZeroThicknessElement(el,"liq") ;
    }
  }

  {
    IntFct_t* intfct = Element_GetIntFct(el) ;
    int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) + 1; //ADDED in frostaco and CO2coal there is  +1

    /** Define the length of tables */
    Element_GetNbOfImplicitTerms(el) = NVI*NbOfIntPoints ;
    Element_GetNbOfExplicitTerms(el) = NVE*NbOfIntPoints ;
    Element_GetNbOfConstantTerms(el) = NV0*NbOfIntPoints ; //ADDED
  }

  return(0) ;
}



int  ComputeLoads(Element_t* el,double t,double dt,Load_t* cg,double* r)
/** Compute the residu (r) due to loads
 *  Return 0 if succeeded and -1 if failed */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int dim = Element_GetDimensionOfSpace(el) ;
  IntFct_t* fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ; //ADDED
  int ndof = nn*NEQ ;
  FEM_t* fem = FEM_GetInstance(el) ;
  int    i ;
  double*  vc0    = Element_GetConstantTerm(el) ; //added
  double*  vc    = Element_GetConstantTerm(el) ;

  {
    double* r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;

    for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;

    /* thermic */  //ADDED
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_the) {
      for(i = 0 ; i < nn ; i++) {

        double tem = Element_GetValueOfNodalUnknown(el,u,i,U_tem) ;

        //R(i,E_the) /= tem_0 ;//!

      }
    } else if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_liq) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;

    /* other */
    } else {
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }

  }

  return(0) ;
#undef R
}


int ComputeInitialState(Element_t* el)
{
  double* vim0  = Element_GetImplicitTerm(el) ;
  double* vex0  = Element_GetExplicitTerm(el) ;
  double* vc0  = Element_GetConstantTerm(el) ;   //ADDED
  double* vc  = Element_GetConstantTerm(el) ;   //ADDED
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  DataFile_t* datafile = Element_GetDataFile(el) ;
  int    p ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ; //ADDED
  FEM_t* fem = FEM_GetInstance(el) ; //ADDED
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  //ADDED calculate average y coordiante for this element through the mean of the nodes
  int nnodes = Element_GetNbOfNodes(el)    ;
  double depth = 0. ; //ADDED average y coordinate of the element
  double* coord ;
  int i ;
  // for(i = 0 ; i < nnodes ; i++ ) {
  //   coord = Element_GetNodeCoordinate(el,i);
  //   depth += coord[1]/nnodes ;
  // }


  /* Pre-initialization */
  for(p = 0 ; p < NbOfIntPoints ; p++) {

    double p_l0 = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;//ADDED
    double tem_0  = FEM_ComputeUnknown(fem,u,intfct,p,U_tem) ; //ADDED


    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;

      /* Initial density */ //ADDED
      PHI = phi0;

      double alpha_ini = alpha_0 ;

      for( i=0; i < t0 ; i++ ){
        alpha_ini += 1./3600. * AlphaDot(alpha_ini,p_l0,273. + tem_0) ;
      }

      ALPHA = alpha_ini ;

      S_F = 1.0 ;

      S_TOT = - 1.e6 * rho_s * alpha_0 / tem_0 ;


      /* Initial stresses */
      //calculate coordinates based on gradients given as material parameters and the depth
      coord = Element_GetNodeCoordinate(el,p);
      depth = coord[1];
       /* How to account for partial initialization? */
      if(DataFile_ContextIsPartialInitialization(datafile)) {
       for(i = 0 ; i < 9 ; i++) SIG0[i] = SIG[i] ;
      } else {
       //for(i = 0 ; i < 9 ; i++) SIG0[i] = sig0[i] ;
       // SIG0[0] = depth*c_a_11+c_b_11 ;
       // SIG0[4] = depth*c_a_22+c_b_22 ;
       // SIG0[8] = depth*c_a_33+c_b_33 ;
       SIG0[0] = sig0 ;
       SIG0[4] = sig0 ;
       SIG0[8] = sig0 ;
       for(i = 0 ; i < 9 ; i++) SIG[i]  = SIG0[i] ;

       HARDV = hardv0 ;

      }

      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = 0. ;

    }



    /* storage in vc */
    {
      double* vc  = vc0 + p*NV0 ;
      P_L0 = p_l0 ;//ADDED
      TEMP0 = tem_0 ; //ADDED

    }

  }

  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,0,0,p) ;

    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;

      M_L = x[I_M_L] ;

      PHI = x[I_PHI] ;

      // S_F = x[I_S_F] ;

      // ALPHA = x[I_ALPHA] ;

      S_TOT = x[I_S_TOT] ;                                /*ADDED entropy */

      for(i = 0 ; i < 3 ; i++) W_L[i] = x[I_W_L + i] ;

      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;

      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;

      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;

      for(i = 0 ; i < 3 ; i++) W_THE[i]  = x[I_W_THE + i] ;   /*ADDED heat flux */

      CRIT = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;
      ALPHA = x[I_ALPHA] ;
    }




    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
      double rho_l = x[I_RHO_L] ;
      double tem = x[I_TEM] ;
      double tem_0 = TEMP0 ;
      double p_l = x[I_P_L] ;


      mu_l = visc_w(tem,p_l) ;   //! nonlinearized

      double k_l[3] = {0,0,0} ;
      rotate1(k_l,0) = rho_l*k_int_1/mu_l ;
      rotate1(k_l,1) = rho_l*k_int_1/mu_l ;
      rotate1(k_l,2) = rho_l*k_int_3/mu_l ;
      //double k_l = rho_l*k_int/mu_l ;

      double k_t[3] = {0,0,0} ;
      rotate1(k_t,0) = lam_1;
      rotate1(k_t,1) = lam_1;
      rotate1(k_t,2) = lam_3 ;
      //double k_t = lam_d;

      for(i = 0 ; i < 3 ; i++)  K_L[i] = k_l[i] ;   //changed to vector
      for(i = 0 ; i < 3 ; i++)  K_T[i] = k_t[i]/TEMP0;      //changed to vector
      //TODO here we take the initial temperature, as this is the initialization
    }


  }
  return(0) ;
}


int  ComputeExplicitTerms(Element_t* el,double t)
/** Compute the explicit terms */
{
  double* vim0 = Element_GetPreviousImplicitTerm(el) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  double* vc   = Element_GetConstantTerm(el) ; //ADDED
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    p ;

  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;

  /*
    Input data
  */
  GetProperties(el) ;

  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {

    /* read initial temperature */   //ADDED
    double tem_0 = TEMP0;

    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,t,0,p) ;

    double tem = x[I_TEM]; //ADDED
    double p_l = x[I_P_L]; //ADDED


    /* fluid mass density */
    double rho_l = x[I_RHO_L] ;
    mu_l = visc_w(tem,p_l) ;   //! nonlinearized

    // printf("rho_l = %e, mu_l = %e, k(S_r) = %e, k(alpha) = %e\n", rho_l, mu_l, Permeability(x[I_S_F]), pow(10,-11-3*x[I_ALPHA]));

    double k_l[3] = {0,0,0} ;
    k_l[0] = rho_l*pow(10,-16-3*x[I_ALPHA])/mu_l * Permeability(x[I_S_F]) ;
    k_l[1] = rho_l*pow(10,-16-3*x[I_ALPHA])/mu_l * Permeability(x[I_S_F]) ;
    k_l[2] = rho_l*pow(10,-16-3*x[I_ALPHA])/mu_l * Permeability(x[I_S_F]) ;

    // rotate1(k_l,0) = rho_l*k_int_1/mu_l ;
    // rotate1(k_l,1) = rho_l*k_int_1/mu_l ;
    // rotate1(k_l,2) = rho_l*k_int_3/mu_l ;
    //double k_l = rho_l*k_int/mu_l ;

    double k_t[3] = {0,0,0} ;
    rotate1(k_t,0) = lam_1;
    rotate1(k_t,1) = lam_1;
    rotate1(k_t,2) = lam_3 ;
    //double k_t = lam_d;

    /* storage in vex */
    {
      double* vex  = vex0 + p*NVE ;
      int i;
      for(i = 0 ; i < 3 ; i++)  K_L[i] = k_l[i] ;   //changed to vector
      for(i = 0 ; i < 3 ; i++)  K_T[i] = k_t[i]/tem_0;      //changed to vector
      int m=0; //! seem to work
      //TODO here the current temperature is used as in frostaco
    }
  }

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
  GetProperties(el) ;

  /* Loop on integration points */
  for(p = 0 ; p < NbOfIntPoints ; p++) {
    /* Variables */

    double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ;

    /* storage in vim */
    {
      double* vim  = vim0 + p*NVI ;
      int    i ;

      M_L = x[I_M_L] ;
      PHI = x[I_PHI] ;
      ALPHA = x[I_ALPHA] ;
      S_F = x[I_S_F];
      S_TOT   = x[I_S_TOT] ;    /* ADDED entropy*/

      for(i = 0 ; i < 3 ; i++) W_L[i] = x[I_W_L + i] ;

      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;

      for(i = 0 ; i < 3 ; i++) F_MASS[i] = x[I_Fmass + i] ;

      for(i = 0 ; i < 9 ; i++) EPS_P[i]  = x[I_EPS_P + i] ;

      for(i = 0 ; i < 3 ; i++) W_THE[i]  = x[I_W_THE + i] ;   /*ADDED heat flux */


      CRIT = x[I_CRIT] ;
      HARDV = x[I_HARDV] ;
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
  GetProperties(el) ;


  /*
  ** Poromechanics matrix
  */
  {
    double c[IntFct_MaxNbOfIntPoints*121] ;   //! changed
    int dec = ComputeTangentCoefficients(fem,t,dt,c) ; //t ADDED
    double* kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,2) ;  //the last int describes fluid mass and temperature
    /* The matrix kp is stored as (u for displacement, p for pressure, t for temperature)
     * | Kuu(9x9) Kup(9x1) Kut(9x1)|
     * | Kpu(1x9) Kpp(1x1) Kpt(1x1)|
     * | Ktu(1x9) Ktp(1x1) Ktt(1x1)|
     * i.e. the displacements u are in the positions 0 to dim-1 and
     * the pressure p is in the position dim, t in dim+1.
     * So we need to store the matrix by accounting for the right indexes.
     *
     * the arguments are FEM_ComputePoroelasticMatrix6(FEM_t* fem,IntFct_t* fi,const double* c,const int dec,const int n_dif,const int idis)
     * Return a pointer on a FE poroelastic matrix (Ndof x Ndof).
     *
     *  Ndof = nb of degrees of freedom (= NN*Neq)
     *  NN   = nb of nodes
     *  Neq  = nb of equations (= Dim + n_dif)
     *  Dim  = dimension of the problem
     *
     *  The inputs are:
     *
     *  n_dif = nb of Biot-like coupling terms (pressure, temperature, etc...)//? in our case this is pressure and temperature
     *
     *  idis = position index of the first displacement in the unknown vector //? this is generally by default 0
     *
     *  c = the entry matrix which should be given in the following order:
     *
     *  K0 to Kn then A0 to An etc.. with n = n_dif: //? in our case we have 9*9
     *
     *  | K0(9x9) K1(9x1) K2(9x1) ... |
     *  | A0(1x9) A1(1x1) A2(1x1) ... |
     *  | B0(1x9) B1(1x1) B2(1x1) ... |
     *  | ........................... |
     *
     *   K0    = Stiffness matrix                   //? in our case we have 9*9
     *   Kn    = Mechanic-hydraulic coupling terms  //? in our case we have 2*(9*1)
     *   A0,B0 = Hydraulic-mechanic coupling terms  //? in our case we have 2*(9*1)
     *   Ai,Bj = Hydraulic terms                    //? in our case we have 4*(1*1)
     *                                              //? total 9*9+2*2*9+4=121
     *
     *  The outputs is provided in an order depending on "idis".
     *  The first displacement unknown U1 will be positionned at "idis".
     *  (1  ... idis ... Neq)
     *   |  ...  |   ...  |
     *   v       v        v
     *  (P1 ...  U1  ...  Pn)
     */

    #define KP(i,j)   (kp[(i)*ndof + (j)])

    if(E_mec == 0) {
      int i ;

      for(i = 0 ; i < ndof*ndof ; i++) {
        k[i] = kp[i] ;
      }
    } else {
      int n ;

      for(n = 0 ; n < nn ; n++) {
        int m ;
        /* 2 nested loops (n and m) over all nodes */
        //in the K and KP matrices, always add n*NEQ and m*NEQ in the indices
        //the matrices K start always at (E ... , U ... )
        //the matrices KP start at (x ... , x ... ), wher x is 0 (meca), dim (p_l) or dim+1 (thermal)
        //then loop over correct indexes i (and j)
        // KP always has a fixed arrangement starting with mechanical matrix, then hydraulic coupling, then thermal coupling
        // KP coefficients are the C coefficients (tangent matrices, e.g. stiffness matrix) multiplied by the jacobian matrices
        // KP is then rearranged into K, where K has the order given by the definition of E_mech, E_hyd and E_term defined by the user

        for(m = 0 ; m < nn ; m++) {

          /* Mechanics */
          {
            int i ;

            /* Stiffness matrix */
            for(i = 0 ; i < dim ; i++) {
              int j ;

              for(j = 0 ; j < dim ; j++) {
                K(E_mec + i + n*NEQ,U_u + j + m*NEQ) = KP(0 + i + n*NEQ,0 + j + m*NEQ) ;
              }
            }

            /* Coupling matrix with p*/
            for(i = 0 ; i < dim ; i++) {
              K(E_mec + i + n*NEQ,U_p_l + m*NEQ) = KP(0 + i + n*NEQ,dim + m*NEQ) ;
            }

            /* Coupling matrix with T*/    //ADDED
            for(i = 0 ; i < dim ; i++) {
              K(E_mec + i + n*NEQ,U_tem + m*NEQ) = KP(0 + i + n*NEQ,dim + 1 + m*NEQ) ;
            }
          }

          /* Hydraulics */
          {
            int j ;

            /* Coupling matrix with u*/
            for(j = 0 ; j < dim ; j++) {
              K(E_liq + n*NEQ,U_u + j + m*NEQ) = KP(dim + n*NEQ,0 + j + m*NEQ) ;
            }

            /* Storage matrix */
            K(E_liq + n*NEQ,U_p_l + m*NEQ) = KP(dim + n*NEQ,dim + m*NEQ) ;

            /* Coupling matrix with T*/   //ADDED
            K(E_liq + n*NEQ,U_tem + m*NEQ) = KP(dim + n*NEQ,dim + 1 + m*NEQ) ;
          }


          /* Thermal ADDED*/
          {
            int j ;
            /* Storage matrix */
            K(E_the + n*NEQ,U_tem + m*NEQ) = KP(dim + 1 + n*NEQ,dim + 1 + m*NEQ) ;

            /* Coupling matrix with u */
            for(j = 0 ; j < dim ; j++) {
              K(E_the + n*NEQ,U_u + j + m*NEQ) = KP(dim + 1 + n*NEQ,0 + j + m*NEQ) ;
            }

            /* Coupling matrix with p */
            K(E_the + n*NEQ,U_p_l + m*NEQ) = KP(dim + 1 + n*NEQ,dim + m*NEQ) ;
          }


        }
      }
    }
    #undef KP
  }

  /*
  ** Conduction Matrix
  */
  {
    int n = 121; //ADDED here we have a (3 fluid + 3 thermal) = 6x6 conduction matrix //! here 6*6 or 11*11?
    double c[IntFct_MaxNbOfIntPoints*n] ;   /* ADDED changed number of terms*/
    int dec = ComputeTransferCoefficients(fem,dt,c) ;
    double* kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
    double* kc_the  = FEM_ComputeConductionMatrix(fem,intfct,c+9*3,dec) ; /* ADDED */


    int    i ;

    for(i = 0 ; i < nn ; i++) {
      int    j ;

      for(j = 0 ; j < nn ; j++) {
        K(E_liq + i*NEQ,U_p_l + j*NEQ) += dt*kc[i*nn + j] ;
        K(E_the  + i*NEQ,U_tem  + j*NEQ) += dt*kc_the[i*nn + j] ; /* ADDED */
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
  double* vim_1 = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
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
    double* vim = vim_1 ;
    double* rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;

    for(i = 0 ; i < nn ; i++) {
      int j ;

      for(j = 0 ; j < dim ; j++) R(i,E_mec + j) -= rw[i*dim + j] ;
    }

  }

  /* 1.2 Body forces */
  {
    double* vim = vim_1 ;
    double* rbf = FEM_ComputeBodyForceResidu(fem,intfct,F_MASS + dim - 1,NVI) ;

    for(i = 0 ; i < nn ; i++) {
      R(i,E_mec + dim - 1) -= -rbf[i] ;
    }

  }


  /* 2. Hydraulics */

  /* 2.1 Accumulation Terms */
  {
    double* vim = vim_1 ;
    double* vim_n = Element_GetPreviousImplicitTerm(el) ;
    double g1[IntFct_MaxNbOfIntPoints] ;

    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) {g1[i] = M_L - M_L_n;}

    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;

      for(i = 0 ; i < nn ; i++) R(i,E_liq) -= ra[i] ;
    }
  }

  /* 2.2 Transport Terms */
  {
    double* vim = vim_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_L,NVI) ;

    for(i = 0 ; i < nn ; i++) R(i,E_liq) -= -dt*rf[i] ;
  }










  /* 3. Entropy balance ADDED*/

  /* 3.1 Accumulation Terms ADDED*/
  {
    double* vim = vim_1 ;
    double* vim_n = Element_GetPreviousImplicitTerm(el) ; // already defined before
    double g1[IntFct_MaxNbOfIntPoints] ;

    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = S_TOT - S_TOTn ;



    {
      double* ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;

      for(i = 0 ; i < nn ; i++) R(i,E_the) -= ra[i] ;
    }
  }

  /* 3.2 Transport Terms  ADDED*/
  {
    double* vim = vim_1 ;
    double* rf = FEM_ComputeFluxResidu(fem,intfct,W_THE,NVI) ;

    for(i = 0 ; i < nn ; i++) R(i,E_the) -= -dt*rf[i] ;
  }





  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 18 ; //ADDED, +3 due to test of initial conditions //! check number !
  double* vex  = Element_GetExplicitTerm(el) ;
  double* vim  = Element_GetCurrentImplicitTerm(el) ;
  double* vim_n  = Element_GetPreviousImplicitTerm(el) ;
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
  plasty = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double modulus = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  //double ACC_M = Plasticity_GetACC_M(plasty);
  double ACC_M = Plasticity_GetParameter(plasty)[1] ;

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int p = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    /* Pressure */
    double p_l = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    /* Displacement */
    double dis[3] = {0,0,0} ;
    /* strains */
    double eps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    //double eps_tot[9] = {0,0,0,0,0,0,0,0,0} ;
    double w_l[3] = {0,0,0} ;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double hardv = 0 ;
    double crit = 0 ;
    double k_l = K_L[0] ;


    double w_the[3] = {0,0,0} ; /* ADDED heat flux */
    double k_t[3] = {0,0,0} ;  /* ADDED thermal conductivity */
    double tem = FEM_ComputeUnknown(fem,u,intfct,p,U_tem) ;  /* ADDED */

    double alpha = ALPHA;
    double sf = S_F;
    double phi = PHI;
    double m_l = M_L;

    int    i ;
    for(i = 0 ; i < dim ; i++) {
      dis[i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }

    double* eps_tot = FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;


    /* Averaging */
    for(i = 0 ; i < np ; i++ , vim += NVI , vex += NVE ) {


      int j ;

      for(j = 0 ; j < 3 ; j++) w_l[j] += W_L[j]/np ;

      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;

      for(j = 0 ; j < 9 ; j++) eps_p[j] += EPS_P[j]/np ;

      for(j = 0 ; j < 3 ; j++) w_the[j]  += W_THE[j]/np ;    /* ADDED heat flux */

      hardv += HARDV/np ;
      crit += CRIT/np ;
    }

    double pp    = (sig[0] + sig[4] + sig[8])/3. + p_l ;
    double q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;

    i = 0 ;
    Result_Store(r + i++,dis      ,"Displacements",3) ;
    Result_Store(r + i++,sig      ,"Stresses",9) ;
    Result_Store(r + i++,eps_tot     ,"Strains",9) ;
    Result_Store(r + i++,eps_p    ,"Plastic_strains",9) ;
    Result_Store(r + i++,&phi     ,"Porosity",1) ;
    Result_Store(r + i++,&p_l     ,"Pore_pressure",1) ;
    Result_Store(r + i++,&sf    ,"Saturation Degree",1) ;
    Result_Store(r + i++,&m_l     ,"Fluid mass",1) ;
    Result_Store(r + i++,&tem     ,"Temperature",1) ;
    Result_Store(r + i++,&alpha     ,"Hydration_Degree",1) ;
    Result_Store(r + i++,&hardv   ,"Hardening_variable",1) ;
    Result_Store(r + i++,&ACC_M   ,"M",1) ;
    Result_Store(r + i++,&crit    ,"Yield_criterion",1) ;
    Result_Store(r + i++,&modulus    ,"Young Modulus",1) ;
    Result_Store(r + i++,&poisson    ,"Poisson's ratio",1) ;
    Result_Store(r + i++,&pp    ,"Pprime",1) ;
    Result_Store(r + i++,&q    ,"q",1) ;
    Result_Store(r + i++,&k_l    ,"Permeability",1) ;
  }

  return(NbOfOutputs) ;
}



int ComputeTangentCoefficients(FEM_t* fem,double t,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
  //#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
//#define T4(a,i,j,k,l)  ((a)[((AXIS(i)*3+AXIS(j))*3+AXIS(k))*3+AXIS(l)])
//#define T4(a,i,j,k,l)  ((a)[((AXIS(i)*3+AXIS(j))*3+AXIS(k))*3+AXIS(l)])
//#define T2(a,i,j)      ((a)[AXIS(i)*3+AXIS(j)])
#define T2(a,i,j)      ((a)[axis_v[i]*3+axis_v[j]])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define B1(i,j)        T2(c1,i,j)
  Element_t* el  = FEM_GetElement(fem) ;
  double*  vim0  = Element_GetCurrentImplicitTerm(el) ;
  double*  vim_n = Element_GetPreviousImplicitTerm(el) ;
  double*  vc    = Element_GetConstantTerm(el) ;              //ADDED
  //  double*  vex0  = Element_GetExplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  double dui[NEQ] ; //ADDED for numerical derivative

  int    dec = 121 ;   //CHANGED
  int    p ;
  double zero = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ; //ADDED

  ObVal_t* obval = Element_GetObjectiveValue(el) ;

  {
    int i ;

    for(i = 0 ; i < NEQ ; i++) {
      dui[i] =  1.e-4*ObVal_GetValue(obval + i) ;
    }
  }



  for(p = 0 ; p < np ; p++) {
    double* vim  = vim0 + p*NVI ;
    double* c0 = c + p*dec ;

    /* Initial state variables */
    double* x = ComputeVariables(el,u,u_n,vim_n,t,dt,p) ; //ADDED
    double p_l0 = P_L0 ; //ADDED
    double tem_0  = TEMP0 ;  //ADDED


    /* Pressure */
    double pl  = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;
    double tem = FEM_ComputeUnknown(fem,u,intfct,p,U_tem) ; //ADDED
    //double pl  = x[I_P_L] ;

    /* Yield and potential function gradients */
    double fcg = 0 ;

    /* Criterion */
    double crit = CRIT ;

    /* initialization */
    {
      int i ;

      for(i = 0 ; i < dec ; i++) c0[i] = zero ;
    }

    /* The derivative of equations w.r.t unknowns */

    /* Derivatives w.r.t strains */
    {
      /* Tangent stiffness matrix: stress w.r.t strains*/
      {
        double sig[9] ;
        int i ;
        for(i = 0 ; i < 9 ; i++) sig[i] = SIG[i] ;

        /* Elastic effective stresses */
        sig[0] += beta*pl ;
        sig[4] += beta*pl ;
        sig[8] += beta*pl ;

        /* Tangent stiffness matrix */
        double* c1 = c0 ;

        //OLD CopyElasticTensor(c1) ;
        for(i = 0 ; i < 81 ; i++) {
          c1[i] = cijkl[i] ;
        }

        {
          /* Criterion */ //this guy activates when we are already in plasticity in the current state
          //replaces the elastic stiffness tensor with the elastoplastic one


          if(crit >= 0.) {
            double hardv[2] = {HARDV,ALPHA} ; //TODO update here for temperature effect???
            double* crit1 = ComputeFunctionGradients(sig,hardv) ;

            fcg = UpdateElastoplasticTensor(c1) ;//!

            if(fcg < 0) return(-1) ;
          }
        }
      }

      /* Coupling matrices: Variables w.r.t strains */ //ADDED all new
      {
        double deps = dui[U_u] ;  //increments strains with this

        int i ;
        for(i = 0 ; i < 3 ; i++) {
          int i_strain = I_EPS+i*3+i ;

          double* dx = ComputeVariableDerivatives(el,t,dt,x,deps,i_strain) ;
          //loop has to go through all strain components due to anisotropy
          //x contains all current variables which has been calculated before
          //deps is the ammount which a selected variable is incremented
          //I_EPS is the variable which is incremented, here we take only normal strains at the diagonal
          //dx is then the numerical derivative dx=delta_all_variables/delta_incremented_variable

          /* Biot's coupling matrix: fluid mass w.r.t strains */
          {
            double* c1 = c0 + 81 + 2*9 ;
            B1(i,i) = dx[I_M_L] ; //assigns only values on the diagonal

          }

          /* Coupling matrix: entropy w.r.t strains */
          {
            double* c1 = c0 + 81 + 9*3 + 1 + 1 ;
            B1(i,i) = dx[I_S_TOT] ; //assigns only values on the diagonal
          }
        }
      }
    }

    /* Derivatives w.r.t pore pressure */
    {
      double dp = dui[U_p_l] ;  //increments pore pressure with this
      double* dx = ComputeVariableDerivatives(el,t,dt,x,dp,I_P_L) ;

      /* Biot's coupling matrix: stress w.r.t pp*/
      {
        double* c1 = c0 + 81 ;
        int i;
        // B1(0,0) = b ;
        // B1(1,1) = b_3 ;
        // B1(2,2) = b ;
        for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_SIG+i*3+i] ;
        //this is only b
      }

      /* Storage term: fluid mass w.r.t pp */
      {
        double* c1 = c0 + 81 + 9 * 3 ;
        c1[0] = dx[I_M_L] ;
      }

      /* Coupling matrix: entropy w.r.t pp */
      {
        double* c1 = c0 + 81 + 9*3 + 1 + 1 + 9 ;
        c1[0] = dx[I_S_TOT] ;
      }
    }

    /* Derivatives w.r.t temperature */
    {
      double dtem = dui[U_tem] ;  //increments temp with this
      double* dx = ComputeVariableDerivatives(el,t,dt,x,dtem,I_TEM) ;

      /* Coupling matrix: stress w.r.t temperature*/
      //  =C_ijkl*alpha_kl
      {
        double* c1 = c0 + 81 + 9;
        double* dsigdtem = dx + I_SIG ;

        int i ;
        //for(i = 0 ; i < 3 ; i++) B1(i,i) = - alpha_s * K_d ;   //alpha_s is volumetric!!!
        for(i = 0 ; i < 9 ; i++) c1[i] = dsigdtem[i] ;
      }

      /* Coupling matrix: fluid mass w.r.t temperature */
      {
        double* c1 = c0 + 81 + 9*3 + 1 ;
        c1[0] = dx[I_M_L] ;
      }

      /* Storage term: entropy w.r.t temperature */
      {
        double* c1 = c0 + 81 + 9*3 + 1 + 1 + 9 + 1;
        c1[0] = dx[I_S_TOT] ;
        //c1[0] = C_d / TEMP0 ; //- s_l * rho_l0 * 3. * alpha_m ;       //corresponds to alpha_m of Coussy
      }
    }
  }

  return(dec) ;
#undef C1
#undef B1
#undef T2
#undef T4
}



int ComputeTransferCoefficients(FEM_t* fem,double dt,double* c)
/*
*  Conduction matrix (c) and shift (dec)
*
*   each submatrix matrix is a 3*3 matrix
*  | K_1(3x3) 0(3x3)    ... |
*  | 0(3x3)   K_2(3x3)  ... |
*
* here K_1 represents the hydraulic conduction matrix, K_2 the thermal conduction matrix
* K_1 and K_2 contain each the hydraulic and thermal conduction coefficient on the diaoganal, else zero
*
*/
{
  Element_t* el = FEM_GetElement(fem) ;
  double* vex0 = Element_GetExplicitTerm(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int    dec = 6*6 ;  //ADDED changed to 6x6 entries of the conduction matrix
  int    p ;
  double zero = 0. ;


  for(p = 0 ; p < np ; p++) {
    int i ;
    double* c1 = c + p*dec ;

    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = zero ;

    {
      double* vex  = vex0 + p*NVE ;

       /* Permeability tensor */   //ADDED changed to c2
      {
        double* c2 = c1 ;

        c2[0] = K_L[0] ;
        c2[4] = K_L[1] ;
        c2[8] = K_L[2] ;
      }

      /* Thermal conductivity tensor */   //ADDED
      {
        double* c2 = c1 + 3 * 9 ; //"jumps over" K_1(3x3), 0(3x3), 0(3x3) to write K_2(3x3)

        c2[0] = K_T[0] ;
        c2[4] = K_T[1] ;
        c2[8] = K_T[2] ;
      }
    }
  }

  return(dec) ;
}





double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double t,double dt,int p)
/** This locally defined function computes the internal variables at
 *  the interpolation point p, from the nodal unknowns.
 *  Return a pointer on the locally defined array of the variables.
 *  u gives current timestep, u_n the previous one.
 *  */
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double*   x      = Variable ;
  double*  x_n = Variable_n ;

  /* Primary Variables */
  /* Load the primary variables in x */
  {
    int    i ;

    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,p,U_u + i) ;
    }

    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }

    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,p,U_u) ;

      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }

      FEM_FreeBufferFrom(fem,eps) ;
    }

    /* Pressure */
    x[I_P_L] = FEM_ComputeUnknown(fem,u,intfct,p,U_p_l) ;

    /* Pressure gradient */
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_p_l) ;

      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_P_L + i] = grd[i] ;
      }

      FEM_FreeBufferFrom(fem,grd) ;
    }

    /* Temperature */    //ADDED
    x[I_TEM] = FEM_ComputeUnknown(fem,u,intfct,p,U_tem) ;

    /* Temperature gradient */           //ADDED
    {
      double* grd = FEM_ComputeUnknownGradient(fem,u,intfct,p,U_tem) ;

      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_TEM + i] = grd[i] ;
      }

      FEM_FreeBufferFrom(fem,grd) ;
    }

  }


  /* Needed variables to compute secondary variables */
  {
    int i ;

    /* Stresses, strains at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,p,U_u) ;
      double* vim_n = f_n + p*NVI ;

      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS_n   + i] = eps_n[i] ;
        x[I_SIG_n   + i] = SIG_n[i] ;
        x[I_EPS_P_n + i] = EPS_P_n[i] ;
      }

      x[I_HARDV_n] = HARDV_n ;

      x[I_PHI_n] = PHI_n ; //ADDED

      x[I_ALPHA_n] = ALPHA_n ;


      FEM_FreeBufferFrom(fem,eps_n) ;
    }

    /* Pressure at previous time step */
    x[I_P_L_n] = FEM_ComputeUnknown(fem,u_n,intfct,p,U_p_l) ;

    /* Temperature at previous time step (u_n) */                    //ADDED
    x[I_TEM_n] = FEM_ComputeUnknown(fem,u_n,intfct,p,U_tem) ;

    /* Transfer coefficient */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + p*NVE ;

      for(i = 0 ; i < 3 ; i++) {
        x[I_K_H + i]  = K_L[i] ; //TODO this is zero in the firs couple runs
        x[I_KTH + i]  = K_T[i] ;  //ADDED
      }
    }

    double* vim_n = f_n + p*NVI ;

  }

  ComputeSecondaryVariables(el,t,dt,x) ;

  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double t,double dt,double* x)
{
  int dim = Element_GetDimensionOfSpace(el) ;
  double* vc   = Element_GetConstantTerm(el) ; //ADDED

  GetProperties(el);

  /* Strains */
  double* eps   =  x + I_EPS ;
  double* eps_n =  x + I_EPS_n ;
  /* Plastic strains */
  double* eps_p  = x + I_EPS_P ;
  double* eps_pn = x + I_EPS_P_n ;
  /* Pressure */
  double  pl   = x[I_P_L] ;
  double  pl_n = x[I_P_L_n] ;
  double  pl_0 = P_L0 ;           //ADDED

  double phi_n = x[I_PHI_n] ; //ADDED
  double rho_l = x[I_RHO_L] ; //ADDED
  double alpha_n = x[I_ALPHA_n] ;

  /* Temperature */  //ADDED
  double tem    = x[I_TEM] ;
  double tem_n  = x[I_TEM_n] ;
  double tem_0  = TEMP0 ;
  double theta_n = tem_n - tem_0 ;
  double theta = tem - tem_0 ;

  double* sig   = x + I_SIG ;
  double* sig_n = x + I_SIG_n ;

  /* Update Hydration Degree */

  double alpha = alpha_n + dt/3600. * AlphaDot(alpha_n,x[I_P_L],273.15+x[I_TEM]) ;
  x[I_ALPHA] = alpha ;

  /* Update Saturation Degree */

  double s_f = SaturationDegree(x[I_P_L],alpha) ;

  x[I_S_F] = s_f ;

  /* Update permeability with Saturation Degree */

  int i ;

  // for(i = 0 ; i < 3 ; i++) {
  //   x[I_K_H + i]  = Permeability(s_f) * x[I_K_H + i]  ;
  // }
  /* Update poroelastic parameters with hydration */

  double biot = BiotCoefficientFromHydration(alpha)*Chi(s_f) ;
  double K_d = ( GetProperty("young") * ModulusFromHydration(alpha) )/ ( 3 * ( 1 - 2*PoissonFromHydration(alpha)) ) ;
  double K_s = K_d / ( 1 - biot )  ;

  /* Update elasticity tensor with hydration */

  double B_Skem = biot/K_d / (biot/K_d + phi_n / 2.2E+9) ;
  double young = GetProperty("young") * ModulusFromHydration(alpha) ;
  double poisson = PoissonFromHydration(alpha) ;

  plasty = Element_FindMaterialData(el,Plasticity_t,"Plasticity") ;
  
  #if 0
  {
    double ACC_M = GetProperty("ACC_M") ;
    double ACC_N = GetProperty("ACC_N") ;
    double h_eps = GetProperty("volumetric_strain_hardening_parameter") ;
    Plasticity_SetParameter(plasty,"ACC M", alpha * ACC_M) ;
    Plasticity_SetParameter(plasty,"ACC N", alpha * ACC_N) ;
    Plasticity_SetParameter(plasty,"ACC vol strain hardening parameter beta", alpha * h_eps ) ;
  }
  #endif
  
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    Elasticity_SetParameters(elasty,young,poisson) ;
    cijkl   = Elasticity_GetStiffnessTensor(elasty) ;
    Elasticity_ComputeStiffnessTensor(elasty,cijkl) ;
  }



  double b_ij[9] = {0,0,0,0,0,0,0,0,0}  ;  //ADDED
  rotate(b_ij,0,0) = biot ;
  rotate(b_ij,1,1) = biot ;
  rotate(b_ij,2,2) = biot ;

  double  alpha_ij[9] = {0,0,0,0,0,0,0,0,0}  ;  //ADDED
  rotate(alpha_ij,0,0) = alpha_1;
  rotate(alpha_ij,1,1) = alpha_1 ;
  rotate(alpha_ij,2,2) = alpha_3 ;

  /* Backup stresses, plastic strains */
  {


    double pc_n = x[I_HARDV_n] ;
    double  dpl   = pl - pl_n ;
    double  dtem  = tem - tem_n ;
    double  dalpha = alpha - alpha_n ;

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
          sig[i] += C(i,j)*deps[j] - C(i,j)*alpha_ij[j]* dtem ;
        }

        sig[i] += - b_ij[i] * dpl ;
      }
      #undef C

      /* Elastic trial effective stresses (with beta coefficient) */
      sig[0] += Chi(s_f)*pl ;
      sig[4] += Chi(s_f)*pl ;
      sig[8] += Chi(s_f)*pl ;


      /* Plastic strains */
      for(i = 0 ; i < 9 ; i++) eps_p[i] = eps_pn[i] ;

      //calculate the vol plastic strain from the previous step //ADDED
      double eps_p_v = eps_p[0] + eps_p[4] + eps_p[8];


      //calculate trial preconsolidation pressure, using the eps_p_v_n from the previous step n and the temperature from the current step n+1
      //theta is the integrated temperature change at n+1 obtained through tem - tem_0
      // we need the parameters for thermo-plasticity:
      //double pc0 = Plasticity_GetACCInitialPreconsolidationPressure(plasty);
      //double beta_TEMP0 = Plasticity_GetThermalHardeningParameter(plasty);
      //double beta_eps = Plasticity_GetVolumetricStrainHardeningParameter(plasty);
      //double beta_alpha = Plasticity_GetHydrationHardeningParameter(plasty);


      // pc_n is the preconsolidation pressure from the previous state, accounting for eps_p_v_n and theta_n
      ////double pc_t = pc0*(beta_TEMP0*theta + 1)*exp(beta_eps*eps_p_v);   //added here to be able to check free thermal hardening
      ////pc_t = MIN(pc_n,pc_t) ; //! deactivates softening

      // incremental hardening form:
      // thermal hardening is added incremetally to pc_n

      double dpc_alpha =  alpha > 0.03 ? beta_alpha * dalpha : 0. ;// beta_alpha * MIN(dalpha,alpha-0.2) : 0. ;
      // double dpc_alpha = beta_alpha * dalpha ;

      double pc_t = pc_n + dpc_alpha ;//pc_n * (1 + beta_TEMP0 / (beta_TEMP0*theta_n + 1)*dtem ) + beta_alpha * dalpha  ; //!


      /* Projection */
      // returns plastic strains in eps_p, the criterion f and the hardening variable
      // the condition criterion f >0  is checked within the function ReturnMapping
      {
        double hardv[2] = {pc_t,x[I_ALPHA]} ;
        double* crit = ReturnMapping(sig,eps_p,hardv);  // changed here directly to ACC !

        x[I_CRIT]  = crit[0] ;
        x[I_HARDV] = hardv[0] ;
      }


      //then the stresses which were removed are added again


      /* Total stresses */
      sig[0] -= Chi(s_f)*pl ;
      sig[4] -= Chi(s_f)*pl ;
      sig[8] -= Chi(s_f)*pl ;
    }
  }


  /* Backup mass flow */
  {
    /* Porosity */
    double tre   = eps[0] + eps[4] + eps[8] ;
    double tre_n   = eps_n[0] + eps_n[4] + eps_n[8] ; //ADDED
    double dtre   = tre - tre_n;  //ADDED

    double tre_p = eps_p[0] + eps_p[4] + eps_p[8] ;
    double tre_pn = eps_pn[0] + eps_pn[4] + eps_pn[8] ;
    double dtre_p = tre_p-tre_pn;
    double phi_p = beta*tre_p ; //!!

    double phi_li = (b_ij[0]*(eps[0]-eps_p[0]) + b_ij[4]*(eps[4]-eps_p[4]) + b_ij[8]*(eps[8]-eps_p[8])) + (biot - phi0)/K_s*(pl - pl_0) - (2*biot*alpha_1 + biot*alpha_3 - phi0*alpha_s)*(tem-tem_0) ; //ADDED thermal
    double phi_alpha = beta_ksi_phi*alpha ;
    double phi = phi0 + phi_p + phi_li + phi_alpha ;

    double rho_l = dens_w(tem,pl);

    double d_S ;
    int i ;
    #define C(i,j)  (cijkl[(i)*9+(j)])
      for(i = 0 ; i < 9 ; i++) {
        int  j ;
        for(j = 0 ; j < 9 ; j++) {
          d_S += C(i,j)*alpha_ij[j]*(eps[i]-eps_p[i]);//ADDED thermal stresses
        }
      }
    #undef C

    d_S = 0.; //! temporaily
    //double s_tot = s_tot0 - s_l*phi0*rho_l0 + s_l*rho_l*phi + d_S - (2*b*alpha_1 + b_3*alpha_3 - phi0*alpha_s)*(tem-tem_0) + (C_d/tem_0)*(tem-tem_0)  ; //ADDED
    double s_tot = s_tot0 + (C_d/tem_0)*(tem-tem_0) - rho_s * 1.e6 * alpha / tem  ; //ADDED

    x[I_S_TOT] = s_tot ;       //ADDED


    /* Fluid mass flow */
    {
      /* Transfer coefficient */
      double* k_l = x + I_K_H ;

      /* Pressure gradient */
      double* gpl = x + I_GRD_P_L ;

      /* Mass flow */
      double* w_l = x + I_W_L ;

      /* temperature */              //ADDED
      double* k_t = x + I_KTH ;
      /* Temperature gradient */   //ADDED
      double* gtem = x + I_GRD_TEM ;
      /* Heat flow */ //ADDED
      double* w_the = x + I_W_THE ;

      for(i = 0 ; i < 3 ; i++) {
        w_l[i] = - k_l[i]*gpl[i] ;      //TODO kl and kt is zero
        w_the[i] = - k_t[i]*gtem[i] ;      //ADDED
      }
      w_l[dim - 1] += k_l[dim - 1]*rho_l0*gravite ; //TODO check this for future incremental form
    }

    /* Liquid mass content, body force */
    {
      double m_l = rho_l * ( phi*s_f + GetProperty("beta_ksi_w") * alpha ) ;
      double* f_mass = x + I_Fmass ;

      x[I_M_L] = m_l ;
      x[I_RHO_L] = rho_l ;
      x[I_PHI] = phi ;

      for(i = 0 ; i < 3 ; i++) f_mass[i] = 0 ;
      f_mass[dim - 1] = (rho_s + m_l)*gravite ;
    }
  }

}

double* ComputeVariableDerivatives(Element_t* el,double t,double dt,double* x,double dxi,int i)
{

  double* dx = dVariable ;
  int j ;

  /* Primary Variables */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] = x[j] ;
  }

  /* We increment the variable as (x + dx) */
  dx[i] += dxi ;

  ComputeSecondaryVariables(el,t,dt,dx) ;

  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}

