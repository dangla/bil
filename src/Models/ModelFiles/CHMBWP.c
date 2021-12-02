#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#include "FEM.h"

#define TITLE   "Chemo-hydro-mechanics of bituminized waste products"
#define AUTHORS "G.MELOT"

#include "PredefinedMethods.h"

/* Nb of equations */
#define NEQ     (2 + dim)

#define NVI     (53)       /* Nb of implicit terms */
#define NVE     (8)        /* Nb of explicit terms */
#define NV0     (0)        /* Nb of constant terms */


/* Indices of equation */
#define E_Mech    (0)
#define E_Mass    (dim)
#define E_Salt    (1 + dim)


/* Indices of unknown (generic indices) */
#define U_Mech    E_Mech
#define U_Mass    E_Mass
#define U_Salt    E_Salt


/* Method chosen at compiling time.
 * Each equation is associated to a specific unknown.
 * Each unknown can deal with specific modelings.
 * Uncomment/comment to let only one unknown per equation */
 
/* Mechanics */
#define U_dis      U_Mech

/* Pressure */
#define U_p_l      U_Mass

/* Mass fraction of salt */
#define U_w_s      U_Salt



/* We define some names for implicit terms (vim must be used as pointer below) */
#define W_total     (vim)             // Liquid flux (water + dissolved salt)
#define M_total     (vim   + 3)[0]    // Total mass
#define M_totaln    (vim_n + 3)[0]

#define W_salt      (vim   + 4)       // Salt flux (dissolved salt + crystal salt)
#define M_salt      (vim   + 7)[0]    // Mass of dissolved salt + crystal salt
#define M_saltn     (vim_n + 7)[0]

#define SIG         (vim   + 8)
#define SIG_n       (vim_n + 8)

#define Phi_c       (vim   + 17)[0]    // Volumetric fraction of crystal salt
#define Phi_c_n     (vim_n + 17)[0]

#define Phi_l       (vim   + 18)[0]    // Porosity
#define Phi_ln      (vim_n + 18)[0]

#define Rho_l       (vim   + 19)[0]   // Density of liquide phases ( not equal to Rho_w)
#define Rho_ln      (vim_n + 19)[0]

#define Rho_w       (vim   + 20)[0]     // Density of water

#define Dm_s        (vim   + 21)[0]     // Quantity of crystal salt dissolved

#define SIG_M       (vim   + 22)[0]
#define SIG_M_n     (vim_n + 22)[0]

#define SIG_D       (vim   + 23)
#define SIG_D_n     (vim_n + 23)

#define EPSI_e      (vim   + 32)

#define EPSI_vp     (vim   + 41)
#define EPSI_vpn    (vim_n + 41)

#define M_W_water   (vim   + 50)[0]   // Mass of water transported through each element in each time step
#define M_W_watern  (vim_n + 50)[0]

#define M_W_salt    (vim   + 51)[0]   // Mass of salt transported through each element in each time step
#define M_W_saltn   (vim_n + 51)[0]

#define Poro_E      (vim   + 52)[0]  // Eulerian porosity
#define Poro_En     (vim_n + 52)[0]  // Eulerian porosity





/* We define some names for explicit terms (vex must be used as pointer below) */
#define K_l        (vex)[0]
#define D          (vex + 1)[0]
#define Tau        (vex + 2)[0]

#define M_Wc_water (vex + 3)[0]   // Mass of water transported through each element from the begining
#define M_Wc_salt  (vex + 4)[0]   // Mass of salt transported through each element from the begining

#define K_Darcy        (vex + 5)[0]
#define K_Osmosis      (vex + 6)[0]
#define D_Fick         (vex + 7)[0]



/* We define some names for constant terms (v0 must be used as pointer below) */
// None



/* Physical constant */
#define RT        (2436.)   /* product of R=8.3143 and T=293 (Pa.m3/mol) */ 


/* Fluid properties
 * ---------------- */
#define WaterDensity(p) \
        (rho_w_i * exp(compressibility_w*((p) - p_i)))

#define V_m_w  (1.e-3) /* Specific volume of water */
        
#define LiquidDensity(p,w_s) \
        (WaterDensity(p) * ((rho_l_dep == 0 ) ? 1 : (exp(-rho_w_i*(V_m_w-v_m_s)*(w_s)))))


/* Material properties 
 * ------------------- */
/* Mechanical properties (Mori-Tanaka) */
#define BulkModulusMoriTanaka(Ks,Gs,n) \
        ((1 - (n))*(Ks)*4*(Gs)/(3*(n)*(Ks) + 4*(Gs)))
        
#define ShearModulusMoriTanaka(Ks,Gs,n) \
        ((1 - (n))*(Gs)*(9*(Ks) + 8*(Gs))/(9*(Ks)*(1 + 2*(n)/3) + 8*(Gs)*(1 + 3*(n)/2)))
        
#define BiotCoefficientMoriTanaka(Ks,Gs,n) \
        (1 - (1 - (n))*4*(Gs)/(3*(n)*(Ks) + 4*(Gs)))

/* Permeability */
#define IntrinsicPermeability(n) \
        ((K_dep == 0) ? k_i : ((K_dep == 1) ? IntrinsicPermeability_UPC(n) : ((K_dep == 2) ? IntrinsicPermeability_Maxwell(n) : 0)))
/* Implementation */
#define IntrinsicPermeability_UPC(n) \
        (k_i * (pow(((n)/phi_l_i),3)*(pow(((1-phi_l_i)/(1-(n))),2))))
#define IntrinsicPermeability_Maxwell(n) \
        (k_i * (1 + 3*(n)/(1 - (n))))
        
/* Diffusion coefficient */
#define DiffusionCoefficient(n,tau) \
        ((D_dep == 0) ? d_i : ((D_dep == 1) ? DiffusionCoefficient_1(n,tau) : ((D_dep == 2) ? DiffusionCoefficient_2(n) : ((D_dep == 3) ? DiffusionCoefficient_3(tau) : ((D_dep == 4) ? DiffusionCoefficient_4(n) : 0)))))
/* Implementation */
#define DiffusionCoefficient_1(n,tau) \
        (DiffusionCoefficient_2(n) * (1 - (tau)))
#define DiffusionCoefficient_2(n) \
        (d_i * (1 + 3*(n)/(1 - (n))))
#define DiffusionCoefficient_3(tau) \
        (d_i * (1 - (tau)))
#define DiffusionCoefficient_4(n) \
        (d_i * (n) / phi_l_i)
        
/* Osmotic efficiency coefficient */
#define OsmoticEfficiencyCoefficient(n) \
        ((Tau_dep == 0) ? tau_i : ((Tau_dep == 2) ? OsmoticEfficiencyCoefficient_UPC(n) : ((Tau_dep == 3) ?OsmoticEfficiencyCoefficient_Maxwell01(n) : ((Tau_dep == 4) ? OsmoticEfficiencyCoefficient_Maxwell05(n) : ((Tau_dep == 5) ? OsmoticEfficiencyCoefficient_5(n) : ((Tau_dep == 6) ? OsmoticEfficiencyCoefficient_Maxwell08(n) : ((Tau_dep == 7) ? OsmoticEfficiencyCoefficient_Sigmoid(n) : 0)))))))
/* Implementation */
#define OsmoticEfficiencyCoefficient_UPC(n) \
        (tau_i * (phi_l_i/(n)))
#define OsmoticEfficiencyCoefficient_Maxwell01(n) \
        (tau_i * (1 - 3*0.1*(n)/(2 + (n))))
#define OsmoticEfficiencyCoefficient_Maxwell05(n) \
        (tau_i * (1 - 3*0.5*(n)/(2 + (n))))
#define OsmoticEfficiencyCoefficient_5(n) \
        (((n) < 0.55) ? tau_i : tau_i*(phi_l_i/((n) - 0.55)))
#define OsmoticEfficiencyCoefficient_Maxwell08(n) \
        (tau_i * (1 - 3*0.8*(n)/(2 + (n))))
#define OsmoticEfficiencyCoefficient_Sigmoid(n) \
        (tau_i * (1 - 1/(1 + exp(-A_tau*((n) - B_tau)))))
#define A_tau  30
#define B_tau  0.3
        
        
/* Crystal properties
 * ------------------ */
#define CrystalVolumeFraction(phi_c_n,w_s,dt) \
        ((phi_c_n) * (1 + (dt) * sigma_c * beta * ((w_s)/w_s_sat - 1) / rho_c))



/* Functions */
static int     pm(const char*) ;

static int     ComputeTangentCoefficients(Element_t*,double,double*) ;
static int     ComputeTransferCoefficients(Element_t*,double,double*) ;

static double* ComputeVariables(Element_t*,double**,double**,double*,double,int) ;
static void    ComputeSecondaryVariables(Element_t*,double,double*,double*) ;
static double* ComputeVariableDerivatives(Element_t*,double,double*,double,int) ;



/************************/
/** Material parameters */
/************************/
/**********************************************************************/
//Transport parameters
static double  k_i;        //Initial intrinsic permeability
static double  d_i;        //Initial diffusion coefficient
static double  tau_i;      //Initial efficiency coefficient
static double  phi_l_i;      //Initial porosity
//Water property
static double  rho_w_i;      //Liquid density
static double  viscosity_w;    //Dynamic viscosity of water    
static double  compressibility_w;  //Compressibility of water
//Salt property
static double  rho_c;      //Crystal salt density
static double  beta;      //Dissolution rate
static double  sigma_c;      //Crystal specific surface  
static double  Molarm_s;    //Solute molar mass
static double  w_s_sat;      //Solute mass fraction at saturation
static double  v_m_s;      //Specific volume of disolve salt
//Elasticity law parameters
static double  young;      //Young modulus
static double  poisson;      //Poisson coefficient
//Creep law parameters 
static double  eta_v;      //Bitumen volumetric viscosity (Maxwell law)
static double  eta_d;      //Bitumen deviatoric viscosity (Maxwell law)
/**********************************************************************/
/** Initial condition */
/**********************************************************************/
static double sig0_11;      //Initial longitudinal stress
static double sig0_22;      //Initial stress
static double sig0_33;      //Initial stress 
static double w_s_i;      //Initial solute mass fraction
static double phi_c_i;      //Initial crystal volumetric fraction
static double p_i;        //Initial pore water pressure
/**********************************************************************/
//static double biot;        //Biot coefficient
/**********************************************************************/
/** Coefficients dependency  */
/**********************************************************************/
static double D_dep;      // 0 -> D = cte / 1 -> D = f(phi_l) 
static double K_dep;      // 0 -> K = cte / 1 -> K = f(phi_l) 
static double Tau_dep;      // 0 -> Tau = cte / 1 -> Tau = f(phi_l)
static double rho_l_dep;      // 0 -> rho_l = cte = p / 1 -> rho_l = f(w_s)
/**********************************************************************/
static double SampleSurface;  //Sample exchange surface [m²]
static double ElementSize;     //Element size [m]  




enum {
I_U = 0,
I_U2 = I_U + 2,

I_P,

I_W_S,
I_W_S2 = I_W_S + 2,

I_EPS,
I_EPS8 = I_EPS + 8,

I_SIG,
I_SIG8 = I_SIG + 8,

I_SIG_M,

I_SIG_D,
I_SIG_D8 = I_SIG_D + 8,

I_PHI_C,
I_PHI_L,

I_RHO_W,
I_RHO_L,

I_M_TOT,
I_M_SALT,

I_K_L,
I_D,
I_TAU,

I_K_Darcy,
I_D_Fick,
I_K_Osmosis,

I_GRD_P,
I_GRD_P2 = I_GRD_P + 2,

I_GRD_W_S,
I_GRD_W_S2 = I_GRD_W_S + 2,

I_W_F,
I_W_F2 = I_W_F + 2,

I_W_D,
I_W_D2 = I_W_D + 2,

I_W_O,
I_W_O2 = I_W_O + 2,

I_W_L,
I_W_L2 = I_W_L + 2,

I_W_SALT,
I_W_SALT2 = I_W_SALT + 2,

I_Poro_E,

I_Last
} ;


#define NbOfVariables     (I_Last)
static double Variable[NbOfVariables] ;
static double Variable_n[NbOfVariables] ;
static double dVariable[NbOfVariables] ;








int pm(const char* s)
{
       if(!strcmp(s,"k_i"))               return(0) ;
  else if(!strcmp(s,"d_i"))               return(1) ;
  else if(!strcmp(s,"tau_i"))             return(2) ;
  else if(!strcmp(s,"phi_l_i"))           return(3) ;
  else if(!strcmp(s,"rho_w_i"))           return(4) ;
  else if(!strcmp(s,"viscosity_w"))       return(5) ;
  else if(!strcmp(s,"compressibility_w")) return(6) ;
  else if(!strcmp(s,"rho_c"))             return(7) ;
  else if(!strcmp(s,"beta"))              return(8) ;
  else if(!strcmp(s,"sigma_c"))           return(9) ;
  else if(!strcmp(s,"Molarm_s"))          return(10);
  else if(!strcmp(s,"w_s_sat"))           return(11);
  else if(!strcmp(s,"v_m_s"))             return(12);
  else if(!strcmp(s,"young"))             return(13);
  else if(!strcmp(s,"poisson"))           return(14);
  else if(!strcmp(s,"eta_v"))             return(15);
  else if(!strcmp(s,"eta_d"))             return(16);
  else if(!strcmp(s,"sig0_11"))           return(17);
  else if(!strcmp(s,"sig0_22"))           return(18);
  else if(!strcmp(s,"sig0_33"))           return(19);
  else if(!strcmp(s,"w_s_i"))             return(20);
  else if(!strcmp(s,"phi_c_i"))           return(21);
  else if(!strcmp(s,"p_i"))               return(22);
  else if(!strcmp(s,"b"))                 return(23);
  else if(!strcmp(s,"D_dep"))             return(24);
  else if(!strcmp(s,"K_dep"))             return(25);
  else if(!strcmp(s,"Tau_dep"))           return(26);
  else if(!strcmp(s,"rho_l_dep"))         return(27); 
  else if(!strcmp(s,"SampleSurface"))     return(28);
  else if(!strcmp(s,"ElementSize"))       return(29); 
  else return(-1);
}




/* They are retrieved automatically by calling the following function */
static void    GetProperties(Element_t*) ;
void GetProperties(Element_t* el)// copy the value of material properties
{

/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
  k_i       = GetProperty("k_i");
  d_i        = GetProperty("d_i");
  tau_i      = GetProperty("tau_i");
  phi_l_i    = GetProperty("phi_l_i");
  rho_w_i    = GetProperty("rho_w_i");
  viscosity_w  = GetProperty("viscosity_w");
  compressibility_w  = GetProperty("compressibility_w");
  rho_c      = GetProperty("rho_c");
  beta      = GetProperty("beta");
  sigma_c    = GetProperty("sigma_c");
  Molarm_s    = GetProperty("Molarm_s");
  w_s_sat    = GetProperty("w_s_sat");
  v_m_s      = GetProperty("v_m_s");
  young      = GetProperty("young");
  poisson    = GetProperty("poisson");
  eta_v      = GetProperty("eta_v");
  eta_d      = GetProperty("eta_d");
  sig0_11    = GetProperty("sig0_11");
  sig0_22    = GetProperty("sig0_22");
  sig0_33    = GetProperty("sig0_33");
  w_s_i      = GetProperty("w_s_i");
  phi_c_i    = GetProperty("phi_c_i");
  p_i        = GetProperty("p_i");
  //biot       = GetProperty("b");
  D_dep         = GetProperty("D_dep");
  K_dep         = GetProperty("K_dep");
  Tau_dep       = GetProperty("Tau_dep");
  rho_l_dep     = GetProperty("rho_l_dep");
  SampleSurface  = GetProperty("SampleSurface");
  ElementSize  = GetProperty("ElementSize");
#undef GetProperty
}



int SetModelProp(Model_t *model)
{
  int dim = Model_GetDimension(model) ;
  char name_eqn[3][7] = {"meca_1","meca_2","meca_3"} ;
  char name_unk[3][4] = {"u_1","u_2","u_3"} ;
  int i ;
  
  /** Number of equations to be solved */
  Model_GetNbOfEquations(model) = NEQ ;
  
  /** Names of these equations */  
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfEquation(model,E_Mech + i,name_eqn[i]) ;
  }
  Model_CopyNameOfEquation(model,E_Mass,"Mass") ;
  Model_CopyNameOfEquation(model,E_Salt,"Salt") ;

  /** Names of the main (nodal) unknowns */
  for(i = 0 ; i < dim ; i++) {
    Model_CopyNameOfUnknown(model,U_dis + i,name_unk[i]) ;
  }
  
#if defined (U_p_l)
  Model_CopyNameOfUnknown(model,U_Mass,"p") ;
#else
  #error "Undefined unknown"
#endif

#if defined (U_w_s)
  Model_CopyNameOfUnknown(model,U_Salt,"w_s") ;
#else
  #error "Undefined unknown"
#endif

  
  return(0) ;
}



int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Reading of material properties in file ficd */
{
  int  NbOfProp = 30;
  Material_ScanProperties(mat,datafile,pm) ;

  return(NbOfProp) ;
}



int PrintModelChar(Model_t* model,FILE* ficd)
/** Print the model characteristics */
{
  printf(TITLE) ;
  
  if(!ficd) return(0) ;
  
  printf("\n") ;
  printf("The set of equations is:\n") ;
  printf("\t- Mechanical Equilibrium    (mec)\n") ;
  printf("\t- Total mass conservation    (mass)\n") ;
  printf("\t- Salt mass conservation    (salt)\n") ;


  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Displacements          (u_1,u_2,u_3) \n") ;
  printf("\t- Pore water pressure      (p) \n") ;
  printf("\t- Solute mass fraction      (w_s) \n") ;
  
  printf("\n") ;
  printf("Example of input data\n") ;

  fprintf(ficd,"k_i = 4e-25       # (m2)    Initial intrinsic permeability\n") ;
  fprintf(ficd,"d_i = 3e-14       # (m2/s)  Initial diffusion coefficient\n") ;
  fprintf(ficd,"tau_i = 0.9996        # (-)     Initial efficiency coefficient\n") ;
  fprintf(ficd,"phi_l_i = 0.01       # (-)     Initial porosity\n") ;
  fprintf(ficd,"rho_w_i = 1000      # (kg/m3) Water density \n") ;
  fprintf(ficd,"viscosity_w = 1e-3     # (Pa.s)  Dynamic viscosity of water\n") ;
  fprintf(ficd,"compressibility_w = 5e-10   # (1/Pa)  Compressibility of water\n") ;
  fprintf(ficd,"rho_c = 2260        # (kg/m3) Crystal salt density\n") ;
  fprintf(ficd,"beta = 1e-5        # (kg/(s.m3)) Dissolution rate\n") ;
  fprintf(ficd,"sigma_c = 148       # (m2/m3) Crystal specific surface\n") ;
  fprintf(ficd,"Molarm_s = 85e-3       # (kg/mol) Solute molar mass\n") ;
  fprintf(ficd,"w_s_sat = 0.47      # (-)     Solute mass fraction at saturation\n") ;
  fprintf(ficd,"v_m_s = 1.689e-3    # (m3/kgl) Specific volume of dissolved salt\n") ;
  fprintf(ficd,"young = 1e8         # (Pa)    Young modulus\n") ;
  fprintf(ficd,"poisson = 0.33      # (-)     Poisson's coefficient\n") ;
  fprintf(ficd,"eta_v = 3.33e14      # (Pa.s)  Bitumen volumetric viscosity (Maxwell law)\n") ;
  fprintf(ficd,"eta_d = 1e15      # (Pa.s)  Bitumen deviatoric viscosity (Maxwell law)\n") ;
  fprintf(ficd,"sig0_11 = 0         # (Pa)    Initial longitudinal stress\n") ;
  fprintf(ficd,"sig0_22 = 0          # (Pa)    Initial stress\n") ;
  fprintf(ficd,"sig0_33 = 0         # (Pa)    Initial stress\n") ;
  fprintf(ficd,"w_s_i = 0.47        # (-)     Initial solute mass fraction\n") ;
  fprintf(ficd,"phi_c_i = 0.16      # (-)     Initial crystal volumetric fraction\n") ;
  fprintf(ficd,"p_i = 1e5         # (Pa)    Initial pore water pressure\n") ;
  fprintf(ficd,"b = 1           # (-)     Coefficient de Biot\n") ;
  fprintf(ficd,"D_dep = 2           # (-)     Diffusion coef. dependency with phi_lt\n") ;
  fprintf(ficd,"K_dep = 2           # (-)     Permeability dependency with phi_l\n") ;
  fprintf(ficd,"Tau_dep = 0         # (-)     Osmosis coef. dependency with phi_l\n") ;
  fprintf(ficd,"rho_l_dep = 1         # (-)     Liquid density dependency with w_s\n") ;
  fprintf(ficd,"SampleSurface = 19.63e-4 # (m2)    Sample exchange surface\n") ;
  fprintf(ficd,"ElementSize = 2e-5     # (m)    Element size\n") ;
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



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/** Compute the residu (r) due to loads */
{
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t *fi = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  {
    double *r1 = FEM_ComputeSurfaceLoadResidu(fem,fi,cg,t,dt) ;
    /* hydraulic */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_Mass) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }   
    /* Salt */
    if(Element_FindEquationPositionIndex(el,Load_GetNameOfEquation(cg)) == E_Salt) {
      for(i = 0 ; i < ndof ; i++) r[i] = -r1[i] ;
  }
    // other //  
     else {
      for(i = 0 ; i < ndof ; i++) r[i] = r1[i] ;
    }
  }
  return(0) ;
}



int ComputeInitialState(Element_t* el)
{
/** Compute the initial state i.e. the constant terms,the explicit terms,the implicit terms. */
  double* vim0 = Element_GetImplicitTerm(el) ;
  double *vex0 = Element_GetExplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    m ;
 
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;

  /* Pre-initialization */
  for(m = 0 ; m < NbOfIntPoints ; m++) {
    double p = FEM_ComputeUnknown(fem,u,intfct,m,U_Mass) ;
    double w_s = FEM_ComputeUnknown(fem,u,intfct,m,U_Salt) ;
    double rho_l = LiquidDensity(p,w_s) ;
    
    /* Storage in vim */
    {
      double* vim = vim0 + m*NVI ;
      int i ;
      
      for(i = 0 ; i < 9 ; i++) SIG[i] = 0 ;
      
      SIG[0] = sig0_11 ;
      SIG[4] = sig0_22 ;
      SIG[8] = sig0_33 ;
      
      Phi_c = phi_c_i ;
      
      {
        double K_s = young/(3 - 6*poisson) ;
        double G_s = young/(2 + 2*poisson) ;
        double n   = phi_l_i ;
        double biot = BiotCoefficientMoriTanaka(K_s,G_s,n) ;
        double sigm = (SIG[0]+SIG[4]+SIG[8])/3 ;
      
        for(i = 0 ; i < 9 ; i++) SIG_D[i] = SIG[i] ;
        
        SIG_D[0] -= sigm ;
        SIG_D[4] -= sigm ;
        SIG_D[8] -= sigm ; 
      
        SIG_M = sigm - biot*(p - p_i) ;
      }
    }
  
    /* storage in vex */
    {
      double * vex = vex0 + m*NVE ;
      /** transfert coefficient */
      double n    = phi_l_i ;
      double k_l  = IntrinsicPermeability(n) ;
      double tau  = OsmoticEfficiencyCoefficient(n) ;
      double d    = DiffusionCoefficient(n,tau) ;
    
      K_l = k_l ;
      D   = d ;
      Tau = tau ;
      
      {
        double k_darcy = rho_l * k_l / viscosity_w ;
        double d_fick = rho_l * d ;
        double k_osmosis = k_darcy*(RT*rho_l/Molarm_s)*tau ;
      
        K_Darcy = k_darcy ;
        D_Fick  = d_fick ;
        K_Osmosis = k_osmosis ;
      }
    }
  }
  
  /* Loop on integration points */
  for(m = 0 ; m < NbOfIntPoints ; m++) {
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0,0,m) ;

    /* Storage in vim */
    {
      double* vim = vim0 + m*NVI ;
      int i ;
    
      M_total = x[I_M_TOT] ;
      M_salt  = x[I_M_SALT] ;
    
      for(i = 0 ; i < 3 ; i++) W_total[i] = x[I_W_L + i] ;
      for(i = 0 ; i < 3 ; i++) W_salt[i]  = x[I_W_SALT + i] ;
    
      Rho_l = x[I_RHO_L] ;    
      Rho_w = x[I_RHO_W] ;
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      for(i = 0 ; i < 9 ; i++) SIG_D[i] = x[I_SIG_D + i] ;
      
      SIG_M = x[I_SIG_M] ;
    
      Phi_l = x[I_PHI_L] ;
      Phi_c = x[I_PHI_C] ;
      Poro_E = x[I_Poro_E] ;
    
      /* For post-treatment of M_salt_leached and m_water_absorbed */
      M_W_water = 0 ;
      M_W_salt = 0 ;  
    }
  }
  
  return(0) ;
}



int  ComputeExplicitTerms(Element_t *el,double t)
/* Explicites terms  */
{
  double *vex0 = Element_GetExplicitTerm(el) ;
  double *vim0_n = Element_GetPreviousImplicitTerm(el) ;   
  double** u = Element_ComputePointerToPreviousNodalUnknowns(el) ;   
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    m ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;
  
  /* Loop on integration points */
  for(m = 0 ; m < NbOfIntPoints ; m++) {
    double * vex = vex0 + m*NVE ;
    double * vim_n = vim0_n + m*NVI ;
    /* Variables */
    double* x = ComputeVariables(el,u,u,vim0_n,0,m) ;
    
    /* densities */
    double rho_l = x[I_RHO_L] ;
    
    /* eulerian porosity */
    double n     = x[I_Poro_E] ;
    
    /* Transfer coefficients */
    double k_l = IntrinsicPermeability(n) ;
    double tau = OsmoticEfficiencyCoefficient(n) ;
    double d   = DiffusionCoefficient(n,tau) ;
    
    /* Storage in vex */
    {
      K_l = k_l ;
      Tau = tau ;
      D   = d ;
    }
      
    {
      double k_darcy = rho_l * k_l / viscosity_w ;
      double d_fick = rho_l * d ;
      double k_osmosis = k_darcy*(RT*rho_l/Molarm_s)*tau ;
      
      K_Darcy = k_darcy ;
      D_Fick  = d_fick ;
      K_Osmosis = k_osmosis ;
    }
    
    
    /* For post-treatment of M_salt_leached and m_water_absorbed */
    M_Wc_salt  += M_W_saltn ;
    M_Wc_water += M_W_watern;
  }

  return(0) ;
}



int  ComputeImplicitTerms(Element_t* el,double t,double dt)
/** Compute the implicit terms */
{
  double* vim0   = Element_GetImplicitTerm(el) ;
  double* vim0_n = Element_GetPreviousImplicitTerm(el) ;
  double** u   = Element_ComputePointerToNodalUnknowns(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int NbOfIntPoints = IntFct_GetNbOfPoints(intfct) ;
  int    m ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;

  /* Loop on integration points */
  for(m = 0 ; m < NbOfIntPoints ; m++) {
    double * vim = vim0 + m*NVI ;
    /* Variables */
    double* x = ComputeVariables(el,u,u_n,vim0_n,dt,m) ;
    
    /* Storage in vim */
    {
      int i ;
    
      M_total = x[I_M_TOT] ;
      M_salt  = x[I_M_SALT] ;
    
      for(i = 0 ; i < 3 ; i++) W_total[i] = x[I_W_L + i] ;
      for(i = 0 ; i < 3 ; i++) W_salt[i] = x[I_W_SALT + i] ;
    
      Rho_l = x[I_RHO_L] ;    
      Rho_w = x[I_RHO_W] ;    
      // manque Dm_s par rapport à CHM21 ... utile ??
    
      for(i = 0 ; i < 9 ; i++) SIG[i] = x[I_SIG + i] ;
      for(i = 0 ; i < 9 ; i++) SIG_D[i] = x[I_SIG_D + i] ;
      SIG_M = x[I_SIG_M] ;
    
      Phi_l = x[I_PHI_L] ;
      Phi_c = x[I_PHI_C] ;
      Poro_E = x[I_Poro_E] ;
    
      /* For post-treatment of M_salt_leached and m_water_absorbed */
      M_W_water = (W_total[0] - W_salt[0])*dt ;
      M_W_salt = W_salt[0]*dt ;
    }
  }
  
  return(0) ;
}



int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Compute the matrix (k) */
{
#define K(i,j)    (k[(i)*ndof + (j)])
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i,j ;
  double zero = 0. ;

  
  /* Initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;
  
  /* We skip if the element is a submanifold */
  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* Input data */
  GetProperties(el) ;
  
  
  /***** Poromechanic Matrix *****/
  {
    int n = 81 + 2*9 + 2*(9 + 2) ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTangentCoefficients(el,dt,c) ;
    double *kp = FEM_ComputePoroelasticMatrix(fem,intfct,c,dec,2,E_Mech) ;
    
    for(i = 0 ; i < ndof*ndof ; i++) {        
      k[i] = kp[i] ;
    }
  }
  
  /*
  ** Conduction matrix
  */
  {
    int n = 4*9 ;
    double c[IntFct_MaxNbOfIntPoints*n] ;
    int dec = ComputeTransferCoefficients(el,dt,c) ;
    
    /***** Conduction Matrix liquid (dW_total/dp) *****/
    {
      double *kc = FEM_ComputeConductionMatrix(fem,intfct,c,dec) ;
  
      for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
        K(E_Mass + i*NEQ,U_Mass + j*NEQ) -= dt*kc[i*nn + j] ;
      }
    }
    
    /***** Conduction Matrix liquid (dW_total/dw_s) *****/
    {
      double *kc = FEM_ComputeConductionMatrix(fem,intfct,c+9,dec) ;
  
      for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
        K(E_Mass + i*NEQ,U_Salt + j*NEQ) -= dt*kc[i*nn + j] ;
      }
    }
    
    /***** Conduction Matrix salt (dW_salt/dp) *****/
    {
      double *kc = FEM_ComputeConductionMatrix(fem,intfct,c+2*9,dec) ;
  
      for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
        K(E_Salt + i*NEQ,U_Mass + j*NEQ) -=  dt*kc[i*nn + j] ;
      }
    }
    
    /***** Conduction Matrix salt (dW_salt/dw_s) *****/
    {
      double *kc = FEM_ComputeConductionMatrix(fem,intfct,c+3*9,dec) ;
  
      for(i = 0 ; i < nn ; i++) for(j = 0 ; j < nn ; j++) {
        K(E_Salt + i*NEQ,U_Salt + j*NEQ) -= dt*kc[i*nn + j] ;
      }
    }
  }
  

  /** Elements P2P1 */
  {
    char* method = Material_GetMethod(Element_GetMaterial(el)) ;
    
    if(strstr(method,"P2P1")) {
      FEM_TransformMatrixFromDegree2IntoDegree1(fem,U_Mass,E_Mass,k) ;
      FEM_TransformMatrixFromDegree2IntoDegree1(fem,U_Salt,E_Salt,k) ;
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
  double *vim_n1 = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int ndof = nn*NEQ ;
  FEM_t *fem = FEM_GetInstance(el) ;
  int    i ;
  double zero = 0. ;

  /* Initialisation */
  for(i = 0 ; i < ndof ; i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;


  /***** 1. Mechanics *****/
  /** 1.1 Stresses */
  {
    double *vim = vim_1 ;
    double *rw = FEM_ComputeStrainWorkResidu(fem,intfct,SIG,NVI) ;
    
    
    for(i = 0 ; i < nn ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) R(i,E_Mech + j) -= rw[i*dim + j] ;
    }
  }

  /***** 2. Total mass conservation *****/
  /** 2.1 Accumulation Terms */
  {
    double *vim = vim_1 ;
    double *vim_n = vim_n1 ;
    double g1[IntFct_MaxNbOfIntPoints] ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_total - M_totaln ;
    
    {
      double *ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
      
      for(i = 0 ; i < nn ; i++) R(i,E_Mass) -= ra[i] ;
    }
  }
  /** 2.2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_total,NVI) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_Mass) -= -dt*rf[i] ; 
  }
  
  /***** 3. Salt mass conservation *****/
  /** 3.1 Accumulation Terms */
  {
    double *vim = vim_1 ;
    double *vim_n = vim_n1;
    double g1[IntFct_MaxNbOfIntPoints] ;
    
    for(i = 0 ; i < np ; i++ , vim += NVI , vim_n += NVI) g1[i] = M_salt - M_saltn ;
    
    {
      double *ra = FEM_ComputeBodyForceResidu(fem,intfct,g1,1) ;
      
      for(i = 0 ; i < nn ; i++) R(i,E_Salt) -= ra[i] ;
    }
  }
  /** 3.2 Transport Terms */
  {
    double *vim = vim_1 ;
    double *rf = FEM_ComputeFluxResidu(fem,intfct,W_salt,NVI) ;
    
    for(i = 0 ; i < nn ; i++) R(i,E_Salt) -= -dt*rf[i] ; 
  }
  
  /** 4. Elements P2P1 */
  {
    char* method = Material_GetMethod(Element_GetMaterial(el)) ;
    
    if(strstr(method,"P2P1")) {
      FEM_TransformResiduFromDegree2IntoDegree1(fem,U_Mass,E_Mass,r) ;
      FEM_TransformResiduFromDegree2IntoDegree1(fem,U_Salt,E_Salt,r) ;
    }
  }
  
  return(0) ;
#undef R  
}



int  ComputeOutputs(Element_t* el,double t,double* s,Result_t* r)
/** Compute the outputs (r) */
{
  int NbOfOutputs = 17 ;
  double** u  = Element_ComputePointerToNodalUnknowns(el) ;
  double *vim0 = Element_GetImplicitTerm(el) ;
  double *vex0 = Element_GetExplicitTerm(el) ;
  IntFct_t  *intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  int dim = Geometry_GetDimension(Element_GetGeometry(el)) ;
  FEM_t *fem = FEM_GetInstance(el) ;
  double zero = 0. ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /* initialization */
  {
    int    i,j ;
    
    for(i = 0 ; i < NbOfOutputs ; i++) for(j = 0 ; j < 9 ; j++) {
      Result_GetValue(r + i)[j] = zero ;
    }
  }

  {
    /* Interpolation functions at s */
    double* a = Element_ComputeCoordinateInReferenceFrame(el,s) ;
    int m = IntFct_ComputeFunctionIndexAtPointOfReferenceFrame(intfct,a) ;
    
    /* unknown */
    double p   = FEM_ComputeUnknown(fem,u,intfct,m,U_Mass) ;
    double w_s = FEM_ComputeUnknown(fem,u,intfct,m,U_Salt) ;
    
    double dis[3] = {0,0,0} ;
    
    /* VI */
    double w_l[3] = {0,0,0} ;
    double w_salt[3] = {0,0,0} ;
    double phil = 0 ;
    double phil_Eulerian = 0 ;
    double phic = 0;
    double sig[9] = {0,0,0,0,0,0,0,0,0} ;
    double tre = 0 ; 
    
    /* For post-treatment of M_salt_leached and m_water_absorbed */   
    double mwcsalt = 0 ;
    double mwcwater = 0 ;
    double msaltleached = 0 ;
    double mwaterabsorbed = 0 ;
    ////////////////////////////////////////////////////////////////////
    double msalt_element = 0 ;
    double mwater_element = 0 ;
    //double indicateur =0 ;
    
    /* VE = k,d and tau */
    double k = 0 ;
    double d = 0 ;
    double tau = 0 ;
    int i ;
    
    
    for(i = 0 ; i < dim ; i++) {                    // Displacement calculation
      dis[i] = FEM_ComputeUnknown(fem,u,intfct,m,U_Mech + i) ;
    }
    
    
    /* Averaging */
    for(i = 0 ; i < np ; i++) {    // np = 2; permet de moyenner sur chaque élément 
      double * vim = vim0 + i*NVI ;
      double * vex = vex0 + i*NVE ;
      int j ;

      for(j = 0 ; j < dim ; j++) w_l[j] += W_total[j]/np ;
      for(j = 0 ; j < dim ; j++) w_salt[j] += W_salt[j]/np ;
      
      for(j = 0 ; j < 9 ; j++) sig[j] += SIG[j]/np ;
      
      phil += Phi_l/np ;
      phic += Phi_c/np ;
      mwcsalt += M_Wc_salt/np ;
      mwcwater += M_Wc_water/np ;
          ///////////////////////////////////////////////////////////////////////////////////////////////////
      //~ printf("\n valeur 1 %.2f | valeur deux  %.2f | valeur 3 %i ",Phi_l,phil,np);
     ///////////////////////////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////////////////////////
      //~ printf("\n valeur i %i | valeur np  %i | valeur 3 %i ",i,np,2);
     ///////////////////////////////////////////////////////////////////////////////////////////////////

      {
        double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,i,U_Mech) ;
        
        tre   += (eps[0] + eps[4] + eps[8])/np ;
      }
      
      k += K_l/np ;
      d += D/np ;
      tau += Tau/np ;
      ////////////////////////////////////////////////
      msalt_element += M_salt* SampleSurface * ElementSize /np ;  // msalt_element [kg] = M_Salt [kg/m3] x SampleSurface [m2] x ElementSize [m] / np [-]
      mwater_element += (M_total - M_salt)* SampleSurface * ElementSize /np ;
      
      phil_Eulerian += Poro_E/np ;
    }
    
    
    msaltleached = -mwcsalt * SampleSurface *1e6 ;
    mwaterabsorbed = mwcwater * SampleSurface *1e6;
      
    i = 0 ;
    Result_Store(r + i++,dis,"Displacement vector",3) ;
    Result_Store(r + i++,&p,"Liquid pressure",1) ;
    Result_Store(r + i++,&w_s,"Solute mass fraction",1) ;
    
    Result_Store(r + i++,w_l,"Liquid flux vector",3) ;
    Result_Store(r + i++,w_salt,"Solute flux vector",3) ;
    Result_Store(r + i++,sig,"Stress",9) ;
    Result_Store(r + i++,&phil,"Lagrangian porosity",1) ;
    Result_Store(r + i++,&phic,"Crystal volumetric fraction",1) ;
    Result_Store(r + i++,&msaltleached,"msaltleached",1) ;
    Result_Store(r + i++,&mwaterabsorbed,"mwaterabsorbed",1) ;
    
    Result_Store(r + i++,&k,"k",1) ;
    Result_Store(r + i++,&d,"d",1) ;
    Result_Store(r + i++,&tau,"tau",1) ;
    
    Result_Store(r + i++,&tre,"tre",1) ;
     ///////////////////////////////////////////////////////////////////////////////////////////////////
     //~ printf("\n juste avant la sauvegarde de la valeur  %.2f ",phil_Eulerian);
     ///////////////////////////////////////////////////////////////////////////////////////////////////

    Result_Store(r + i++,&phil_Eulerian,"Eulerian porosity",1) ;   
    /////////////////////////
    Result_Store(r + i++,&msalt_element,"msalt_element",1) ;
    Result_Store(r + i++,&mwater_element,"mwater_element",1) ;
           
    if(i != NbOfOutputs) arret("ComputeOutputs") ;
  }

  return (NbOfOutputs) ; 
}



int ComputeTangentCoefficients(Element_t* el,double dt,double* c)
/*
**  Tangent matrix (c), return the shift (dec).
*/
{
#define T4(a,i,j,k,l)  ((a)[(((i)*3+(j))*3+(k))*3+(l)])
#define T2(a,i,j)      ((a)[(i)*3+(j)])
#define C1(i,j,k,l)    T4(c1,i,j,k,l)
#define B1(i,j)        T2(c1,i,j)
  double*  vim0_n = Element_GetPreviousImplicitTerm(el) ;
  double** u     = Element_ComputePointerToCurrentNodalUnknowns(el) ;
  double** u_n   = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  double dui[NEQ] ;
  int    dec = 81 + 2*9 + 2*(9 + 2) ;
  
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  {
    ObVal_t* obval = Element_GetObjectiveValue(el) ;
    int i ;
    
    for(i = 0 ; i < NEQ ; i++) {
      dui[i] =  1.e-2*ObVal_GetValue(obval + i) ;
    }
  }


  {
    IntFct_t*  intfct = Element_GetIntFct(el) ;
    int np = IntFct_GetNbOfPoints(intfct) ;
    int    p ;
    
    for(p = 0 ; p < np ; p++) {
      double* vim_n  = vim0_n + p*NVI ;
      /* Variables */
      double* x = ComputeVariables(el,u,u_n,vim0_n,dt,p) ;
      double* c0 = c + p*dec ;


      /* initialization */
      {
        int i ;
      
        for(i = 0 ; i < dec ; i++) c0[i] = 0. ;
      }
      
      
      /* The derivative of equations w.r.t unknowns */
      

      /* Derivatives w.r.t strains */
      {
        double deps = 1.e-6 ;
        double* dx = ComputeVariableDerivatives(el,dt,x,deps,I_EPS) ;
        
        /* Tangent stiffness matrix */
        {
          double* c1 = c0 ;
          int i ;
    
          double K_s = young/(3 - 6*poisson) ;
          double G_s = young/(2 + 2*poisson) ;
          double n   = Poro_En ;
    
          double K   = BulkModulusMoriTanaka(K_s,G_s,n) ;
          double G   = ShearModulusMoriTanaka(K_s,G_s,n) ;

          double eta_v_mt = BulkModulusMoriTanaka(eta_v,eta_d,n) ;
          double eta_d_mt = ShearModulusMoriTanaka(eta_v,eta_d,n);
          
          double Kt  = K / (1 + (K*dt)/eta_v_mt) ;
          double Gt  = G / (1 + (2*G*dt)/eta_d_mt) ;
          double Lambdat = Kt - 2./3.*Gt ;

          for(i = 0 ; i < 3 ; i++) {
            int j ;
            
            for(j = 0 ; j < 3 ; j++) {
              C1(i,i,j,j) += Lambdat ;
              C1(i,j,i,j) += Gt ;
              C1(i,j,j,i) += Gt ;
            }
          }
        }
        
        /* Coupling matrices */
        {
          double* c00 = c0 + 81 + 2*9 ;
          
          /* assuming to be the same for the derivatives wrt I_EPS+4 and I_EPS+8
           * and zero for the derivatives w.r.t others */
          {
            double* c1 = c00 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 11 ;
            int i ;

            for(i = 0 ; i < 3 ; i++) B1(i,i) = dx[I_M_SALT] ;
          }
        }
      }
      
      /* Derivatives w.r.t pressure */
      {
        double  dp = dui[U_Mass] ;
        double* dx = ComputeVariableDerivatives(el,dt,x,dp,I_P) ;
        
        /* Tangent Biot's coefficient */
        {
          double* dsigdp = dx + I_SIG ;
          double* c1 = c0 + 81 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdp[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 2*9 + 9 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 11 ;
        
            c1[0] = dx[I_M_SALT] ;
          }
        }
      }
      
      /* Derivatives w.r.t salt mass fraction */
      {
        double  dc = dui[U_Salt] ;
        double* dx = ComputeVariableDerivatives(el,dt,x,dc,I_W_S) ;
        
        /* Tangent Biot's coefficient */
        {
          double* dsigdc = dx + I_SIG ;
          double* c1 = c0 + 81 + 9 ;
          int i ;

          for(i = 0 ; i < 9 ; i++) c1[i] = dsigdc[i] ;
        }
        
    
        /* General storage matrix */
        {
          double* c00 = c0 + 81 + 2*9 + 9 + 1 ;
          
          {
            double* c1 = c00 ;
        
            c1[0] = dx[I_M_TOT] ;
          }
          {
            double* c1 = c00 + 11 ;
        
            c1[0] = dx[I_M_SALT] ;
          }
        }
      }
    }
  }

  return(dec) ;
#undef C1
#undef B1
#undef T2
#undef T4
}



int ComputeTransferCoefficients(Element_t* el,double dt,double* c)
/*
**  Conduction matrix (c) and shift (dec)
*/
{
  double* ve0 = Element_GetExplicitTerm(el) ;
  double** u_n = Element_ComputePointerToPreviousNodalUnknowns(el) ;
  IntFct_t*  intfct = Element_GetIntFct(el) ;
  int np = IntFct_GetNbOfPoints(intfct) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  int    dec = 4 * 9 ;
  int    p ;
  

  for(p = 0 ; p < np ; p++) {
    double* c1 = c + p*dec ;
    double w_s_n = FEM_ComputeUnknown(fem,u_n,intfct,p,U_Salt) ;
    int i ;
    
    /* initialization */
    for(i = 0 ; i < dec ; i++) c1[i] = 0. ;
    
    {
      double* vex  = ve0 + p*NVE ;
      
      /* 1. Transport of liquid
       * ---------------------- */
      /* Darcy permeability tensor */
      {
        double* c2 = c1 ;
        
        c2[0] = -K_Darcy ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      
      /* Osmotic permeability tensor */
      {
        double* c2 = c1 + 9 ;
        
        c2[0] = K_Osmosis ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      
      /* 2. Transport of salt 
       * -------------------- */
      /* Advective term */
      {
        double* c2 = c1 + 2 * 9 ;
        
        c2[0] = -(1 - Tau) * w_s_n * K_Darcy ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
      
      /* Fick diffusion + advective term */
      {
        double* c2 = c1 + 3 * 9 ;
        
        c2[0] = -D_Fick + (1 - Tau) * w_s_n * K_Osmosis ;
        c2[4] = c2[0] ;
        c2[8] = c2[0] ;
      }
    }
  }

  return(dec) ;
}






double* ComputeVariables(Element_t* el,double** u,double** u_n,double* f_n,double dt,int m)
{
  IntFct_t* intfct = Element_GetIntFct(el) ;
  FEM_t*    fem    = FEM_GetInstance(el) ;
  //Model_t*  model  = Element_GetModel(el) ;
  int dim = Element_GetDimensionOfSpace(el) ;
  //double*   x      = Model_GetVariable(model,p) ;
  double*  x   = Variable ;
  double*  x_n = Variable_n ;
  
    
  /* Primary Variables */
  {
    int    i ;
    
    /* Displacements */
    for(i = 0 ; i < dim ; i++) {
      x[I_U + i] = FEM_ComputeUnknown(fem,u,intfct,m,U_Mech + i) ;
    }
    
    for(i = dim ; i < 3 ; i++) {
      x[I_U + i] = 0 ;
    }
    
    /* Strains */
    {
      double* eps =  FEM_ComputeLinearStrainTensor(fem,u,intfct,m,U_Mech) ;
    
      for(i = 0 ; i < 9 ; i++) {
        x[I_EPS + i] = eps[i] ;
      }
      
      FEM_FreeBufferFrom(fem,eps) ;
    }
    
    /* Pressure */
    x[I_P] = FEM_ComputeUnknown(fem,u,intfct,m,U_Mass) ;
    
    /* Pressure gradient */
    {
      double* grd_p = FEM_ComputeUnknownGradient(fem,u,intfct,m,U_Mass) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_P + i] = grd_p[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd_p) ;
    }
    
    /* Salt mass fraction */
    x[I_W_S] = FEM_ComputeUnknown(fem,u,intfct,m,U_Salt) ;
    
    /* Salt mass fraction gradient */
    {
      double* grd_w_s = FEM_ComputeUnknownGradient(fem,u,intfct,m,U_Salt) ;
    
      for(i = 0 ; i < 3 ; i++) {
        x[I_GRD_W_S + i] = grd_w_s[i] ;
      }
      
      FEM_FreeBufferFrom(fem,grd_w_s) ;
    }
  }
  
  
  /* Needed variables to compute secondary variables */
  {
    int    i ;
    
    /* Implicit terms at previous time step */
    {
      double* eps_n =  FEM_ComputeLinearStrainTensor(fem,u_n,intfct,m,U_Mech) ;
      double* vim_n = f_n + m*NVI ;
    
      for(i = 0 ; i < 9 ; i++) {
        x_n[I_EPS   + i] = eps_n[i] ;
        x_n[I_SIG   + i] = SIG_n[i] ;
        x_n[I_SIG_D + i] = SIG_D_n[i] ;
      }
      
      x_n[I_SIG_M] = SIG_M_n ;
      
      FEM_FreeBufferFrom(fem,eps_n) ;
    
      /* Pressure at previous time step */
      x_n[I_P] = FEM_ComputeUnknown(fem,u_n,intfct,m,U_Mass) ;
    
      /* Salt mass fraction at previous time step */
      x_n[I_W_S] = FEM_ComputeUnknown(fem,u_n,intfct,m,U_Salt) ;
    
      /* Crystal salt volumique fraction at previous time step */
      x_n[I_PHI_C] = Phi_c_n ;
    
      /* liquide density at previous time step */
      x_n[I_RHO_L] = Rho_ln ;
    
      /* Eulerian porosity at previous time step */
      x_n[I_Poro_E] = Poro_En ;
    }
    
    /* Explicit terms (Transfer coefficients) */
    {
      double* vex0 = Element_GetExplicitTerm(el) ;
      double* vex  = vex0 + m*NVE ;
      
      x[I_K_L]  = K_l ;
      x[I_D]    = D ;
      x[I_TAU]  = Tau ;
      
      x[I_K_Darcy]   = K_Darcy ;
      x[I_K_Osmosis] = K_Osmosis ;
      x[I_D_Fick]    = D_Fick ;
    }
  }
    
  ComputeSecondaryVariables(el,dt,x_n,x) ;
  
  return(x) ;
}



void  ComputeSecondaryVariables(Element_t* el,double dt,double* x_n,double* x)
{
  /* Strains */
  double* eps   =  x   + I_EPS ;
  double* eps_n =  x_n + I_EPS ;
  /* Pressure */
  double  p   = x[I_P] ;
  double  p_n = x_n[I_P] ;
  /* Solute mass fraction */
  double w_s   = x[I_W_S] ;
  double w_s_n = x_n[I_W_S] ;
  
  /* Backup */
  
  /* Fluid densities */
  {
    double rho_w = WaterDensity(p) ;
    double rho_l = LiquidDensity(p,w_s) ;
    
    x[I_RHO_W] = rho_w ;
    x[I_RHO_L] = rho_l ;
  }
  
  /* Porosities */
  {
    double K_s = young/(3 - 6*poisson) ;
    double G_s = young/(2 + 2*poisson) ;
    double n_n = x_n[I_Poro_E] ;
    double biot = BiotCoefficientMoriTanaka(K_s,G_s,n_n) ;
    /** Strain */
    double tre   = eps[0] + eps[4] + eps[8] ;
    double dphi_l_epsi = biot*tre ;
    /** crystal salt dissolution */
    double phi_c_n = x_n[I_PHI_C] ;
    double phi_c = CrystalVolumeFraction(phi_c_n,w_s,dt) ;
    double dphi_l_salt = phi_c_i - phi_c ;
    /** Porosity evolution */
    double phi_l = phi_l_i + dphi_l_epsi + dphi_l_salt ;
    double n   = phi_l / (1 + tre) ;
    
    if (phi_l < 0) {
      Message_Direct("\nNegative porosity: %e\n",phi_l) ;
      Exception_Interrupt ;
    }
    
    x[I_PHI_C] = phi_c ;
    x[I_PHI_L] = phi_l ;
    x[I_Poro_E] = n ;
  }
  
  /* Mass contents */
  {
    double rho_l = x[I_RHO_L] ;
    double phi_c = x[I_PHI_C] ;
    double phi_l = x[I_PHI_L] ;
    double m_total = rho_l*phi_l     + rho_c*phi_c ;
    double m_salt  = rho_l*phi_l*w_s + rho_c*phi_c ;
    
    x[I_M_TOT]  = m_total ;
    x[I_M_SALT] = m_salt ;
  }
  
    
  /* Fluxes */
  {
    /* Transfer coefficient */
    double tau = x[I_TAU] ;
    double k_darcy = x[I_K_Darcy] ;
    double k_osmosis = x[I_K_Osmosis] ;
    double d_fick = x[I_D_Fick] ;
    
    /* Gradients */
    double* grd_p   = x + I_GRD_P ;
    double* grd_w_s = x + I_GRD_W_S ;
    
    /* Fluxes */
    double* w_f = x + I_W_F ;
    double* w_d = x + I_W_D ;
    double* w_o = x + I_W_O ;
    double* w_l = x + I_W_L ;
    double* w_salt = x + I_W_SALT ;
    int i ;
    
    for(i = 0 ; i < 3 ; i++){
      w_d[i]    = - k_darcy * grd_p[i] ;
      w_o[i]    =   k_osmosis * grd_w_s[i] ;
      w_l[i]    =   w_d[i] + w_o[i] ;
      w_f[i]    = - d_fick * grd_w_s[i] ;
      w_salt[i] =   w_f[i] + (1 - tau)*w_s_n*w_l[i] ;
    }
  }
  
  /* Stresses */
  {
    double* sig     = x   + I_SIG ;
    double* sig_d   = x   + I_SIG_D ;
    
    /* Mori Tanaka  */
    double K_s = young/(3 - 6*poisson) ;
    double G_s = young/(2 + 2*poisson) ;
    double n   = x_n[I_Poro_E] ; /* Explicit in porosity */
    
    double K   = BulkModulusMoriTanaka(K_s,G_s,n) ;
    double G   = ShearModulusMoriTanaka(K_s,G_s,n) ;
  
    double eta_v_mt = BulkModulusMoriTanaka(eta_v,eta_d,n) ;
    double eta_d_mt = ShearModulusMoriTanaka(eta_v,eta_d,n) ;
    
    double biot = BiotCoefficientMoriTanaka(K_s,G_s,n) ;
    
    double  deps[9] ;
    int    i ;
      
    /* Incremental deformations */
    for(i = 0 ; i < 9 ; i++) {
      deps[i] =  eps[i] - eps_n[i] ;
    }
    
    {     
      double  deps_d[9] ;
      double trde = deps[0] + deps[4] + deps[8] ;
      
      /* Deviatoric strains */
      for(i = 0 ; i < 9 ; i++) deps_d[i] = deps[i] ;
      
      deps_d[0] -= trde / 3 ;
      deps_d[4] -= trde / 3 ;
      deps_d[8] -= trde / 3 ;
      
      /* Deviatoric stresses */
      {
        double* sig_d_n = x_n + I_SIG_D ;
        
        for(i = 0 ; i < 9 ; i++) {
          sig_d[i] = (sig_d_n[i] + 2*G*deps_d[i])/(1 + (2*G*dt)/eta_d_mt) ;
        }
      }
  
      /* Mean effective stress */
      {
        double sig_m_n = x_n[I_SIG_M] ;
        double sig_m   = (sig_m_n + K*trde)/(1 + (K*dt)/eta_v_mt) ;
      
        x[I_SIG_M] = sig_m ;
      }
  
      /* Total stresses */
      {
        double sig_m = x[I_SIG_M] - biot*(p - p_i) ;
        
        for(i = 0 ; i < 9 ; i++) sig[i] = sig_d[i] ;

        sig[0] += sig_m ;
        sig[4] += sig_m ;
        sig[8] += sig_m ;
      }
    }
  }
}



double* ComputeVariableDerivatives(Element_t* el,double dt,double* x,double dxi,int i)
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
  
  ComputeSecondaryVariables(el,dt,x_n,dx) ;
  
  /* The numerical derivative as (f(x + dx) - f(x))/dx */
  for(j = 0 ; j < NbOfVariables ; j++) {
    dx[j] -= x[j] ;
    dx[j] /= dxi ;
  }

  return(dx) ;
}

