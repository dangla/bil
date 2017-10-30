#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "Common.h"
#include "FVM.h"


#define TITLE   "Drying-Wetting with Salt (1D case)"
#define AUTHORS "Dangla"

#include "PredefinedMethods.h"


/* Nb of equations */
#define NEQ 	  (2)
/* Nb of (im/ex)plicit terms and constant terms */
#define NVI     (12)
#define NVE     (9)
#define NV0     (0)

/* Equation index */
#define E_eau	  (0)
#define E_salt	(1)

/* Unknown index */
#define U_H_r   (0)
#define U_C_s   (1)


#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])
#define UNKNOWNn(n,i)    (u_n[n][Element_GetNodalUnknownPosition(el,n,i)])


/* We define some names for primary unknowns */
#define H_r(n)      (UNKNOWN(n,U_H_r))
#define C_s(n)      (UNKNOWN(n,U_C_s))

#define M_W(n)      (f[(n)])
#define N_S(n)      (f[(n+2)])
#define W_W         (f[4])
#define W_S         (f[5])
#define P_l(n)      (f[(8+n)])
#define P_c(n)      (f[(10+n)])

#define M_Wn(n)     (f_n[(n)])
#define N_Sn(n)     (f_n[(n+2)])

#define KD_W        (va[(0)])
#define KF_V        (va[(1)])

#define KD_A        (va[(2)])
#define KD_C        (va[(3)])
#define KF_A        (va[(4)])
#define KF_C        (va[(5)])

#define KD_S        (va[(6)])
#define KF_S        (va[(7)])
#define KM_S        (va[(8)])

/* valences */
#define z_cl      (-1.)
#define z_na      (1.)
#define z_so4     (-2.)

/* coefficients stoechiometriques */
#define nu_h2o_nacl   (0)
#define nu_na_nacl    (1)
#define nu_cl_nacl    (1)
#define nu_h2o_na2so4 (10)
#define nu_na_na2so4  (2)
#define nu_so4_na2so4 (1)

/* volumes molaires partiels liquides (m3/mole) */
#define v_na      (1.87e-5)
#define v_cl      (2.52e-6)
#define v_so4     (-8.334e-6)
#define V_H2O     (1.80e-5)

/* volumes molaires solides (m3/mole) */
#define v_nacl    (24.5e-6)
#define v_na2so4  (220.e-6)

/* coefficients de diffusion moleculaire (m2/s) */
#define do_cl     (2.032e-9)
#define do_na     (1.334e-9)
#define do_so4    (1.065e-9)

#define do_va     (2.42e-5)

/* constante d'equilibre */
#define K_nacl    (6.e3)      /* Solubilite de NaCl (moles/m3) */
#define K_na2so4  (1.24166e3) /* Solubilite de Na2SO4.10H2O (moles/m3) */

/* tensions superficielles (N/m) */
#define Gamma_cl_nacl   (0.1)
#define Gamma_cl_na2so4 (0.1)

#define GAMMA_LG        (0.07)


/* constantes physiques */
#define FARADAY   (9.64846e4) /* Faraday (C/mole) */
#define T         (293.)      /* Temperature (K) */
#define RT        (2436.)     /* Produit R = 8.3143 et T = 293. (J/mole) */

/* viscosites (Pa.s) */
#define mu_g      (1.8e-5)
#define mu_l      (1.002e-3)

/* Masses molaires (kg/mole) */
#define M_h2o     (1.8e-2)
#define M_air     (2.896e-2)

/* autres */
#define p_atm     (1.01325e5) /* Pression atmospherique (Pa) */
#define p_g       (0.)        /* Pression du gaz (Pa) */

/* Pression de vapeur (Pa) */
#define P_VS(T)	  (609.14*pow(10.,7.45*(T - 273.)/(T - 38.)))


/* Types de sel */
#define NaCl    (0)
#define Na2SO4  (1)
#define KNO3    (2)

/* Choix du sel */
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
#define D_A       do_so4
#define D_C       do_na
#define K_SALT    K_na2so4
#define GAMMA_CL  Gamma_cl_na2so4
#else
#error "Type de sel non prevu"
#endif

#define NU_AC     (NU_A + NU_C)
#define V_AC      (NU_A*V_A + NU_C*V_C)
#define V_ACH     (NU_H2O*V_H2O + V_AC)




/* Material properties 
 * ------------------- */
#define SATURATION_CURVE           (Element_GetCurve(el) + 0)
#define RELATIVEPERMLIQ_CURVE      (Element_GetCurve(el) + 1)
#define TORTUOSITYLIQ_CURVE        (Element_GetCurve(el) + 2)

#define RelativePermeabilityToLiquid(x)  (Curve_ComputeValue(RELATIVEPERMLIQ_CURVE,x))
#define RelativePermeabilityToGas(x)     (Curve_ComputeValue(RELATIVEPERMGAS_CURVE,x))
#define TortuosityToGas(f,sg)            ((sg > 0) ? pow(f,2.67)*pow(sg,4.67) : 0)
#define TortuosityToLiquid(f,sl)         (tortuosite_l(phi)*Curve_ComputeValue(TORTUOSITYLIQ_CURVE,sl))
#define aa                        (1.74)  /* 1/3 Millington, Thiery 1.74 */
#define bb                        (3.2)   /* 7/3 Millington, Thiery 3.2 */


/* To retrieve the material properties */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])

/* Fonctions */
static int     pm(const char *s) ;
static void    GetProperties(Element_t*) ;

static double* ComputeComponents(Element_t*,double**,double*,double,int) ;
static void    ComputeSecondaryComponents(Element_t*,double,double*) ;
static double* ComputeComponentDerivatives(Element_t*,double,double*,double,int) ;

static void    ComputeTransferCoefficients(Element_t*,double**,double*) ;
static double* ComputeFluxes(Element_t*,double**) ;
static double* Fluxes(Element_t*,double*) ;

extern double lna_i(double,double,double,double,double,double) ;
extern double lng_TQN(double,double,double,double,double,double,double,double) ;

static double tortuosite_l(double) ;

static double saturation_l(double,double,Curve_t*) ;
static double saturation_c(double,double,Curve_t*) ;

static double dsaturation_ll(double,double,Curve_t*) ;
static double dsaturation_lc(double,double,Curve_t*) ;
static double dsaturation_cl(double,double,Curve_t*) ;
static double dsaturation_cc(double,double,Curve_t*) ;

static double dsaturation(double,double,Curve_t*,double (*)(double,double,Curve_t*),char) ;

static double activite_w(double,double) ;
static double activite_s(double,double) ;
static double activite_w_ideal(double,double) ;
static double activite_s_ideal(double,double) ;

static double (*xactivite_w[])(double,double) = {activite_w,activite_w_ideal} ;
static double (*xactivite_s[])(double,double) = {activite_s,activite_s_ideal} ;

#define SATURATION_L(x,y)   saturation_l(x,y,SATURATION_CURVE)
#define SATURATION_C(x,y)   saturation_c(x,y,SATURATION_CURVE)

#define DSATURATION_LL(x,y) dsaturation_ll(x,y,SATURATION_CURVE)
#define DSATURATION_LC(x,y) dsaturation_lc(x,y,SATURATION_CURVE)
#define DSATURATION_CL(x,y) dsaturation_cl(x,y,SATURATION_CURVE)
#define DSATURATION_CC(x,y) dsaturation_cc(x,y,SATURATION_CURVE)


/*
#define DSATURATION_LL(x,y) (-dsaturation(x,y,el.mat->cb[0],saturation_l,'x'))
#define DSATURATION_LC(x,y) dsaturation(x,y,el.mat->cb[0],saturation_l,'y')
#define DSATURATION_CL(x,y) (-dsaturation(x,y,el.mat->cb[0],saturation_c,'x'))
#define DSATURATION_CC(x,y) dsaturation(x,y,el.mat->cb[0],saturation_c,'y')
*/

/*
#define ACTIVITE_W(a)   courbe(a,el.mat->cb[3])
#define ACTIVITE_S(a)   courbe(a,el.mat->cb[4])
*/
#define ACTIVITE_W(a)     xactivite_w[1](a,T)
#define ACTIVITE_S(a)     xactivite_s[1](a,T)

/* Parametres */
static double phi0,r_d,k_int ;
static double d_cl ;
static double lna_w0 ;
static double lna_s0 ;
#define NN                (2)
#define NbOfComponents    (16)
static double Components[NN][NbOfComponents] ;
static double dComponents[NbOfComponents] ;


#define I_M_W          (2)
#define I_N_S          (3)

#define I_P_L          (4)
#define I_P_C          (5)
#define I_S_L          (6)
#define I_S_C          (7)
#define I_S_G          (8)
#define I_LNA_W        (9)
#define I_LNA_S        (10)
#define I_RHO_V        (11)
#define I_RHO_W        (12)
#define I_C_W          (13)
#define I_C_S          (14)
#define I_H_R          (15)
  
  

#define NbOfComponentFluxes    (6)
static double ComponentFluxes[NbOfComponentFluxes] ;

#define I_W_W           (0)
#define I_W_S           (1)

#define I_W_L           (2)
#define I_W_V           (3)
#define I_W_Ani         (4)
#define I_W_Cat         (5)



int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0) return (0) ;
  else if(strcmp(s,"D_Cl") == 0) return (1);  
  else if(strcmp(s,"r_d") == 0) return (2);
  else if(strcmp(s,"k_int") == 0) return (3);
  else if(strcmp(s,"lna_w0") == 0) return (4);
  else if(strcmp(s,"lna_s0") == 0) return (5);
  else if(strcmp(s,"courbes") == 0) return (6);
  else return (-1) ;
}


void GetProperties(Element_t *el)
{
  phi0    = GetProperty("porosite") ;
  r_d     = GetProperty("r_d") ;
  k_int   = GetProperty("k_int") ;
  d_cl    = GetProperty("D_Cl") ;
  lna_w0  = GetProperty("lna_w0") ;
  lna_s0  = GetProperty("lna_s0") ;
}


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_eau,"mass") ;
  Model_CopyNameOfEquation(model,E_salt,"air") ;

  Model_CopyNameOfUnknown(model,U_H_r,"h_r") ;
  Model_CopyNameOfUnknown(model,U_C_s,"c_s") ;
  
  return(0) ;
}


int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    NbOfProp = 7 ;

  {
    int dim = Material_GetDimension(mat) ;
    
    if(dim > 1) arret("M28: dimension > 1 not available") ;
  }

  Material_ScanProperties(mat,datafile,pm) ;


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
          
        y_lnaw[i] = ACTIVITE_W(c_s) ;
        y_lnas[i] = ACTIVITE_S(c_s) ;
      }
      free(x) ;
    }

    /* on met a jour le nb de proprietes */
    Material_GetNbOfCurves(mat) += 2 ;
  }
#endif

  {
    double c_s = K_SALT ;
    
    Material_GetProperty(mat)[pm("lna_w0")] = ACTIVITE_W(c_s) ;
    Material_GetProperty(mat)[pm("lna_s0")] = ACTIVITE_S(c_s) ;
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
  printf("The equations to be solved are:\n") ;
  printf("\t- Water Mass balance     (liq)\n") ;
  printf("\t- Salt Mass Balance      (sel)\n") ;
  
  printf("\n") ;
  printf("The primary unknowns are:\n") ;
  printf("\t- Relative humidity        (h_r)\n") ;
  printf("\t- Salt concentration       (c_s)\n") ;

  printf("\n") ;
  printf("Some other informations\n") ;
  printf("Example of input data\n") ;
  printf("\n") ;

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
  int    i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;

  
  /*
    Input data
  */
  GetProperties(el) ;
  

  /* Contenus molaires */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;

    /* contenus molaires en cations et anions */
    N_S(i)   = x[I_N_S] ;
    /* masse d eau */
    M_W(i)   = x[I_M_W] ;

    /* sauvegarde */
    P_l(i)    = x[I_P_L] ;
    P_c(i)    = x[I_P_C] ;
  }

  /* Transfer coefficients */
  ComputeTransferCoefficients(el,u,f) ;

  /* Fluxes */
  {
    double* w = ComputeFluxes(el,u) ;

    W_W     = w[I_W_W] ;
    W_S     = w[I_W_S] ;
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
   
  /* Contenus molaires */

  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;

    /* Backup */
    N_S(i)   = x[I_N_S] ;
    M_W(i)   = x[I_M_W] ;

    /* sauvegarde */
    P_l(i)    = x[I_P_L] ;
    P_c(i)    = x[I_P_C] ;
    
    {
      double h_r = x[I_H_R] ;
      double c_s = x[I_C_S] ;
      double c_w = x[I_C_W] ;

      if(h_r <= 0 || c_s < 0. || c_w < 0.) {
        double xx = Element_ComputePointerToNodalCoordinates(el)[i][0] ;
        double s_l = x[I_S_L] ;
        double s_c = x[I_S_C] ;
      
        printf("\n") ;
        printf("x       = %e\n",xx) ;
        printf("h_r     = %e\n",h_r) ;
        printf("c_s     = %e\n",c_s) ;
        printf("c_w     = %e\n",c_w) ;
        printf("s_l     = %e\n",s_l) ;
        printf("s_c     = %e\n",s_c) ;
        return(-1) ;
      }
    }
    /* printf("\nS_L = %e ; S_C = %e\n",s_l,s_c) ; */
    /* printf("\nc_s = %e\n",c_s) ; */
  }
  

  /* Fluxes */
  {
    double* w = ComputeFluxes(el,u) ;

    W_W     = w[I_W_W] ;
    W_S     = w[I_W_S] ;
  }
  
  return(0) ;
}


int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrix (k) */
{
#define K(i,j)    (k[(i)*2*NEQ + (j)])
  double* f  = Element_GetCurrentImplicitTerm(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  int    nn = Element_GetNbOfNodes(el) ;
  int    ndof = nn*NEQ ;
  FVM_t* fvm = FVM_GetInstance(el) ;
  double* volume = FVM_ComputeCellVolumes(fvm) ;
  double surf ;
  int    i ;
  double zero = 0. ;
  
  /* initialisation */
  for(i = 0 ; i < ndof*ndof ; i++) k[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  
  /*
    Input data
  */
  GetProperties(el) ;
  
  /* Boundary Surface Area */
  {
    double *area = FVM_ComputeCellSurfaceAreas(fvm) ;
    surf = area[1] ;
  }

  double Dp_lSDh_r[2] ;
  double Dp_lSDc_s[2] ;

  /* termes d'accumulation */
  for(i = 0 ; i < 2 ; i++) {
    /* Components */
    double *x = ComputeComponents(el,u,f,0,i) ;
    
    double c_s    = x[I_C_S] ; 
    double h_r    = x[I_H_R] ;
    double p_vs   = P_VS(T) ; 
    /* activite de l'eau */
    double lna_w  = x[I_LNA_W] ;
    /* activite du sel */
    double lna_s  = x[I_LNA_S] ;
    /* pressions */
    double p_l = x[I_P_L] ;
    double p_c = x[I_P_C] ;
    /* masses volumiques */
    double rho_v = x[I_RHO_V] ;
    double rho_w = x[I_RHO_W] ;
    /* saturations */
    double s_l  = x[I_S_L] ;
    double s_g  = x[I_S_G] ;

    /* derivees */
    /* ... par rapport aux pressions */	
    double ds_lsdp_l = DSATURATION_LL(p_g - p_l,p_c - p_g) ;
    double ds_lsdp_c = DSATURATION_LC(p_g - p_l,p_c - p_g) ;
    double ds_csdp_l = DSATURATION_CL(p_g - p_l,p_c - p_g) ;
    double ds_csdp_c = DSATURATION_CC(p_g - p_l,p_c - p_g) ;
    double ds_gsdp_l = - ds_lsdp_l - ds_csdp_l ;
    double ds_gsdp_c = - ds_lsdp_c - ds_csdp_c ;
    /* ... par rapport a h_r */
    double dp_lsdh_r    = RT/V_H2O/h_r ;
    double dp_csdh_r    = V_ACH*dp_lsdh_r/V_SALT ;
    double ds_lsdh_r    = ds_lsdp_l*dp_lsdh_r + ds_lsdp_c*dp_csdh_r ;
    double ds_csdh_r    = ds_csdp_l*dp_lsdh_r + ds_csdp_c*dp_csdh_r ;
    double ds_gsdh_r    = ds_gsdp_l*dp_lsdh_r + ds_gsdp_c*dp_csdh_r ;
    double drho_vsdh_r  = M_h2o/RT*p_vs ;
    double dn_wsdh_r    = phi0*(drho_vsdh_r*s_g + rho_v*ds_gsdh_r	\
				+ rho_w*ds_lsdh_r) ;
    double dn_ssdh_r    = phi0*(ds_lsdh_r*c_s + ds_csdh_r/V_SALT) ;
    /* ... par rapport a c_s */
    double dc_s        = K_nacl*1.e-4 ;
    double c_s2        = c_s + dc_s ;
    double lna_w2      = ACTIVITE_W(c_s2) ;
    double lna_s2      = ACTIVITE_S(c_s2) ;
    double dlna_wsdc_s = (lna_w2 - lna_w)/dc_s ;
    double dlna_ssdc_s = (lna_s2 - lna_s)/dc_s ;
    double dp_lsdc_s   = -RT/V_H2O*dlna_wsdc_s ;
    double dp_csdc_s   = V_ACH*dp_lsdc_s/V_SALT				\
                       + RT*(NU_H2O*dlna_wsdc_s + dlna_ssdc_s)/V_SALT ;
    double ds_lsdc_s   = ds_lsdp_l*dp_lsdc_s + ds_lsdp_c*dp_csdc_s ;
    double ds_csdc_s   = ds_csdp_l*dp_lsdc_s + ds_csdp_c*dp_csdc_s ;
    double ds_gsdc_s   = ds_gsdp_l*dp_lsdc_s + ds_gsdp_c*dp_csdc_s ;
    double drho_wsdc_s = -M_h2o*V_AC/V_H2O ;
    double dn_wsdc_s   = phi0*(drho_wsdc_s*s_l + rho_w*ds_lsdc_s	\
			       + rho_v*ds_gsdc_s) ;
    double dn_ssdc_s   = phi0*(s_l + ds_lsdc_s*c_s + ds_csdc_s/V_SALT) ;
    /*
      Conservation du sel : (n_s1 - n_sn) + dt * div(w_s) = 0
    */
    K(i*NEQ+E_salt,i*NEQ+U_C_s)  += volume[i]*NU_AC*dn_ssdc_s ;
    K(i*NEQ+E_salt,i*NEQ+U_H_r)  += volume[i]*NU_AC*dn_ssdh_r ;

    /*
      Conservation de l'eau : (n_w1 - n_wn) + dt * div(w_w) = 0
    */
    K(i*NEQ+E_eau,i*NEQ+U_H_r)  += volume[i]*dn_wsdh_r ;
    K(i*NEQ+E_eau,i*NEQ+U_C_s)  += volume[i]*dn_wsdc_s ;

    /* sauvegardes pour les termes de transport */
    Dp_lSDh_r[i] = dp_lsdh_r ;
    Dp_lSDc_s[i] = dp_lsdc_s ;
  }

  /* termes d'ecoulement */
  {
    double** coor = Element_ComputePointerToNodalCoordinates(el) ;
    double dx = coor[1][0] - coor[0][0] ;
    double tr  = dt*surf/dx ;
    double dz2 = KF_A*Z_A*Z_A + KF_C*Z_C*Z_C ;
    double dpsisdc_s = -(Z_A*KF_A*NU_A + Z_C*KF_C*NU_C)/dz2 ;
    
    double trd_s = tr*KD_S ;
    double trf_s = tr*KF_S ;
    double tre_s = tr*KM_S ;

    double trd_w  = tr*KD_W ;
    double trf_v  = tr*KF_V ;

    double c[2] ;

  /*
    Conservation du sel : (n_s1 - n_sn) + dt * div(w_s) = 0
  */
  for(i=0;i<2;i++) {
    c[i] = trd_s*Dp_lSDc_s[i] + (trf_s + tre_s*dpsisdc_s) ;
  }
  K(E_salt,U_C_s)            += + c[0] ;
  K(E_salt,U_C_s+NEQ)        += - c[1] ;
  K(E_salt+NEQ,U_C_s)        += - c[0] ;
  K(E_salt+NEQ,U_C_s+NEQ)    += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = trd_s*Dp_lSDh_r[i] ;
  }
  K(E_salt,U_H_r)             += + c[0] ;
  K(E_salt,U_H_r+NEQ)         += - c[1] ;
  K(E_salt+NEQ,U_H_r)         += - c[0] ;
  K(E_salt+NEQ,U_H_r+NEQ)     += + c[1] ;
  
  /*
    Conservation de l'eau : (n_w1 - n_wn) + dt * div(w_w) = 0
  */
  for(i=0;i<2;i++) {
    c[i] = trd_w*Dp_lSDc_s[i] ;
  }
  K(E_eau,U_C_s)           += + c[0] ;
  K(E_eau,U_C_s+NEQ)       += - c[1] ;
  K(E_eau+NEQ,U_C_s)       += - c[0] ;
  K(E_eau+NEQ,U_C_s+NEQ)   += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = trd_w*Dp_lSDh_r[i] + trf_v ;
  }
  K(E_eau,U_H_r)            += + c[0] ;
  K(E_eau,U_H_r+NEQ)        += - c[1] ;
  K(E_eau+NEQ,U_H_r)        += - c[0] ;
  K(E_eau+NEQ,U_H_r+NEQ)    += + c[1] ;
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
    Conservation of salt mass  : (n_s1 - n_sn) + dt * div(w_s) = 0
  */
  R(0,E_salt) -= volume[0]*(N_S(0) - N_Sn(0)) + dt*surf*W_S ;
  R(1,E_salt) -= volume[1]*(N_S(1) - N_Sn(1)) - dt*surf*W_S ;

  /*
    Conservation of water mass  : (n_w1 - n_wn) + dt * div(w_V) = 0
  */
  R(0,E_eau) -= volume[0]*(M_W(0) - M_Wn(0)) + dt*surf*W_W ;
  R(1,E_eau) -= volume[1]*(M_W(1) - M_Wn(1)) - dt*surf*W_W ;
  
  return(0) ;

#undef R
}


int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  FVM_t *fvm = FVM_GetInstance(el) ;
  double** u = Element_ComputePointerToNodalUnknowns(el) ;
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
  {
    int    j = FVM_FindLocalCellIndex(fvm,s) ;
    /* Components */
    double* x    = ComputeComponents(el,u,f,0,j) ;
    /* concentrations */
    double c_s    = x[I_C_S] ;
    double h_r    = x[I_H_R] ;
        /* activite de l'eau */
    double lna_w  = x[I_LNA_W] ;
    /* activite du sel */
    double lna_s  = x[I_LNA_S] ;
    double dlna_s = lna_s - lna_s0 ;
    /* pressions */
    double p_l = x[I_P_L] ;
    double p_c = x[I_P_C] ;
    /* saturations */
    double s_l  = x[I_S_L] ;
    double s_c  = x[I_S_C] ;	

    /* Fluxes */
    double* w    = ComputeFluxes(el,u) ;
    double  w_l  = w[I_W_L] ;
    double  w_v  = w[I_W_V] ;
    
    double w_ani =  w[I_W_Ani] ;
    double w_cat =  w[I_W_Cat] ;
    
    /* quantites exploitees */
    i = 0 ;
    Result_Store(r + i++,&h_r,"h_r",1) ;
    Result_Store(r + i++,&c_s,"Concentration sel libre",1) ;
    Result_Store(r + i++,&s_l,"Saturation_liquide",1) ;
    Result_Store(r + i++,&s_c,"Saturation en sel solide",1) ;
    Result_Store(r + i++,&p_l,"pression liquide",1) ;
    Result_Store(r + i++,&p_c,"pression cristal",1) ;
    Result_Store(r + i++,&lna_w,"Log(a_w)",1) ;
    Result_Store(r + i++,&dlna_s,"Log(a_s)",1) ;
    {
      double lna_w_ideal  = activite_w_ideal(c_s,T) ;
      Result_Store(r + i++,&lna_w_ideal,"Log(a_w) ideal",1) ;
    }
    {
      double c_s0 = K_SALT ;
      double lna_s_ideal  = activite_s_ideal(c_s,T) ;
      double lna_s_ideal0 = activite_s_ideal(c_s0,T) ;
      double dlna_s_ideal = lna_s_ideal - lna_s_ideal0 ;
      Result_Store(r + i++,&dlna_s_ideal,"Log(a_s) ideal",1) ;
    }
    {
      double n_sel = phi0*NU_AC*(s_l*c_s + s_c/V_SALT)  ;
      Result_Store(r + i++,&n_sel,"Sel total",1) ;
    }
    Result_Store(r + i++,&w_l,"flux eau liquide",1) ;
    Result_Store(r + i++,&w_v,"flux vapeur",1) ;
    Result_Store(r + i++,&w_ani,"Flux anions",1) ;
    Result_Store(r + i++,&w_cat,"Flux cations",1) ;
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
    
    double c_s    = C_s(i) ;
    double p_vs   = P_VS(T) ;
    /* concentrations */
    double c_w    = (1. - V_AC*c_s)/V_H2O ;
    /* saturations */
    double s_l  = x[I_S_L] ;
    double s_g  = x[I_S_G] ;
    /* permeabilite */
    double k_rl    = RelativePermeabilityToLiquid(s_l) ;
    double k_h     = k_int/mu_l*k_rl ;
    /* tortuosites gaz et liquide*/
    double phi = phi0 ;
    double tau_g   = TortuosityToGas(phi,s_g) ;
    double tau_l   = TortuosityToLiquid(phi,s_l) ;
    /* tortuosites anions et cations */
    double tau_ani = tau_l*d_cl/(tortuosite_l(phi0)*do_cl) ;
    double tau_cat = tau_ani/r_d ;
    
    double kd_w = M_h2o*c_w*k_h ;
    double kd_a = NU_A*c_s*k_h ;
    double kd_c = NU_C*c_s*k_h ;
    double kd_s = kd_a + kd_c ;
    
    double kf_v = tau_g*do_va*M_h2o/RT*p_vs ;
    double kf_a = tau_ani*D_A ;
    double kf_c = tau_cat*D_C ;
    double kf_s = NU_A*kf_a + NU_C*kf_c ;
    
    double km_s = Z_A*kf_a  + Z_C*kf_c ;
  
    /* sauvegarde */
    KD_W   += kd_w ;
    KD_A   += kd_a ;
    KD_C   += kd_c ;
    KD_S   += kd_s ;

    KF_V   += kf_v ;
    KF_A   += kf_a ;
    KF_C   += kf_c ;
    KF_S   += kf_s ;

    KM_S   += km_s ;
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
  double grd_c_s   = grd[I_C_S] ;
  double grd_p_l   = grd[I_P_L] ;
  double grd_h_r   = grd[I_H_R] ;
  double grd_c_ani = NU_A*grd_c_s ;
  double grd_c_cat = NU_C*grd_c_s ;
  double dz2       = KF_A*Z_A*Z_A + KF_C*Z_C*Z_C ;  
  double grd_psi   = -(Z_A*KF_A*grd_c_ani + Z_C*KF_C*grd_c_cat)/dz2 ;
 
  /* Fluxes */
  double w_l = - KD_W*grd_p_l ;
  double w_v = - KF_V*grd_h_r ;
  double w_w =   w_l + w_v ;
  double j_s = - KF_S*grd_c_s - KM_S*grd_psi ;
  double w_s = - KD_S*grd_p_l + j_s ;
  double w_ani =  - KD_A*grd_p_l - KF_A*grd_c_ani - KF_A*Z_A*grd_psi ;
  double w_cat =  - KD_C*grd_p_l - KF_C*grd_c_cat - KF_C*Z_C*grd_psi ;
  
  w[I_W_L ]  = w_l ;
  w[I_W_V ]  = w_v ;
  w[I_W_Ani] = w_ani ;
  w[I_W_Cat] = w_cat ;
  
  w[I_W_W ]  = w_w ;
  w[I_W_S ]  = w_s ;
   
  return(w) ;
}


double tortuosite_l(double phi)
{
  double phic = 0.18 ;
  double dphi = (phi > phic) ? phi - phic : 0. ;
  return(phi*(0.001 + 0.07*phi*phi + 1.8*dphi*dphi)) ;
}

double saturation_l(double p_gl,double p_cg,Curve_t* cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double s_l ;

  if(p_gl > p_cg*y) {
    s_l = Curve_ComputeValue(cb,p_gl) ;
  } else {
    double p_cl = p_cg + p_gl ;
    s_l = Curve_ComputeValue(cb,p_cl*x) ;
  }

  return(s_l) ;
}

double saturation_c(double p_gl,double p_cg,Curve_t* cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double s_c ;

  if(p_gl > p_cg*y) {
    s_c = 1. - Curve_ComputeValue(cb,p_cg*y) ;
  } else {
    double p_cl = p_cg + p_gl ;
    s_c = 1. - Curve_ComputeValue(cb,p_cl*x) ;
  }

  return(s_c) ;
}


double dsaturation_ll(double p_gl,double p_cg,Curve_t* cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_l ;

  if(p_gl > p_cg*y) {
    ds_l = - Curve_ComputeDerivative(cb,p_gl) ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_l = - x*Curve_ComputeDerivative(cb,p_cl*x) ;
  }

  return(ds_l) ;
}


double dsaturation_lc(double p_gl,double p_cg,Curve_t* cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_l ;

  if(p_gl > p_cg*y) {
    ds_l = 0. ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_l = x*Curve_ComputeDerivative(cb,p_cl*x) ;
  }

  return(ds_l) ;
}

double dsaturation_cc(double p_gl,double p_cg,Curve_t* cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_c ;

  if(p_gl > p_cg*y) {
    ds_c = - y*Curve_ComputeDerivative(cb,p_cg*y) ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_c = - x*Curve_ComputeDerivative(cb,p_cl*x) ;
  }

  return(ds_c) ;
}

double dsaturation_cl(double p_gl,double p_cg,Curve_t* cb)
{
  double x = GAMMA_LG/GAMMA_CL ;
  double y = x/(1 - x) ;
  double ds_c ;

  if(p_gl > p_cg*y) {
    ds_c = 0. ;
  } else {
    double p_cl = p_cg + p_gl ;
    ds_c = x*Curve_ComputeDerivative(cb,p_cl*x) ;
  }

  return(ds_c) ;
}

double dsaturation(double x,double y,Curve_t* cb,double (*saturation)(double,double,Curve_t*),char v)
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  char scale = Curve_GetScaleType(cb) ;
  double* pa  = Curve_GetXRange(cb) ;
  double a1 = pa[0] ;
  double a2 = pa[1] ;
  double da ;
  double s1,s2,ds ;

  if(scale == 'n') {
    da = (a2 - a1)/n_i ;
  } else if(scale == 'l') {
    double loga1 = log10(a1),loga2 = log10(a2) ;
    double dloga = (loga2 - loga1)/n_i ;
    double a = (v == 'x') ? x : y ;
    double ada = a*pow(10.,dloga) ;
    da = ada - a ;
  } else {
    arret("dsaturation_1 : option non prevue") ;
  }

  if(v == 'x') {
    s1 = (*saturation)(x - da,y,cb) ;
    s2 = (*saturation)(x + da,y,cb) ;
  } else {
    s1 = (*saturation)(x,y - da,cb) ;
    s2 = (*saturation)(x,y + da,cb) ;
  }

  ds = (s2 - s1)*0.5/da ;

  return(ds) ;
}

double activite_w(double c_s,double Ta)
/* L'activite chimique de l'eau d'une solution saline */
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
  double b0   = sqrt(M_h2o),S0   = pow(M_h2o,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_cl_nacl = 1.827/b0,S_cl_nacl = 19.245/S0 ;
  double b_na_nacl = 4.352/b0,S_na_nacl = 26.448/S0 ;
  /* Na2SO4 (d'apres Lin et Lee) */
  double b_so4_na2so4 = 1.662/b0,S_so4_na2so4 = 0.022/S0 ;
  double b_na_na2so4 = 1.552/b0,S_na_na2so4 = 3.464/S0 ;

  double epsi = 0.0007*(Ta - T0)*(Ta - T0) - 0.3918*(Ta - T0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*Ta,1.5)/b0 ;

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
    lna_w = m_ani*lna_i(Ta,I,Z_A,b_ani,S_ani,A)
          + m_cat*lna_i(Ta,I,Z_C,b_cat,S_cat,A) ;
    return(lna_w) ;
  } else {
    return(0.) ;
  }
}

double activite_s(double c_s,double Ta)
/* L'activite chimique du sel d'une solution saline */
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
  double b0   = sqrt(M_h2o),S0   = pow(M_h2o,1.29) ;
  /* NaCl (d'apres Lin et Lee) */
  double b_cl_nacl = 1.827/b0,S_cl_nacl = 19.245/S0 ;
  double b_na_nacl = 4.352/b0,S_na_nacl = 26.448/S0 ;
  /* Na2SO4 (d'apres Lin et Lee) */
  double b_so4_na2so4 = 1.662/b0,S_so4_na2so4 = 0.022/S0 ;
  double b_na_na2so4 = 1.552/b0,S_na_na2so4 = 3.464/S0 ;

  double epsi = 0.0007*(Ta - T0)*(Ta - T0) - 0.3918*(Ta - T0) + 87.663 ;
  double A    = 1398779.816/pow(epsi*Ta,1.5)/b0 ;
  
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
    double lna_w   = activite_w(c_s,Ta) ;
    double lng_ani = lng_TQN(Ta,I,Z_A,b_ani,S_ani,A,lna_w,m_t) ;
    double lng_cat = lng_TQN(Ta,I,Z_C,b_cat,S_cat,A,lna_w,m_t) ;
    double lna_s   = (NU_A*lng_ani + NU_C*lng_cat) + NU_AC*log(m_s) ;
    return(lna_s) ;
  }
}

double activite_w_ideal(double c_s,double Ta)
/* L'activite chimique de l'eau d'une solution ideale */
{
  double c_w   = (1. - V_AC*c_s)/V_H2O ;
  double c_cat = NU_C*c_s ;
  double c_ani = NU_A*c_s ;
  double c_t   = c_w + c_cat + c_ani ;
  double lna_w = log(c_w/c_t) ;
  return(lna_w) ;
}

double activite_s_ideal(double c_s1,double Ta)
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







double* ComputeComponents(Element_t *el,double **u,double *f_n,double dt,int n)
{
  double *x = Components[n] ;
  
  /* Primary Variables */
  x[U_H_r ] = H_r(n) ;
  x[U_C_s ] = C_s(n) ;
  
  /* Needed variables to compute secondary components */
    
  ComputeSecondaryComponents(el,dt,x) ;
  return(x) ;
}


double* ComputeComponentDerivatives(Element_t *el,double dt,double *x,double dxi,int i)
{
  double *dx = dComponents ;
  int j ;
  
  /* Primary Variables */
  dx[U_H_r ] = x[U_H_r ] ;
  dx[U_C_s ] = x[U_C_s ] ;
  
  /* Needed variables to compute secondary components */
  
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
  double c_s    = x[U_C_s] ; 
  double h_r    = x[U_H_r] ;
    
  double p_vs   = P_VS(T) ; 
  double p_v    = h_r*p_vs ;
    
  /* concentration en eau liquide */
  double c_w    = (1. - V_AC*c_s)/V_H2O ;
  /* activite de l'eau */
  double lna_w  = ACTIVITE_W(c_s) ;
  double dlna_w = lna_w - lna_w0 ;
  /* activite du sel */
  double lna_s  = ACTIVITE_S(c_s) ;
  double dlna_s = lna_s - lna_s0 ;
  /* pressions */
  double p_l = RT/V_H2O*(log(h_r) - lna_w) + p_g ;
  double p_c = p_g + V_ACH*(p_l - p_g)/V_SALT		\
             + RT*(NU_H2O*dlna_w + dlna_s)/V_SALT ;
  /* masses volumiques */
  double rho_v = M_h2o/RT*p_v ;
  double rho_w = M_h2o*c_w ;
  /* saturations */
  double s_l  = SATURATION_L(p_g - p_l,p_c - p_g) ;
  double s_c  = SATURATION_C(p_g - p_l,p_c - p_g) ;
  double s_g  = 1. - s_l - s_c ;

  /* contenus molaires en cations et anions */
  double n_s   = phi0*NU_AC*(s_l*c_s + s_c/V_SALT) ;
  /* masse d eau */
  double m_w   = phi0*(rho_w*s_l + rho_v*s_g) ;

  /* backup */
  x[I_P_L      ] = p_l ;
  x[I_P_C      ] = p_c ;
  x[I_S_L      ] = s_l ;
  x[I_S_C      ] = s_c ;
  x[I_S_G      ] = s_g ;
  x[I_LNA_W    ] = lna_w ;
  x[I_LNA_S    ] = lna_s ;
  x[I_RHO_V    ] = rho_v ;
  x[I_RHO_W    ] = rho_w ;
  x[I_C_W      ] = c_w ;
  x[I_C_S      ] = c_s ;
  x[I_H_R      ] = h_r ;
  
  /* Fluid contents */
  x[I_N_S      ] = n_s ;
  x[I_M_W      ] = m_w ;
}
