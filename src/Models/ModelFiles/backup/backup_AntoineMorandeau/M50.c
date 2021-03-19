/*
Deux reactions avec cinetique :
1. acido-basique : HCO3- + H20 <-> H2CO3 + OH-    (loi de Danckwerts)
2. dissolution de portlandite : Ca(OH)2 <-> Ca2+ + 2OH-    (loi en log)
Cristaux de portlandite spherique, 
Diffusion gaz MT,
Transferts avec diffusion ionique
Electroneutralite
Ajout des alcalins (sodium et potassium) : Na+, K+, NaOH, KOH, NaHCO3, NaCO3-
Ajout d'espèces CaOH+, CaHCO3+, CaCO3aq, CaOH2aq
Prise en compte de la quantitée totale en alcalins
Modele de solution solide pour les C-S-H
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"
#include "FVM.h"

#define TITLE   "Atmospheric Carbonation Of CBM (2012)"
#define AUTHORS "Morandeau"

#include "PredefinedMethods.h"

/* Macros */

#define NEQ    	(9)
#define NVE    	(58)
#define NVI     (30)
#define NV0     (0)

#define E_C       (0)
#define E_q       (1)
#define E_mass    (2)
#define E_Ca      (3)
#define E_k       (4)
#define E_el      (5)
#define E_Na      (6)
#define E_K       (7)
#define E_Si      (8)

#define I_CO2     (0)
#define I_OH      (5)
#define I_P_l     (2)
#define I_CC      (3)
#define I_HCO3    (4)
#define I_psi     (1)
#define I_Na      (6)
#define I_K       (7)
#define I_Si_S    (8)

#define RHO     1
#define LOG_RHO 2
#define Ln10    2.302585093
#define U_CO2   LOG_RHO

#define UNKNOWN(n,i)     (u[n][Element_GetNodalUnknownPosition(el,n,i)])

#if (U_CO2 == LOG_RHO)
#define X_CO2(n)	(exp(Ln10*UNKNOWN(n,I_CO2)))
#else
#define X_CO2(n)	(UNKNOWN(n,I_CO2))
#endif

#define N_Si_S(n)		(UNKNOWN(n,I_Si_S))
#define X_OH(n)     (UNKNOWN(n,I_OH))
#define N_CC(n)     (UNKNOWN(n,I_CC))
#define P_l(n)      (UNKNOWN(n,I_P_l))
#define X_HCO3(n)   (UNKNOWN(n,I_HCO3))
#define PSI(n)      (UNKNOWN(n,I_psi))
#define X_Na(n)     (UNKNOWN(n,I_Na))
#define X_K(n)      (UNKNOWN(n,I_K))


#define N_C(n)     		(f[(n)])
#define N_q(n)     		(f[(2+n)])
#define Mass(n)       (f[(4+n)])
#define N_Ca(n)    		(f[(6+n)])
#define N_k(n)     		(f[(8+n)])
#define N_Na(n)    		(f[(10+n)])
#define N_K(n)     		(f[(12+n)])
#define N_Si(n)    		(f[(16+n)])
#define XI(n)      		(f[(18+n)])
#define W_C        		(f[20])
#define W_q        		(f[21])
#define W_m        		(f[22])
#define W_Ca       		(f[23])
#define W_k        		(f[24])
#define W_Na       		(f[25])
#define W_K        		(f[26])
#define W_Si       		(f[27])
#define N_CH(n) 	    (f[(28+n)])

#define N_Cn(n)    		(f_n[(n)])
#define N_qn(n)    		(f_n[(2+n)])
#define Mass_n(n)     (f_n[(4+n)])
#define N_Can(n)   	  (f_n[(6+n)])
#define N_kn(n)    		(f_n[(8+n)])
#define N_Nan(n)   	  (f_n[(10+n)])
#define N_Kn(n)    		(f_n[(12+n)])
#define N_Sin(n)   		(f_n[(16+n)])
#define N_CHn(n) 	    (f_n[(28+n)])

#define KD_Ca      		  (va[(0)])
#define KD_OH      	    (va[(1)])
#define KD_H       		  (va[(2)])
#define KD_H2CO3   	    (va[(3)])
#define KD_HCO3    	    (va[(4)])
#define KD_CO3     	    (va[(5)])
#define KD_Na      		  (va[(6)])
#define KD_NaOH    	    (va[(7)])
#define KD_NaHCO3  	    (va[(8)])
#define KD_NaCO3   	    (va[(9)])
#define KD_m       		  (va[(10)])
#define KD_K       		  (va[(11)])
#define KD_KOH     	    (va[(12)])
#define KD_CaOH        	(va[(13)])
#define KD_CaHCO3  	    (va[(14)])
#define KD_CaCO3aq     	(va[(15)])
#define KD_CaOH2aq 	    (va[(16)])
#define KD_H3SiO4  	    (va[(17)])
#define KD_H2SiO4  	    (va[(18)])
#define KD_H4SiO4  	    (va[(19)])
#define KD_CaH2SiO4 	  (va[(20)])
#define KD_CaH3SiO4 	  (va[(21)])

#define KF_CO2     	    (va[(22)])
#define KF_Ca      		  (va[(23)])
#define KF_OH      		  (va[(24)])
#define KF_H       	  	(va[(25)])
#define KF_H2CO3   	    (va[(26)])
#define KF_HCO3    	    (va[(27)])
#define KF_CO3     	    (va[(28)])
#define KF_Na      		  (va[(29)])
#define KF_NaOH    	    (va[(30)])
#define KF_NaHCO3     	(va[(31)])
#define KF_NaCO3      	(va[(32)])
#define KF_K       	  	(va[(33)])
#define KF_KOH     	    (va[(34)])
#define KF_CaOH    	    (va[(35)])
#define KF_CaHCO3  	    (va[(36)])
#define KF_CaCO3aq 	    (va[(37)])
#define KF_CaOH2aq 	    (va[(38)])
#define KF_H3SiO4  	    (va[(39)])
#define KF_H2SiO4  	    (va[(40)])
#define KF_H4SiO4  	    (va[(41)])
#define KF_CaH2SiO4 	  (va[(42)])
#define KF_CaH3SiO4 	  (va[(43)])

#define Kpsi_Ca     		(va[(44)])
#define Kpsi_OH     		(va[(45)])
#define Kpsi_H      		(va[(46)])
#define Kpsi_HCO3   		(va[(47)])
#define Kpsi_CO3    		(va[(48)])
#define Kpsi_Na     		(va[(49)])
#define Kpsi_NaCO3  		(va[(50)])
#define Kpsi_q      		(va[(51)])
#define Kpsi_K      		(va[(52)])
#define Kpsi_CaOH   		(va[(53)])
#define Kpsi_CaHCO3 		(va[(54)])
#define Kpsi_H3SiO4 		(va[(55)])
#define Kpsi_H2SiO4 		(va[(56)])
#define Kpsi_CaH3SiO4	  (va[(57)])


/* les valences */
#define z_ca    	 (2.)
#define z_h     	 (1.)
#define z_oh    	 (-1.)
#define z_hco3  	 (-1.)
#define z_co3   	 (-2.)
#define z_na    	 (1.)
#define z_naco3 	 (-1.)
#define z_k      	 (1.)
#define z_caoh   	 (1.)
#define z_cahco3 	 (1.)
#define z_h3sio4 	 (-1.)
#define z_h2sio4 	 (-2.)
#define z_cah3sio4 (1.)


/* Stoechimetrie des C-S-H */
#define x_jen 		  (1.5)
#define y_jen 		  (0.9)
#define z_jen 		  (0.9)
#define zz_jen 		  (1.83)
#define x_tobII 		(1.5)
#define y_tobII 		(1.8)
#define z_tobII 		(1.2)
#define zz_tobII		(2.7)
#define x_tobI 		  (2)
#define y_tobI 		  (2.4)
#define z_tobI 		  (2)
#define zz_tobI 		(3.6)
#define x_sil			  (0.)
#define y_sil			  (1.)
#define z_sil			  (2.)


/* volumes molaires partiels des ions (dm3/mole) */
#define v_h       		(-5.5e-3)     /* (-5.50e-3)  d'apres TQN */
#define v_oh       		(23.50e-3)    /* (23.50e-3)  d'apres TQN */
#define v_h2o      		(18.e-3)
#define v_h2co3    	  (50.e-3)
#define v_hco3     		(50.e-3)
#define v_co3      		(-2.3e-3)     /* (-2.3e-3) d'apres 3w */
#define v_ca       		(-18.7e-3)    /* (-18.7e-3 d'apres */
#define v_na       		(22.4e-3)
#define v_naoh     		(22.4e-3)
#define v_nahco3  	  (22.4e-3) 
#define v_naco3    	  (22.4e-3) 
#define v_k        		(43.93e-3)
#define v_koh      		(27.44e-3)
#define v_caoh     		(26.20e-3)  /*a modifier*/
#define v_cahco3   	  (26.20e-3)  /*a modifier*/
#define v_caco3aq  	  (26.20e-3)  /*a modifier*/
#define v_caoh2aq  	  (26.20e-3)  /*a modifier*/
#define v_cah2sio4 	  (26.20e-3)  /*a modifier*/
#define v_cah3sio4 	  (26.20e-3)  /*a modifier*/
#define v_h3sio4   	  (50.e-3)    /*a modifier*/
#define v_h2sio4   	  (50.e-3)    /*a modifier*/
#define v_h4sio4   	  (50.e-3)    /*a modifier*/

/* volumes molaires solides (dm3/mole) */
#define V_CH	    	  (33.e-3)
#define V_CC	      	(37.e-3)
#define v_jen       	(73.6e-3)    /* Jennite (CSH sain) C1.5-S0.9-H1.83 */
#define v_tobII     	(98.6e-3)    /* Tobermorite  C1.5-S1.8-H2.7 */
#define v_tobI      	(131.47e-3)  /* Tobermorite  C2-S2.4-H3.6 */
#define v_sil       	(43.e-3)     /* Silice Amorphe S-H2 */
#define v_caco3csh    V_CC

/* Masses molaires (unite arbitraire = M_H) */
#define M_Ca    		 (40.1)
#define M_H2CO3  	   (62.)
#define M_HCO3   		 (61.)
#define M_CO3    		 (60.)
#define M_OH     		 (17.)
#define M_H      		 (1.)
#define M_H2O    		 (18.)
#define M_Na	 		   (23.)
#define M_NaOH	  	 (40.)
#define M_NaHCO3  	 (84.)
#define M_NaCO3	  	 (83.)
#define M_CO2     	 (44.)
#define M_CaOH2   	 (74.)
#define M_CaCO3   	 (100.)
#define M_CaCO3CSH 	 (100.)
#define M_CaO 		   (56.)
#define M_SiO2 		   (60.)
#define M_K      	 	 (39.)  
#define M_KOH    		 (56.)
#define M_CaOH    	 (57.)
#define M_CaHCO3  	 (101.)
#define M_CaCO3aq 	 (100.)
#define M_CaOH2aq  	 (74.)
#define M_CaH2SiO4 	 (134.)
#define M_CaH3SiO4   (135.)
#define M_H3SiO4  	 (95)
#define M_H2SiO4  	 (94)
#define M_H4SiO4  	 (96)

#define M_CSH      	 (x_csh*M_CaO + M_SiO2 + z_csh*M_H2O)

/* coefficients de diffusion moleculaire (dm2/s) */
#define d_oh       		(5.273e-7)    /* (5.273e-7) d'apres TQN */
#define d_h       		(9.310e-7)
#define d_ca       		(7.92e-8)
#define d_h2co3    	  (7.2e-8)
#define d_hco3     		(11.8e-8)
#define d_co3      		(9.55e-8)
#define d_na       		(1.33e-7)
#define d_naoh     		(1.33e-7)
#define d_nahco3  	  (1.33e-7)
#define d_naco3    	  (1.33e-7)
#define d_k        		(1.957e-7)
#define d_koh      		(1.957e-7)
#define d_caoh     		(7.92e-8) /*a modifier */
#define d_cahco3  	  (7.92e-8) /*a modifier */
#define d_caco3aq 	  (7.92e-8) /*a modifier */
#define d_caoh2aq  	  (7.92e-8) /*a modifier */
#define d_cah2sio4   	(7.92e-8) /*a modifier */
#define d_cah3sio4   	(7.92e-8) /*a modifier */
#define d_h3sio4   	  (11.8e-8) /*a modifier */
#define d_h2sio4   	  (11.8e-8) /*a modifier */
#define d_h4sio4     	(11.8e-8) /*a modifier */

#define d_co2      		(1.6e-3)

/* constantes d'equilibre (ref = 1 mole/L) */
#define k_e       		(1.e-14)                /* autoprotolyse de l'eau */
#define k_h        		(1.)                    /* cste de Henry */
#define k_co3      		(4.570881896148751e3)   /* Equilibre de HCO3  <-> CO3 */
#define k_ca       		(3.890451449942805e-9)  /* Equilibre de CaCO3 */
#define k_1        		(2.187761623949552e-8)  /* Equilibre de H2CO3 <-> HCO3 */
#define k_2        		(6.456542290346550e-6)  /* Equilibre de Ca(OH)2 */
#define k_naoh     		(6.60693448e-15)        /* Equilibre de Na <-> NaOH */
#define k_nahco3     	(0.5623413252)          /* Equilibre de Na <-> NaHCO3  */
#define k_naco3    	  (8.72971368e-10)        /* Equilibre de Na <-> NaCO3 */
#define k_caoh     		(1.65958691e-13)        /* Equilibre de Ca <-> CaOH */
#define k_cahco3   	  (12.76438809)           /* Equilibre de Ca <-> CaHCO3 */
#define k_caco3      	(7.852356346e-8)        /* Equilibre de Ca <-> CaCO3aqueux */
#define k_koh      		(3.4673685e-15)         /* Equilibre de K + H20 <-> KOH + H */
#define k_caoh2    	  (1.)                    /* Equilibre de Ca(OH)2 en solution */
#define k_h2sio4     	(2.13796209e+13)        /* Equilibre de H2SiO4-- + H+ <-> H3SiO4- */
#define k_h4sio4     	(1.380384265e+23)       /* Equilibre de H2SiO4 -- + 2H+ <-> H4SiO4 */
#define k_jen      		(1.412537545e-12)       /* Solubilité des CSH-jennite */
#define k_tobII    		(2.238721139e-14)       /* Solubilité des CSH-tobermoriteII */
#define k_tobI     		(6.30957344e-19)        /* Solubilité des CSH-tobermoriteII */
#define k_sil      		(1.936421964e-3)        /* Solubilité de la silice amorphe */
#define k_cah2sio4 	  (39810.71)              /* Ca <-> CaH2SiO4*/
#define k_cah3sio4 	  (15.84)                 /* Ca <-> CaH3SiO4*/

/* viscosite (Pa.s) */
#define mu_l       	(1.e-3)

/* constantes physiques */
#define FARADAY   (9.64846e7) /* Faraday (C/mole) */
#define RT        (2436.e3)   /* produit de R=8.3143 et T=293 (J/mole) */

/*puissance de la loi cinétique de dissolution de la portlandite*/
#define Xp   	   	 (0.5)
#define Xcinetique (1)

/*Loi de Davies pour la prise en compte de l'activite ionique*/
#define A_DAVIES  (0.5)
#define B_DAVIES  (0.24)

/* Curves */
#define SATURATION(p)    (Curve_ComputeValue(Element_GetCurve(el),p))
#define DSATURATION(p)   (Curve_ComputeDerivative(Element_GetCurve(el),p))
#define RELATIVEPERM(p)  (Curve_ComputeValue(Element_GetCurve(el) + 1,p))
#define X_CSH(q)         (Curve_ComputeValue(Element_GetCurve(el) + 2,q))
#define DX_CSH(q)        (Curve_ComputeDerivative(Element_GetCurve(el) + 2,q))
#define Z_CSH(q)         (Curve_ComputeValue(Element_GetCurve(el) + 3,q))
#define DZ_CSH(q)        (Curve_ComputeDerivative(Element_GetCurve(el) + 3,q))
#define V_CSH(q)         (Curve_ComputeValue(Element_GetCurve(el) + 4,q))
#define DV_CSH(q)        (Curve_ComputeDerivative(Element_GetCurve(el) + 4,q))
#define Q_SH(q)          (Curve_ComputeValue(Element_GetCurve(el) + 5,q))
#define DQ_SH(q)         (Curve_ComputeDerivative(Element_GetCurve(el) + 5,q))


/* Ion Activity Product of amorpheous silica */
#define IAP_SH(q)                 (k_sil*Q_SH(q))
#define DIAP_SH(q)                (k_sil*DQ_SH(q))

/* Fonctions */
#define GetProperty(a)   (Element_GetProperty(el)[pm(a)])
static int    pm(char *s) ;

static double* ComputeAqueousParameters(Element_t*,double**,int) ;

static void concentration1(double x_co2, double x_na_tot, double x_k_tot, double *pointeur_x_oh, double *pointeur_x_h3sio4, double *pointeur_x_na, double *pointeur_x_k) ;

static double force_ionique(double x_ca, double x_cahco3, double x_caoh, double x_k, double x_na, double x_h, double x_h3sio4, double x_naco3, double x_hco3, double x_oh, double x_h2sio4, double x_co3) ;

static void concentration_vers_activite(double x_ca, double x_cahco3, double x_caoh, double x_k, double x_na, double x_h, double x_h3sio4, double naco3, double x_hco3, double x_oh, double x_h2sio4, double x_co3,double *pointeur_x_ca,double *pointeur_x_cahco3, double *pointeur_x_caoh, double *pointeur_x_k, double *pointeur_x_na, double *pointeur_x_h, double *pointeur_x_h3sio4, double *pointeur_x_naco3, double *pointeur_x_hco3, double *pointeur_x_oh, double *pointeur_x_h2sio4, double *pointeur_x_co3) ;

static double dn1_caoh2sdt(double,double) ;

static void   transfert(Element_t*,double**,double*) ;
static void   flux(Element_t*,double**) ;


int SetModelProp(Model_t *model)
{
  Model_GetNbOfEquations(model) = NEQ ;
  
  Model_CopyNameOfEquation(model,E_C,"carbone") ;
  Model_CopyNameOfEquation(model,E_q,"charge") ;
  Model_CopyNameOfEquation(model,E_mass,"masse") ;
  Model_CopyNameOfEquation(model,E_Ca,"calcium") ;
  Model_CopyNameOfEquation(model,E_k,"E_k") ;
  Model_CopyNameOfEquation(model,E_el,"E_el") ;
  Model_CopyNameOfEquation(model,E_Na,"sodium") ;
  Model_CopyNameOfEquation(model,E_K,"potassium") ;
  Model_CopyNameOfEquation(model,E_Si,"silice") ;
  
  
#if (U_CO2 == LOG_RHO)
  Model_CopyNameOfUnknown(model,I_CO2,"logc_co2") ;
#else
  Model_CopyNameOfUnknown(model,I_CO2,"c_co2") ;
#endif
  Model_CopyNameOfUnknown(model,I_Si_S,"n_si_s") ;
  Model_CopyNameOfUnknown(model,I_OH,"c_oh") ;
  Model_CopyNameOfUnknown(model,I_P_l,"p_l") ;
  Model_CopyNameOfUnknown(model,I_CC,"c_caco3") ;
  Model_CopyNameOfUnknown(model,I_HCO3,"c_hco3") ;
  Model_CopyNameOfUnknown(model,I_psi,"psi") ;
  Model_CopyNameOfUnknown(model,I_Na,"c_na") ;
  Model_CopyNameOfUnknown(model,I_K,"c_k") ;
  
  return(0) ;
}


/* Parametres */
static double phii,k_int,a_1,a_2,c_2,n_ch0,x_na0,x_k0,n_si_s0 ;
static double p_g = 0. ;
static double var[36] ;

#define I_H          (9)
#define I_H2O        (10)
#define I_H2CO3      (11)
#define I_CO3        (12)
#define I_Ca         (13)
#define I_CaOH       (14)
#define I_CaHCO3     (15)
#define I_CaCO3aq    (16)
#define I_CaOH2aq    (17)
#define I_H2SiO4     (18)
#define I_H3SiO4     (19)
#define I_H4SiO4     (20)
#define I_CaH2SiO4   (21)
#define I_CaH3SiO4   (22)
#define I_NaOH       (23)
#define I_NaHCO3     (24)
#define I_NaCO3      (25)
#define I_KOH        (26)
#define I_Q_CH       (27)
#define I_RHO_L      (28)
#define I_N_Q        (29)
#define I_C_L        (30)
#define I_Ca_L       (31)
#define I_Na_L       (32)
#define I_K_L        (33)
#define I_Si_L       (35)

int pm(char *s)
{
  if(strcmp(s,"porosite") == 0)     return (0) ;
  else if(strcmp(s,"k_int") == 0)   return (1) ;
  else if(strcmp(s,"N_CaOH2") == 0) return (2) ;
  else if(strcmp(s,"N_Si") == 0)    return (4) ;
  else if(strcmp(s,"X_K") == 0)     return (5) ;
  else if(strcmp(s,"X_Na") == 0)    return (6) ;
  else if(strcmp(s,"A_1") == 0)     return (7) ;
  else if(strcmp(s,"A_2") == 0)     return (8) ;
  else if(strcmp(s,"C_2") == 0)     return (9) ;
  else if(strcmp(s,"R_CaOH2") == 0) return (10) ;
  else if(strcmp(s,"D") == 0) 	    return (11) ;
  else return(-1) ;
}

int ReadMatProp(Material_t *mat,DataFile_t *datafile)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  FILE *ficd = DataFile_GetFileStream(datafile) ;
  int  n_donnees = 12 ;
  /* Initialisation automatique */
  double h     = 5.6e-6 ;  /* (moles/dm2/s) */
  double R_0 ;
  double D  ;
  double n_ch0 ;
  double a_2 ;
  double c_2 ;

  Material_ScanProperties(mat,ficd,pm) ;
    
  R_0 = Material_GetProperty(mat)[pm("R_CaOH2")] ;
  D = Material_GetProperty(mat)[pm("D")] ;
  n_ch0 = Material_GetProperty(mat)[pm("N_CaOH2")] ; /* contenu molaire initial en CaOH2 */
  a_2 = 3*h/R_0*n_ch0*V_CH ; /* (dm/mole/s) these MT p 227 */
  c_2 = h*R_0/D ;     /* (1/dm) these MT p 228 */
  
  Material_GetProperty(mat)[pm("A_2")] = a_2 ;
  Material_GetProperty(mat)[pm("C_2")] = c_2 ;
  
  return(n_donnees) ;
}


int PrintModelChar(Model_t *model,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 7 equations:\n\
\t- la conservation de la masse de C      (c_co2)\n\
\t- la conservation de la charge          (psi)\n\
\t- la conservation de la masse totale    (p_l)\n\
\t- la conservation de la masse de Ca     (c_caco3)\n\
\t- 1 equation de cinetique               (c_hco3)\n\
\t- Electroneutralite                     (c_oh)\n\
\t- la conservation de la masse de Na     (c_na)\n\
\t- la conservation de la masse de K      (c_k)\n") ;

  printf("\n\
ATTENTION aux unites : \n\
\t longueur : dm !\n\
\t temps    : s !\n\
\t pression : Pa !\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.38   # La porosite\n") ;
  fprintf(ficd,"k_int = 5e-19     # Permeabilite intrinseque (dm2)\n") ;
  fprintf(ficd,"N_CaOH2 = 6.1     # Contenu initial en Ca(OH)2 (moles/L)\n") ;
  fprintf(ficd,"N_Si = 2.4       # contenu initial en CSH (moles/L)\n") ;
  fprintf(ficd,"T_csh = 1.8e-4    # k_h/T_csh = temps caracteristique de carbo des CSH (1/s)\n") ;
  fprintf(ficd,"A_1 = 150         # Coef de la cinetique 1 (dm/mole/s)\n") ;
  fprintf(ficd,"A_2 = 1e-2        # Coef de la cinetique 2 (dm/mole/s)\n") ;
  fprintf(ficd,"C_2 = 0.14e6      # Coef de la cinetique 2 (1/dm)\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l k_rl\n") ;  

  return(NEQ) ;
}


int DefineElementProp(Element_t *el,IntFcts_t *intfcts)
{
  Element_GetNbOfImplicitTerms(el) = NVI ;
  Element_GetNbOfExplicitTerms(el) = (Element_IsSubmanifold(el)) ? 0 : NVE ;
  Element_GetNbOfConstantTerms(el) = NV0 ;
  return(0) ;
}



int  ComputeLoads(Element_t *el,double t,double dt,Load_t *cg,double *r)
/* Residu du aux chargements (r) */
{
  int nn = Element_GetNbOfNodes(el) ;
  int    i ;

  FVM_ComputeSurfaceLoadResidu(el,cg,t,dt,r) ;
  for(i = 0 ; i < NEQ*nn ; i++) r[i] = -r[i] ;
  return(0) ;
}


int ComputeInitialState(Element_t *el)
/* Initialise les variables du systeme (f,va) */ 
{
  double *f = Element_GetImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Donnees
  */
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  n_si_s0  = GetProperty("N_Si") ;
  x_na0    = GetProperty("X_Na") ;
  x_k0     = GetProperty("X_K") ;
  
  {
    double x_na_tot 	= x_na0 ;
    double x_k_tot  	= x_k0 ;
    double x_co2   	= 2.884031968e-15 ;
    double x_na,x_k,x_h3sio4,x_oh ;
    
    concentration1(x_co2,x_na_tot,x_k_tot,&x_oh,&x_h3sio4,&x_na,&x_k) ;

    for(i = 0 ; i < 2 ; i++) {
#if (U_CO2 == LOG_RHO)
      UNKNOWN(i,I_CO2) = log(x_co2)/Ln10 ;
#else
      X_CO2(i)    = x_co2 ;
#endif
      X_Na(i)     = x_na ;
      X_K(i)      = x_k ;
      X_OH(i)     = x_oh ;
      X_HCO3(i)   = k_h*x_co2*x_oh/k_1 ;
    }
  }
  

  for(i = 0 ; i < 2 ; i++) {
    double p_l     	= P_l(i) ;
    double p_c     	= p_g - p_l ;
    double s_l     	= SATURATION(p_c) ;
    double s_g     	= 1 - s_l ;
    
    double *x = ComputeAqueousParameters(el,u,i) ;
    
    double x_oh    	  = x[I_OH] ;
    
    double x_hco3  	  = x[I_HCO3] ;
    double x_h2co3 	  = x[I_H2CO3] ;
    double x_co3   	  = x[I_CO3] ;
    
    double x_cahco3  	= x[I_CaHCO3] ;
    double x_caco3aq 	= x[I_CaCO3aq] ;
    
    double x_nahco3  	= x[I_NaHCO3] ;
    double x_naco3 	  = x[I_NaCO3] ;
    
    /* masse volumique liquide */
    double rho_l      = x[I_RHO_L] ;
    double rho_g      = M_CO2*x[I_CO2] ;

    /* solides */
    double q_ch       = x[I_Q_CH] ;
    double x_csh      = X_CSH(q_ch) ;
    double z_csh      = Z_CSH(q_ch) ;
    double v_csh      = V_CSH(q_ch) ;

    double n_cc       = N_CC(i) ;
    double n_ch       = n_ch0 ;
    double n_si_s     = N_Si_S(i) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
    double n_c_s      = n_cc ;

    /* porosite */
    double v_csh0   = V_CSH(1) ;
    double phi = phii + V_CH*(n_ch0 - n_ch) + v_csh0*n_si_s0 - v_csh*n_si_s + V_CC*(-n_cc) ;
    double phi_l = phi*s_l ;
    double phi_g = phi*s_g ;

    /* Molar contents */
    double n_c_g  = phi_g*x[I_CO2] ;
    double n_c_l  = phi_l*x[I_C_L] ;
    double n_ca_l = phi_l*x[I_Ca_L] ;
    double n_na_l = phi_l*x[I_Na_L] ;
    double n_k_l  = phi_l*x[I_K_L] ;
    double n_si_l = phi_l*x[I_Si_L] ;
    
    /* contenus atomiques */
    N_C(i)  = n_c_l  + n_c_s  + n_c_g ;
    N_Ca(i) = n_ca_l + n_ca_s ;
    N_Na(i) = n_na_l ; 
    N_K(i)  = n_k_l  ;
    N_Si(i) = n_si_l + n_si_s ;
   
    /* masse totale */
    Mass(i) = phi_g*rho_g + phi_l*rho_l + M_CaOH2*n_ch + M_CaCO3*n_cc + M_CSH*n_si_s ;

    /* cinetique */
    N_k(i)  = phi_l*(x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq) + n_cc ;
    XI(i)   = phi*pow(s_l,Xcinetique)*a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

    /* densite de charge */
    N_q(i)  = x[I_N_Q] ;


    /* contenus solides */
    N_CH(i) = n_ch ;
  }

  /* Coefficient de transfert */
  transfert(el,u,f) ;

  /* Flux */
  flux(el,u) ;
  return(0) ;
}


int  ComputeExplicitTerms(Element_t *el,double t)
/* Thermes explicites (va)  */
{
  double *f = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetPreviousNodalUnknown(el,i) ;
  }
  
  /*
    Coefficients de transfert
  */
  transfert(el,u,f) ;

  return(0) ;
}


int  ComputeImplicitTerms(Element_t *el,double t,double dt)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double *u[MAX_NOEUDS] ;
  int i ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetCurrentNodalUnknown(el,i) ;
  }
  
  
  /*
    Donnees
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  n_si_s0  = GetProperty("N_Si") ;
  
  
  /* Contenus molaires */
  for(i = 0 ; i < nn ; i++) {
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    double s_l     = SATURATION(p_c) ;
    double s_g     = 1 - s_l ;
    
    /* Aqueous components */
    double *x = ComputeAqueousParameters(el,u,i) ;
    
    double x_oh    	  = x[I_OH] ;
    double x_hco3  	  = x[I_HCO3] ;
    double x_h2co3 	  = x[I_H2CO3] ;
    double x_co3   	  = x[I_CO3] ;
    double x_nahco3  	= x[I_NaHCO3] ;
    double x_naco3 	  = x[I_NaCO3] ;
    double x_cahco3  	= x[I_CaHCO3] ;
    double x_caco3aq 	= x[I_CaCO3aq] ;
    
    /* masse volumique liquide */
    double rho_l      = x[I_RHO_L] ;
    double rho_g      = M_CO2*x[I_CO2] ;


    /* Solids */
    double q_ch       = x[I_Q_CH] ;
    double x_csh      = X_CSH(q_ch) ;
    double z_csh      = Z_CSH(q_ch) ;
    double v_csh      = V_CSH(q_ch) ;

    double n_cc       = N_CC(i) ;
    double n_si_s     = N_Si_S(i) ;

    /* double av = n_cc/n_ch0 ; */

    double av = 1 - N_CHn(i)/n_ch0 ;
    double dn1sdt = a_2*dn1_caoh2sdt(av,c_2) ;
    /*cinétique de dissolution de la portlandite en loi log*/
    double dn_chsdt = dn1sdt*log((k_ca*x_oh/(k_co3*k_2*x_hco3))) ;    
    /*cinétique de dissolution de la portlandite en loi puissance*/ 
    /* double dn_chsdt = - dn1sdt*(pow(fabs((k_ca*x_oh/(k_co3*k_2*x_hco3)-1)),Xp)) ; */
      
    double n_ch    = MAX(N_CHn(i) + dt*dn_chsdt , 0.) ;
    
    double n_ca_s  = n_ch + n_cc + x_csh*n_si_s ;
    
    double n_c_s   = n_cc ;
    

    /* porosite */
    double v_csh0   = V_CSH(1) ;
    double phi = phii + V_CH*(n_ch0 - n_ch) + v_csh0*n_si_s0 - v_csh*n_si_s + V_CC*(-n_cc) ;
    double phi_l = phi*s_l ;
    double phi_g = phi*s_g ;

    /* Molar contents */
    double n_c_g  = phi_g*x[I_CO2] ;
    double n_c_l  = phi_l*x[I_C_L] ;
    double n_ca_l = phi_l*x[I_Ca_L] ;
    double n_na_l = phi_l*x[I_Na_L] ;
    double n_k_l  = phi_l*x[I_K_L] ;
    double n_si_l = phi_l*x[I_Si_L] ;
    
    /* contenus atomiques */
    N_C(i)  = n_c_l  + n_c_s  + n_c_g ;
    N_Ca(i) = n_ca_l + n_ca_s ;
    N_Na(i) = n_na_l ; 
    N_K(i)  = n_k_l  ;
    N_Si(i) = n_si_l + n_si_s ;
   
    /* masse totale */
    Mass(i) = phi_g*rho_g + phi_l*rho_l + M_CaOH2*n_ch + M_CaCO3*n_cc + M_CSH*n_si_s ;

    /* cinetique */
    N_k(i)  = phi_l*(x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq) + n_cc ;
    XI(i)   = phi*pow(s_l,Xcinetique)*a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;

    /* densite de charge */
    N_q(i)  = x[I_N_Q] ;

    /* contenus solides */
    N_CH(i) = n_ch ;
    
    {
      double x_co2      = x[I_CO2] ;
      double x_oh    	  = x[I_OH] ;
      double x_hco3  	  = x[I_HCO3] ;
      double x_ca    	  = x[I_Ca] ;
      double x_na    	  = x[I_Na] ;
      double x_k     	  = x[I_K] ;
      double x_naoh    	= x[I_NaOH] ;
      double x_nahco3  	= x[I_NaHCO3] ;
      double x_naco3 	  = x[I_NaCO3] ;
      double x_h2o      = x[I_H2O] ;
      if(x_co2 < 0. || x_oh <= 0. || x_h2o <= 0. || x_hco3 <= 0. || x_na < 0. || x_k < 0. || x_ca < 0. || n_si_s < 0.) {
        double x = Element_GetNodeCoordinate(el,i)[0] ;
        printf("\n") ;
        printf("en x     = %e\n",x) ;
        printf("x_co2    = %e\n",x_co2) ;
        printf("x_oh     = %e\n",x_oh) ;
        printf("x_h2o    = %e\n",x_h2o) ;
        printf("x_hco3   = %e\n",x_hco3) ;
        printf("n_cc     = %e\n",n_cc) ;
        printf("x_na     = %e\n",x_na) ;
        printf("x_k      = %e\n",x_k) ;
        printf("x_naoh   = %e\n",x_naoh) ;
        printf("x_nahco3 = %e\n",x_nahco3) ;
        printf("x_naco3  = %e\n",x_naco3) ;
        return(1) ;
      }
    }
  }
  
  if(Element_IsSubmanifold(el)) return(0) ;
  /* Flux */
  flux(el,u) ;

  return(0) ;
}


int  ComputeMatrix(Element_t *el,double t,double dt,double *k)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
  Symmetry_t sym = Element_GetSymmetry(el) ;
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  double *u[MAX_NOEUDS] ;
  int nn = Element_GetNbOfNodes(el) ;

  double dx ;
  double volume[2],surf ;
  int    i ;
  double c[2] ;
  
  double Dx_h2co3SDx_co2[2] ;
  
  double Dx_hSDx_oh[2] ;
  double Dx_co3SDx_oh[2] ;
  double Dx_caSDx_oh[2] ;
  double Dx_naohSDx_oh[2] ;
  double Dx_naco3SDx_oh[2] ;
  double Dx_kohSDx_oh[2] ;
  double Dx_cahco3SDx_oh[2] ;
  double Dx_caoh2aqSDx_oh[2] ;
  double Dx_h2sio4SDx_oh[2] ;
  double Dx_h3sio4SDx_oh[2] ;
  double Dx_h4sio4SDx_oh[2] ;
  double Dx_cah2sio4SDx_oh[2] ;
  double Dx_cah3sio4SDx_oh[2] ;
  
  double Dx_co3SDx_hco3[2] ;
  double Dx_caSDx_hco3[2] ;
  double Dx_nahco3SDx_hco3[2] ;
  double Dx_naco3SDx_hco3[2] ;
  double Dx_caohSDx_hco3[2] ;
  double Dx_caoh2aqSDx_hco3[2] ;
  double Dx_h2sio4SDx_hco3[2] ;
  double Dx_h3sio4SDx_hco3[2] ;
  double Dx_h4sio4SDx_hco3[2] ;
  double Dx_cah2sio4SDx_hco3[2] ;
  double Dx_cah3sio4SDx_hco3[2] ;
  
  double Dx_naohSDx_na[2] ;
  double Dx_nahco3SDx_na[2] ;
  double Dx_naco3SDx_na[2] ;
  
  double Dx_kohSDx_k[2] ;

  
  /*
    Initialisation 
  */
  for(i = 0 ; i < 4*NEQ*NEQ ; i++) k[i] = 0. ;

  
  if(Element_IsSubmanifold(el)) return(0) ;
  
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Donnees 
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  n_si_s0  = GetProperty("N_Si") ;

  /*
    CALCUL DE volume ET DE surf 
  */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double xm = (x1 + x0)*0.5 ;
    dx = x1 - x0 ;

    for(i=0;i<nn;i++) {
      double x = Element_GetNodeCoordinate(el,i)[0] ;
      volume[i] = fabs(dx)*0.5 ; 
      if(sym == AXIS) volume[i] *= M_PI*(x + xm) ; 
    }
    if(sym == AXIS) surf = 2*M_PI*xm ; else surf = 1 ;
  }
  
  /*
    termes d'accumulation
  */
  for(i=0;i<2;i++) {
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    double s_l     = SATURATION(p_c) ;
    double s_g     = 1 - s_l ;
    
    /* Aqueous components */
    double *x = ComputeAqueousParameters(el,u,i) ;
    
    double x_co2      = x[I_CO2] ;
    double x_oh    	  = x[I_OH] ;
    double x_hco3  	  = x[I_HCO3] ;
    double x_h2co3 	  = x[I_H2CO3] ;
    double x_co3   	  = x[I_CO3] ;
    double x_h     	  = x[I_H] ;
    double x_ca    	  = x[I_Ca] ;
    double x_caoh  	  = x[I_CaOH] ;
    double x_cahco3  	= x[I_CaHCO3] ;
    double x_caco3aq 	= x[I_CaCO3aq] ;
    double x_caoh2aq 	= x[I_CaOH2aq] ;
    double x_h4sio4   = x[I_H4SiO4] ;
    double x_h3sio4   = x[I_H3SiO4] ;
    double x_h2sio4   = x[I_H2SiO4] ;
    double x_cah2sio4 = x[I_CaH2SiO4];
    double x_cah3sio4 = x[I_CaH3SiO4] ;
    double x_na    	  = x[I_Na] ;
    double x_nahco3  	= x[I_NaHCO3] ;
    double x_naco3 	  = x[I_NaCO3] ;
    double x_k     	  = x[I_K] ;
    
    /* Atom Concentrations in the liquid phase */
    double x_c_l  = x[I_C_L] ;
    double x_ca_l = x[I_Ca_L] ;
    double x_na_l = x[I_Na_L] ;
    double x_k_l  = x[I_K_L] ;
    double x_si_l = x[I_Si_L] ;
    
    /* masse volumique liquide */
    double q_ch       = x[I_Q_CH] ;
    double rho_l      = x[I_RHO_L] ;
    double rho_g      = M_CO2*x_co2 ;

    /* solides */
    double x_csh      = X_CSH(q_ch) ;
    double z_csh      = Z_CSH(q_ch) ;
    double v_csh      = V_CSH(q_ch) ;
    
    double n_cc       = N_CC(i) ;
    double n_ch       = N_CH(i) ;
    double n_si_s     = N_Si_S(i) ;

    double xi      = a_1*(x_oh*x_h2co3 - k_1*x_hco3) ;
    
    double av   = 1 - N_CHn(i)/n_ch0 ;
    double dn1  = (n_ch > 0) ? dt*a_2*dn1_caoh2sdt(av,c_2) : 0 ;
    

    /* porosite */
    double v_csh0   = V_CSH(1) ;
    double phi = phii + V_CH*(n_ch0 - n_ch) + v_csh0*n_si_s0 - v_csh*n_si_s + V_CC*(-n_cc) ;
    double phi_l = phi*s_l ;
    double phi_g = phi*s_g ;



    /* derivees ... */
    /* ... par rapport a p_l */
    double ds_lsdp_c = DSATURATION(p_c) ;
    double ds_lsdp_l = -ds_lsdp_c ;

    /* ... par rapport a x_co2 */
    double dx_h2co3sdx_co2    = k_h ;
    
    double dx_h2osdx_co2      = -(v_h2co3*dx_h2co3sdx_co2)/v_h2o ;

    double drho_lsdx_co2      = M_H2CO3*dx_h2co3sdx_co2 + M_H2O*dx_h2osdx_co2 ;
    
    double dxisdx_co2         =  a_1*x_oh*dx_h2co3sdx_co2 ;
    
    double dx_c_lsdx_co2  = dx_h2co3sdx_co2 ;
    
    /* ... par rapport a x_oh */
    /* liquid */
    double dx_hsdx_oh        = - x_h/x_oh ;
    
    double dx_co3sdx_oh      =   k_co3*x_hco3 ;
    
    double dx_casdx_oh       = - x_ca/x_oh ;
    double dx_cahco3sdx_oh   = - x_cahco3/x_oh ;
    
    double dq_chsdx_oh       =   q_ch/x_oh ;
    double dx_h4sio4sdx_oh   =   DIAP_SH(q_ch)*dq_chsdx_oh ;
    double dx_h3sio4sdx_oh   =   x_h3sio4*(dx_h4sio4sdx_oh/x_h4sio4 + 1/x_oh) ;
    double dx_h2sio4sdx_oh   =   x_h2sio4*(1/x_oh + dx_h3sio4sdx_oh/x_h3sio4) ;
    double dx_caoh2aqsdx_oh  =   x_caoh2aq/x_oh ;
    double dx_cah2sio4sdx_oh =   x_cah2sio4*(dx_h2sio4sdx_oh/x_h2sio4 + dx_casdx_oh/x_ca) ;
    double dx_cah3sio4sdx_oh =   x_cah3sio4*(dx_h3sio4sdx_oh/x_h3sio4 + dx_casdx_oh/x_ca) ;
    
    double dx_naohsdx_oh     =   k_naoh/k_e*x_na ;
    double dx_naco3sdx_oh    =   k_naco3/k_e*x_na*x_hco3 ;
    double dx_kohsdx_oh      =   x_k*k_koh/k_e ;
    
    double dx_h2osdx_oh      = -(v_co3*dx_co3sdx_oh + v_oh + v_h*dx_hsdx_oh + v_ca*dx_casdx_oh + v_naoh*dx_naohsdx_oh + v_naco3*dx_naco3sdx_oh + v_koh*dx_kohsdx_oh + v_cahco3*dx_cahco3sdx_oh + v_caoh2aq*dx_caoh2aqsdx_oh + v_h2sio4*dx_h2sio4sdx_oh + v_h3sio4*dx_h3sio4sdx_oh + v_h4sio4*dx_h4sio4sdx_oh + v_cah2sio4*dx_cah2sio4sdx_oh + v_cah3sio4*dx_cah3sio4sdx_oh )/v_h2o ;
    
    double drho_lsdx_oh      = M_H*dx_hsdx_oh + M_OH + M_H2O*dx_h2osdx_oh + M_CO3*dx_co3sdx_oh + M_Ca*dx_casdx_oh  + M_NaOH*dx_naohsdx_oh + M_NaCO3*dx_naco3sdx_oh + M_KOH*dx_kohsdx_oh + M_CaHCO3*dx_cahco3sdx_oh + M_CaOH2aq*dx_caoh2aqsdx_oh + M_H2SiO4*dx_h2sio4sdx_oh + M_H3SiO4*dx_h3sio4sdx_oh + M_H4SiO4*dx_h4sio4sdx_oh + M_CaH2SiO4*dx_cah2sio4sdx_oh + M_CaH3SiO4*dx_cah3sio4sdx_oh ;
    
    double dxisdx_oh          =  a_1*x_h2co3 ;
         
    double dx_c_lsdx_oh  = dx_co3sdx_oh + dx_cahco3sdx_oh + dx_naco3sdx_oh ;
    double dx_ca_lsdx_oh = dx_casdx_oh + dx_cahco3sdx_oh + dx_caoh2aqsdx_oh + dx_cah2sio4sdx_oh + dx_cah3sio4sdx_oh ;
    double dx_na_lsdx_oh = dx_naohsdx_oh + dx_naco3sdx_oh ;
    double dx_k_lsdx_oh  = dx_kohsdx_oh ;
    double dx_si_lsdx_oh = dx_h2sio4sdx_oh + dx_h3sio4sdx_oh + dx_h4sio4sdx_oh + dx_cah2sio4sdx_oh + dx_cah3sio4sdx_oh ;
    double dx_qsdx_oh    = z_h*dx_hsdx_oh + z_oh \
       + z_co3*dx_co3sdx_oh \
       + z_ca*dx_casdx_oh + z_cahco3*dx_cahco3sdx_oh \
       + z_h3sio4*dx_h3sio4sdx_oh + z_h2sio4*dx_h2sio4sdx_oh \
       + z_cah3sio4*dx_cah3sio4sdx_oh \
       + z_naco3*dx_naco3sdx_oh ;

    /* solid */
    double dx_cshsdx_oh      = DX_CSH(q_ch)*dq_chsdx_oh ;
    double dz_cshsdx_oh      = DZ_CSH(q_ch)*dq_chsdx_oh ;
    double dv_cshsdx_oh      = DV_CSH(q_ch)*dq_chsdx_oh ;
    double dM_CSHsdx_oh      = dx_cshsdx_oh*M_CaO + dz_cshsdx_oh*M_H2O ;
    
    double dn_chsdx_oh        =  dn1/x_oh ;
    double dn_ca_ssdx_oh      = dn_chsdx_oh + dx_cshsdx_oh*n_si_s ;
    
    double dphisdx_oh         = - V_CH*dn_chsdx_oh - dv_cshsdx_oh*n_si_s  ;
    
    /* ... par rapport a x_hco3 */
    /* liquid */
    double dx_co3sdx_hco3      = k_co3*x_oh ;
    
    double dx_casdx_hco3       = - x_ca/x_hco3 ;
    double dx_caohsdx_hco3     = - x_caoh/x_hco3 ;
    double dx_caoh2aqsdx_hco3  = - x_caoh2aq/x_hco3 ;
    
    double dq_chsdx_hco3       = - q_ch/x_hco3 ;
    double dx_h4sio4sdx_hco3   =   DIAP_SH(q_ch)*dq_chsdx_hco3 ;
    double dx_h3sio4sdx_hco3   =   x_h3sio4*(dx_h4sio4sdx_hco3/x_h4sio4) ;
    double dx_h2sio4sdx_hco3   =   x_h2sio4*(dx_h3sio4sdx_hco3/x_h3sio4) ;
    double dx_cah2sio4sdx_hco3 =   x_cah2sio4*(dx_h2sio4sdx_hco3/x_h2sio4 + dx_casdx_hco3/x_ca) ;
    double dx_cah3sio4sdx_hco3 =   x_cah3sio4*(dx_h3sio4sdx_hco3/x_h3sio4 + dx_casdx_hco3/x_ca) ;
    
    double dx_nahco3sdx_hco3   =   k_nahco3*x_na  ;
    double dx_naco3sdx_hco3    =   k_naco3/k_e*x_na*x_oh ;
    
    double dx_h2osdx_hco3     = -(v_hco3 + v_co3*dx_co3sdx_hco3 + v_ca*dx_casdx_hco3 + v_nahco3*dx_nahco3sdx_hco3 + v_naco3*dx_naco3sdx_hco3 + v_caoh*dx_caohsdx_hco3 + v_caoh2aq*dx_caoh2aqsdx_hco3 + v_h2sio4*dx_h2sio4sdx_hco3 + v_h3sio4*dx_h3sio4sdx_hco3 + v_h4sio4*dx_h4sio4sdx_hco3 + v_cah2sio4*dx_cah2sio4sdx_hco3 + v_cah3sio4*dx_cah3sio4sdx_hco3)/v_h2o ;
    
    double drho_lsdx_hco3     = M_H2O*dx_h2osdx_hco3 + M_HCO3 + M_CO3*dx_co3sdx_hco3 + M_Ca*dx_casdx_hco3 + M_NaHCO3*dx_nahco3sdx_hco3 + M_NaCO3*dx_naco3sdx_hco3 + M_CaOH*dx_caohsdx_hco3 + M_CaOH2aq*dx_caoh2aqsdx_hco3 + M_H2SiO4*dx_h2sio4sdx_hco3 + M_H3SiO4*dx_h3sio4sdx_hco3 + M_H4SiO4*dx_h4sio4sdx_hco3 + M_CaH2SiO4*dx_cah2sio4sdx_hco3 + M_CaH3SiO4*dx_cah3sio4sdx_hco3 ;
    
    double dxisdx_hco3        = -a_1*k_1 ;
    
    double dx_c_lsdx_hco3  = 1 + dx_co3sdx_hco3 + dx_nahco3sdx_hco3 + dx_naco3sdx_hco3 ;
    double dx_ca_lsdx_hco3 = dx_casdx_hco3 + dx_caohsdx_hco3 + dx_caoh2aqsdx_hco3 + dx_cah2sio4sdx_hco3 + dx_cah3sio4sdx_hco3 ;
    double dx_na_lsdx_hco3 = dx_nahco3sdx_hco3 + dx_naco3sdx_hco3 ;
    double dx_si_lsdx_hco3 = dx_h2sio4sdx_hco3 + dx_h3sio4sdx_hco3 + dx_h4sio4sdx_hco3 + dx_cah2sio4sdx_hco3 + dx_cah3sio4sdx_hco3 ;
    double dx_qsdx_hco3    = z_hco3 + z_co3*dx_co3sdx_hco3 \
       + z_ca*dx_casdx_hco3 + z_caoh*dx_caohsdx_hco3 \
       + z_h3sio4*dx_h3sio4sdx_hco3 + z_h2sio4*dx_h2sio4sdx_hco3 \
       + z_cah3sio4*dx_cah3sio4sdx_hco3 \
       + z_naco3*dx_naco3sdx_hco3 ;
    
    /* solid */
    double dx_cshsdx_hco3      = DX_CSH(q_ch)*dq_chsdx_hco3 ;
    double dz_cshsdx_hco3      = DZ_CSH(q_ch)*dq_chsdx_hco3 ;
    double dv_cshsdx_hco3      = DV_CSH(q_ch)*dq_chsdx_hco3 ;
    double dM_CSHsdx_hco3      = dx_cshsdx_hco3*M_CaO + dz_cshsdx_hco3*M_H2O ;
    
    double dn_chsdx_hco3      = -dn1/x_hco3 ;
    double dn_ca_ssdx_hco3    = dn_chsdx_hco3 + dx_cshsdx_hco3*n_si_s ;
    
    double dphisdx_hco3       = - V_CH*dn_chsdx_hco3 - dv_cshsdx_hco3*n_si_s ;
    
    /* ... par rapport a x_na */
    double dx_naohsdx_na     =  k_naoh/k_e*x_oh ;
    double dx_nahco3sdx_na   =  k_nahco3*x_hco3 ;
    double dx_naco3sdx_na    =  k_naco3/k_e*x_oh*x_hco3 ;
    
    double dx_h2osdx_na       = -(v_na + v_naoh*dx_naohsdx_na + v_nahco3*dx_nahco3sdx_na + v_naco3*dx_naco3sdx_na)/v_h2o ;
    
    double drho_lsdx_na       = M_H2O*dx_h2osdx_na + M_Na + M_NaOH*dx_naohsdx_na + M_NaHCO3*dx_nahco3sdx_na + M_NaCO3*dx_naco3sdx_na ;
    
    double dx_c_lsdx_na  = dx_nahco3sdx_na + dx_naco3sdx_na ;
    double dx_na_lsdx_na = 1 + dx_naohsdx_na + dx_nahco3sdx_na + dx_naco3sdx_na ;
    double dx_qsdx_na    = z_na + z_naco3*dx_naco3sdx_na ;
    
    /* ... par rapport a x_k */
    double dx_kohsdx_k        = x_oh*k_koh/k_e ;
    
    double dx_h2osdx_k        = -(v_k + v_koh*dx_kohsdx_k)/v_h2o ;
    
    double drho_lsdx_k        = M_H2O*dx_h2osdx_k + M_K + M_KOH*dx_kohsdx_k ;

    double dx_k_lsdx_k  = 1 + dx_kohsdx_k ;
    double dx_qsdx_k    = z_k ;

    
    /* ... par rapport a n_cc */
    double dn_chsdn_cc     =  0. ;
    double dn_ca_ssdn_cc   = 1 + dn_chsdn_cc ;
      
    double dphisdn_cc      = - V_CC - V_CH*dn_chsdn_cc ;
    
    
    /* ... par rapport a n_si_s */
    double dn_ca_ssdn_si_s   = x_csh ;
    
    double dphisdn_si_s      = - v_csh ;
      
      
      
    int j = i*NEQ ;
    /*
      Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
    */
    K(E_C+j,I_P_l+j)   += volume[i]*(phi*ds_lsdp_l*(-x_co2 + x_c_l)) ;
    K(E_C+j,I_CO2+j)   += volume[i]*(phi_g + phi_l*dx_c_lsdx_co2) ;
    K(E_C+j,I_OH+j)    += volume[i]*(phi_l*dx_c_lsdx_oh  + dphisdx_oh*(s_g*x_co2 + s_l*x_c_l)) ;
    K(E_C+j,I_HCO3+j)  += volume[i]*(phi_l*dx_c_lsdx_hco3 + dphisdx_hco3*(s_g*x_co2 + s_l*x_c_l)) ;
    K(E_C+j,I_CC+j)    += volume[i]*(1 + dphisdn_cc*(s_g*x_co2 + s_l*x_c_l)) ;
    K(E_C+j,I_Na+j)    += volume[i]*(phi_l*dx_c_lsdx_na) ;
    K(E_C+j,I_Si_S+j)  += volume[i]*(dphisdn_si_s*(s_g*x_co2 + s_l*x_c_l)) ;

    /*
      Conservation de la charge  : div(w_q) = 0
    */

    /*
      Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
    */
    K(E_mass+j,I_P_l+j)   += volume[i]*phi*ds_lsdp_l*(-rho_g + rho_l) ;
    K(E_mass+j,I_CO2+j)   += volume[i]*phi*(s_g*M_CO2 + s_l*drho_lsdx_co2) ;
    K(E_mass+j,I_OH+j)    += volume[i]*(phi_l*drho_lsdx_oh + dphisdx_oh*(s_l*rho_l + s_g*rho_g) + M_CaOH2*dn_chsdx_oh + dM_CSHsdx_oh*n_si_s) ;
    K(E_mass+j,I_HCO3+j)  += volume[i]*(phi_l*drho_lsdx_hco3 + dphisdx_hco3*(s_l*rho_l + s_g*rho_g) + M_CaOH2*dn_chsdx_hco3  + dM_CSHsdx_hco3*n_si_s) ;
    K(E_mass+j,I_CC+j)    += volume[i]*(dphisdn_cc*(s_l*rho_l + s_g*rho_g) + M_CaOH2*dn_chsdn_cc + M_CaCO3) ;
    K(E_mass+j,I_Na+j)    += volume[i]*(phi_l*drho_lsdx_na) ;
    K(E_mass+j,I_K+j)     += volume[i]*(phi_l*drho_lsdx_k) ;
    K(E_mass+j,I_Si_S+j)  += volume[i]*(dphisdn_si_s*(s_l*rho_l + s_g*rho_g) + M_CSH) ;


    /*
      Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
    */
    K(E_Ca+j,I_P_l+j)    += volume[i]*(phi*ds_lsdp_l*x_ca_l) ;
    K(E_Ca+j,I_OH+j)     += volume[i]*(phi_l*dx_ca_lsdx_oh + dn_ca_ssdx_oh + dphisdx_oh*s_l*x_ca_l) ;
    K(E_Ca+j,I_HCO3+j)   += volume[i]*(phi_l*dx_ca_lsdx_hco3 + dn_ca_ssdx_hco3 + dphisdx_hco3*s_l*x_ca_l) ;
    K(E_Ca+j,I_CC+j)     += volume[i]*(dn_ca_ssdn_cc + dphisdn_cc*s_l*x_ca_l) ;
    K(E_Ca+j,I_Si_S+j)   += volume[i]*(dn_ca_ssdn_si_s + dphisdn_si_s*s_l*x_ca_l) ;


    /*
      Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
    */
    K(E_k+j,I_P_l+j)   += volume[i]*(phi*ds_lsdp_l*(x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 )- dt*phi*Xcinetique*ds_lsdp_l*xi*pow(s_l,Xcinetique-1)) ;
    K(E_k+j,I_CO2+j)   += volume[i]*phi*pow(s_l,Xcinetique)*(-dt*dxisdx_co2) ;
    K(E_k+j,I_OH+j)    += volume[i]*(phi_l*(dx_co3sdx_oh + dx_naco3sdx_oh + dx_cahco3sdx_oh ) + dphisdx_oh*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3) - dphisdx_oh*pow(s_l,Xcinetique)*dt*xi - phi*pow(s_l,Xcinetique)*dt*dxisdx_oh );
    K(E_k+j,I_HCO3+j)  += volume[i]*(phi_l*(1 + dx_co3sdx_hco3 + dx_nahco3sdx_hco3 + dx_naco3sdx_hco3  ) + dphisdx_hco3*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 ) - dphisdx_hco3*pow(s_l,Xcinetique)*dt*xi - phi*pow(s_l,Xcinetique)*dt*dxisdx_hco3 );
    K(E_k+j,I_CC+j)    += volume[i]*(1 + dphisdn_cc*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 ) - dphisdn_cc*pow(s_l,Xcinetique)*dt*xi) ;
    K(E_k+j,I_Na+j)    += volume[i]*phi_l*(dx_nahco3sdx_na + dx_naco3sdx_na) ;
    K(E_k+j,I_Si_S+j)  += volume[i]*(dphisdn_si_s*s_l*(x_co3 + x_hco3 + x_nahco3 + x_naco3 + x_caco3aq + x_cahco3 ) - dphisdn_si_s*pow(s_l,Xcinetique)*dt*xi) ;
    
    /*
      Electroneutralite : q = 0
    */
    K(E_el+j,I_OH+j)    += volume[i]*(dx_qsdx_oh) ;
    K(E_el+j,I_HCO3+j)  += volume[i]*(dx_qsdx_hco3) ;
    K(E_el+j,I_Na+j)    += volume[i]*(dx_qsdx_na) ;
    K(E_el+j,I_K+j)     += volume[i]*(dx_qsdx_k) ;


    /*
      Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
    */
    K(E_Na+j,I_P_l+j)   += volume[i]*(phi*ds_lsdp_l*x_na_l) ;
    K(E_Na+j,I_OH+j)    += volume[i]*(phi_l*dx_na_lsdx_oh + dphisdx_oh*s_l*x_na_l) ;
    K(E_Na+j,I_HCO3+j)  += volume[i]*(phi_l*dx_na_lsdx_hco3 + dphisdx_hco3*s_l*x_na_l) ;
    K(E_Na+j,I_CC+j)    += volume[i]*(dphisdn_cc*s_l*x_na_l);
    K(E_Na+j,I_Na+j)    += volume[i]*(phi_l*dx_na_lsdx_na) ;
    K(E_Na+j,I_Si_S+j)  += volume[i]*(dphisdn_si_s*s_l*x_na_l) ;
    
    /*
      Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
    K(E_K+j,I_P_l+j)    += volume[i]*(phi*ds_lsdp_l*x_k_l) ;
    K(E_K+j,I_OH+j)     += volume[i]*(phi_l*dx_k_lsdx_oh + dphisdx_oh*s_l*x_k_l) ;
    K(E_K+j,I_HCO3+j)   += volume[i]*(dphisdx_hco3*s_l*x_k_l) ;
    K(E_K+j,I_CC+j)     += volume[i]*(dphisdn_cc*s_l*x_k_l) ;
    K(E_K+j,I_K+j)      += volume[i]*(phi_l*dx_k_lsdx_k) ;
    K(E_K+j,I_Si_S+j)   += volume[i]*(dphisdn_si_s*s_l*x_k_l) ;

    /*
      Conservation de Si (silice) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
    */
    K(E_Si+j,I_P_l+j)      += volume[i]*(phi*ds_lsdp_l*x_si_l) ;
    K(E_Si+j,I_OH+j)       += volume[i]*(phi_l*dx_si_lsdx_oh + dphisdx_oh*s_l*x_si_l) ;
    K(E_Si+j,I_HCO3+j)     += volume[i]*(phi_l*dx_si_lsdx_hco3  + dphisdx_hco3*s_l*x_si_l) ;
    K(E_Si+j,I_CC+j)       += volume[i]*(dphisdn_cc*s_l*x_si_l) ;
    K(E_Si+j,I_Si_S+j)     += volume[i]*(dphisdn_si_s*s_l*x_si_l + 1) ;
    
    /* sauvegarde pour termes de transport */
    Dx_h2co3SDx_co2[i]   = dx_h2co3sdx_co2;
    
    Dx_hSDx_oh[i]        = dx_hsdx_oh ;
    Dx_co3SDx_oh[i]      = dx_co3sdx_oh ;
    Dx_caSDx_oh[i]       = dx_casdx_oh;
    Dx_naohSDx_oh[i]     = dx_naohsdx_oh;
    Dx_naco3SDx_oh[i]    = dx_naco3sdx_oh;
    Dx_kohSDx_oh[i]      = dx_kohsdx_oh;
    Dx_cahco3SDx_oh[i]   = dx_cahco3sdx_oh;
    Dx_caoh2aqSDx_oh[i]  = dx_caoh2aqsdx_oh;
    Dx_h2sio4SDx_oh[i]   = dx_h2sio4sdx_oh;
    Dx_h3sio4SDx_oh[i]   = dx_h3sio4sdx_oh;
    Dx_h4sio4SDx_oh[i]   = dx_h4sio4sdx_oh;
    Dx_cah2sio4SDx_oh[i] = dx_cah2sio4sdx_oh;
    Dx_cah3sio4SDx_oh[i] = dx_cah3sio4sdx_oh;
    
    Dx_caSDx_hco3[i]       = dx_casdx_hco3;
    Dx_co3SDx_hco3[i]      = dx_co3sdx_hco3;
    Dx_nahco3SDx_hco3[i]   = dx_nahco3sdx_hco3;
    Dx_naco3SDx_hco3[i]    = dx_naco3sdx_hco3;
    Dx_caoh2aqSDx_hco3[i]  = dx_caoh2aqsdx_hco3;
    Dx_caohSDx_hco3[i]     = dx_caohsdx_hco3;
    Dx_h2sio4SDx_hco3[i]   = dx_h2sio4sdx_hco3;
    Dx_h3sio4SDx_hco3[i]   = dx_h3sio4sdx_hco3;
    Dx_h4sio4SDx_hco3[i]   = dx_h4sio4sdx_hco3;
    Dx_cah2sio4SDx_hco3[i] = dx_cah2sio4sdx_hco3;
    Dx_cah3sio4SDx_hco3[i] = dx_cah3sio4sdx_hco3;
    
    Dx_naohSDx_na[i]   = dx_naohsdx_na;
    Dx_nahco3SDx_na[i] = dx_nahco3sdx_na;
    Dx_naco3SDx_na[i]  = dx_naco3sdx_na;
    
    Dx_kohSDx_k[i] = dx_kohsdx_k ;
    
  }

  /* termes d'ecoulement */
  {
  double tr        = dt*surf/dx ;

  double trd_h2co3 = tr*KD_H2CO3 ;
  double trd_hco3  = tr*KD_HCO3 ;
  double trd_co3   = tr*KD_CO3 ;
  /*
  double trd_oh    = tr*KD_OH ;
  double trd_h     = tr*KD_H ;
  */
  double trd_ca    = tr*KD_Ca ;
  double trd_na    = tr*KD_Na ;
  double trd_naoh  = tr*KD_NaOH ;
  double trd_nahco3= tr*KD_NaHCO3 ;
  double trd_naco3 = tr*KD_NaCO3 ;
  double trd_m     = tr*KD_m ;
  double trd_k     = tr*KD_K ;
  double trd_koh   = tr*KD_KOH ;
  double trd_caoh  = tr*KD_CaOH ;
  double trd_cahco3= tr*KD_CaHCO3 ;
  double trd_caco3aq= tr*KD_CaCO3aq ;
  double trd_caoh2aq= tr*KD_CaOH2aq ;
  double trd_h3sio4 = tr*KD_H3SiO4 ;
  double trd_h2sio4 = tr*KD_H2SiO4 ;
  double trd_h4sio4 = tr*KD_H4SiO4 ;
  double trd_cah2sio4=tr*KD_CaH2SiO4 ;
  double trd_cah3sio4=tr*KD_CaH3SiO4 ;

  double trf_co2   = tr*KF_CO2 ;
  double trf_h2co3 = tr*KF_H2CO3 ;
  double trf_hco3  = tr*KF_HCO3 ;
  double trf_co3   = tr*KF_CO3 ;
  double trf_ca    = tr*KF_Ca ;
  double trf_oh    = tr*KF_OH ;
  double trf_h     = tr*KF_H ;
  double trf_na    = tr*KF_Na ;
  double trf_naoh  = tr*KF_NaOH ;
  double trf_nahco3= tr*KF_NaHCO3 ;
  double trf_naco3 = tr*KF_NaCO3 ;
  double trf_k     = tr*KF_K ;
  double trf_koh   = tr*KF_KOH ;
  double trf_caoh  = tr*KF_CaOH ;
  double trf_cahco3= tr*KF_CaHCO3 ;
  /*
  double trf_caco3aq= tr*KF_CaCO3aq ;
  */
  double trf_caoh2aq= tr*KF_CaOH2aq ; 
  double trf_h3sio4 = tr*KF_H3SiO4 ;
  double trf_h2sio4 = tr*KF_H2SiO4 ;
  double trf_h4sio4 = tr*KF_H4SiO4 ;
  double trf_cah2sio4=tr*KF_CaH2SiO4 ;
  double trf_cah3sio4=tr*KF_CaH3SiO4 ;

  double tre_hco3  = tr*Kpsi_HCO3 ;
  double tre_co3   = tr*Kpsi_CO3 ;
  double tre_ca    = tr*Kpsi_Ca ;
  /*
  double tre_oh    = tr*Kpsi_OH ;
  double tre_h     = tr*Kpsi_H ;
  */
  double tre_na    = tr*Kpsi_Na ;
  double tre_naco3 = tr*Kpsi_NaCO3 ;
  double tre_k     = tr*Kpsi_K ;
  double tre_caoh  = tr*Kpsi_CaOH ;
  double tre_cahco3= tr*Kpsi_CaHCO3 ;
  double tre_h3sio4 = tr*Kpsi_H3SiO4 ;
  double tre_h2sio4 = tr*Kpsi_H2SiO4 ;
  double tre_cah3sio4=tr*Kpsi_CaH3SiO4 ;
  double tre_q      = tr*Kpsi_q ;
  
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_h2co3 + trd_hco3 + trd_co3 + trd_nahco3 + trd_naco3  + trd_cahco3 + trd_caco3aq ;
  }
  K(E_C,I_P_l)          += + c[0] ;
  K(E_C,I_P_l+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_P_l)      += - c[0] ;
  K(E_C+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co2 + trf_h2co3*Dx_h2co3SDx_co2[i] ;
  }
  K(E_C,I_CO2)          += + c[0] ;
  K(E_C,I_CO2+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_CO2)      += - c[0] ;
  K(E_C+NEQ,I_CO2+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3 + trf_co3*Dx_co3SDx_hco3[i] + trf_nahco3*Dx_nahco3SDx_hco3[i] + trf_naco3*Dx_naco3SDx_hco3[i] ;
  }
  K(E_C,I_HCO3)         += + c[0] ;
  K(E_C,I_HCO3+NEQ)     += - c[1] ;
  K(E_C+NEQ,I_HCO3)     += - c[0] ;
  K(E_C+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co3*Dx_co3SDx_oh[i] + trf_naco3*Dx_naco3SDx_oh[i] + trf_cahco3*Dx_cahco3SDx_oh[i];
  }
  K(E_C,I_OH)           += + c[0] ;
  K(E_C,I_OH+NEQ)       += - c[1] ;
  K(E_C+NEQ,I_OH)       += - c[0] ;
  K(E_C+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 + tre_naco3 + tre_cahco3;
  }
  K(E_C,I_psi)          += + c[0] ;
  K(E_C,I_psi+NEQ)      += - c[1] ;
  K(E_C+NEQ,I_psi)      += - c[0] ;
  K(E_C+NEQ,I_psi+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
  c[i] = trf_nahco3*Dx_nahco3SDx_na[i] + trf_naco3*Dx_naco3SDx_na[i] ;
  }
  K(E_C,I_Na)           += + c[0] ;
  K(E_C,I_Na+NEQ)       += - c[1] ;
  K(E_C+NEQ,I_Na)       += - c[0] ;
  K(E_C+NEQ,I_Na+NEQ)   += + c[1] ;

  /*
    Conservation de la charge  : div(w_q) = 0
  */
  for(i=0;i<2;i++){
    c[i] = z_hco3*trf_hco3 + z_co3*trf_co3*Dx_co3SDx_hco3[i] + z_ca*trf_ca*Dx_caSDx_hco3[i] + z_naco3*trf_naco3*Dx_naco3SDx_hco3[i] + z_caoh*trf_caoh*Dx_caohSDx_hco3[i] + z_h2sio4*trf_h2sio4*Dx_h2sio4SDx_hco3[i] + z_h3sio4*trf_h3sio4*Dx_h3sio4SDx_hco3[i] + z_cah3sio4*trf_cah3sio4*Dx_cah3sio4SDx_hco3[i]  ;
  }
  K(E_q,I_HCO3)         += + c[0] ;
  K(E_q,I_HCO3+NEQ)     += - c[1] ;
  K(E_q+NEQ,I_HCO3)     += - c[0] ;
  K(E_q+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_h*trf_h*Dx_hSDx_oh[i] + z_oh*trf_oh + z_co3*trf_co3*Dx_co3SDx_oh[i] + z_ca*trf_ca*Dx_caSDx_oh[i] + z_naco3*trf_naco3*Dx_naco3SDx_oh[i] + z_cahco3*trf_cahco3*Dx_cahco3SDx_oh[i] + z_h2sio4*trf_h2sio4*Dx_h2sio4SDx_oh[i] + z_h3sio4*trf_h3sio4*Dx_h3sio4SDx_oh[i] + z_cah3sio4*trf_cah3sio4*Dx_cah3sio4SDx_oh[i] ;
  }
  K(E_q,I_OH)           += + c[0] ;
  K(E_q,I_OH+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_OH)       += - c[0] ;
  K(E_q+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = z_na*trf_na + z_naco3*trf_naco3*Dx_naco3SDx_na[i] ;
  }
  K(E_q,I_Na)           += + c[0] ;
  K(E_q,I_Na+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_Na)       += - c[0] ;
  K(E_q+NEQ,I_Na+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = z_k*trf_k  ;
  }
  K(E_q,I_K)           += + c[0] ;
  K(E_q,I_K+NEQ)       += - c[1] ;
  K(E_q+NEQ,I_K)       += - c[0] ;
  K(E_q+NEQ,I_K+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_q ;
  }
  K(E_q,I_psi)          += + c[0] ;
  K(E_q,I_psi+NEQ)      += - c[1] ;
  K(E_q+NEQ,I_psi)      += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)  += + c[1] ;
  
  

  /*
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  K(E_mass,I_P_l)          += + trd_m ;
  K(E_mass,I_P_l+NEQ)      += - trd_m ;
  K(E_mass+NEQ,I_P_l)      += - trd_m ;
  K(E_mass+NEQ,I_P_l+NEQ)  += + trd_m ;

  for(i=0;i<2;i++){
    c[i] = M_CO2*trf_co2 ;
  }
  K(E_mass,I_CO2)          += + c[0] ;
  K(E_mass,I_CO2+NEQ)      += - c[1] ;
  K(E_mass+NEQ,I_CO2)      += - c[0] ;
  K(E_mass+NEQ,I_CO2+NEQ)  += + c[1] ;

  
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_ca + trd_caoh + trd_cahco3 + trd_caco3aq + trd_caoh2aq + trd_cah2sio4 + trd_cah3sio4;
  }
  K(E_Ca,I_P_l)          += + c[0] ;
  K(E_Ca,I_P_l+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_P_l)      += - c[0] ;
  K(E_Ca+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dx_caSDx_hco3[i] + trf_caoh*Dx_caohSDx_hco3[i] + trf_caoh2aq*Dx_caoh2aqSDx_hco3[i] +trf_cah2sio4*Dx_cah2sio4SDx_hco3[i] + trf_cah3sio4*Dx_cah3sio4SDx_hco3[i] ;
  }
  K(E_Ca,I_HCO3)         += + c[0] ;
  K(E_Ca,I_HCO3+NEQ)     += - c[1] ;
  K(E_Ca+NEQ,I_HCO3)     += - c[0] ;
  K(E_Ca+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_ca*Dx_caSDx_oh[i] + trf_cahco3*Dx_cahco3SDx_oh[i] + trf_caoh2aq*Dx_caoh2aqSDx_oh[i] + trf_cah2sio4*Dx_cah2sio4SDx_oh[i] + trf_cah3sio4*Dx_cah3sio4SDx_oh[i];
  }
  K(E_Ca,I_OH)           += + c[0] ;
  K(E_Ca,I_OH+NEQ)       += - c[1] ;
  K(E_Ca+NEQ,I_OH)       += - c[0] ;
  K(E_Ca+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_ca + tre_caoh + tre_cahco3 + tre_cah3sio4;
  }
  K(E_Ca,I_psi)          += + c[0] ;
  K(E_Ca,I_psi+NEQ)      += - c[1] ;
  K(E_Ca+NEQ,I_psi)      += - c[0] ;
  K(E_Ca+NEQ,I_psi+NEQ)  += + c[1] ;
  
  
  /*
    Cinetique 1 : (n_k1 - n_kn) + dt * div(W_k) - dt * XI = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_hco3 + trd_co3 + trd_nahco3 + trd_naco3 + trd_caco3aq + trd_cahco3 ;
  }
  K(E_k,I_P_l)          += + c[0] ;
  K(E_k,I_P_l+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_P_l)      += - c[0] ;
  K(E_k+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_hco3 + trf_co3*Dx_co3SDx_hco3[i] + trd_nahco3*Dx_nahco3SDx_hco3[i] + trd_naco3*Dx_naco3SDx_hco3[i]  ;
  }
  K(E_k,I_HCO3)          += + c[0] ;
  K(E_k,I_HCO3+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_HCO3)      += - c[0] ;
  K(E_k+NEQ,I_HCO3+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_co3*Dx_co3SDx_oh[i] + trd_naco3*Dx_naco3SDx_oh[i] + trd_cahco3*Dx_cahco3SDx_oh[i] ;
  }
  K(E_k,I_OH)            += + c[0] ;
  K(E_k,I_OH+NEQ)        += - c[1] ;
  K(E_k+NEQ,I_OH)        += - c[0] ;
  K(E_k+NEQ,I_OH+NEQ)    += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_hco3 + tre_co3 + tre_naco3 + tre_cahco3 ;
  }
  K(E_k,I_psi)           += + c[0] ;
  K(E_k,I_psi+NEQ)       += - c[1] ;
  K(E_k+NEQ,I_psi)       += - c[0] ;
  K(E_k+NEQ,I_psi+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] =  trd_nahco3*Dx_nahco3SDx_na[i] + trd_naco3*Dx_naco3SDx_na[i] ;
  }
  K(E_k,I_Na)          += + c[0] ;
  K(E_k,I_Na+NEQ)      += - c[1] ;
  K(E_k+NEQ,I_Na)      += - c[0] ;
  K(E_k+NEQ,I_Na+NEQ)  += + c[1] ;

  /*
    Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_na + trd_naoh + trd_nahco3 + trd_naco3 ;
  }
  K(E_Na,I_P_l)          += + c[0] ;
  K(E_Na,I_P_l+NEQ)      += - c[1] ;
  K(E_Na+NEQ,I_P_l)      += - c[0] ;
  K(E_Na+NEQ,I_P_l+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_nahco3*Dx_nahco3SDx_hco3[i] + trf_naco3*Dx_naco3SDx_hco3[i] ;
  }
  K(E_Na,I_HCO3)         += + c[0] ;
  K(E_Na,I_HCO3+NEQ)     += - c[1] ;
  K(E_Na+NEQ,I_HCO3)     += - c[0] ;
  K(E_Na+NEQ,I_HCO3+NEQ) += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_naoh*Dx_naohSDx_oh[i] + trf_naco3*Dx_naco3SDx_oh[i];
  }
  K(E_Na,I_OH)           += + c[0] ;
  K(E_Na,I_OH+NEQ)       += - c[1] ;
  K(E_Na+NEQ,I_OH)       += - c[0] ;
  K(E_Na+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_na + tre_naco3 ;
  }
  K(E_Na,I_psi)          += + c[0] ;
  K(E_Na,I_psi+NEQ)      += - c[1] ;
  K(E_Na+NEQ,I_psi)      += - c[0] ;
  K(E_Na+NEQ,I_psi+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_na + trf_naoh*Dx_naohSDx_na[i] + trf_nahco3*Dx_nahco3SDx_na[i] + trf_naco3*Dx_naco3SDx_na[i] ;
  }
  K(E_Na,I_Na)          += + c[0] ;
  K(E_Na,I_Na+NEQ)      += - c[1] ;
  K(E_Na+NEQ,I_Na)      += - c[0] ;
  K(E_Na+NEQ,I_Na+NEQ)  += + c[1] ;
  
  /*
    Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_k + trd_koh ;
  }
  K(E_K,I_P_l)          += + c[0] ;
  K(E_K,I_P_l+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_P_l)      += - c[0] ;
  K(E_K+NEQ,I_P_l+NEQ)  += + c[1] ;


  for(i=0;i<2;i++){
    c[i] = trf_koh*Dx_kohSDx_oh[i] ;
  }
  K(E_K,I_OH)           += + c[0] ;
  K(E_K,I_OH+NEQ)       += - c[1] ;
  K(E_K+NEQ,I_OH)       += - c[0] ;
  K(E_K+NEQ,I_OH+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_k  ;
  }
  K(E_K,I_psi)          += + c[0] ;
  K(E_K,I_psi+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_psi)      += - c[0] ;
  K(E_K+NEQ,I_psi+NEQ)  += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = trf_k + trf_koh*Dx_kohSDx_k[i]  ;
  }
  K(E_K,I_K)          += + c[0] ;
  K(E_K,I_K+NEQ)      += - c[1] ;
  K(E_K+NEQ,I_K)      += - c[0] ;
  K(E_K+NEQ,I_K+NEQ)  += + c[1] ;

  /*
    Conservation de Si (silice) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  for(i=0;i<2;i++){
    c[i] = trd_h3sio4 + trd_h2sio4 + trd_h4sio4 + trd_cah2sio4 + trd_cah3sio4;
  }
  K(E_Si,I_P_l)          += + c[0] ;
  K(E_Si,I_P_l+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_P_l)      += - c[0] ;
  K(E_Si+NEQ,I_P_l+NEQ)  += + c[1] ;


  for(i=0;i<2;i++){
    c[i] = trf_h2sio4*Dx_h2sio4SDx_oh[i] + trf_h3sio4*Dx_h3sio4SDx_oh[i] + trf_h4sio4*Dx_h4sio4SDx_oh[i] + trf_cah2sio4*Dx_cah2sio4SDx_oh[i] + trf_cah3sio4*Dx_cah3sio4SDx_oh[i];
  }
  K(E_Si,I_OH)           += + c[0] ;
  K(E_Si,I_OH+NEQ)       += - c[1] ;
  K(E_Si+NEQ,I_OH)       += - c[0] ;
  K(E_Si+NEQ,I_OH+NEQ)   += + c[1] ;
  
  for(i=0;i<2;i++){
    c[i] = trf_h2sio4*Dx_h2sio4SDx_hco3[i] + trf_h3sio4*Dx_h3sio4SDx_hco3[i] + trf_h4sio4*Dx_h4sio4SDx_hco3[i] + trf_cah2sio4*Dx_cah2sio4SDx_hco3[i] + trf_cah3sio4*Dx_cah3sio4SDx_hco3[i];
  }
  K(E_Si,I_HCO3)           += + c[0] ;
  K(E_Si,I_HCO3+NEQ)       += - c[1] ;
  K(E_Si+NEQ,I_HCO3)       += - c[0] ;
  K(E_Si+NEQ,I_HCO3+NEQ)   += + c[1] ;

  for(i=0;i<2;i++){
    c[i] = tre_h3sio4 + tre_h2sio4  + tre_cah3sio4;
  }
  K(E_Si,I_psi)          += + c[0] ;
  K(E_Si,I_psi+NEQ)      += - c[1] ;
  K(E_Si+NEQ,I_psi)      += - c[0] ;
  K(E_Si+NEQ,I_psi+NEQ)  += + c[1] ;
  }


#if (U_CO2 == LOG_RHO)
  for(i=0;i<2*NEQ;i++){
    K(i,I_CO2)     *= Ln10*X_CO2(0) ;
    K(i,I_CO2+NEQ) *= Ln10*X_CO2(1) ;
  }
#endif


  return(0) ;

#undef K
}


int  ComputeResidu(Element_t *el,double t,double dt,double *r)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  Symmetry_t sym = Element_GetSymmetry(el) ;
  double *f = Element_GetCurrentImplicitTerm(el) ;
  double *f_n = Element_GetPreviousImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0. ;
  /*
    INITIALISATION DU RESIDU
  */
  for(i=0;i<NEQ*2;i++) r[i] = zero ;

  if(Element_IsSubmanifold(el)) return(0) ;
  /*
    CALCUL DE volume ET DE surf
  */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    dx = x1 - x0 ;
    xm = (x1 + x0)*0.5 ;
  }
  for(i=0;i<nn;i++) {
    double x = Element_GetNodeCoordinate(el,i)[0] ;
    volume[i] = fabs(dx)*0.5 ; 
    if(sym == AXIS) volume[i] *= M_PI*(x + xm) ; 
  }
  if(sym == AXIS) surf = 2*M_PI*xm ; else surf = 1. ;
  /*
    Conservation de C (carbone) : (n_C1 - n_Cn) + dt * div(w_C) = 0
  */
  R(0,E_C) -= volume[0]*(N_C(0) - N_Cn(0)) + dt*surf*W_C ;
  R(1,E_C) -= volume[1]*(N_C(1) - N_Cn(1)) - dt*surf*W_C ;
  /*
    Conservation de la charge  : div(w_q) = 0
  */
  R(0,E_q) -= + dt*surf*W_q ;
  R(1,E_q) -= - dt*surf*W_q ;
  /*
    Conservation de la masse totale : (m_1 - m_n) + dt * div(w) = 0
  */
  R(0,E_mass) -= volume[0]*(Mass(0) - Mass_n(0)) + dt*surf*W_m ;
  R(1,E_mass) -= volume[1]*(Mass(1) - Mass_n(1)) - dt*surf*W_m ;
  /*
    Conservation de Ca (calcium) : (n_Ca1 - n_Can) + dt * div(w_Ca) = 0
  */
  R(0,E_Ca) -= volume[0]*(N_Ca(0) - N_Can(0)) + dt*surf*W_Ca ;
  R(1,E_Ca) -= volume[1]*(N_Ca(1) - N_Can(1)) - dt*surf*W_Ca ;
  /*
    Cinetique 1 : (n_k11 - n_k1n) + dt * div(W_k) - dt * XI = 0
  */
  R(0,E_k) -= volume[0]*(N_k(0) - N_kn(0) - dt*XI(0)) + dt*surf*W_k ;
  R(1,E_k) -= volume[1]*(N_k(1) - N_kn(1) - dt*XI(1)) - dt*surf*W_k ;
  /*
    Electroneutralite : q = 0
  */
  R(0,E_el) -= volume[0]*N_q(0) ;
  R(1,E_el) -= volume[1]*N_q(1) ;
  /*
    Conservation de Na (sodium) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  R(0,E_Na) -= volume[0]*(N_Na(0) - N_Nan(0)) + dt*surf*W_Na ;
  R(1,E_Na) -= volume[1]*(N_Na(1) - N_Nan(1)) - dt*surf*W_Na ;
  /*
    Conservation de K (potassium) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ;

  /*
    Conservation de Si (silice) : (n_Si1 - n_Sin) + dt * div(w_Si) = 0
  */
  R(0,E_Si) -= volume[0]*(N_Si(0) - N_Sin(0)) + dt*surf*W_Si ;
  R(1,E_Si) -= volume[1]*(N_Si(1) - N_Sin(1)) - dt*surf*W_Si ;

  return(0) ;
#undef R
}



int  ComputeOutputs(Element_t *el,double t,double *s,Result_t *r)
/* Les valeurs exploitees (s) */
{
  double *f = Element_GetCurrentImplicitTerm(el) ;
  int nn = Element_GetNbOfNodes(el) ;
  int    i,nso ;
  double zero = 0. ;
  double *u[MAX_NOEUDS] ;
#define UNKNOWN_s(i)   UNKNOWN(j,i)


  if(Element_IsSubmanifold(el)) return(0) ;
  
  for(i = 0 ; i < nn ; i++) {
    u[i] = Element_GetNodalUnknown(el,i) ;
  }
  
  /*
    Donnees
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  n_si_s0  = GetProperty("N_Si") ;

  /* initialisation */
  nso = 36 ;
  for(i = 0 ; i < nso ; i++) {
    int j ;
    for(j = 0 ; j < 9 ; j++) Result_GetValue(r+i)[j] = zero ;
  }

  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double xm = (x0 + x1)*0.5 ;
    int j = (Element_IsSubmanifold(el)) ? 0 : ((s[0] < xm) ? 0 : 1) ;
    /* pression */
    double p_l     =  UNKNOWN_s(I_P_l) ;
    /* saturation */
    double p_c     = p_g - p_l ;
    double s_l     = SATURATION(p_c) ;
    /* molarites */
    double *x = ComputeAqueousParameters(el,u,j) ;
    
    double x_co2      = x[I_CO2] ;
    double x_oh    	  = x[I_OH] ;
    double x_hco3  	  = x[I_HCO3] ;
    double x_na    	  = x[I_Na] ;
    double x_k     	  = x[I_K] ;
    double x_co3   	  = x[I_CO3] ;
    double x_h     	  = x[I_H] ;
    double x_ca    	  = x[I_Ca] ;
    double x_naoh    	= x[I_NaOH] ;
    double x_nahco3  	= x[I_NaHCO3] ;
    double x_naco3 	  = x[I_NaCO3] ;
    double x_koh   	  = x[I_KOH] ;
    double x_caoh  	  = x[I_CaOH] ;
    double x_cahco3  	= x[I_CaHCO3] ;
    double x_caco3aq 	= x[I_CaCO3aq] ;
    double x_caoh2aq 	= x[I_CaOH2aq] ;
    double q_ch       = x[I_Q_CH] ;
    double x_h4sio4   = x[I_H4SiO4] ;
    double x_h3sio4   = x[I_H3SiO4] ;
    double x_h2sio4   = x[I_H2SiO4] ;
    double x_cah2sio4 = x[I_CaH2SiO4];
    double x_cah3sio4 = x[I_CaH3SiO4] ;

    /* Force Ionique */
    double I = 0.5*(z_ca*z_ca*x_ca + z_h2sio4*z_h2sio4*x_h2sio4 + z_co3*z_co3*x_co3   + z_cahco3*z_cahco3*x_cahco3 + z_caoh*z_caoh*x_caoh + z_k*z_k*x_k + z_na*z_na*x_na + z_h*z_h*x_h + z_h3sio4*z_h3sio4*x_h3sio4 + z_naco3*z_naco3*x_naco3 + z_hco3*z_hco3*x_hco3 + z_oh*z_oh*x_oh ) ;

    /* densite de charge */
    double x_q = x[I_N_Q] ;

    /* contenus solides */
    double n_si_s     = UNKNOWN_s(I_Si_S) ;
    double n_cc 	    = UNKNOWN_s(I_CC) ;
    double n_ch       = (N_CH(0) + N_CH(1))*0.5 ;
    double x_csh      = X_CSH(q_ch) ;
    double v_csh      = V_CSH(q_ch) ;
    double n_ca_s     = n_ch + n_cc + x_csh*n_si_s ;
  
    double CsurS      = (x_csh + n_ch/n_si_s) ;
    

    /* porosite */
    double v_csh0     = V_CSH(1) ;
    double phi = phii + V_CH*(n_ch0 - n_ch) + v_csh0*n_si_s0 - v_csh*n_si_s + V_CC*(-n_cc) ;

    double ph = 14 + log(x_oh)/log(10.) ;
    double n_Na = 0.5*(N_Na(0) + N_Na(1)) ;
    double n_Ca = 0.5*(N_Ca(0) + N_Ca(1)) ;
    double n_Si = 0.5*(N_Si(0) + N_Si(1)) ; 


    /* Transferts */
    double dx        = x1 - x0 ;
    double grd_psi   = (PSI(1) - PSI(0))/dx ;

    /* quantites exploitees */
    i = 0 ;
    Result_Store(r + i++,&x_co2,"x_co2",1) ;
    Result_Store(r + i++,&ph,"ph",1) ;
    Result_Store(r + i++,&n_si_s,"n_Si_s",1) ;
    Result_Store(r + i++,&phi,"porosite",1) ;
    Result_Store(r + i++,&n_ch,"n_CH",1) ;
    Result_Store(r + i++,&x_ca,"x_ca",1) ;
    Result_Store(r + i++,&x_co3,"x_co3",1) ;
    Result_Store(r + i++,&x_hco3,"x_hco3",1) ;
    Result_Store(r + i++,&n_cc,"n_CC",1) ;
    Result_Store(r + i++,&x_h,"x_h",1) ;
    Result_Store(r + i++,&x_oh,"x_oh",1) ;
    Result_Store(r + i++,&s_l,"saturation",1) ;
    Result_Store(r + i++,&grd_psi,"grad_psi",1) ;
    Result_Store(r + i++,&x_q,"charge",1) ;
    Result_Store(r + i++,&x_na,"x_na",1) ;
    Result_Store(r + i++,&x_naoh,"x_naoh",1) ;
    Result_Store(r + i++,&x_nahco3,"x_nahco3",1) ;
    Result_Store(r + i++,&x_naco3,"x_naco3",1) ;
    Result_Store(r + i++,&x_k,"x_k",1) ;
    Result_Store(r + i++,&x_koh,"x_koh",1) ;
    Result_Store(r + i++,&x_caoh,"x_caoh",1) ;
    Result_Store(r + i++,&x_cahco3,"x_cahco3",1) ;
    Result_Store(r + i++,&x_caco3aq,"x_caco3aq",1) ;
    Result_Store(r + i++,&x_caoh2aq,"x_caoh2aq",1) ;
    Result_Store(r + i++,&p_l,"p_l",1) ;
    Result_Store(r + i++,&x_h3sio4,"x_h3sio4",1) ;
    Result_Store(r + i++,&n_Na,"n_Na",1) ;
    Result_Store(r + i++,&n_Ca,"n_Ca",1) ;
    Result_Store(r + i++,&n_Si,"n_Si",1) ;
    Result_Store(r + i++,&n_ca_s,"n_Ca_s",1) ;
    Result_Store(r + i++,&x_cah2sio4,"x_cah2sio4",1) ;
    Result_Store(r + i++,&x_cah3sio4,"x_cah3sio4",1) ;
    Result_Store(r + i++,&CsurS,"CsurS",1) ;
    Result_Store(r + i++,&x_h4sio4,"x_h4sio4",1) ;
    Result_Store(r + i++,&x_h2sio4,"x_h2sio4",1) ;
    Result_Store(r + i++,&I,"I",1) ;
  }
  
  
  if(i != nso) arret("ComputeOutputs") ;
  return(nso) ;
}


void transfert(Element_t *el,double **u,double *f)
/* Termes explicites (va)  */
{
  double *va = Element_GetExplicitTerm(el) ;
  int    i ; 
  /*
    Donnees
  */
  
  phii     = GetProperty("porosite") ;
  k_int    = GetProperty("k_int") ;
  a_1      = GetProperty("A_1") ;
  a_2      = GetProperty("A_2") ;
  c_2      = GetProperty("C_2") ;
  n_ch0    = GetProperty("N_CaOH2") ;
  n_si_s0  = GetProperty("N_Si") ;

  /* initialisation */
  for(i = 0 ; i < NVE ; i++) va[i] = 0. ;
  /*
    Coefficients de transfert
  */
  for(i = 0 ; i < 2 ; i++) {
    /* pressions */
    double p_l     = P_l(i) ;
    double p_c     = p_g - p_l ;
    /* saturation */
    double s_l     = SATURATION(p_c) ;
    
    double *x = ComputeAqueousParameters(el,u,i) ;
    
    double x_h    	  = x[I_H] ;
    double x_oh    	  = x[I_OH] ;
    double x_hco3  	  = x[I_HCO3] ;
    double x_na    	  = x[I_Na] ;
    double x_k     	  = x[I_K] ;
    double x_h2co3 	  = x[I_H2CO3] ;
    double x_co3   	  = x[I_CO3] ;
    double x_ca    	  = x[I_Ca] ;
    double x_naoh    	= x[I_NaOH] ;
    double x_nahco3  	= x[I_NaHCO3] ;
    double x_naco3 	  = x[I_NaCO3] ;
    double x_koh   	  = x[I_KOH] ;
    double x_caoh  	  = x[I_CaOH] ;
    double x_cahco3  	= x[I_CaHCO3] ;
    double x_caco3aq 	= x[I_CaCO3aq] ;
    double x_caoh2aq 	= x[I_CaOH2aq] ;
    double q_ch       = x[I_Q_CH] ;
    double x_h4sio4   = x[I_H4SiO4] ;
    double x_h3sio4   = x[I_H3SiO4] ;
    double x_h2sio4   = x[I_H2SiO4] ;
    double x_cah2sio4 = x[I_CaH2SiO4];
    double x_cah3sio4 = x[I_CaH3SiO4] ;
    
    /* masse volumique liquide */
    double rho_l      = x[I_RHO_L] ;
    
    /* solides */
    double n_ch    = N_CH(i) ;
    double n_cc    = N_CC(i) ;
    double n_si_s  = N_Si_S(i) ;
  
    /* autres solides */
    double v_csh      = V_CSH(q_ch) ;
    
    /* porosite */
    double v_csh0   = V_CSH(1) ;
    double phi = phii + V_CH*(n_ch0 - n_ch) + v_csh0*n_si_s0 - v_csh*n_si_s + V_CC*(-n_cc) ;


    /* permeabilite */
    double k_l  = (k_int/mu_l)*RELATIVEPERM(p_c)*pow(phi/phii,3.)*pow(((1-phii)/(1-phi)),2.) ;
    /* tortuosite gaz */
    /*double tau  = pow(phi,1/3)*pow(s_g,7/3) ;*/
    /* double tau  = pow(phi,1.74)*pow(s_g,3.20) ; */
    /* tortuosite liquide */
    double iff    = 0.00029*exp(9.95*phi)/(1+625*pow((1-s_l),4)) ;

    /* Humidité relative */
    double hr = exp(-p_c*M_H2O/(RT*rho_l)) ;
 
    KD_Ca     	+= x_ca*k_l ;
    KD_H2CO3  	+= x_h2co3*k_l ;
    KD_HCO3   	+= x_hco3*k_l ;
    KD_CO3    	+= x_co3*k_l ;
    KD_OH     	+= x_oh*k_l ;
    KD_H      	+= x_h*k_l ;
    KD_Na     	+= x_na*k_l ;
    KD_NaOH   	+= x_naoh*k_l ;
    KD_NaHCO3 	+= x_nahco3*k_l ;
    KD_NaCO3  	+= x_naco3*k_l ;
    KD_K      	+= x_k*k_l ;
    KD_KOH    	+= x_koh*k_l ;
    KD_CaOH   	+= x_caoh*k_l ;
    KD_CaHCO3 	+= x_cahco3*k_l ;
    KD_CaCO3aq	+= x_caco3aq*k_l ;
    KD_CaOH2aq	+= x_caoh2aq*k_l ;
    KD_H3SiO4 	+= x_h3sio4*k_l ;
    KD_H2SiO4 	+= x_h2sio4*k_l ;
    KD_H4SiO4 	+= x_h4sio4*k_l ;
    KD_CaH3SiO4 += x_cah3sio4*k_l ;
    KD_CaH2SiO4 += x_cah2sio4*k_l ;
    KD_m      	+= rho_l*k_l ;
    
    /* KF_CO2    += phi_g*tau*d_co2 ; */
    KF_CO2    	+= (1.6e-3)*pow(phi,1.8)*pow(1-hr,2.2) ;
    KF_Ca     	+= d_ca*iff ;
    KF_OH     	+= d_oh*iff ;
    KF_H      	+= d_h*iff ;
    KF_H2CO3  	+= d_h2co3*iff ;
    KF_HCO3   	+= d_hco3*iff ;
    KF_CO3    	+= d_co3*iff ;
    KF_Na     	+= d_na*iff ;
    KF_NaOH   	+= d_naoh*iff ;
    KF_NaHCO3 	+= d_nahco3*iff ;
    KF_NaCO3  	+= d_naco3*iff ;
    KF_K      	+= d_k*iff ;
    KF_KOH    	+= d_koh*iff ;
    KF_CaOH   	+= d_caoh*iff ;
    KF_CaHCO3 	+= d_cahco3*iff ;
    KF_CaCO3aq	+= d_caco3aq*iff ;
    KF_CaOH2aq	+= d_caoh2aq*iff ;
    KF_H3SiO4 	+= d_h3sio4*iff ;
    KF_H2SiO4 	+= d_h2sio4*iff ;
    KF_H4SiO4 	+= d_h4sio4*iff ;
    KF_CaH3SiO4 += d_cah3sio4*iff ;
    KF_CaH2SiO4 += d_cah2sio4*iff ;
    
    Kpsi_Ca   	 += FARADAY/RT*KF_Ca*z_ca*x_ca ;
    Kpsi_HCO3 	 += FARADAY/RT*KF_HCO3*z_hco3*x_hco3 ;
    Kpsi_CO3  	 += FARADAY/RT*KF_CO3*z_co3*x_co3 ;
    Kpsi_OH   	 += FARADAY/RT*KF_OH*z_oh*x_oh ;
    Kpsi_H    	 += FARADAY/RT*KF_H*z_h*x_h ;
    Kpsi_Na    	 += FARADAY/RT*KF_Na*z_na*x_na ;
    Kpsi_NaCO3 	 += FARADAY/RT*KF_NaCO3*z_naco3*x_naco3 ;
    Kpsi_K     	 += FARADAY/RT*KF_K*z_k*x_k ;
    Kpsi_CaOH 	 += FARADAY/RT*KF_CaOH*z_caoh*x_caoh ;
    Kpsi_CaHCO3	 +=FARADAY/RT*KF_CaHCO3*z_cahco3*x_cahco3 ;
    Kpsi_H3SiO4	 +=FARADAY/RT*KF_H3SiO4*z_h3sio4*x_h3sio4 ;
    Kpsi_H2SiO4	 +=FARADAY/RT*KF_H2SiO4*z_h2sio4*x_h2sio4 ;
    Kpsi_CaH3SiO4+=FARADAY/RT*KF_CaH3SiO4*z_cah3sio4*x_cah3sio4 ;
    Kpsi_q    	 += z_h*Kpsi_H + z_oh*Kpsi_OH + z_hco3*Kpsi_HCO3 + z_co3*Kpsi_CO3 + z_ca*Kpsi_Ca + z_na*Kpsi_Na + z_naco3*Kpsi_NaCO3 + z_k*Kpsi_K + z_caoh*Kpsi_CaOH + z_cahco3*Kpsi_CaHCO3 + z_h3sio4*Kpsi_H3SiO4 + z_h2sio4*Kpsi_H2SiO4 + z_cah3sio4*Kpsi_CaH3SiO4 ;
  }
  
  /* moyenne */
  for(i = 0 ; i < NVE ; i++) va[i] *= 0.5 ;
}



void flux(Element_t *el,double **u)
/* Les flux (f) */
{
  double *f = Element_GetImplicitTerm(el) ;
  double *va = Element_GetExplicitTerm(el) ;
  
  double r_h[2] ;
  double r_h2co3[2] ;
  double r_co3[2] ;
  double r_ca[2] ;
  double r_naoh[2] ;
  double r_nahco3[2] ;
  double r_naco3[2] ;
  double r_koh[2] ;
  double r_caoh[2] ;
  double r_cahco3[2] ;
  double r_caco3aq[2] ;
  double r_caoh2aq[2] ;
  double r_h2sio4[2] ;
  double r_h3sio4[2] ;
  double r_h4sio4[2] ;
  double r_cah2sio4[2] ;
  double r_cah3sio4[2];
  
  int    i ;

  for(i = 0 ; i < 2 ; i++) {
    double *x = ComputeAqueousParameters(el,u,i) ;
    
    r_h[i]        = x[I_H] ;
    r_h2co3[i]    = x[I_H2CO3] ;
    r_co3[i]      = x[I_CO3] ;
    r_ca[i]       = x[I_Ca] ;
    r_naoh[i]     = x[I_NaOH];
    r_nahco3[i]   = x[I_NaHCO3] ;
    r_naco3[i]    = x[I_NaCO3] ;
    r_koh[i]      = x[I_KOH] ;
    r_caoh[i]     = x[I_CaOH] ;
    r_cahco3[i]   = x[I_CaHCO3] ;
    r_caco3aq[i]  = x[I_CaCO3aq] ;
    r_caoh2aq[i]  = x[I_CaOH2aq] ;
    r_h2sio4[i]   = x[I_H2SiO4] ;
    r_h3sio4[i]   = x[I_H3SiO4] ;
    r_h4sio4[i]   = x[I_H4SiO4] ;
    r_cah2sio4[i] = x[I_CaH2SiO4] ;
    r_cah3sio4[i] = x[I_CaH3SiO4] ;
  }

  /* Gradients */
  {
    double x1 = Element_GetNodeCoordinate(el,1)[0] ;
    double x0 = Element_GetNodeCoordinate(el,0)[0] ;
    double dx = x1 - x0 ;
    
    double grd_p_l      = (P_l(1)      - P_l(0)     )/dx ;
    double grd_co2      = (X_CO2(1)    - X_CO2(0)   )/dx ;
    double grd_oh       = (X_OH(1)     - X_OH(0)    )/dx ;
    double grd_hco3     = (X_HCO3(1)   - X_HCO3(0)  )/dx ;
    double grd_na       = (X_Na(1)     - X_Na(0)    )/dx ;
    double grd_k        = (X_K(1)      - X_K(0)     )/dx ;
    double grd_ca       = (r_ca[1]     - r_ca[0]    )/dx ;
    double grd_h        = (r_h[1]      - r_h[0]     )/dx ;
    double grd_h2co3    = (r_h2co3[1]  - r_h2co3[0] )/dx ;
    double grd_co3      = (r_co3[1]    - r_co3[0]   )/dx ;
    double grd_naoh     = (r_naoh[1]   - r_naoh[0]  )/dx ;
    double grd_nahco3   = (r_nahco3[1] - r_nahco3[0])/dx ;
    double grd_naco3    = (r_naco3[1]  - r_naco3[0] )/dx ;
    double grd_koh      = (r_koh[1]    - r_koh[0]   )/dx ;
    double grd_caoh     = (r_caoh[1]   - r_caoh[0]  )/dx ;
    double grd_cahco3   = (r_cahco3[1] - r_cahco3[0])/dx ;
    double grd_caco3aq  = (r_caco3aq[1]- r_caco3aq[0])/dx ;
    double grd_caoh2aq  = (r_caoh2aq[1]- r_caoh2aq[0])/dx ;
    double grd_h2sio4   = (r_h2sio4[1] - r_h2sio4[0])/dx ;
    double grd_h3sio4   = (r_h3sio4[1] - r_h3sio4[0])/dx ;
    double grd_h4sio4   = (r_h4sio4[1] - r_h4sio4[0])/dx ;
    double grd_cah2sio4 = (r_cah2sio4[1] - r_cah2sio4[0])/dx ;
    double grd_cah3sio4 = (r_cah3sio4[1] - r_cah3sio4[0])/dx ;
    double grd_psi      = (PSI(1)        - PSI(0)     )/dx ;

    /* Flux */
    double w_ca    = - KD_Ca*grd_p_l    - KF_Ca*grd_ca       - Kpsi_Ca*grd_psi    ;
    double w_hco3  = - KD_HCO3*grd_p_l  - KF_HCO3*grd_hco3   - Kpsi_HCO3*grd_psi  ;
    double w_h3sio4= - KD_H3SiO4*grd_p_l- KF_H3SiO4*grd_h3sio4 - Kpsi_H3SiO4*grd_psi ;
    double w_h2sio4= - KD_H2SiO4*grd_p_l- KF_H2SiO4*grd_h2sio4 - Kpsi_H2SiO4*grd_psi ;
    double w_h4sio4= - KD_H4SiO4*grd_p_l- KF_H4SiO4*grd_h4sio4  ;
    double w_co3   = - KD_CO3*grd_p_l   - KF_CO3*grd_co3     - Kpsi_CO3*grd_psi   ; 
    double w_h2co3 = - KD_H2CO3*grd_p_l - KF_H2CO3*grd_h2co3                      ;
    double w_na    = - KD_Na*grd_p_l    - KF_Na*grd_na       - Kpsi_Na*grd_psi    ;
    double w_naoh  = - KD_NaOH*grd_p_l  - KF_NaOH*grd_naoh                        ;
    double w_nahco3= - KD_NaHCO3*grd_p_l- KF_NaHCO3*grd_nahco3                    ;
    double w_naco3 = - KD_NaCO3*grd_p_l - KF_NaCO3*grd_naco3 - Kpsi_NaCO3*grd_psi ;
    double w_k     = - KD_K*grd_p_l     - KF_K*grd_k         - Kpsi_K*grd_psi     ;
    double w_koh   = - KD_KOH*grd_p_l   - KF_KOH*grd_koh                          ;
    double w_caoh  = - KD_CaOH*grd_p_l  - KF_CaOH*grd_caoh   - Kpsi_CaOH*grd_psi  ;
    double w_cahco3= - KD_CaHCO3*grd_p_l- KF_CaHCO3*grd_cahco3-Kpsi_CaHCO3*grd_psi;
    double w_caco3aq=- KD_CaCO3aq*grd_p_l-KF_CaCO3aq*grd_caco3aq                  ;
    double w_caoh2aq=- KD_CaOH2aq*grd_p_l-KF_CaOH2aq*grd_caoh2aq                  ;
    double w_cah3sio4= - KD_CaH3SiO4*grd_p_l- KF_CaH3SiO4*grd_cah3sio4 - Kpsi_CaH3SiO4*grd_psi  ;
    double w_cah2sio4= - KD_CaH2SiO4*grd_p_l- KF_CaH2SiO4*grd_cah2sio4            ;
    double w_co2   =                    - KF_CO2*grd_co2                          ;
    double w_m     = - KD_m*grd_p_l + M_CO2*w_co2 ;
    double w_q     =                    - z_h*KF_H*grd_h             \
                               - z_oh*KF_OH*grd_oh          \
                               - z_hco3*KF_HCO3*grd_hco3    \
                               - z_co3*KF_CO3*grd_co3       \
                               - z_ca*KF_Ca*grd_ca          \
                               - z_na*KF_Na*grd_na          \
                               - z_naco3*KF_NaCO3*grd_naco3 \
                               - z_k*KF_K*grd_k             \
                               - z_caoh*KF_CaOH*grd_caoh    \
                               - z_cahco3*KF_CaHCO3*grd_cahco3 \
                               - z_h3sio4*KF_H3SiO4*grd_h3sio4 \
                               - z_h2sio4*KF_H2SiO4*grd_h2sio4 \
                               - z_cah3sio4*KF_CaH3SiO4*grd_cah3sio4 \
                                                    - Kpsi_q*grd_psi ;


    W_C     = w_co2 + w_h2co3 + w_hco3 + w_co3 + w_nahco3 + w_naco3 + w_cahco3 + w_caco3aq ;
    W_Ca    = w_ca + w_caoh + w_cahco3 + w_caco3aq + w_caoh2aq + w_cah2sio4 + w_cah3sio4;
    W_Na    = w_na + w_naoh + w_nahco3 + w_naco3 ;
    W_m     = w_m ;
    W_k     = w_hco3 + w_co3 + w_nahco3 + w_naco3 + w_caco3aq + w_cahco3 ;
    W_Si    = w_h3sio4 + w_h4sio4 + w_h2sio4 + w_cah2sio4 + w_cah3sio4 ;
    W_q     = w_q ;
    W_K     = w_k + w_koh ;
 }
}


void concentration1(double x_co2, double x_na_tot, double x_k_tot, double *pointeur_x_oh, double *pointeur_x_h3sio4, double *pointeur_x_na, double *pointeur_x_k)
{
/*
Resolution du systeme a quatre equations, quatre inconnues
avec f(x_oh,x_h3sio4,x_na,x_k)=0 electroneutralitée
     g(x_oh,x_h3sio4)=0 phases du CSH
     h_na(x_oh,x_na)=0 conservation Na
     k_k(x_oh,x_k)=0 conservation K
*/
  double x_oh_nouveau, x_h3sio4_nouveau, x_na_nouveau, x_k_nouveau ;
  double x_na = x_na_tot ;
  double x_k   = x_k_tot ;
  double x_oh = x_na + x_k +0.02 ;
  double x_co2_eq = x_co2 ;

  double x_h3sio4 = x_oh*pow(k_jen,1/y_jen)/pow(k_2,x_jen/y_jen) ;


  /* fonction f(x_oh,x_h3sio4,x_na,x_k) = A*x_oh^4+B*x_oh^3+ C*x_oh^2 + D*x_oh^1 + E*x_oh^0 
  + F*x_na*x_oh^2 + Gx_k*x_oh^2
  + K*x_na*x_oh^4
  + I*x_h3sio4*xoh^2 + H*x_h3sio4*x_oh^3 + J*x_h3sio4
 */
  double A = z_co3*x_co2_eq*k_co3*k_h/k_1 ;
  double B = z_oh + z_hco3*x_co2_eq*k_h/k_1 ;
  double C = 0 ;
  double D = z_h*k_e + z_caoh*k_caoh*k_2/k_e + z_cahco3*k_cahco3*k_h*k_2*x_co2_eq/k_1 ;
  double E = z_ca*k_2 ;
  double F = z_na ;
  double G = z_k ;
  double I   = z_h3sio4 ;
  double J  = z_cah3sio4*k_cah3sio4*k_2 ;
  double H  = z_h2sio4/(k_e*k_h2sio4) ;
  double K  = z_naco3*k_naco3*k_h*x_co2_eq/(k_e*k_1) ;

  double x_oh2= x_oh*x_oh ;
  double x_oh3= x_oh2*x_oh ;
  double x_oh4= x_oh2*x_oh2 ;
  double f ;
  double dfsdx_oh ;
  double dfsdx_h3sio4 ;
  double dfsdx_na ;
  double dfsdx_k ;

 
  /* fonction g(x_oh,x_h3sio4) = jen*x_oh^-y_jen*x_h3sio4^y_jen 
                               + tobII*x_oh^-y_tobII*x_h3sio4^y_tobII
                               + tobI*x_oh^-y_tobI*x_h3sio4^y_tobI
                               + sil*x_oh^-y_sil*x_h3sio4^y_sil - 1 */
  double jen = pow(k_2,x_jen)/k_jen ;
  double tobII = pow(k_2,x_tobII)/k_tobII  ;
  double tobI = pow(k_2,x_tobI)/k_tobI  ;
  double sil = k_h4sio4*k_e/(k_sil*k_h2sio4) ;
  double g  ;
  double dgsdx_oh ;
  double dgsdx_h3sio4  ;

  /* fonction h_na(x_oh,x_na) =  A_Na*x_na + B_Na*x_na*x_oh + C_Na*x_na*x_oh^2 -x_na_tot */
  double A_Na = 1 ;
  double B_Na = k_naoh/k_e + k_nahco3*k_h*x_co2_eq/k_1 ;
  double C_Na = k_naco3*k_h*x_co2_eq/(k_1*k_e) ;
  double h_na ;
  double dh_nasdx_oh ;
  double dh_nasdx_na ;

  /* fonction h_k(x_oh,x_k) =  A_K*x_k + B_K*x_k*x_oh - x_k_tot */
  double A_K = 1 ;
  double B_K = k_koh/k_e  ;
  double h_k ;
  double dh_ksdx_oh ;
  double dh_ksdx_k ;

  double a11,a12,a13,a14,a21,a22,a31,a33,a41,a44 ;
  
  double det ;
  
  double tolerance = 1e-3 ;
  double delta_oh  = 1 ;
  double delta_h3sio4 = 1 ;
  double delta_na  = 1 ;
  double delta_k  = 1 ;

  double f_jen,f_tobI,f_tobII,f_sil ;
 

  do
  {
  x_oh2= x_oh*x_oh ;
  x_oh3= x_oh2*x_oh ;
  x_oh4= x_oh2*x_oh2 ;

  f = A*x_oh4 + B*x_oh3+ C*x_oh2 + D*x_oh + E + F*x_na*x_oh2 + G*x_k*x_oh2 + K*x_na*x_oh4 + I*x_h3sio4*x_oh2 + H*x_h3sio4*x_oh3 + J*x_h3sio4 ;
  dfsdx_oh     	= 4*A*x_oh3 + 3*B*x_oh2 + 2*C*x_oh + D + 2*(F*x_na + G*x_k)*x_oh + 4*K*x_na*x_oh3 + 2*I*x_h3sio4*x_oh + 3*H*x_h3sio4*x_oh2 ;
  dfsdx_na     	= (F*x_oh2 + K*x_oh4) ;
  dfsdx_k     	= G*x_oh2 ;
  dfsdx_h3sio4 = I*x_oh2 + H*x_oh3 + J ;
   
    double x_hco3 = k_h*x_co2_eq*x_oh/k_1 ;
    double x_co3   = k_co3*x_oh*x_hco3 ;
    double x_h     = k_e/x_oh ;
    double x_ca    = k_ca/x_co3 ;
    double x_naco3 = k_naco3/k_e*x_na*x_oh*x_hco3 ;
    double x_caoh  = k_caoh/k_e*x_oh*x_ca ;
    double x_cahco3  = k_cahco3*x_hco3*x_ca ;
    double x_caoh2aq = k_caoh2*x_ca*x_oh*x_oh ;
    double x_h2sio4  = x_oh*x_h3sio4/(k_e*k_h2sio4) ;
    double x_cah3sio4 = k_cah3sio4*x_h3sio4*x_ca ;

    double f_bis = (z_h*x_h + z_oh*x_oh + z_ca*x_ca + z_hco3*x_hco3 + z_co3*x_co3 + z_na*x_na + z_naco3*x_naco3 + z_k*x_k + z_caoh*x_caoh + z_cahco3*x_cahco3 + z_cah3sio4*x_cah3sio4 + z_h2sio4*x_h2sio4 + z_h3sio4*x_h3sio4)*x_oh2 ;


  g = jen*pow(x_oh,-y_jen)*pow(x_h3sio4,y_jen) + tobII*pow(x_oh,-y_tobII)*pow(x_h3sio4,y_tobII) + tobI*pow(x_oh,-y_tobI)*pow(x_h3sio4,y_tobI) + sil*pow(x_oh,-y_sil)*pow(x_h3sio4,y_sil) - 1 ;
  dgsdx_oh     = - y_jen*jen*pow(x_oh,-y_jen-1)*pow(x_h3sio4,y_jen) - y_tobII*tobII*pow(x_oh,-y_tobII-1)*pow(x_h3sio4,y_tobII) -y_tobI*tobI*pow(x_oh,-y_tobI-1)*pow(x_h3sio4,y_tobI) - y_sil*sil*pow(x_oh,-y_sil-1)*pow(x_h3sio4,y_sil) ;
  dgsdx_h3sio4 = y_jen*jen*pow(x_oh,-y_jen)*pow(x_h3sio4,y_jen-1) + y_tobII*tobII*pow(x_oh,-y_tobII)*pow(x_h3sio4,y_tobII-1) + y_tobI*tobI*pow(x_oh,-y_tobI)*pow(x_h3sio4,y_tobI-1) + y_sil*sil*pow(x_oh,-y_sil)*pow(x_h3sio4,y_sil-1) ;

  f_jen = jen*pow(x_oh,-y_jen)*pow(x_h3sio4,y_jen) ;
  f_tobII = tobII*pow(x_oh,-y_tobII)*pow(x_h3sio4,y_tobII) ;
  f_tobI = tobI*pow(x_oh,-y_tobI)*pow(x_h3sio4,y_tobI) ;
  f_sil = sil*pow(x_oh,-y_sil)*pow(x_h3sio4,y_sil) ;

   h_na = A_Na*x_na + B_Na*x_na*x_oh + C_Na*x_na*x_oh2 - x_na_tot ;
  dh_nasdx_na = A_Na + B_Na*x_oh + C_Na*x_oh2 ;
  dh_nasdx_oh = B_Na*x_na + 2*C_Na*x_na*x_oh ;

   h_k = A_K*x_k + B_K*x_k*x_oh - x_k_tot ;
  dh_ksdx_k = A_K + B_K*x_oh ;
  dh_ksdx_oh = B_K*x_k ;

  a11 = dfsdx_oh ;
  a12 = dfsdx_h3sio4 ;
  a13 = dfsdx_na ;
  a14 = dfsdx_k ;

  a21 = dgsdx_oh ;
  a22 = dgsdx_h3sio4 ;

  a31 = dh_nasdx_oh ;
  a33 = dh_nasdx_na ;

  a41 = dh_ksdx_oh ;
  a44 = dh_ksdx_k ;
  
  det =  a11*a22*a33*a44-a12*a21*a33*a44-a13*a22*a31*a44-a14*a22*a33*a41 ;




  x_oh_nouveau     = x_oh - (-a13*a22*a44*h_na-a14*a22*a33*h_k-a12*a33*a44*g+a22*a33*a44*f)/det ;
  x_h3sio4_nouveau = x_h3sio4 - (a13*a21*a44*h_na+a14*a21*a33*h_k+((a11*a33-a13*a31)*a44-a14*a33*a41)*g-a21*a33*a44*f)/det ;
  x_na_nouveau = x_na - (((a11*a22-a12*a21)*a44-a14*a22*a41)*h_na+a14*a22*a31*h_k+a12*a31*a44*g-a22*a31*a44*f)/det ;
  x_k_nouveau = x_k - (a13*a22*a41*h_na+((a11*a22-a12*a21)*a33-a13*a22*a31)*h_k+a12*a33*a41*g-a22*a33*a41*f)/det ; 

  delta_oh     = fabs((x_oh_nouveau     - x_oh)/x_oh) ;
  delta_h3sio4 = fabs((x_h3sio4_nouveau - x_h3sio4)/x_h3sio4) ;
  delta_na     = fabs((x_na_nouveau     - x_na)/x_na) ;
  delta_k     = fabs((x_k_nouveau     - x_k)/x_k) ;
  
  x_oh = x_oh_nouveau ;
  x_h3sio4 = x_h3sio4_nouveau ;
  x_na = x_na_nouveau ;
  x_k = x_k_nouveau ;


  
  } while(delta_oh > tolerance && delta_h3sio4 > tolerance && delta_na > tolerance && delta_k > tolerance) ;
  
printf("\n electro= %e \n", f) ;
printf("CSH= %e \n", g) ;
printf("concer Na= %e \n", h_na) ;
printf("concer K= %e \n", h_k) ;
printf("\n x_oh= %e \n", x_oh) ;
printf("x_h3sio4= %e \n", x_h3sio4) ;
printf("x_na= %e \n", x_na) ;
printf("x_k= %e \n", x_k) ;
printf("f_jen= %e \n", f_jen) ;
printf("f_tobII= %e \n", f_tobII) ;
printf("f_tobI= %e \n", f_tobI) ;
printf("f_sil= %e \n", f_sil) ;
printf("Tolerances \n") ;
printf("delta_oh= %e \n", delta_oh) ;
printf("delta_h3sio4= %e \n", delta_h3sio4) ;
printf("delta_na= %e \n", delta_na) ;
printf("delta_k= %e \n", delta_k) ;


  *pointeur_x_oh = x_oh ;
  *pointeur_x_h3sio4 = x_h3sio4 ;
  *pointeur_x_na = x_na ;
  *pointeur_x_k = x_k ;
 
}


double force_ionique(double x_ca, double x_cahco3, double x_caoh, double x_k, double x_na, double x_h, double x_h3sio4, double x_naco3, double x_hco3, double x_oh, double x_h2sio4, double x_co3)
{
  double I=0.5*( z_ca*z_ca*x_ca    + z_h2sio4*z_h2sio4*x_h2sio4 + z_co3*z_co3*x_co3   + z_cahco3*z_cahco3*x_cahco3 + z_caoh*z_caoh*x_caoh + z_k*z_k*x_k + z_na*z_na*x_na + z_h*z_h*x_h    + z_h3sio4*z_h3sio4*x_h3sio4 + z_naco3*z_naco3*x_naco3 + z_hco3*z_hco3*x_hco3 + z_oh*z_oh*x_oh ) ;
  return(I);
}




void concentration_vers_activite(double x_ca, double x_cahco3, double x_caoh, double x_k, double x_na, double x_h, double x_h3sio4, double x_naco3, double x_hco3, double x_oh, double x_h2sio4, double x_co3,double *pointeur_x_ca,double *pointeur_x_cahco3, double *pointeur_x_caoh, double *pointeur_x_k, double *pointeur_x_na, double *pointeur_x_h, double *pointeur_x_h3sio4, double *pointeur_x_naco3, double *pointeur_x_hco3, double *pointeur_x_oh, double *pointeur_x_h2sio4, double *pointeur_x_co3)
{
     
  /*Force Ionique*/
 double I=force_ionique(x_ca,x_cahco3,x_caoh,x_k,  x_na, x_h,  x_h3sio4, x_naco3,  x_hco3, x_oh, x_h2sio4, x_co3);
 /*Loi de Davies*/
 double f_I=A_DAVIES*(sqrt(I)/(1+sqrt(I))-B_DAVIES*I); /*a 25degreC A_DAVIES=0,5 pour l'eau et B_DAVIES est pris entre 0.2 et 0.3*/
 /* Transformation des concentrations en activite*/
 x_ca=pow(10,(-f_I*z_ca*z_ca))*x_ca;
 x_cahco3=pow(10,(-f_I*z_cahco3*z_cahco3))*x_cahco3;
 x_caoh=pow(10,(-f_I*z_caoh*z_caoh))*x_caoh;
 x_k=pow(10,(-f_I*z_k*z_k))*x_k;
 x_na=pow(10,(-f_I*z_na*z_na))*x_na;
 x_h=pow(10,(-f_I*z_h*z_h))*x_h;
 x_h3sio4=pow(10,(-f_I*z_h3sio4*z_h3sio4))*x_h3sio4;
 x_naco3=pow(10,(-f_I*z_naco3*z_naco3))*x_naco3;
 x_hco3=pow(10,(-f_I*z_hco3*z_hco3))*x_hco3;
 x_oh=pow(10,(-f_I*z_oh*z_oh))*x_oh;
 x_h2sio4=pow(10,(-f_I*z_h2sio4*z_h2sio4))*x_h2sio4;
 x_co3=pow(10,(-f_I*z_co3*z_co3))*x_co3;
 
 
 /*Assignation des valeurs par le biais de pointeurs*/
 *pointeur_x_ca = x_ca;
 *pointeur_x_cahco3= x_cahco3;
 *pointeur_x_caoh=x_caoh;
 *pointeur_x_k=x_k;
 *pointeur_x_na=x_na;
 *pointeur_x_h=x_h;
 *pointeur_x_h3sio4=x_h3sio4;
 *pointeur_x_naco3=x_naco3;
 *pointeur_x_hco3=x_hco3;
 *pointeur_x_oh=x_oh;
 *pointeur_x_h2sio4=x_h2sio4;
 *pointeur_x_co3=x_co3;
 
  
  
}

double dn1_caoh2sdt(double av,double c_2)
{
#define AV         ((av < 1.) ? av : 1.)
  double rp = (av < 1.) ? pow(1 - av,1./3) : 0. ;
  double rc = pow(1 - AV + V_CC/V_CH*AV,1./3) ;
  /*
  double alpha2 = -1./3*av*av - 2./3*av + 1 ;
  double alpha  = -5.29478*av*av*av*av + 8.6069*av*av*av - 4.2444*av*av + 0.9325*av ;
  */

  return((rc > 0.) ? rp*rp/(1 + c_2*rp*(1 - rp/rc)) : 0.) ;
  /* return(alpha2/(1 + c_2*alpha)) ; */
#undef AV
}



double* ComputeAqueousParameters(Element_t *el,double **u,int n)
{
  double *x = var ;
  double x_co2      = X_CO2(n) ;
  double x_oh       = X_OH(n) ;
  double x_hco3     = X_HCO3(n) ;
  double x_na       = X_Na(n) ;
  double x_k        = X_K(n) ;
  
  double x_h        = k_e/x_oh ;
  
  double x_h2co3    = k_h*x_co2 ;
  double x_co3      = k_co3*x_oh*x_hco3 ;
  
  double x_ca       = k_ca/x_co3 ;
  double x_caoh     = k_caoh/k_e*x_oh*x_ca ;
  double x_cahco3   = k_cahco3*x_hco3*x_ca ;
  double x_caco3aq  = k_caco3*k_ca/(k_e*k_co3);
  double x_caoh2aq  = k_caoh2*x_ca*x_oh*x_oh ;
  
  double q_ch       = x_ca*x_oh*x_oh/k_2 ;
  double x_h4sio4   = IAP_SH(q_ch) ;
  double x_h3sio4   = (x_h4sio4*x_oh)/(k_e*k_h4sio4/k_h2sio4) ;
  double x_h2sio4   = x_oh*x_h3sio4/(k_e*k_h2sio4) ;
  double x_cah2sio4 = k_cah2sio4*x_h2sio4*x_ca ;
  double x_cah3sio4 = k_cah3sio4*x_h3sio4*x_ca ;
  
  double x_naoh     = k_naoh/k_e*x_na*x_oh ;
  double x_nahco3   = k_nahco3*x_na*x_hco3 ;
  double x_naco3    = k_naco3/k_e*x_na*x_oh*x_hco3 ;
  double x_koh      = k_koh/k_e*x_k*x_oh ;
    
  double x_h2o   = (1 - (x_h*v_h + x_oh*v_oh \
       + x_h2co3*v_h2co3 + x_hco3*v_hco3 + x_co3*v_co3 \
       + x_ca*v_ca + x_caoh*v_caoh + x_cahco3*v_cahco3 \
       + x_caco3aq*v_caco3aq + x_caoh2aq*v_caoh2aq \
       + x_h3sio4*v_h3sio4 + x_h4sio4*v_h4sio4 + x_h2sio4*v_h2sio4 \
       + x_cah2sio4*v_cah2sio4 + x_cah3sio4*v_cah3sio4 \
       + x_na*v_na + x_naoh*v_naoh + x_nahco3*v_nahco3 + x_naco3*v_naco3 \
       + x_k*v_k + x_koh*v_koh))/v_h2o ;
  
  double rho_l   = M_H*x_h + M_OH*x_oh + M_H2O*x_h2o \
       + M_H2CO3*x_h2co3 + M_HCO3*x_hco3 + M_CO3*x_co3 \
       + M_Ca*x_ca + M_CaOH*x_caoh + M_CaHCO3*x_cahco3 \
       + M_CaCO3aq*x_caco3aq + M_CaOH2aq*x_caoh2aq \
       + M_H3SiO4*x_h3sio4 + M_H4SiO4*x_h4sio4 + M_H2SiO4*x_h2sio4 \
       + M_CaH2SiO4*x_cah2sio4 + M_CaH3SiO4*x_cah3sio4 \
       + M_Na*x_na + M_NaOH*x_naoh + M_NaHCO3*x_nahco3 + M_NaCO3*x_naco3 \
       + M_K*x_k + M_KOH*x_koh ;
  
  double x_q = z_h*x_h + z_oh*x_oh \
       + z_hco3*x_hco3 + z_co3*x_co3 \
       + z_ca*x_ca + z_caoh*x_caoh + z_cahco3*x_cahco3 \
       + z_h3sio4*x_h3sio4 + z_h2sio4*x_h2sio4 \
       + z_cah3sio4*x_cah3sio4 \
       + z_na*x_na + z_naco3*x_naco3 \
       + z_k*x_k ;
  
  x[I_H       ] = x_h ;
  x[I_OH      ] = x_oh ;
  x[I_H2O     ] = x_h2o ;
  
  x[I_CO2     ] = x_co2 ;
  x[I_HCO3    ] = x_hco3 ;
  x[I_H2CO3   ] = x_h2co3 ;
  x[I_CO3     ] = x_co3 ;
  
  x[I_Ca      ] = x_ca ;
  x[I_CaOH    ] = x_caoh ;
  x[I_CaHCO3  ] = x_cahco3 ;
  x[I_CaCO3aq ] = x_caco3aq ;
  x[I_CaOH2aq ] = x_caoh2aq ;
  
  x[I_H4SiO4  ] = x_h4sio4 ;
  x[I_H3SiO4  ] = x_h3sio4 ;
  x[I_H2SiO4  ] = x_h2sio4 ;
  x[I_CaH2SiO4] = x_cah2sio4 ;
  x[I_CaH3SiO4] = x_cah3sio4 ;
  
  x[I_Na      ] = x_na ;
  x[I_NaOH    ] = x_naoh ;
  x[I_NaHCO3  ] = x_nahco3 ;
  x[I_NaCO3   ] = x_naco3 ;
  
  x[I_K       ] = x_k ;
  x[I_KOH     ] = x_koh ;
  
  x[I_Q_CH    ] = q_ch ;
  
  x[I_RHO_L   ] = rho_l ;
  x[I_N_Q     ] = x_q ;
  
  x[I_C_L     ] = x_h2co3 + x_hco3 + x_co3 + x_nahco3 + x_naco3 + x_cahco3 + x_caco3aq ;
  x[I_Ca_L    ] = x_ca + x_caoh + x_cahco3 + x_caco3aq + x_caoh2aq + x_cah2sio4 + x_cah3sio4 ;
  x[I_Na_L    ] = x_na + x_naoh + x_nahco3 + x_naco3 ;
  x[I_K_L     ] = x_k + x_koh ;
  x[I_Si_L    ] = x_h2sio4 + x_h3sio4 + x_h4sio4 + x_cah2sio4 + x_cah3sio4 ;
  
  return(x) ;
}

