/*
  (ancien modele 66)
  4 ions : Cl-, Na+, K+, OH-, H+
  Isotherme d'interaction Chlorure - CSH de type Langmuir
  Une autre partie des chlores est fixee instantanement
  dans les sels de Friedel proportionnellement au C3A.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"

#define MODELINDEX  43
#define TITLE "Chlorures dans les betons satures (version simplifiee)"
#define AUTHORS "Nguyen"

#include "OldMethods.h"

/* Macros */

#define NEQ     (4)

#define NVI     (12)
#define NVE     (10)

#define E_Cl    (0)
#define E_q     (1)
#define E_Na    (2)
#define E_K     (3)

#define I_c_cl  (0)
#define I_c_k   (3)
#define I_c_na  (2)
#define I_psi   (1)


#define C_Cl(n)     (u[(n)][I_c_cl])
#define C_K(n)      (u[(n)][I_c_k])
#define C_Na(n)     (u[(n)][I_c_na])
#define PSI(n)      (u[(n)][I_psi])

#define N_Cl(n)     (f[(n)])
#define N_Na(n)     (f[(n+2)])
#define N_K(n)      (f[(n+4)])

#define N_Cln(n)    (f_n[(n)])
#define N_Nan(n)    (f_n[(n+2)])
#define N_Kn(n)     (f_n[(n+4)])

#define W_Cl        (f[8])
#define W_q         (f[9])
#define W_Na        (f[10])
#define W_K         (f[11])


#define KF_Cl       (va[(0)])
#define KF_OH       (va[(1)])
#define KF_Na       (va[(2)])
#define KF_K        (va[(3)])
#define KF_H        (va[(4)])

#define KM_Cl       (va[(5)])
#define KM_OH       (va[(6)])
#define KM_Na       (va[(7)])
#define KM_K        (va[(8)])
#define KM_H        (va[(9)])


/* valences */
#define z_cl    (-1.)
#define z_oh    (-1.)
#define z_na    (1.)
#define z_k     (1.)
#define z_h     (1.)

/* coefficients de diffusion dans l'eau dilue*/
#define do_cl   (2.032e-9)
#define do_oh   (5.273e-9)
#define do_na   (1.334e-9)
#define do_k    (1.957e-9)
#define do_h    (9.310e-9)

/* constantes d'equilibre */
#define k_e     (1.e-8)

/* constantes physiques */
#define F   	(9.64846e4)  /* Faraday (C/mole) */
#define RT      (2436.)      /* produit de R et T */
                             /* R = 8.3143 J/K/mol cste des gaz parfait */
                             /* T = 293 K */

/* Densite initiale de radicaux hydroxyles (valeur arbitraire) */
#define s_oh0   (0.)


/* Fonctions */
static int    pm(const char *s) ;
static double concentration_oh(double,double) ;
static double dconcentration_oh(double,double) ;

/* Parametres */
static double phi,d_cl,r_d ;
static double s_csh,s_c3a,alpha,beta ;

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0) return (0) ;
  else if(strcmp(s,"D_Cl") == 0) return (1);
  else if(strcmp(s,"r_d") == 0) return (2);
  else if(strcmp(s,"s_c3aeq") == 0) return (3) ;
  else if(strcmp(s,"s_csh") == 0) return (4) ;
  else if(strcmp(s,"alpha") == 0) return (5);
  else if(strcmp(s,"beta") == 0) return (6); 
  else {
    printf("donnee \"%s\" non connue (pm43)\n",s) ; exit(0) ;
  }
}

int dm43(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 7 ;

  if(dim > 1) arret("dm43 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_Cl],"E_Cl") ;
  strcpy(mat->eqn[E_q],"E_q") ;
  strcpy(mat->eqn[E_Na],"E_Na") ;
  strcpy(mat->eqn[E_K],"E_K") ;

  strcpy(mat->inc[I_c_cl],"c_cl") ;
  strcpy(mat->inc[I_c_k],"c_k") ;
  strcpy(mat->inc[I_c_na],"c_na") ;
  strcpy(mat->inc[I_psi],"psi") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}


int qm43(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 4 equations :\n\
\t 1. Conservation de la masse de Cl (c_cl)\n\
\t 3. Conservation de la masse de Na (c_na)\n\
\t 4. Conservation de la masse de K  (c_k)\n\
\t 2. Conservation de la charge      (psi)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = 0.121 # Porosite\n") ;
  fprintf(ficd,"D_Cl = 2.6e-12   # Diffusion effective de Cl\n") ;
  fprintf(ficd,"r_d = 1.         # Rapport des tortuosites des anions et des cations\n") ;
  fprintf(ficd,"s_c3aeq = 8.43   # Contenu en C3A equivalents\n") ;
  fprintf(ficd,"s_csh = 635.     # Contenu en CSH\n") ;
  fprintf(ficd,"alpha = 0.12     # Coef. de l\'isotherme\n") ;
  fprintf(ficd,"beta = 2.66      # Coef. de l\'isotherme\n") ;

  return(NEQ) ;
}

void tb43(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ; /* implicite */
  el->n_ve = NVE ; /* explicite */
}

void ch43(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
#define R(n,i)     (r[(n)*NEQ+(i)])

  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;

#undef R
}

void in43(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double c_cl,c_oh,c_na,c_k,c_h,C_OH[2] ;
  double grd_cl,grd_oh,grd_na,grd_k,grd_h,grd_psi ;
  double w_oh,w_h ;
  double s_cl ;
  double dx,ton,top ;
  int    i ;
  double deux = 2. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ;

  s_c3a   = el.mat->pr[pm("s_c3aeq")] ;
  s_csh   = el.mat->pr[pm("s_csh")] ;
  alpha   = el.mat->pr[pm("alpha")] ;
  beta    = el.mat->pr[pm("beta")] ;
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    c_cl    = C_Cl(i) ;
    c_k     = C_K(i) ;
    c_na    = C_Na(i) ;
    c_oh    = concentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
    c_h	    = k_e/c_oh ;
    C_OH[i] = c_oh ;
    
    s_cl     = alpha*s_csh*c_cl*beta/(c_oh+beta*c_cl) + 2*s_c3a ;

    N_Cl(i)  = phi*c_cl + s_cl ;
    N_Na(i)  = phi*c_na ;
    N_K(i)   = phi*c_k ;
  }

  /* Coefficient de transfert */
  c_cl     = (C_Cl(0) + C_Cl(1))/deux ;
  c_k      = (C_K(0)  + C_K(1))/deux ;
  c_na     = (C_Na(0) + C_Na(1))/deux ;
  c_oh     = concentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
  c_h      =  k_e/c_oh ;
  
  ton = d_cl/do_cl ;
  top = ton/r_d ;
  
  KF_Cl  = ton*do_cl;
  KF_OH  = ton*do_oh ;
  KF_Na  = top*do_na;
  KF_K   = top*do_k;
  KF_H   = top*do_h;
  
  KM_Cl  = z_cl*KF_Cl*c_cl*F/RT;	
  KM_OH  = z_oh*KF_OH*c_oh*F/RT;
  KM_Na  = z_na*KF_Na*c_na*F/RT;
  KM_K   = z_k*KF_K*c_k*F/RT;	
  KM_H   = z_h*KF_H*c_h*F/RT;	

  /* Gradients */
  dx      = x[1][0] - x[0][0] ;
  
  grd_cl  = (C_Cl(1) - C_Cl(0))/dx ;
  grd_k   = (C_K(1)  - C_K(0))/dx ;
  grd_na  = (C_Na(1) - C_Na(0))/dx ;
  grd_oh  = (C_OH[1] - C_OH[0])/dx ;
  grd_h   = (k_e/C_OH[1] - k_e/C_OH[0])/dx ;
  grd_psi = (PSI(1) - PSI(0))/dx ;
  
  /* Flux */
  w_oh   = - KF_OH*grd_oh - KM_OH*grd_psi ;
  w_h    = - KF_H*grd_h   - KM_H*grd_psi ;

  W_Cl   = - KF_Cl*grd_cl - KM_Cl*grd_psi ;
  W_Na   = - KF_Na*grd_na - KM_Na*grd_psi ;
  W_K    = - KF_K*grd_k   - KM_K*grd_psi ;
  W_q    = z_cl*W_Cl + z_na*W_Na + z_k*W_K + z_oh*w_oh + z_h*w_h ;
}


int ex43(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double c_cl,c_oh,c_na,c_k,c_h ;
  double ton,top ;
  double deux = 2. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ;

  s_c3a   = el.mat->pr[pm("s_c3aeq")] ;
  s_csh   = el.mat->pr[pm("s_csh")] ;
  alpha   = el.mat->pr[pm("alpha")] ;
  beta    = el.mat->pr[pm("beta")] ;
  
  c_cl     = (C_Cl(0) + C_Cl(1))/deux ;
  c_k      = (C_K(0)  + C_K(1))/deux ;
  c_na     = (C_Na(0) + C_Na(1))/deux ;
  c_oh     = concentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
  c_h      =  k_e/c_oh ;
  
  ton = d_cl/do_cl ;
  top = ton/r_d ;
  
  KF_Cl  = ton*do_cl ;
  KF_OH  = ton*do_oh ;
  KF_Na  = top*do_na ;
  KF_K   = top*do_k ;
  KF_H   = top*do_h ;
  
  KM_Cl  = z_cl*KF_Cl*c_cl*F/RT ;	
  KM_OH  = z_oh*KF_OH*c_oh*F/RT ;
  KM_Na  = z_na*KF_Na*c_na*F/RT ;
  KM_K   = z_k*KF_K*c_k*F/RT ;	
  KM_H   = z_h*KF_H*c_h*F/RT ;

  return(0) ;
}

int ct43(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double c_cl, c_oh, c_na, c_k, c_h, C_OH[2] ;
  double grd_cl,grd_oh, grd_na, grd_k,grd_h,grd_psi ;
  double w_oh,w_h ;
  double s_cl ;
  double dx ;
  int    i ;

  
  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ;

  s_c3a    = el.mat->pr[pm("s_c3aeq")] ;
  s_csh    = el.mat->pr[pm("s_csh")] ;
  alpha    = el.mat->pr[pm("alpha")] ;
  beta     = el.mat->pr[pm("beta")] ;
   
  /* Contenus molaires */

 for(i=0;i<2;i++) {
    c_cl     = C_Cl(i) ;
    c_k      = C_K(i) ;
    c_na     = C_Na(i) ;
    c_oh     = concentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
    c_h	     = k_e/c_oh ;
    C_OH[i]  = c_oh ;
    /*
      if(c_cl <0. || c_oh <= 0. || c_na <0. || c_k < 0.|| c_h <= 0.) {
      printf("\n\en x  = %e\n\
      c_cl    = %e\n\
      c_oh    = %e\n\
      c_na    = %e\n\
      c_k     = %e\n",x[i][0],c_cl,c_oh,c_na,c_k) ;
      return(-1) ;
      }
    */
    
    /* clorures fixes sur les CSH et les sels de Friedel */
    s_cl     = alpha*s_csh*c_cl*beta/(c_oh+beta*c_cl) + 2*s_c3a ;
   
    N_Cl(i)  = phi*c_cl + s_cl ;
    N_Na(i)  = phi*c_na ;
    N_K(i)   = phi*c_k ;
 }

 /* Gradients */	
 dx      = x[1][0] - x[0][0] ;
 
 grd_cl  = (C_Cl(1) - C_Cl(0))/dx ;
 grd_k   = (C_K(1)  - C_K(0))/dx ;
 grd_na  = (C_Na(1) - C_Na(0))/dx ;
 grd_oh  = (C_OH[1] - C_OH[0])/dx ;
 grd_h   = (k_e/C_OH[1] - k_e/C_OH[0])/dx ;
 grd_psi = (PSI(1) - PSI(0))/dx ;

 /* Flux */
 w_oh   = - KF_OH*grd_oh - KM_OH*grd_psi ;
 w_h    = - KF_H*grd_h   - KM_H*grd_psi ;

 W_Cl   = - KF_Cl*grd_cl - KM_Cl*grd_psi ;
 W_Na   = - KF_Na*grd_na - KM_Na*grd_psi ;
 W_K    = - KF_K*grd_k   - KM_K*grd_psi ;
 W_q    = z_cl*W_Cl + z_na*W_Na + z_k*W_K + z_oh*w_oh + z_h*w_h ;
 
 return(0);
}

int mx43(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])
  
  double c_cl,c_oh,c_h,c_na,c_k,s_cl,s_clcsh,s_clcshsc_cl ;
  double dc_ohsdq,dc_ohsdc_cl[2],dc_ohsdc_na[2],dc_ohsdc_k[2] ;
  double dc_hsdc_cl[2],dc_hsdc_na[2],dc_hsdc_k[2] ;
  double ds_clsdc_cl,ds_clsdc_na,ds_clsdc_k,ds_clsdc_oh ;
  double tr ;
  double trf_cl,trf_oh,trf_h,trf_na,trf_k ;
  double tre_cl,tre_oh,tre_h,tre_na,tre_k ;
  double c[2] ;
  double dx,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
 
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;

  s_c3a    = el.mat->pr[pm("s_c3aeq")] ;
  s_csh    = el.mat->pr[pm("s_csh")] ;
  alpha    = el.mat->pr[pm("alpha")] ;
  beta     = el.mat->pr[pm("beta")] ;
  
  /*
    CALCUL DE volume ET DE surf 
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;

  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    c_cl     = C_Cl(i) ;
    c_k      = C_K(i) ;
    c_na     = C_Na(i) ;
    c_oh     = concentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
    c_h      = k_e/c_oh ;

    /* clorures fixes sur les CSH et les sels de Friedel */
    s_cl     = alpha*s_csh*c_cl*beta/(c_oh + beta*c_cl) + 2*s_c3a ;
    s_clcshsc_cl = alpha*s_csh*beta/(c_oh + beta*c_cl) ;
    s_clcsh      = s_clcshsc_cl*c_cl ;
    ds_clsdc_oh  = - s_clcsh/(c_oh + beta*c_cl) ;

    /* derivees */
    dc_ohsdq       = dconcentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
    dc_ohsdc_cl[i] = z_cl*dc_ohsdq ;
    dc_ohsdc_na[i] = z_na*dc_ohsdq ;
    dc_ohsdc_k[i]  = z_k*dc_ohsdq ;

    dc_hsdc_cl[i]  = -c_h/c_oh*dc_ohsdc_cl[i] ;
    dc_hsdc_na[i]  = -c_h/c_oh*dc_ohsdc_na[i] ;
    dc_hsdc_k[i]   = -c_h/c_oh*dc_ohsdc_k[i] ;

    ds_clsdc_cl =  s_clcshsc_cl/(c_oh + beta*c_cl)*c_oh + ds_clsdc_oh*dc_ohsdc_cl[i] ;
    ds_clsdc_na = ds_clsdc_oh*dc_ohsdc_na[i] ;
    ds_clsdc_k  = ds_clsdc_oh*dc_ohsdc_k[i] ;

    /*
      Conservation de Cl : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
    */
    K(i*NEQ+E_Cl,i*NEQ+I_c_cl) += volume[i]*(phi + ds_clsdc_cl) ;
    K(i*NEQ+E_Cl,i*NEQ+I_c_na) += volume[i]*(ds_clsdc_na) ;
    K(i*NEQ+E_Cl,i*NEQ+I_c_k)  += volume[i]*(ds_clsdc_k) ;
    /*
      Conservation de la charge
    */
    /*
      Conservation de Na : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
    */
    K(i*NEQ+E_Na,i*NEQ+I_c_na) +=  volume[i]*(phi) ;
    /*
      Conservation de K : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
    K(i*NEQ+E_K,i*NEQ+I_c_k)   += volume[i]*(phi) ;
    }

  /* termes d'ecoulement */
  tr     = dt*surf/dx ;

  trf_cl    = tr*KF_Cl ;
  trf_oh    = tr*KF_OH ;
  trf_h     = tr*KF_H ;
  trf_na    = tr*KF_Na ;
  trf_k     = tr*KF_K ;

  tre_cl    = tr*KM_Cl ;
  tre_oh    = tr*KM_OH ;
  tre_h     = tr*KM_H ;
  tre_na    = tr*KM_Na ;
  tre_k     = tr*KM_K ;
   
 /*
    Conservation de Cl (chlore) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
      
  K(E_Cl,I_c_cl)             += + trf_cl ;
  K(E_Cl,I_c_cl+NEQ)         += - trf_cl ;
  K(E_Cl+NEQ,I_c_cl)         += - trf_cl ;
  K(E_Cl+NEQ,I_c_cl+NEQ)     += + trf_cl ;

  K(E_Cl,I_psi)              += + tre_cl ;
  K(E_Cl,I_psi+NEQ)          += - tre_cl ;
  K(E_Cl+NEQ,I_psi)          += - tre_cl ;
  K(E_Cl+NEQ,I_psi+NEQ)      += + tre_cl ;

  /*
    Conservation de la charge
  */
  for(i=0;i<2;i++) {
    c[i] = z_cl*trf_cl + z_oh*trf_oh*dc_ohsdc_cl[i] + z_h*trf_h*dc_hsdc_cl[i] ;
  }
  K(E_q,I_c_cl)              += + c[0] ;
  K(E_q,I_c_cl+NEQ)          += - c[1] ;
  K(E_q+NEQ,I_c_cl)          += - c[0] ;
  K(E_q+NEQ,I_c_cl+NEQ)      += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = z_na*trf_na + z_oh*trf_oh*dc_ohsdc_na[i] + z_h*trf_h*dc_hsdc_na[i] ;
  }
  K(E_q,I_c_na)              += + c[0] ;
  K(E_q,I_c_na+NEQ)          += - c[1] ;
  K(E_q+NEQ,I_c_na)          += - c[0] ;
  K(E_q+NEQ,I_c_na+NEQ)      += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = z_k*trf_k + z_oh*trf_oh*dc_ohsdc_k[i] + z_h*trf_h*dc_hsdc_k[i] ;
  }
  K(E_q,I_c_k)               += + c[0] ;
  K(E_q,I_c_k+NEQ)           += - c[1] ;
  K(E_q+NEQ,I_c_k)           += - c[0] ;
  K(E_q+NEQ,I_c_k+NEQ)       += + c[1] ;

  for(i=0;i<2;i++) {
    c[i] = z_cl*tre_cl + z_na*tre_na + z_k*tre_k + z_oh*tre_oh + z_h*tre_h ;
  }
  K(E_q,I_psi)               += + c[0] ;
  K(E_q,I_psi+NEQ)           += - c[1] ;
  K(E_q+NEQ,I_psi)           += - c[0] ;
  K(E_q+NEQ,I_psi+NEQ)       += + c[1] ;

  /*
    Conservation de Na : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  K(E_Na,I_c_na)             += + trf_na ;
  K(E_Na,I_c_na+NEQ)         += - trf_na ;
  K(E_Na+NEQ,I_c_na)         += - trf_na ;
  K(E_Na+NEQ,I_c_na+NEQ)     += + trf_na ;

  K(E_Na,I_psi)              += + tre_na ;
  K(E_Na,I_psi+NEQ)          += - tre_na ;
  K(E_Na+NEQ,I_psi)          += - tre_na ;
  K(E_Na+NEQ,I_psi+NEQ)      += + tre_na ; 

  /*
    Conservation de K : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  K(E_K,I_c_k)               += + trf_k ;
  K(E_K,I_c_k+NEQ)           += - trf_k ;
  K(E_K+NEQ,I_c_k)           += - trf_k ;
  K(E_K+NEQ,I_c_k+NEQ)       += + trf_k ;

  K(E_K,I_psi)               += + tre_k ;
  K(E_K,I_psi+NEQ)           += - tre_k ;
  K(E_K+NEQ,I_psi)           += - tre_k ;
  K(E_K+NEQ,I_psi+NEQ)       += + tre_k ;

  return(0) ;
}

void rs43(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])

  double dx ,xm ;
  double volume[2],surf ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;

  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++) {
    volume[i] = fabs(dx)/deux ; 
    if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
  }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  /*
    Conservation de Cl (chlore) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*W_Cl ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*W_Cl ;
  /*
      Conservation de la charge :  dt * div(w_q) = 0
    */
  R(0,E_q)  -= + dt*surf*W_q ;
  R(1,E_q)  -= - dt*surf*W_q ;
  /*
    Conservation de Na (natrie) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  R(0,E_Na) -= volume[0]*(N_Na(0) - N_Nan(0)) + dt*surf*W_Na ;
  R(1,E_Na) -= volume[1]*(N_Na(1) - N_Nan(1)) - dt*surf*W_Na ;
 /*
    Conservation de K (kali) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K)  -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*W_K ;
  R(1,E_K)  -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*W_K ;

#undef R
}

int so43(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;
  
  /*
    Donnees
  */
  phi      = el.mat->pr[pm("porosite")] ;
  d_cl     = el.mat->pr[pm("D_Cl")] ;
  r_d	   = el.mat->pr[pm("r_d")] ;

  s_c3a    = el.mat->pr[pm("s_c3aeq")] ;
  s_csh    = el.mat->pr[pm("s_csh")] ;
  alpha    = el.mat->pr[pm("alpha")] ;
  beta     = el.mat->pr[pm("beta")] ;

  
  /* initialisation */
  nso = 10;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  {
    /* concentration */
    double c_cl    = param(u,h_s,el.nn,I_c_cl) ;
    double c_na    = param(u,h_s,el.nn,I_c_na) ;
    double c_k     = param(u,h_s,el.nn,I_c_k) ;  
    
    double psi     = param(u,h_s,el.nn,I_psi) ;
    
    double c_oh    = concentration_oh(z_cl*c_cl + z_na*c_na + z_k*c_k,k_e) ;
    double c_h     = k_e/c_oh ;
    
    double s_cl    = alpha*s_csh*c_cl*beta/(c_oh + beta*c_cl) + 2*s_c3a ;
    
    double dx      = x[1][0] - x[0][0] ;
    double grd_psi = (PSI(1) - PSI(0))/dx ;
    
    /* quantites exploitees */
    strcpy(r[0].text,"Cl libre") ; r[0].n = 1 ;
    r[0].v[0] = c_cl ;
    strcpy(r[1].text,"Cl totaux") ; r[1].n = 1 ;
    r[1].v[0] = s_cl + phi*c_cl; 
    strcpy(r[2].text,"c_OH") ; r[2].n = 1 ;
    r[2].v[0] = c_oh ;
    strcpy(r[3].text,"c_Na") ; r[3].n = 1 ;
    r[3].v[0] = c_na ;
    strcpy(r[4].text,"c_K") ; r[4].n = 1 ;
    r[4].v[0] = c_k ;
    strcpy(r[5].text,"c_H") ; r[5].n = 1 ;
    r[5].v[0] = c_h ;
    strcpy(r[6].text,"potentiel electrique") ; r[6].n = 1 ;
    r[6].v[0] = psi ;
    strcpy(r[7].text,"Flux de Cl") ; r[7].n = 1 ;
    r[7].v[0] = W_Cl ;
    strcpy(r[8].text,"Flux de diffusion") ; r[8].n = 1 ;
    r[8].v[0] =  W_Cl + KM_Cl*grd_psi ;
    strcpy(r[9].text,"Flux de Migration") ; r[9].n = 1 ;
    r[9].v[0] = -KM_Cl*grd_psi ;
  }

  return(nso) ;
}

double concentration_oh(double q,double k_eau)
{
  return (0.5*(q + sqrt(q*q + 4*k_eau))) ;
}

double dconcentration_oh(double q,double k_eau)
{
  return (0.5*(1. + q/sqrt(q*q + 4*k_eau))) ;
}
