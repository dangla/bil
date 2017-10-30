/* Carbonatation d'un lit de concassés de béton */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  48
#define TITLE "Carbonatation d'un lit de béton concasse"
#define AUTHORS "Thiery"

#include "OldMethods.h"

/* Macros */
#define NEQ      (2)
#define E_co2    (0)
#define E_chimie (1)
#define I_x_co2  (0)
#define I_av     (1)
#define NVE      (0)
#define NVI      (7)
#define Puiss    (0.)
/* Fonctions */
static int    pm(char *s) ;

/* Parametres */
static double phi,S,R_0,a,b,phi_be,S_be,a_be,b_be,s_CH,k_henry,D_0 ;

int pm(char *s)
{
  if(strcmp(s,"phi") == 0) return (0) ;
  else if(strcmp(s,"S") == 0) return (1) ;
  else if(strcmp(s,"R_0") == 0) return (2) ;
  else if(strcmp(s,"a") == 0) return (3) ;
  else if(strcmp(s,"b") == 0) return (4) ;
  else if(strcmp(s,"phi_be") == 0) return (5) ;
  else if(strcmp(s,"S_be") == 0) return (6) ;
  else if(strcmp(s,"a_be") == 0) return (7) ;
  else if(strcmp(s,"b_be") == 0) return (8) ;
  else if(strcmp(s,"s_CH") == 0) return (9) ;
  else if(strcmp(s,"k_henry") == 0) return (10) ;
  else if(strcmp(s,"D_0") == 0) return (11) ;
  else
    { printf("donnee \"%s\" non connue (pm48)\n",s) ; exit(0) ; }
}

int dm48(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 12 ;

  if(dim > 1) arret("dm48 : dimension > 1 non prevue") ;
  
  mat->neq    = NEQ ;
  strcpy(mat->eqn[E_co2],"co2") ;
  strcpy(mat->eqn[E_chimie],"chimie") ;
  strcpy(mat->inc[I_x_co2],"x_co2") ;
  strcpy(mat->inc[I_av],"av") ;

  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm48(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
 Le systeme d\'equations est forme de 2 equations :\n\
\t 1. Conservation de la quantite de matiere de CO2 (x_co2) \n\
\t 2. Avancement de la reaction chimique de carbonatation (av) \n\
Exemple de donnees :\n") ;

  fprintf(ficd,"phi    = 0.5        # Porosite du milieu granulaire\n") ;
  fprintf(ficd,"S      = 0.5        # Taux de saturation du milieu granulaire\n") ;
  fprintf(ficd,"R_0    = 1e-2       # Rayon des granulats [m]\n") ;
  fprintf(ficd,"a      = 4/3        # Coeff. Millington sur la porosite (milieu granulaire)\n") ;
  fprintf(ficd,"b      = 10/3       # Coeff. Millington sur la saturation (milieu granulaire)\n") ;
  fprintf(ficd,"phi_be = 0.12       # Porosite du beton\n") ;
  fprintf(ficd,"S_be   = 0.5        # Taux de saturation du beton\n") ;
  fprintf(ficd,"a_be   = 2.74       # Coeff. Millington sur la porosite (beton)\n") ;
  fprintf(ficd,"b_be   = 4.20       # Coeff. Millington sur la saturation (beton)\n") ;
  fprintf(ficd,"s_CH   = 1.5        # Teneur en CH (mol.m-3 de beton)\n") ;
  fprintf(ficd,"k_henry= 1          # Constante de Henry\n") ;
  fprintf(ficd,"D_0    = 1e-5       # Coefficient de diffusion hors milieu poreux\n") ;
  
  return(NEQ) ;
}

void tb48(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ;
  el->n_ve = NVE ;
}

void ch48(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int   i ;
  
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}

void in48(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define X_co2(i)     (u[(i)][I_x_co2])
#define Av(i)        (u[(i)][I_av])
#define N_co2(i)     (f[(i)])
#define N_chimie(i)  (f[(i+2)])
#define Xsi(i)       (f[(i+4)])
#define W_co2        (f[(6)])

  double x_co2, n_co2_g, n_co2_l, n_co2_f, d_co2, D_co2, x_R, av, xsi, dx ;
  int    i ;


  if(el.dim < dim) return ;
  
  /*
    Donnees
  */
  phi      = el.mat->pr[pm("phi")] ;
  phi_be   = el.mat->pr[pm("phi_be")] ;
  S        = el.mat->pr[pm("S")] ;
  S_be     = el.mat->pr[pm("S_be")] ;
  a        = el.mat->pr[pm("a")] ;
  b        = el.mat->pr[pm("b")] ;
  a_be     = el.mat->pr[pm("a_be")] ;
  b_be     = el.mat->pr[pm("b_be")] ;
  R_0      = el.mat->pr[pm("R_0")] ;
  k_henry  = el.mat->pr[pm("k_henry")] ;
  D_0      = el.mat->pr[pm("D_0")] ;
  s_CH     = el.mat->pr[pm("s_CH")] ;
  
  d_co2       = D_0*pow(phi_be,a_be)*pow(1 - S_be,b_be) ; /* coeff. diff. effectif dans le beton */	
 	    for(i=0;i<2;i++)
	    {
	      x_co2       = X_co2(i) ;
	      av          = Av(i) ;
	      n_co2_g     = phi*(1 - S)*x_co2 ;   /* CO2 dans la phase gazeuse */
	      n_co2_l     = phi*S*k_henry*x_co2 ; /* CO2 en solution à la surface des grains concasses */
	      n_co2_f     = (1 - phi)*s_CH*av ;   /* CO2 fixé dans les concasses */
	      x_R         = pow(1 - av,1./3.) ;
	      xsi         = 3*(d_co2*k_henry*S*x_co2/R_0*x_R/(1 - x_R))*(1 - phi)/R_0*pow(1 - av,Puiss) ; /* Cinetique */ 
	      Xsi(i)      = xsi ;
	      N_co2(i)    = n_co2_g + n_co2_l + n_co2_f ;
	      N_chimie(i) = n_co2_f ;
	    }
  dx      = x[1][0] - x[0][0] ;
  D_co2   = D_0*pow(phi,a)*pow(1 - S,b) ; /* coeff. diff. effectif à travers le milieu granulaire */
  W_co2   = -D_co2*(X_co2(1) - X_co2(0))/dx ;

  
#undef X_co2
#undef Av
#undef N_co2
#undef N_chimie
#undef W_co2
}

int ex48(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  return(0) ;
}

int ct48(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define X_co21(i)     (u_1[(i)][I_x_co2])
#define Av1(i)        (u_1[(i)][I_av])
#define Avn(i)        (u_n[(i)][I_av])
#define N_co21(i)     (f_1[(i)])
#define N_chimie1(i)  (f_1[(i+2)])
#define Xsi1(i)       (f_1[(i+4)])
#define W_co21        (f_1[(6)])
  
  double x_co2, n_co2_g, n_co2_l, n_co2_f, d_co2, D_co2, x_R, av, xsi, dx ;
  int    i ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi      = el.mat->pr[pm("phi")] ;
  phi_be   = el.mat->pr[pm("phi_be")] ;
  S        = el.mat->pr[pm("S")] ;
  S_be     = el.mat->pr[pm("S_be")] ;
  a        = el.mat->pr[pm("a")] ;
  b        = el.mat->pr[pm("b")] ;
  a_be     = el.mat->pr[pm("a_be")] ;
  b_be     = el.mat->pr[pm("b_be")] ;
  R_0      = el.mat->pr[pm("R_0")] ;
  k_henry  = el.mat->pr[pm("k_henry")] ;
  D_0      = el.mat->pr[pm("D_0")] ;
  s_CH     = el.mat->pr[pm("s_CH")] ;
  
  d_co2       = D_0*pow(phi_be,a_be)*pow(1 - S_be,b_be) ;	
  for(i=0;i<2;i++)
    {
      x_co2       = X_co21(i) ;
      av          = Avn(i) ;
      
      /* if(x_co2 < 0. || av <= 0. || av >= 1.) {
	 printf("\n\
	 en x    = %e\n				\
	 x_co2   = %e\n				\
	 av      = %e\n",x[i][0],x_co2,av) ;
	 return(1) ;
	 }*/
      
      n_co2_g     = phi*(1 - S)*x_co2 ;
      n_co2_l     = phi*S*k_henry*x_co2 ;
      n_co2_f     = (1 - phi)*s_CH*Av1(i) ;
      
      x_R         = pow(1 - av,1./3.) ;
      xsi         = 3*((d_co2*k_henry*S*x_co2/R_0)*x_R/(1 - x_R))*(1 - phi)/R_0*pow(1 - av,Puiss) ;
      if(av >= 1.) {xsi = 0. ;} ;
      Xsi1(i)     = xsi ;
      N_co21(i)   = n_co2_g + n_co2_l + n_co2_f ;
      N_chimie1(i)= n_co2_f ;
    }
  dx       = x[1][0] - x[0][0] ;
  D_co2    = D_0*pow(phi,a)*pow(1 - S,b) ;
  W_co21   = - D_co2*(X_co21(1) - X_co21(0))/dx ;
  
  return(0) ;
  
#undef X_co21
#undef Av1
#undef Avn
#undef N_co21
#undef N_chimie1
#undef W_co21
}

int  mx48(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])
#define X_co21(i) (u_1[(i)][I_x_co2])
#define Avn(i)    (u_n[(i)][I_av])

  double x_co2, av, xsi, x_R, d_co2, D_co2, tr_co2, dx_Rsdav, dxsisdx_co2, dxsisdav ;
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
  phi      = el.mat->pr[pm("phi")] ;
  phi_be   = el.mat->pr[pm("phi_be")] ;
  S        = el.mat->pr[pm("S")] ;
  S_be     = el.mat->pr[pm("S_be")] ;
  a        = el.mat->pr[pm("a")] ;
  b        = el.mat->pr[pm("b")] ;
  a_be     = el.mat->pr[pm("a_be")] ;
  b_be     = el.mat->pr[pm("b_be")] ;
  R_0      = el.mat->pr[pm("R_0")] ;
  k_henry  = el.mat->pr[pm("k_henry")] ;
  D_0      = el.mat->pr[pm("D_0")] ;
  s_CH     = el.mat->pr[pm("s_CH")] ;

  /*
    CALCUL DE volume ET DE surf
  */
  dx = x[1][0] - x[0][0] ;
  xm = (x[1][0] + x[0][0])/deux ;
  for(i=0;i<2;i++)
    {
      volume[i] = fabs(dx)/deux ; 
      if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
    }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  
  /* termes d'accumulation */
  d_co2       = D_0*pow(phi_be,a_be)*pow(1 - S_be,b_be) ;

  for(i=0;i<2;i++)
    {
	x_co2       = X_co21(i) ;
	av          = Avn(i) ;
	x_R         = pow(1 - av,1./3.) ;
	xsi         = 3*((d_co2*k_henry*S*x_co2/R_0)*x_R/(1 - x_R))*(1 - phi)/R_0*pow(1-av,Puiss) ;
	dx_Rsdav    = -(1/3)*pow(1 - av, -2./3.) ;
	dxsisdx_co2 = 3*((d_co2/R_0)*k_henry*S*x_R/(1 - x_R))*(1 - phi)/R_0*pow(1-av,Puiss) ;
	dxsisdav    = 3*((d_co2*k_henry*S*x_co2/R_0)*(1 - phi)/R_0)*pow(1 - x_R,-2.)*dx_Rsdav*pow(1-av,Puiss) - Puiss*xsi/(1 - av) ;
        if(av >= 1.) {dxsisdx_co2 = 0. ; dxsisdav = 0. ; }

      /*
        CONSERVATION DU CO2 : (N_co21 - N_co2n) + dt * div(w_co21) = 0
      */
      K(i*NEQ+E_co2,i*NEQ+I_x_co2) += volume[i]*phi*(1 - S + k_henry*S) ;
      K(i*NEQ+E_co2,i*NEQ+I_av)    += volume[i]*(1 - phi)*s_CH ;
      /*
        CINETIQUE CHIMIQUE : (N_chimie1 - N_chimien) - dt * Xsi1 = 0
      */ 
      K(i*NEQ+E_chimie,i*NEQ+I_x_co2) += -dt*dxsisdx_co2 ;
      K(i*NEQ+E_chimie,i*NEQ+I_av) += (1 - phi)*s_CH ;

    }
  /* termes d'ecoulement */
	D_co2  = D_0*pow(phi,a)*pow(1 - S,b) ;
	tr_co2 = dt*surf*D_co2/dx ;

  /*
        CONSERVATION DU CO2 : (N_co21 - N_co2n) + dt * div(w_co21) = 0
  */
  K(E_co2,I_x_co2)         += + tr_co2 ;
  K(E_co2,NEQ+I_x_co2)     += - tr_co2 ;
  K(NEQ+E_co2,I_x_co2)     += - tr_co2 ;
  K(NEQ+E_co2,NEQ+I_x_co2) += + tr_co2 ;

  return(0) ;
  
#undef Avn
#undef X_co21
#undef K
}

void rs48(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{

#define N_co21(i)     (f_1[(i)])
#define N_co2n(i)     (f_n[(i)])
#define N_chimie1(i)  (f_1[(i+2)])
#define N_chimien(i)  (f_n[(i+2)])
#define Xsi1(i)       (f_1[(i+4)])
#define W_co21        (f_1[(6)])
#define R(n,i)        (r[(n)*NEQ+(i)])

  double dx,xm ;
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
  for(i=0;i<2;i++)
    {
      volume[i] = fabs(dx)/deux ; 
      if(geom == AXIS) volume[i] *= M_PI*(x[i][0] + xm) ; 
    }
  if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
  
  /*
    CONSERVATION DU CO2 : (N_co21 - N_co2n) + dt * div(w_co21) = 0
  */
  R(0,E_co2)    -= volume[0]*(N_co21(0) - N_co2n(0)) + dt*surf*W_co21 ;
  R(1,E_co2)    -= volume[1]*(N_co21(1) - N_co2n(1)) - dt*surf*W_co21 ;
  /*
    CINETIQUE CHIMIQUE : (N_chimie1 - N_chimien) - dt * Xsi1 = 0
  */ 
  R(0,E_chimie) -= (N_chimie1(0) - N_chimien(0)) - dt*Xsi1(0) ;
  R(1,E_chimie) -= (N_chimie1(1) - N_chimien(1)) - dt*Xsi1(1) ;
  
#undef N_co21
#undef N_co2n
#undef N_chimie1
#undef N_chimien
#undef W_co21
#undef Xsi1
#undef R
}

int so48(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define X_co2(i)     (u[(i)][I_x_co2])
#define Av(i)        (u[(i)][I_av])
#define Xsi(i)       (f[(i+4)])
  double dx,x_co2,av,d_co2,x_R,xsi,n_co2_f,D_co2,w_co2 ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi      = el.mat->pr[pm("phi")] ;
  phi_be   = el.mat->pr[pm("phi_be")] ;
  S        = el.mat->pr[pm("S")] ;
  S_be     = el.mat->pr[pm("S_be")] ;
  a        = el.mat->pr[pm("a")] ;
  b        = el.mat->pr[pm("b")] ;
  a_be     = el.mat->pr[pm("a_be")] ;
  b_be     = el.mat->pr[pm("b_be")] ;
  R_0      = el.mat->pr[pm("R_0")] ;
  k_henry  = el.mat->pr[pm("k_henry")] ;
  D_0      = el.mat->pr[pm("D_0")] ;
  s_CH     = el.mat->pr[pm("s_CH")] ;

  /* initialisation */
  nso = 5 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
 
      x_co2       = param(u,h_s,el.nn,I_x_co2) ;
      av          = param(u,h_s,el.nn,I_av) ;

      d_co2       = D_0*pow(phi_be,a_be)*pow(1 - S_be,b_be) ;	
      x_R         = pow(1 - av,1./3.) ;
      xsi         = (Xsi(0) + Xsi(1))/2 ;
      
      n_co2_f     = (1 - phi)*s_CH*av ;

      dx          = x[1][0] - x[0][0] ;
      D_co2       = D_0*pow(phi,a)*pow(1 - S,b) ;
      w_co2       = -D_co2*(X_co2(1) - X_co2(0))/dx ;

  /* quantites exploitees */
  strcpy(r[0].text,"x_co2") ; r[0].n = 1 ;
  r[0].v[0] = x_co2 ;
  strcpy(r[1].text,"av") ; r[1].n = 1 ;
  r[1].v[0] = av ;
  strcpy(r[2].text,"n_co2_f") ; r[2].n = 1 ;
  r[2].v[0] = n_co2_f ;
  strcpy(r[3].text,"w_co2") ; r[3].n = 1 ;
  r[3].v[0] = w_co2 ;
  strcpy(r[4].text,"xsi") ; r[4].n = 1 ;
  r[4].v[0] = xsi ;
  return (nso) ;
  
#undef X_co2
#undef Av
#undef Xsi

}


#undef NEQ
#undef E_co2
#undef E_chimie
#undef I_x_co2
#undef I_av
#undef NVE
#undef NVI
#undef Puiss
