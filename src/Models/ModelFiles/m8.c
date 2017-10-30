#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  8
#define TITLE "Stockage d\'eau de pluie en surface"
#define AUTHORS "Berthier-Dangla"

#include "OldMethods.h"

/* Fichier */
static FILE  *fpluie ;
/* Fonctions */
static int    pm(const char *s) ;
static double h_pluie(double t,FILE *fpluie) ;
/* #define  H_PLUIE(t)      h_pluie(t,fpluie) */
#define  H_PLUIE(t)      courbe(t,el.mat->cb[0])
/* Parametres */
static double gravite,rho_l,mu_l,p_g,L ; /* L=longueur caracteristique */

int pm(const char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"rho_l") == 0) return (1) ;
  else if(strcmp(s,"mu_l") == 0) return (2) ;
  else if(strcmp(s,"p_g") == 0) return (3) ;
  else if(strcmp(s,"L") == 0) return (4) ;
  else if(strcmp(s,"courbes") == 0) return (5) ;
  else return(-1) ;
}

int dm8(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  nd = 6 ;

  if(dim > 1) arret("dm8 : dimension > 1 non prevue") ;

  mat->neq    = 1 ;
  strcpy(mat->eqn[0],"liq") ;
  strcpy(mat->inc[0],"p_l") ;

  dmat(mat,ficd,pm,nd) ;
  /*
      fpluie = fopen("tombee_h","r") ;
      if(!fpluie) {
        printf("erreur a l'ouverture du fichier %s \n",nom) ;
        exit(0) ;
      }
  */
  return(mat->n) ;
}

int qm8(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(1) ;
  
  printf("\n\n\
Modelisation du stockage de l\'eau de pluie en surface.\n\
L\'inconnue est la pression de liquide.\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 9.81   # Gravite\n") ;
  fprintf(ficd,"rho_l = 1.e3     # Masse volumique de liquide\n") ;
  fprintf(ficd,"mu_l = 1.e-3     # Viscosite du liquide\n") ;
  fprintf(ficd,"p_g = 1.e5       # Pression de gaz\n") ; 
  fprintf(ficd,"L = 100          # Longueur caracteristique de ruisselement\n") ;
  fprintf(ficd,"courbe = my_file # Nom du fichier : t h\n") ;
  
  return(1) ;
}

void tb8(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 0 ;
  el->n_ve = 0 ;
}


void ch8(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
}

void in8(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
}

int ex8(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  return(0) ;
}

int ct8(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  return(0) ;
}


int mx8(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
  double pc ;
  double h_1,h_n,dhsdpl ;
  double u_0,alpha = 0.05,zeta = 0.5 ;
  double zero = 0.,un = 1. ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  L       = el.mat->pr[pm("L")] ;
  u_0     = zeta*rho_l*fabs(gravite)*sin(alpha)/(12.*mu_l)/L ;
  
  /*
    CONSERVATION DE L'EAU DE PLUIE : rho_l*(h_1-h_n)-dt*(w_pluie-w_l-w_rui)=0
  */
  pc     = p_g - u_1[0][0] ;
  h_1    = (pc>zero) ? zero : -pc/(rho_l*fabs(gravite)) ;
  dhsdpl = (pc>zero) ? zero : un/(rho_l*fabs(gravite)) ;
  pc     = p_g - u_n[0][0] ;
  h_n    = (pc>zero) ? zero : -pc/(rho_l*fabs(gravite)) ;
  
  k[0]   = rho_l*(un+u_0*dt*(3*h_1*h_1+2*h_1*h_n+h_n*h_n))*dhsdpl ;

  return(0) ;
}

void rs8(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  double w_pluie,w_rui ;
  double pc ;
  double h_1,h_n,dhsdpl ;
  double u_0,alpha = 0.05,zeta = 0.5 ;
  double zero = 0.,un = 1. ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  L       = el.mat->pr[pm("L")] ;
  u_0     = zeta*rho_l*fabs(gravite)*sin(alpha)/(12.*mu_l)/L ;
  
  /*
    CONSERVATION DE L'EAU DE PLUIE : rho_l*(h_1-h_n)-dt*(w_pluie-w_l-w_rui)=0
  */
  pc     = p_g - u_1[0][0] ;
  h_1    = (pc>zero) ? zero : -pc/(rho_l*fabs(gravite)) ;
  dhsdpl = (pc>zero) ? zero : un/(rho_l*fabs(gravite)) ;
  pc     = p_g - u_n[0][0] ;
  h_n    = (pc>zero) ? zero : -pc/(rho_l*fabs(gravite)) ;
  
  w_pluie = rho_l*(H_PLUIE(t) - H_PLUIE(t-dt)) ;
  w_rui   = rho_l*u_0*dt*(h_1*h_1+h_n*h_n)*(h_1+h_n) ;
  r[0]    = zero ;
  r[0]   -= rho_l*(h_1-h_n)-w_pluie+w_rui ;
}

int so8(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (f) */
{
  double w_rui ;
  double pc,h ;
  double u_0,alpha = 0.05,zeta = 0.5 ;
  double zero = 0. ;
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  L       = el.mat->pr[pm("L")] ;
  
  pc   = p_g - u[0][0] ;
  h    = (pc>zero) ? zero : -pc/(rho_l*fabs(gravite)) ;
  u_0  = zeta*rho_l*fabs(gravite)*sin(alpha)/(12.*mu_l)/L ;
  w_rui = rho_l*u_0*4.*h*h*h ;
  /* quantites exploitees par element */
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] = u[0][0] ;
  strcpy(r[1].text,"w_rui") ; r[1].n = 1 ;
  r[1].v[0] = w_rui ;
  strcpy(r[2].text,"hauteur") ; r[2].n = 1 ;
  r[2].v[0] = h ;
  return(3) ;
}

double h_pluie(double t,FILE *fstrm)
/* hauteur d'eau tombee a l'instant t */
{
  long int ligne,decal ;
  double   h,h1,h2,t1,t2,dt=3600. ;
  ligne = (long int) t/dt ;
  decal = 26L*ligne ;
  fseek(fstrm,decal,SEEK_SET) ;
  fscanf(fstrm,"%le %le",&t1,&h1) ;
  fscanf(fstrm,"%le %le",&t2,&h2) ;
  h  = h1 + (h2-h1)*(t-t1)/dt ;
  return (h) ;
}
