#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  12
#define TITLE "Corrosion en surface"
#define AUTHORS "Dridi"

#include "OldMethods.h"

/* Macros */
#define NEQ   (2)
#define E_liq (0)
#define E_hyd (1)
#define I_p_l (0)
#define I_p_h (1)
/* Fonctions */
static int    pm(char *s) ;
/* Parametres */
static double M_h2,rho_l,M_l,RT_0,k_hen,k_ox0,k_ox1,DG_0,p_l0,C_g0,A_1,C_fer ;

int pm(char *s)
{
if(strcmp(s,"M_h2") == 0) return (0) ;
else if(strcmp(s,"rho_l") == 0) return (1) ;
else if(strcmp(s,"M_l") == 0) return (2) ;
else if(strcmp(s,"RT_0") == 0) return (3) ;
else if(strcmp(s,"k_hen") == 0) return (4) ;
else if(strcmp(s,"k_ox0") == 0) return (5) ;
else if(strcmp(s,"k_ox1") == 0) return (6) ;
else if(strcmp(s,"DG_0") == 0) return (7) ;
else if(strcmp(s,"p_l0") == 0) return (8) ;
else if(strcmp(s,"C_g0") == 0) return (9) ;
else if(strcmp(s,"A_1") == 0) return (10) ;
else if(strcmp(s,"C_fer") == 0) return (11) ;
else
  { printf("donnee \"%s\" non connue (pm12)\n",s) ; exit(0) ; }
}

int dm12(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  n_donnees = 12 ;
  
  if(dim > 1) {
    arret("Ce modele n\'est pas prevu en dimension > 1 (dm12)") ;
    exit(0) ;
  }
  
  mat->neq        = NEQ ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->eqn[E_hyd],"hyd") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  strcpy(mat->inc[I_p_h],"p_h") ;

  dmat(mat,ficd,pm,n_donnees) ;
  return(mat->n) ;
}

int qm12(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 2 CL portant sur :\n\
\t 1. le flux d\'H20 (p_l)\n\
\t 2. le flux de H2  (p_h)\n") ;
  
  printf("\n\
Exemple de donnees :\n\n") ;
  
  fprintf(ficd,"M_h2 = 0.002      # Masse molaire de H2\n") ;
  fprintf(ficd,"rho_l = 1000      # Masse volumique du liquide\n") ;
  fprintf(ficd,"M_l = 0.018       # Masse molaire du liquide\n") ;
  fprintf(ficd,"RT_0 = 2479.      # Constante des gaz parfaits fois la temperature\n") ;
  fprintf(ficd,"k_hen = 7.6e-06   # Constante de Henry\n") ;
  fprintf(ficd,"k_ox0 = 1.e17     # Constante de transport dans la couche d\'oxyde\n") ;
  fprintf(ficd,"k_ox1 = 7.e22     # Coefficient de transport dans la couche d\'oxyde\n") ;
  fprintf(ficd,"DG_0 = -21.8125e3 # Enthalpie libre standard de la reaction de corrosion\n") ;
  fprintf(ficd,"p_l0 = 100000.    # Pression de liquide de reference\n") ;
  fprintf(ficd,"C_g0 = 35270.     # Concentration de H2 dissous de reference\n") ;
  fprintf(ficd,"A_1 = 1.e-11      # Constante de vitesse de corrosion\n") ;
  fprintf(ficd,"C_fer = 141000.   # Concentration molaire du fer\n") ;
  
  return(NEQ) ;
}

void tb12(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 2 ;
  el->n_ve = 0 ;
}


void ch12(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
}


void in12(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
double cg,Vl ;
double zero = 0. ;

#define P_l   (u[0][I_p_l])
#define P_h   (u[0][I_p_h])
#define V_c   (f[0])
#define Ep    (f[1])
/*
  Donnees
*/
M_h2  = el.mat->pr[pm("M_h2")] ;
rho_l = el.mat->pr[pm("rho_l")] ;
M_l   = el.mat->pr[pm("M_l")] ;
RT_0  = el.mat->pr[pm("RT_0")] ;
k_hen = el.mat->pr[pm("k_hen")] ;
k_ox0 = el.mat->pr[pm("k_ox0")] ;
k_ox1 = el.mat->pr[pm("k_ox1")] ;
DG_0  = el.mat->pr[pm("DG_0")] ;
p_l0  = el.mat->pr[pm("p_l0")] ;
C_g0  = el.mat->pr[pm("C_g0")] ;
A_1   = el.mat->pr[pm("A_1")] ;
C_fer = el.mat->pr[pm("C_fer")] ;

/* epaisseur */
Ep   = zero ;
/* vitesse de corrosion */
cg   = k_hen*P_h ;
Vl   = M_l/rho_l ;
V_c  = -A_1*(DG_0 + 4./3.*RT_0*(cg*Vl + log(cg/C_g0))-4./3.*Vl*(P_l-p_l0)) ;

#undef P_l
#undef P_h
#undef V_c
#undef Ep
}

int ex12(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom
,double t) 
/* Termes explicites (va)  */
{
  return(0) ;
}


int ct12(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
double cg,Vl,dr ;
double deux = 2. ;

#define P_l1   (u_1[0][I_p_l])
#define P_h1   (u_1[0][I_p_h])
#define V_c1   (f_1[0])
#define Ep1    (f_1[1])
#define V_cn   (f_n[0])
#define Epn    (f_n[1])
/*
Donnees
*/
M_h2  = el.mat->pr[pm("M_h2")] ;
rho_l = el.mat->pr[pm("rho_l")] ;
M_l   = el.mat->pr[pm("M_l")] ;
RT_0  = el.mat->pr[pm("RT_0")] ;
k_hen = el.mat->pr[pm("k_hen")] ;
k_ox0 = el.mat->pr[pm("k_ox0")] ;
k_ox1 = el.mat->pr[pm("k_ox1")] ;
DG_0  = el.mat->pr[pm("DG_0")] ;
p_l0  = el.mat->pr[pm("p_l0")] ;
C_g0  = el.mat->pr[pm("C_g0")] ;
A_1   = el.mat->pr[pm("A_1")] ;
C_fer = el.mat->pr[pm("C_fer")] ;

/* vitesse de corrosion */
cg    = k_hen*P_h1 ;
Vl    = M_l/rho_l ;
dr    = 1/A_1 + (k_ox0+k_ox1/P_h1)*Epn;
V_c1  = -(DG_0+ 4./3.*RT_0*(cg*Vl + log(cg/C_g0))-4./3.*Vl*(P_l1-p_l0))/dr ;
/* epaisseur */
Ep1   = Epn + 2.3*dt*(V_c1+V_cn)/deux/C_fer ;

 return(0) ;

#undef P_l1
#undef P_h1
#undef V_c1
#undef Ep1
#undef V_cn
#undef Epn
}

int mx12(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
double cg,Vl,dvcsdpl,dvcsdpg,dr,ddrsdpg;
double xm,surf ;
double zero = 0.,un = 1.,deux = 2. ;
int    i ;
#define K(i,j)   (k[(i)*NEQ+(j)])
#define P_l1     (u_1[0][I_p_l])
#define P_h1     (u_1[0][I_p_h])
#define Epn      (f_n[1])

/*
Donnees
*/
M_h2  = el.mat->pr[pm("M_h2")] ;
rho_l = el.mat->pr[pm("rho_l")] ;
M_l   = el.mat->pr[pm("M_l")] ;
RT_0  = el.mat->pr[pm("RT_0")] ;
k_hen = el.mat->pr[pm("k_hen")] ;
k_ox0 = el.mat->pr[pm("k_ox0")] ;
k_ox1 = el.mat->pr[pm("k_ox1")] ;
DG_0  = el.mat->pr[pm("DG_0")] ;
p_l0  = el.mat->pr[pm("p_l0")] ;
C_g0  = el.mat->pr[pm("C_g0")] ;
A_1   = el.mat->pr[pm("A_1")] ;
C_fer = el.mat->pr[pm("C_fer")] ;

for(i=0;i<NEQ*NEQ;i++) k[i] = zero ;

/*
CALCUL DE surf
*/
xm = x[0][0] ;
if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
/*
Derivees de la vitesse de corrosion
*/
Vl      = M_l/rho_l ;
cg    = k_hen*P_h1 ;
dr    = 1/A_1 + (k_ox0+k_ox1/P_h1)*Epn ;
ddrsdpg = -k_ox1*Epn/P_h1/P_h1 ;
dvcsdpl = 4./3.*Vl/dr ;
dvcsdpg = -4./3.*RT_0*(1/P_h1+k_hen*Vl)/dr
         -(DG_0+ 4./3.*RT_0*(cg*Vl + log(cg/C_g0))-4./3.*Vl*(P_l1-p_l0))*ddrsdpg/(dr*dr) ;
/*
CONDITION SUR LE LIQUIDE : f(pl,pg)-dt*w_l.n = 0
*/
K(E_liq,I_p_l) = 4./3.*dt*surf*(M_l*dvcsdpl) ;
K(E_liq,I_p_h) = 4./3.*dt*surf*(M_l*dvcsdpg) ;

/*
CONDITION SUR H2 : f(pl,pg)-dt*w_g.n = 0
*/
K(E_hyd,I_p_l) = 4./3.*dt*surf*(-M_h2*dvcsdpl);
K(E_hyd,I_p_h) = 4./3.*dt*surf*(-M_h2*dvcsdpg) ;

 return(0) ;

#undef K
#undef P_l1
#undef P_h1
#undef Epn
}

void rs12(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
double xm,surf ;
int    i ;
double zero = 0.,un = 1.,deux = 2. ;

#define V_c1      (f_1[0])

/*
Donnees
*/
M_h2  = el.mat->pr[pm("M_h2")] ;
rho_l = el.mat->pr[pm("rho_l")] ;
M_l   = el.mat->pr[pm("M_l")] ;
RT_0  = el.mat->pr[pm("RT_0")] ;
k_hen = el.mat->pr[pm("k_hen")] ;
k_ox0 = el.mat->pr[pm("k_ox0")] ;
k_ox1 = el.mat->pr[pm("k_ox1")] ;
DG_0  = el.mat->pr[pm("DG_0")] ;
p_l0  = el.mat->pr[pm("p_l0")] ;
C_g0  = el.mat->pr[pm("C_g0")] ;
A_1   = el.mat->pr[pm("A_1")] ;
C_fer = el.mat->pr[pm("C_fer")] ;
/*
CALCUL DE volume ET DE surf
*/
xm = x[0][0] ;
if(geom == AXIS) surf = deux*M_PI*xm ; else surf = un ;
/* Initialisation */
for(i=0;i<NEQ;i++) r[i] = zero ;
/*
CONDITION SUR LE LIQUIDE : f(pl,pg)-dt*w_l.n = 0
*/
r[E_liq]   -= 4./3.*dt*surf*(M_l*V_c1) ;
/*
CONDITION SUR H2 : f(pl,pg)-dt*w_g.n = 0
*/
r[E_hyd]   -= 4./3.*dt*surf*(-M_h2*V_c1) ;
#undef V_c1
}

int so12(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define P_l     (u[0][I_p_l])
#define P_h     (u[0][I_p_h])
#define V_c     (f[0])
#define Ep      (f[1])
  double Am, cg, Vl ;
  
  /*
    Donnees
  */
  M_h2  = el.mat->pr[pm("M_h2")] ;
  rho_l = el.mat->pr[pm("rho_l")] ;
  M_l   = el.mat->pr[pm("M_l")] ;
  RT_0  = el.mat->pr[pm("RT_0")] ;
  k_hen = el.mat->pr[pm("k_hen")] ;
  k_ox0 = el.mat->pr[pm("k_ox0")] ;
  k_ox1 = el.mat->pr[pm("k_ox1")] ;
  DG_0  = el.mat->pr[pm("DG_0")] ;
  p_l0  = el.mat->pr[pm("p_l0")] ;
  C_g0  = el.mat->pr[pm("C_g0")] ;
  A_1   = el.mat->pr[pm("A_1")] ;
  C_fer = el.mat->pr[pm("C_fer")] ;
  
  /* dr */
  cg    = k_hen*P_h ;
  Vl    = M_l/rho_l ;
  Am  = -(DG_0 + 4./3.*RT_0*(cg*Vl + log(cg/C_g0))-4./3.*Vl*(P_l-p_l0));
  /* quantites exploitees par element */
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] = P_l ;
  strcpy(r[1].text,"pression-H2") ; r[1].n = 1 ;
  r[1].v[0] = P_h ;
  strcpy(r[2].text,"epaisseur-corrodee") ; r[2].n = 1 ;
  r[2].v[0] = Ep/2.3 ;
  strcpy(r[3].text,"vitesse-corrosion") ; r[3].n = 1 ;
  r[3].v[0] = V_c/C_fer ;
  strcpy(r[4].text,"Affinite") ; r[4].n = 1 ;
  r[4].v[0] = Am ;
  return(5) ;
#undef P_l
#undef P_h
#undef V_c
#undef Ep
}
