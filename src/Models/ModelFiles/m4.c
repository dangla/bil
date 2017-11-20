#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  4
#define TITLE "Poromecanique des argiles gonflantes (Total)"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Constantes */
#define NEQ   (3)
#define E_eau (0)
#define E_sel (1)
#define E_mec (2)
#define I_p_e (0)
#define I_p_s (1)
#define I_dep (2)
/* Fonctions */
static int    pm(const char *s) ;
static double pi_g(double ps,double phi,crbe_t *) ;
static void   courbes(double h,double v,char *nom) ;
/* Parametres */
static double phi0,p_e0,p_s0,k_ee,k_ss,Ms,T_0,omega,young,poisson,Vs,R,Me,Ve ;

int pm(const char *s)
{
if(strcmp(s,"phi") == 0) return (0) ;
else if(strcmp(s,"M_e") == 0) return (1) ;
else if(strcmp(s,"V_e") == 0) return (2) ;
else if(strcmp(s,"M_s") == 0) return (3) ;
else if(strcmp(s,"k_ee") == 0) return (4) ;
else if(strcmp(s,"k_ss") == 0) return (5) ;
else if(strcmp(s,"omega") == 0) return (6) ;
else if(strcmp(s,"T_0") == 0) return (7) ;
else if(strcmp(s,"young") == 0) return (8) ;
else if(strcmp(s,"poisson") == 0) return (9) ;
else if(strcmp(s,"p_e0") == 0) return (10) ;
else if(strcmp(s,"p_s0") == 0) return (11) ;
else if(strcmp(s,"h") == 0) return (12) ;
else if(strcmp(s,"v") == 0) return (13) ;
else if(strcmp(s,"V_s") == 0) return (14) ;
else if(strcmp(s,"R") == 0) return (15) ;
else if(strcmp(s,"n_p") == 0) return (16) ;
else if(strcmp(s,"p_s1") == 0) return (17) ;
else if(strcmp(s,"p_s2") == 0) return (18) ;
else if(strcmp(s,"courbes") == 0) return (19) ;
else
  { printf("donnee \"%s\" non connue (pm4)\n",s) ; exit(0) ; }
}


int dm4(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int i,nd,n_donnees = 20 ;
  double ps,dp,h,dh ;
  char mot[80],nom[80] ;
  int long pos ;
  FILE *fict ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_eau],"eau") ;
  strcpy(mat->eqn[E_sel],"sel") ;
  strcpy(mat->eqn[E_mec],"meca") ;
  strcpy(mat->inc[I_p_e],"p_e") ;
  strcpy(mat->inc[I_p_s],"p_s") ;
  strcpy(mat->inc[I_dep],"u_x") ;

  /* position dans le fichier */
  pos = ftell(ficd) ;
  
  nd = 19 ;
  for(i=0;i<nd;i++) {
    fscanf(ficd," %[^= ] =",mot) ; fscanf(ficd,"%lf",mat->pr+pm(mot)) ;
  }
  fscanf(ficd," %[^= ] = %s",mot,nom) ;
  fict = fopen(nom,"r") ;
  /* fabrication des courbes */
  if(fict == NULL) {
    h  = mat->pr[pm("h")] ;
    dh = h/100. ;
    courbes(h,mat->pr[pm("v")],"dce1") ;
    courbes(h+dh,mat->pr[pm("v")],"dce2") ;
    ps = 2.*8.315*mat->pr[pm("T_0")]*1.e3 ;
    dp = h/mat->pr[pm("phi")]/dh ;
    printf("Suivre les etapes suivantes pour creer le fichier %s (sous unix):\n\
\t1. Executez les fichiers dce1 et dce2\n\
\t2. Tapez la commande \"awk \'FNR<=3 {print $0};FNR>3 {print $5}\' dce1.p1 > %s\"\n\
\t3. Supprimez le fichier dce1.p1\n\
\t4. Tapez la commande \"paste dce2.p1 %s > dce1.p1\"\n\
\t5. Supprimez le fichier %s\n\
\t6. Tapez la commande \"awk \'FNR>3 {print %g*$1,$8,%g*($5-$8)}\' dce1.p1 > %s\"\n\
%s est maintenant cree. Relancez !\n",nom,nom,nom,nom,ps,dp,nom,nom) ;
    exit(0) ;
  }

  /* on retourne en debut de lecture */
  fseek(ficd,pos,SEEK_SET) ;

  dmat(mat,ficd,pm,n_donnees) ;
  
  return(mat->n) ;
}


int qm4(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est formee de : deux equations de conservation.\n\
\t 1. l\'equation de conservation de la masse d\'eau (p_e)\n\
\t 2. l\'equation de conservation de la masse de sel (p_s)\n\
\t 3. l\'equation d\'equilibre mecanique (u_x)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;


  fprintf(ficd,"h = 1.e-9         # Demi-distance interparticule\n") ;
  fprintf(ficd,"v = 2             # Valence des ions\n") ;
  fprintf(ficd,"young = 1.e9      # Module d'Young\n") ;
  fprintf(ficd,"poisson = 0.      # Coefficient de Poisson\n") ;
  fprintf(ficd,"phi = 0.3         # Porosite\n") ;
  fprintf(ficd,"M_e = 18.e-3      # Masse molaire de l'eau\n") ;
  fprintf(ficd,"V_se = 18.e-6     # Volume molaire partiel de l'eau\n") ;
  fprintf(ficd,"M_s = 111.e-3     # Masse molaire du sel\n") ;
  fprintf(ficd,"V_s = 25.e-6      # Volume molaire partiel du sel\n") ;
  fprintf(ficd,"R = 8.315         # Constante des gaz parfaits\n") ;
  fprintf(ficd,"p_e0 = 1.31e7     # Pression initiale de reference de l'eau\n") ;
  fprintf(ficd,"p_s0 = 0.19e7     # Pression initiale de reference du sel\n") ;
  fprintf(ficd,"k_ee = 4.e-18     # Coefficient k_ee\n") ;
  fprintf(ficd,"k_ss = 6.e-18     # Coefficient k_ss\n") ;
  fprintf(ficd,"omega = 0.4       # Omega\n") ;
  fprintf(ficd,"T_0 = 300         # Temperature\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_s l dl/dphi\n") ;
  
  return(NEQ) ;
}

void tb4(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = 8 ;
  el->n_ve = 4 ;
}

void ch4(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int   i ;
  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
 
  if(isdigit(cg.eqn[0]) && (atoi(cg.eqn) - 1) == E_eau) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  } else if(!strcmp(cg.eqn,el.mat->eqn[E_eau])) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  }

  if(isdigit(cg.eqn[0]) && (atoi(cg.eqn) - 1) == E_sel) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  } else if(!strcmp(cg.eqn,el.mat->eqn[E_sel])) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  }
}


void in4(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define P_e(i)   (u[(i)][I_p_e])
#define P_s(i)   (u[(i)][I_p_s])
#define U_x(i)   (u[(i)][I_dep])
#define M_e(i)   (f[(i)])
#define M_s(i)   (f[(i+2)])
#define W_e      (f[(4)])
#define W_s      (f[(5)])
#define S_xx     (f[(6)])
#define S_tt     (f[(7)])
#define K_ee     (va[(0)])
#define K_es     (va[(1)])
#define K_se     (va[(2)])
#define K_ss     (va[(3)])
  double rho_e,rho_s,Mss2RT,MssVs,MesVe ;
  double dx,p_sij,p_eij,lambda,b_s,e_xx,e_tt,dmu,lame,pie,pie0,phi,dphi ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  crbe_t *cb = el.mat->cb ;
  
  if(el.dim < dim) return ;
  /*
    Donnees
  */
  young = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0  = el.mat->pr[pm("phi")] ;
  Me    = el.mat->pr[pm("M_e")] ;
  Ve    = el.mat->pr[pm("V_e")] ;
  k_ee  = el.mat->pr[pm("k_ee")] ;
  k_ss  = el.mat->pr[pm("k_ss")] ;
  R     = el.mat->pr[pm("R")] ;
  Ms    = el.mat->pr[pm("M_s")] ;
  Vs    = el.mat->pr[pm("V_s")] ;
  omega = el.mat->pr[pm("omega")] ;
  T_0   = el.mat->pr[pm("T_0")] ;
  p_e0  = el.mat->pr[pm("p_e0")] ;
  p_s0  = el.mat->pr[pm("p_s0")] ;
  dmu   = young/(un+poisson) ;
  lame  = dmu*poisson/(un-deux*poisson) ;
  Mss2RT = Ms/deux/R/T_0 ;
  MesVe = Me/Ve ;
  MssVs = Ms/Vs ;
  
  dx = x[1][0] - x[0][0] ;
  /* deformations */
  e_xx = (U_x(1) - U_x(0))/dx ;
  if(geom == AXIS) e_tt = (U_x(1)+U_x(0))/(x[1][0]+x[0][0]) ; else e_tt = zero ;
  dphi  = e_xx + e_tt ;
  phi   = phi0 + dphi ;
  /* pressions */
  p_eij = (P_e(0)+P_e(1))/deux ;
  p_sij = (P_s(0)+P_s(1))/deux ;
  /* contraintes */
  pie0 = p_e0+p_s0+pi_g(p_s0,phi0,cb) ;
  pie  = p_eij+p_sij+pi_g(p_sij,phi0,cb) ;
  S_xx = lame*(e_xx+e_tt) + dmu*e_xx - pie + pie0 ;
  S_tt = lame*(e_xx+e_tt) + dmu*e_tt - pie + pie0 ;
  /* masses fluides */
  for(i=0;i<2;i++) {
    lambda = courbe(P_s(i),cb[0]) + dphi*courbe(P_s(i),cb[1]) ;
    rho_s  = Mss2RT*P_s(i)*lambda ;
    rho_e  = MesVe*(un-rho_s/MssVs) ;
    M_e(i) = rho_e*phi ;
    M_s(i) = rho_s*phi ;
  }
  /* coefficients de transfert */
  lambda = courbe(p_sij,cb[0]) + dphi*courbe(p_sij,cb[1]) ;
  b_s    = courbe(p_sij,cb[0]) + phi0*courbe(p_sij,cb[1]) ;
  rho_s  = Mss2RT*p_sij*lambda ;
  rho_e  = MesVe*(un-rho_s/MssVs) ;
  
  K_ee   = rho_e*k_ee ;
  K_es   = rho_e*k_ee*(un - omega) ;
  if(b_s > zero) K_se  = rho_s*k_ee*(un-omega)/b_s ; else K_se = zero ;
  K_ss   = rho_s*k_ss*b_s ;
  
  /* autres tentatives */
  /*
    K_ee   = rho_e*k_ee ;
    K_es   = rho_e*k_ee*lambda ;
    K_se   = rho_s*k_ee ;
    K_ss   = rho_s*k_ss*lambda ;
  */
  /* flux */
  W_e   = - K_ee*(P_e(1) - P_e(0))/dx - K_es*(P_s(1) - P_s(0))/dx ;
  W_s   = - K_se*(P_e(1) - P_e(0))/dx - K_ss*(P_s(1) - P_s(0))/dx ;
  
#undef P_e
#undef P_s
#undef U_x
#undef M_e
#undef M_s
#undef W_e
#undef W_s
#undef S_xx
#undef S_tt
#undef K_ee
#undef K_es
#undef K_se
#undef K_ss
}

int ex4(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
#define P_e(i)   (u[(i)][I_p_e])
#define P_s(i)   (u[(i)][I_p_s])
#define U_x(i)   (u[(i)][I_dep])
#define K_ee     (va[(0)])
#define K_es     (va[(1)])
#define K_se     (va[(2)])
#define K_ss     (va[(3)])
  double rho_e,rho_s,Mss2RT,MssVs,MesVe ;
  double dx,p_sij,p_eij,lambda,b_s,e_xx,e_tt,phi,dphi ;
  double zero = 0.,un = 1.,deux = 2. ;
  crbe_t *cb = el.mat->cb ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  phi0  = el.mat->pr[pm("phi")] ;
  Me    = el.mat->pr[pm("M_e")] ;
  Ve    = el.mat->pr[pm("V_e")] ;
  k_ee  = el.mat->pr[pm("k_ee")] ;
  k_ss  = el.mat->pr[pm("k_ss")] ;
  Ms    = el.mat->pr[pm("M_s")] ;
  Vs    = el.mat->pr[pm("V_s")] ;
  R     = el.mat->pr[pm("R")] ;
  omega = el.mat->pr[pm("omega")] ;
  T_0   = el.mat->pr[pm("T_0")] ;
  Mss2RT = Ms/R/deux/T_0 ;
  MesVe = Me/Ve ;
  MssVs = Ms/Vs ;
  
  dx = x[1][0] - x[0][0] ;
  /* deformations */
  e_xx = (U_x(1) - U_x(0))/dx ;
  if(geom == AXIS) e_tt = (U_x(1)+U_x(0))/(x[1][0]+x[0][0]) ; else e_tt = zero ;
  dphi  = e_xx + e_tt ;
  phi   = phi0 + dphi ;
  /* pressions */
  p_eij = (P_e(0)+P_e(1))/deux ;
  p_sij = (P_s(0)+P_s(1))/deux ;
  /* COEFFICIENTS DE TRANSFERT */
  lambda = courbe(p_sij,cb[0]) + dphi*courbe(p_sij,cb[1]) ;
  b_s    = courbe(p_sij,cb[0]) + phi0*courbe(p_sij,cb[1]) ;
  rho_s  = Mss2RT*p_sij*lambda ;
  rho_e  = MesVe*(un-rho_s/MssVs) ;
  
  /* rho_s  = Mss2RT*el.mat->pr[pm("p_s0")]*lambda ; */ /* cas lineaire */
  
  K_ee   = rho_e*k_ee ;
  K_es   = rho_e*k_ee*(un - omega) ;
  if(b_s > zero) K_se  = rho_s*k_ee*(un-omega)/b_s ; else K_se = zero ;
  K_ss   = rho_s*k_ss*b_s ;
  
  /* autres tentatives */
  /*
    K_ee   = rho_e*k_ee ;
    K_es   = rho_e*k_ee*lambda ;
    K_se   = rho_s*k_ee ;
    K_ss   = rho_s*k_ss*lambda ;
  */
  return(0) ;

#undef P_e
#undef P_s
#undef U_x
#undef K_ee
#undef K_es
#undef K_se
#undef K_ss
}

int ct4(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
#define P_e1(i)  (u_1[(i)][I_p_e])
#define P_s1(i)  (u_1[(i)][I_p_s])
#define U_x1(i)  (u_1[(i)][I_dep])
#define M_e1(i)  (f_1[(i)])
#define M_s1(i)  (f_1[(i+2)])
#define W_e1     (f_1[(4)])
#define W_s1     (f_1[(5)])
#define S_xx1    (f_1[(6)])
#define S_tt1    (f_1[(7)])
#define K_ee     (va[(0)])
#define K_es     (va[(1)])
#define K_se     (va[(2)])
#define K_ss     (va[(3)])
  double rho_e,rho_s,Mss2RT,MssVs,MesVe ;
  double dx,p_sij,p_eij,lambda,e_xx,e_tt,dmu,lame,pie,pie0,phi,dphi ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  crbe_t *cb = el.mat->cb ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  young = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0  = el.mat->pr[pm("phi")] ;
  Me    = el.mat->pr[pm("M_e")] ;
  Ve    = el.mat->pr[pm("V_e")] ;
  k_ee  = el.mat->pr[pm("k_ee")] ;
  k_ss  = el.mat->pr[pm("k_ss")] ;
  Ms    = el.mat->pr[pm("M_s")] ;
  Vs    = el.mat->pr[pm("V_s")] ;
  R     = el.mat->pr[pm("R")] ;
  omega = el.mat->pr[pm("omega")] ;
  T_0   = el.mat->pr[pm("T_0")] ;
  p_e0  = el.mat->pr[pm("p_e0")] ;
  p_s0  = el.mat->pr[pm("p_s0")] ;
  dmu   = young/(un+poisson) ;
  lame  = dmu*poisson/(un-deux*poisson) ;
  Mss2RT = Ms/R/deux/T_0 ;
  MesVe = Me/Ve ;
  MssVs = Ms/Vs ;
  
  dx = x[1][0] - x[0][0] ;
  /* deformations */
  e_xx = (U_x1(1) - U_x1(0))/dx ;
  if(geom == AXIS) e_tt = (U_x1(1)+U_x1(0))/(x[1][0]+x[0][0]) ; else e_tt = zero ;
  dphi  = e_xx + e_tt ;
  phi   = phi0 + dphi ;
  /* pressions */
  p_eij = (P_e1(0)+P_e1(1))/deux ;
  p_sij = (P_s1(0)+P_s1(1))/deux ;
  /* contraintes */
  pie0  = p_e0+p_s0+pi_g(p_s0,phi0,cb) ;
  pie   = p_eij+p_sij+pi_g(p_sij,phi0,cb) ;
  S_xx1 = lame*(e_xx+e_tt) + dmu*e_xx - pie + pie0 ;
  S_tt1 = lame*(e_xx+e_tt) + dmu*e_tt - pie + pie0 ;
  /* masses */
  for(i=0;i<2;i++) {
    lambda  = courbe(P_s1(i),cb[0]) + dphi*courbe(P_s1(i),cb[1]) ;
    rho_s   = Mss2RT*P_s1(i)*lambda ;
    rho_e   = MesVe*(un-rho_s/MssVs) ;
    M_e1(i) = rho_e*phi ;
    M_s1(i) = rho_s*phi ;
  }
  /* flux */
  W_e1 = - K_ee*(P_e1(1) - P_e1(0))/dx - K_es*(P_s1(1) - P_s1(0))/dx ;
  W_s1 = - K_se*(P_e1(1) - P_e1(0))/dx - K_ss*(P_s1(1) - P_s1(0))/dx ;

  return(0) ;
  
#undef P_e1
#undef P_s1
#undef U_x1
#undef M_e1
#undef M_s1
#undef W_e1
#undef W_s1
#undef S_xx1
#undef S_tt1
#undef K_ee
#undef K_es
#undef K_se
#undef K_ss
}

int  mx4(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*2*NEQ+(j)])

#define P_e1(i)  (u_1[(i)][I_p_e])
#define P_s1(i)  (u_1[(i)][I_p_s])
#define U_x1(i)  (u_1[(i)][I_dep])
#define K_ee     (va[(0)])
#define K_es     (va[(1)])
#define K_se     (va[(2)])
#define K_ss     (va[(3)])
  double rho_e,rho_s,Mss2RT,MssVs,MesVe,MeVss2VeRT ;
  double p_sij,p_eij,e_xx,e_tt,lambda,dlambda,b_s,dmu,lame,phi,dphi ;
  double dx,xm,drayon,tr_ee,tr_es,tr_se,tr_ss ;
  double volume[2],surf,volumt,vsdx ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  crbe_t *cb = el.mat->cb ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  young = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0  = el.mat->pr[pm("phi")] ;
  Me    = el.mat->pr[pm("M_e")] ;
  Ve    = el.mat->pr[pm("V_e")] ;
  k_ee  = el.mat->pr[pm("k_ee")] ;
  k_ss  = el.mat->pr[pm("k_ss")] ;
  Ms    = el.mat->pr[pm("M_s")] ;
  Vs    = el.mat->pr[pm("V_s")] ;
  R     = el.mat->pr[pm("R")] ;
  omega = el.mat->pr[pm("omega")] ;
  T_0   = el.mat->pr[pm("T_0")] ;
  dmu   = young/(un+poisson) ;
  lame  = dmu*poisson/(un-deux*poisson) ;
  Mss2RT = Ms/R/deux/T_0 ;
  MesVe = Me/Ve ;
  MssVs = Ms/Vs ;
  MeVss2VeRT = Me*Vs/Ve/R/deux/T_0 ;

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
  volumt = volume[0] + volume[1] ;
  drayon = x[1][0] + x[0][0] ;
  
  /* deformations */
  e_xx = (U_x1(1) - U_x1(0))/dx ;
  if(geom == AXIS) e_tt = (U_x1(1)+U_x1(0))/(x[1][0]+x[0][0]) ; else e_tt = zero ;
  dphi  = e_xx + e_tt ;
  phi   = phi0 + dphi ;
  /* pressions */
  p_eij = (P_e1(0)+P_e1(1))/deux ;
  p_sij = (P_s1(0)+P_s1(1))/deux ;
  /* 
     MECANIQUE 
  */
  b_s   = courbe(p_sij,cb[0]) + phi0*courbe(p_sij,cb[1]) ;
  vsdx  = volumt/dx ;
  K(E_mec,I_dep)         += -vsdx*(-(lame+dmu)/dx) ;
  K(E_mec,I_p_e)         += -vsdx*(-un/deux) ;
  K(E_mec,I_p_s)         += -vsdx*(-b_s/deux) ;
  K(E_mec,NEQ+I_dep)     += -vsdx*((lame+dmu)/dx) ;
  K(E_mec,NEQ+I_p_e)     += -vsdx*(-un/deux) ;
  K(E_mec,NEQ+I_p_s)     += -vsdx*(-b_s/deux) ;
  K(NEQ+E_mec,I_dep)     += +vsdx*(-(lame+dmu)/dx) ;
  K(NEQ+E_mec,I_p_e)     += +vsdx*(-un/deux) ;
  K(NEQ+E_mec,I_p_s)     += +vsdx*(-b_s/deux) ;
  K(NEQ+E_mec,NEQ+I_dep) += +vsdx*((lame+dmu)/dx) ;
  K(NEQ+E_mec,NEQ+I_p_e) += +vsdx*(-un/deux) ;
  K(NEQ+E_mec,NEQ+I_p_s) += +vsdx*(-b_s/deux) ;
  if(geom == AXIS) {
    /* effet de la deformation e_tt dans S_rr */
    K(E_mec,I_dep)         += -vsdx*(lame/drayon) ;
    K(E_mec,NEQ+I_dep)     += -vsdx*(lame/drayon) ;
    K(NEQ+E_mec,I_dep)     += +vsdx*(lame/drayon) ;
    K(NEQ+E_mec,NEQ+I_dep) += +vsdx*(lame/drayon) ;
    /* travail de la contrainte S_tt dans e_tt */
    vsdx = volumt/drayon ;
    K(E_mec,I_dep)         += +vsdx*(+(lame+dmu)/drayon-lame/dx) ;
    K(E_mec,I_p_e)         += +vsdx*(-un/deux) ;
    K(E_mec,I_p_s)         += +vsdx*(-b_s/deux) ;
    K(E_mec,NEQ+I_dep)     += +vsdx*(+(lame+dmu)/drayon+lame/dx) ;
    K(E_mec,NEQ+I_p_e)     += +vsdx*(-un/deux) ;
    K(E_mec,NEQ+I_p_s)     += +vsdx*(-b_s/deux) ;
    K(NEQ+E_mec,I_dep)     += +vsdx*(+(lame+dmu)/drayon-lame/dx) ;
    K(NEQ+E_mec,I_p_e)     += +vsdx*(-un/deux) ;
    K(NEQ+E_mec,I_p_s)     += +vsdx*(-b_s/deux) ;
    K(NEQ+E_mec,NEQ+I_dep) += +vsdx*(+(lame+dmu)/drayon+lame/dx) ;
    K(NEQ+E_mec,NEQ+I_p_e) += +vsdx*(-un/deux) ;
    K(NEQ+E_mec,NEQ+I_p_s) += +vsdx*(-b_s/deux) ;
  }
  /*
    CONSERVATION DE L'EAU : (m_e1 - m_en) + dt * div(w_e1) = 0
  */
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    lambda   = courbe(P_s1(i),cb[0]) + dphi*courbe(P_s1(i),cb[1]) ;
    dlambda  = dcourbe(P_s1(i),cb[0]) + dphi*dcourbe(P_s1(i),cb[1]) ;
    K(i*NEQ+E_eau,i*NEQ+I_p_s) += volume[i]*MeVss2VeRT*phi*(-lambda - P_s1(i)*dlambda) ;
  }
  lambda  = courbe(P_s1(0),cb[0]) + dphi*courbe(P_s1(0),cb[1]) ;
  dlambda = dcourbe(P_s1(0),cb[0]) + dphi*dcourbe(P_s1(0),cb[1]) ;
  rho_s  = Mss2RT*P_s1(0)*lambda ;
  rho_e  = MesVe*(un-rho_s/MssVs) ;
  K(E_eau,I_dep) += volume[0]*(rho_e 
			       - MeVss2VeRT*P_s1(0)*phi*dlambda)*(-un/dx) ;
  K(E_eau,NEQ+I_dep) += volume[0]*(rho_e 
				   - MeVss2VeRT*P_s1(0)*phi*dlambda)*(un/dx) ;
  if(geom == AXIS) {
    K(E_eau,I_dep) += volume[0]*(rho_e
				 - MeVss2VeRT*P_s1(0)*phi*dlambda)*(un/drayon) ;
    K(E_eau,NEQ+I_dep) += volume[0]*(rho_e
				     - MeVss2VeRT*P_s1(0)*phi*dlambda)*(un/drayon) ;
  }
  lambda  = courbe(P_s1(1),cb[0]) + dphi*courbe(P_s1(1),cb[1]) ;
  dlambda = dcourbe(P_s1(1),cb[0]) + dphi*dcourbe(P_s1(1),cb[1]) ;
  rho_s  = Mss2RT*P_s1(1)*lambda ;
  rho_e  = MesVe*(un-rho_s/MssVs) ;
  K(NEQ+E_eau,NEQ+I_dep) += volume[1]*(rho_e 
				       - MeVss2VeRT*P_s1(1)*phi*dlambda)*(un/dx) ;
  K(NEQ+E_eau,I_dep) += volume[1]*(rho_e 
				   - MeVss2VeRT*P_s1(1)*phi*dlambda)*(-un/dx) ;
  if(geom == AXIS) {
    K(NEQ+E_eau,I_dep) += volume[1]*(rho_e
				     - MeVss2VeRT*P_s1(1)*phi*dlambda)*(un/drayon) ;
    K(NEQ+E_eau,NEQ+I_dep) += volume[1]*(rho_e
					 - MeVss2VeRT*P_s1(1)*phi*dlambda)*(un/drayon) ;
  }
  /* termes d'ecoulement */
  tr_ee = dt*surf*K_ee/dx ;
  tr_es = dt*surf*K_es/dx ;
  K(E_eau,I_p_e) += + tr_ee ;
  K(E_eau,NEQ+I_p_e) += - tr_ee ;
  K(E_eau,I_p_s) += + tr_es ;
  K(E_eau,NEQ+I_p_s) += - tr_es ;
  K(NEQ+E_eau,NEQ+I_p_e) += + tr_ee ;
  K(NEQ+E_eau,I_p_e) += - tr_ee ;
  K(NEQ+E_eau,NEQ+I_p_s) += + tr_es ;
  K(NEQ+E_eau,I_p_s) += - tr_es ;
  /*
    CONSERVATION DU SEL : (m_s1 - m_sn) + dt * div(w_s1) = 0
  */
  /* termes d'accumulation */
  for(i=0;i<2;i++) {
    lambda   = courbe(P_s1(i),cb[0]) + dphi*courbe(P_s1(i),cb[1]) ;
    dlambda  = dcourbe(P_s1(i),cb[0]) + dphi*dcourbe(P_s1(i),cb[1]) ;
    K(i*NEQ+E_sel,i*NEQ+I_p_s) += volume[i]*Mss2RT*phi*(lambda + P_s1(i)*dlambda) ;
  }
  lambda  = courbe(P_s1(0),cb[0]) + dphi*courbe(P_s1(0),cb[1]) ;
  dlambda = dcourbe(P_s1(0),cb[0]) + dphi*dcourbe(P_s1(0),cb[1]) ;
  rho_s   = Mss2RT*P_s1(0)*lambda ;
  K(E_sel,I_dep) += volume[0]*(rho_s + Mss2RT*P_s1(0)*phi*dlambda)*(-un/dx) ;
  K(E_sel,NEQ+I_dep) += volume[0]*(rho_s + Mss2RT*P_s1(0)*phi*dlambda)*(un/dx) ;
  if(geom == AXIS) {
    K(E_sel,I_dep) += volume[0]*(rho_s
				 +Mss2RT*P_s1(0)*phi*dlambda)*(un/drayon) ;
    K(E_sel,NEQ+I_dep) += volume[0]*(rho_s
				     +Mss2RT*P_s1(0)*phi*dlambda)*(un/drayon) ;
  }
  lambda  = courbe(P_s1(1),cb[0]) + dphi*courbe(P_s1(1),cb[1]) ;
  dlambda = dcourbe(P_s1(1),cb[0]) + dphi*dcourbe(P_s1(1),cb[1]) ;
  rho_s   = Mss2RT*P_s1(1)*lambda ;
  K(NEQ+E_sel,NEQ+I_dep) += volume[1]*(rho_s + Mss2RT*P_s1(1)*phi*dlambda)*(un/dx) ;
  K(NEQ+E_sel,I_dep) += volume[1]*(rho_s + Mss2RT*P_s1(1)*phi*dlambda)*(-un/dx) ;
  if(geom == AXIS) {
    K(NEQ+E_sel,NEQ+I_dep) += volume[1]*(rho_s
					 +Mss2RT*P_s1(1)*phi*dlambda)*(un/drayon) ;
    K(NEQ+E_sel,I_dep) += volume[1]*(rho_s
				     +Mss2RT*P_s1(1)*phi*dlambda)*(un/drayon) ;
  }
  /* termes d'ecoulement */
  tr_se = dt*surf*K_se/dx ;
  tr_ss = dt*surf*K_ss/dx ;
  K(E_sel,I_p_e) += + tr_se ;
  K(E_sel,NEQ+I_p_e) += - tr_se ;
  K(E_sel,I_p_s) += + tr_ss ;
  K(E_sel,NEQ+I_p_s) += - tr_ss ;
  K(NEQ+E_sel,NEQ+I_p_e) += + tr_se ;
  K(NEQ+E_sel,I_p_e) += - tr_se ;
  K(NEQ+E_sel,NEQ+I_p_s) += + tr_ss ;
  K(NEQ+E_sel,I_p_s) += - tr_ss ;

  return(0) ;
  
#undef K
#undef P_e1
#undef P_s1
#undef U_x1
#undef K_ee
#undef K_es
#undef K_se
#undef K_ss
}

void rs4(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])

#define M_e1(i)  (f_1[(i)])
#define M_en(i)  (f_n[(i)])
#define M_s1(i)  (f_1[(i+2)])
#define M_sn(i)  (f_n[(i+2)])
#define W_e1     (f_1[(4)])
#define W_s1     (f_1[(5)])
#define S_xx1    (f_1[(6)])
#define S_tt1    (f_1[(7)])
  double dx,xm,drayon ;
  double volume[2],surf,volumt,vsdx ;
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
  volumt = volume[0] + volume[1] ;
  drayon = x[1][0] + x[0][0] ;
  
  /* 
     MECANIQUE 
  */
  vsdx  = volumt/dx ;
  R(0,E_mec)        -= -vsdx*S_xx1 ;
  R(1,E_mec)        -= +vsdx*S_xx1 ;
  if(geom == AXIS) {
    /* travail de la contrainte S_tt dans e_tt */
    vsdx = volumt/drayon ;
    R(0,E_mec)        -= +vsdx*S_tt1 ;
    R(1,E_mec)        -= +vsdx*S_tt1 ;
  }
  /*
    CONSERVATION DE L'EAU : (m_e1 - m_en) + dt * div(w_e1) = 0
  */
  R(0,E_eau)        -= volume[0]*(M_e1(0) - M_en(0)) + dt*surf*W_e1 ;
  R(1,E_eau)        -= volume[1]*(M_e1(1) - M_en(1)) - dt*surf*W_e1 ;
  /*
    CONSERVATION DU SEL : (m_s1 - m_sn) + dt * div(w_s1) = 0
  */
  R(0,E_sel)        -= volume[0]*(M_s1(0) - M_sn(0)) + dt*surf*W_s1 ;
  R(1,E_sel)        -= volume[1]*(M_s1(1) - M_sn(1)) - dt*surf*W_s1 ;

#undef R
#undef M_e1
#undef M_s1
#undef W_e1
#undef W_s1
#undef S_xx1
#undef S_tt1
}

int so4(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
#define P_e(i)   (u[(i)][I_p_e])
#define P_s(i)   (u[(i)][I_p_s])
#define U_x(i)   (u[(i)][I_dep])
#define K_ee     (va[(0)])
#define K_es     (va[(1)])
#define K_se     (va[(2)])
#define K_ss     (va[(3)])
  double pe,ps,ux,we,ws,sxx,stt,lam,pig ;
  double dx,exx,ett,dmu,lame,pie,pie0,phi,dphi ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0.,un = 1.,deux = 2. ;
  crbe_t *cb = el.mat->cb ;
  
  if(el.dim < dim) return(0) ;
  
  /*
    Donnees
  */
  young = el.mat->pr[pm("young")] ;
  poisson = el.mat->pr[pm("poisson")] ;
  phi0  = el.mat->pr[pm("phi")] ;
  Me    = el.mat->pr[pm("M_e")] ;
  Ve    = el.mat->pr[pm("V_e")] ;
  k_ee  = el.mat->pr[pm("k_ee")] ;
  k_ss  = el.mat->pr[pm("k_ss")] ;
  Ms    = el.mat->pr[pm("M_s")] ;
  Vs    = el.mat->pr[pm("V_s")] ;
  R     = el.mat->pr[pm("R")] ;
  omega = el.mat->pr[pm("omega")] ;
  T_0   = el.mat->pr[pm("T_0")] ;
  p_e0  = el.mat->pr[pm("p_e0")] ;
  p_s0  = el.mat->pr[pm("p_s0")] ;
  dmu   = young/(un+poisson) ;
  lame  = dmu*poisson/(un-deux*poisson) ;

  /* initialisation */
  nso = 9 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  dx   = x[1][0] - x[0][0] ;
  /* deformations */
  exx  = (U_x(1) - U_x(0))/dx ;
  if(geom == AXIS) ett = (U_x(1)+U_x(0))/(x[1][0]+x[0][0]) ; else ett = zero ;
  dphi = exx + ett ;
  phi  = phi0 + dphi ;
  /* pressions */
  pe   =  param(u,h_s,el.nn,I_p_e) ;
  ps   =  param(u,h_s,el.nn,I_p_s) ;
  /* deplacement */
  ux   =  param(u,h_s,el.nn,I_dep) ;
  /* pression interstitielle equivalente */
  pig  = pi_g(ps,phi0,cb) ;
  /* contraintes */
  pie0 = p_e0+p_s0+pi_g(p_s0,phi0,cb) ;
  pie  = pe+ps+pig ;
  sxx  = lame*(exx+ett) + dmu*exx - pie + pie0 ;
  stt  = lame*(exx+ett) + dmu*ett - pie + pie0 ;
  /* activite */
  lam  = courbe(ps,cb[0]) + dphi*courbe(ps,cb[1]) ;
  /* flux */
  we   = - K_ee*(P_e(1) - P_e(0))/dx - K_es*(P_s(1) - P_s(0))/dx ;
  ws   = - K_se*(P_e(1) - P_e(0))/dx - K_ss*(P_s(1) - P_s(0))/dx ;
  /* quantites exploitees */
  strcpy(r[0].text,"p_e") ; r[0].n = 1 ;
  r[0].v[0] = pe ;
  strcpy(r[1].text,"p_s") ; r[1].n = 1 ;
  r[1].v[0]= ps ;
  strcpy(r[2].text,"deplacement") ; r[2].n = 1 ;
  r[2].v[0] = ux ;
  strcpy(r[3].text,"w_e") ; r[3].n = 1 ;
  r[3].v[0] = we ;
  strcpy(r[4].text,"w_s") ; r[4].n = 1 ;
  r[4].v[0] = ws ;
  strcpy(r[5].text,"s_xx") ; r[5].n = 1 ;
  r[5].v[0] = sxx ;
  strcpy(r[6].text,"s_tt") ; r[6].n = 1 ;
  r[6].v[0] = stt ;
  strcpy(r[7].text,"lambda") ; r[7].n = 1 ;
  r[7].v[0] = lam ;
  strcpy(r[8].text,"pi_g") ; r[8].n = 1 ;
  r[8].v[0] = pig ;
  /*
  strcpy(r[9].text,"?") ; r[9].n = 1 ;
  r[9].v[0] = Me/Ve*(un-Vs/(deux*R*T_0)*ps*lam) ;
  strcpy(r[10].text,"?") ; r[10].n = 1 ;
  r[10].v[0] = Ms/(deux*R*T_0)*ps*lam ;
  */
  return (nso) ;
  
#undef P_e
#undef P_s
#undef U_x
#undef K_ee
#undef K_es
#undef K_se
#undef K_ss
}


double pi_g(double ps,double phi,crbe_t *cb)
{
  double p_s1,p_s2 ;
  double p1,p2,b1,b2,pig ;
  int    i ;
  double zero = 0.,un = 1.,deux = 2. ;
  double a[3] = {0.93246951420,0.66120938646,0.23861918608} ;
  double w[3] = {0.17132449237,0.36076157304,0.46791393457} ;
  
  p_s1 = cb[0].a[0] ;
  p_s2 = cb[0].a[1] ;
  pig  = zero ;
  for(i=0;i<3;i++) {
    p1   = p_s2*(un + a[i])/deux + ps*(un - a[i])/deux ;
    p2   = p_s2*(un - a[i])/deux + ps*(un + a[i])/deux ;
    b1   = courbe(p1,cb[0]) + phi*courbe(p1,cb[1]) ;
    b2   = courbe(p2,cb[0]) + phi*courbe(p2,cb[1]) ;
    pig += w[i]*(deux - b1 - b2) ;
  }
  pig *= (p_s2 - ps)/deux ;
  return(pig) ;
}


void   courbes(double h,double v,char *nom)
{
  FILE *fic ;
  fic = fopen(nom,"w") ;
  fprintf(fic,"DIME\n		                \
 1 plan\n") ;
  fprintf(fic,"MAIL\n		                \
 3 0 %g %g\n	                                \
 %g\n		                                \
 100 1\n	                                \
 1 1\n",h,h,h/100) ;
  fprintf(fic,"MATE\n		                \
 3 \n						\
 sigma = -0.2\n					\
 e = 1.6e-19\n					\
 v = %g\n					\
 epsilon = 7.0832e-10\n				\
 k_B = 1.38e-23\n				\
 T_0 = 300\n					\
 n_0 = 0\n					\
 n_1 = 6.022e26\n				\
 ne = 100\n",v) ;
  fprintf(fic,"CHMP\n				\
 1\n						\
 V = 1. G = 0. P = 0.\n") ;
  fprintf(fic,"INIT\n				\
 0\n") ;
  fprintf(fic,"FONC\n				\
 0\n") ;
  fprintf(fic,"COND\n				\
 0\n") ;
  fprintf(fic,"CHAR\n				\
 1\n						\
 R = 2 E = eq_PB T = force C = 1 F = 0\n") ;
  fprintf(fic,"POIN\n				\
 1\n						\
 0.\n") ;
  fprintf(fic,"TEMP\n				\
 3\n						\
 0. 0.1 4.\n") ;
  fprintf(fic,"OBJE\n			      	\
 psi = 1.\n") ;
  fprintf(fic,"ALGO\n			       	\
 Iter = 1000 Tol = 1.e-10 Rec = 0\n	        \
 Dtini = 0.01 Dtmax = 0.01") ;
  fclose(fic) ;
}
