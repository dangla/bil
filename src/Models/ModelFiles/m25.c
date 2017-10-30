#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"

#define MODELINDEX  25
#define TITLE "Corrosion en surface a l\'oxygene (3D)"
#define AUTHORS "Dridi"

#include "OldMethods.h"

/* Macros */
#define NEQ      (5)

#define NVI      (2)
#define NVE      (0)

#define E_h      (0)
#define E_i      (2)
#define E_o2     (1)
#define E_m      (3)
#define E_fe     (4)

#define I_p_l    (0)
#define I_psi    (2)
#define I_c_o2   (1)
#define I_psi_m  (3)
#define I_c_fe   (4)

/* Fonctions */
static int    pm(const char *s) ;
static int    c25(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c) ;
/* Parametres */
static double i0_a,i0_c,tafel_a,tafel_c,phi,c0_o2,p_g=0.  ;

int pm(const char *s)
{
  if(strcmp(s,"i_a") == 0) return (0) ;
  else if(strcmp(s,"i_c") == 0) return (1) ;
  else if(strcmp(s,"tafel_a") == 0) return (2) ;
  else if(strcmp(s,"tafel_c") == 0) return (3) ;
  else if(strcmp(s,"phi") == 0) return (4) ;
  else if(strcmp(s,"c_o2") == 0) return (5) ;
  else if(strcmp(s,"courbes") == 0) return (6) ;
  else
    { printf("donnee \"%s\" non connue (pm25)\n",s) ; exit(0) ; }
}

int dm25(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 7 ;

  mat->neq = NEQ ;

  strcpy(mat->eqn[E_h]    ,"E_h") ;
  strcpy(mat->eqn[E_i]    ,"E_i") ;
  strcpy(mat->eqn[E_o2]   ,"E_o2") ;
  strcpy(mat->eqn[E_m]    ,"E_m") ;
  strcpy(mat->eqn[E_fe]   ,"E_fe") ;

  strcpy(mat->inc[I_p_l]  ,"p_l") ;
  strcpy(mat->inc[I_psi]  ,"psi") ;
  strcpy(mat->inc[I_c_o2] ,"c_o2") ;
  strcpy(mat->inc[I_c_fe] ,"c_fe") ;
  strcpy(mat->inc[I_psi_m],"psi_m") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}


int qm25(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 5 CL portant sur :\n\
\t 1. Le flux de H (p_l)\n\
\t 2. Le flux cathodique ou consommation de O2 (c_o2)\n\
\t 3. Le flux anodique ou production de Fe (c_fe)\n\
\t 4. Le courant (psi)\n\
\t 5. Egalite des courants anodique et cathodique (psi_m)\n") ;
  
  printf("\n\
Exemple de donnees :\n\n") ;
  
  fprintf(ficd,"i_a = 1.e-6       # Constante cinetique de la reaction anodique\n") ;
  fprintf(ficd,"i_c = 2.5e-6      # Constante cinetique de la reaction cathodique\n") ;
  fprintf(ficd,"tafel_a = 60.e-3  # Coefficient de Tafel anodique\n") ;
  fprintf(ficd,"tafel_c = 160.e-3 # Coefficient de Tafel cathodique\n") ;
  fprintf(ficd,"c_o2 = 0.25       # Concentration en oxygene de reference\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier : p_c S_l\n") ;
  fprintf(ficd,"phi = 0.28        # La porosite\n") ;
  
  return(NEQ) ;
}

void tb25(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI*el->fi->np ;
  el->n_ve = NVE*el->fi->np ;
}

void ch25(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
}

void in25(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double v_a,v_c;
  double pl,pc,sl,c_o2,psi,psi_m ;
  int    p ;
  double *h ;

  if(el.dim != dim - 1) return ;

  /*
  Donnees
  */
  i0_a    = el.mat->pr[pm("i_a")] ;
  i0_c    = el.mat->pr[pm("i_c")] ;
  tafel_a = el.mat->pr[pm("tafel_a")] ;
  tafel_c = el.mat->pr[pm("tafel_c")] ;
  phi     = el.mat->pr[pm("phi")] ;
  c0_o2   = el.mat->pr[pm("c_o2")] ;


  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pression */
    pl  = param(u,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* concentrations */
    c_o2    = param(u,h,el.nn,I_c_o2) ;
    /* potentiel electrique */
    psi   = param(u,h,el.nn,I_psi) ;
    psi_m = param(u,h,el.nn,I_psi_m) ;  
    /* vitesses des reactions partielles */
    v_a  = i0_a*phi*sl*exp((psi_m - psi)/tafel_a) ;
    v_c  = i0_c*c_o2/c0_o2*phi*sl*exp(-(psi_m - psi)/tafel_c) ;
    /* rangement dans f */
    f[p*NVI]   = v_a ;   /* vitesse anodique */
    f[p*NVI+1] = v_c ;   /* vitesse cathodique */
  }
}

int ex25(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
  return(0) ;
}

int ct25(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double v_a,v_c;
  double pl,pc,c_o2,psi,psi_m,sl ;
  int    p ;
  double *h ;

  if(el.dim != dim - 1) return(0) ;

  /*
  Donnees
  */
  i0_a    = el.mat->pr[pm("i_a")] ;
  i0_c    = el.mat->pr[pm("i_c")] ;
  tafel_a = el.mat->pr[pm("tafel_a")] ;
  tafel_c = el.mat->pr[pm("tafel_c")] ;
  phi     = el.mat->pr[pm("phi")] ;
  c0_o2   = el.mat->pr[pm("c_o2")] ;


  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pression */
    pl  = param(u_1,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* concentrations */
    c_o2    = param(u_1,h,el.nn,I_c_o2) ;
    /* potentiel electrique */
    psi   = param(u_1,h,el.nn,I_psi) ;
    psi_m = param(u_1,h,el.nn,I_psi_m) ;  
    /* vitesses des reactions partielles */
    if(c_o2 < 0.) {
      printf("\n\
      c_o2    = %e\n",c_o2) ;
      return(-1) ;
    }
    v_a  = i0_a*phi*sl*exp((psi_m - psi)/tafel_a) ;
    v_c  = i0_c*c_o2/c0_o2*phi*sl*exp(-(psi_m - psi)/tafel_c) ;
    /* rangement dans f_1 */
    f_1[p*NVI]   = v_a ;   /* vitesse anodique */
    f_1[p*NVI+1] = v_c ;   /* vitesse cathodique */
  }
  return(0) ;
}

int mx25(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define NN        (el.nn)
  int    i,dec ;
  double c[MAX_PGAUSS*NEQ*NEQ] ;
  double kb[NEQ*NEQ*MAX_NOEUDS*MAX_NOEUDS] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<NN*NN*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim != dim - 1) return(0) ;

  /* Matrice de conduction */
  dec = c25(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxcmss(kb,x,*el.fi,c,dim,dec,geom,NEQ) ;
  for(i=0;i<NN*NN*NEQ*NEQ;i++) k[i] += dt*kb[i] ;
  return(0) ;
#undef NN
}

void rs25(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  int    i ;
  double rb[MAX_NOEUDS],g1[MAX_PGAUSS] ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0. ;

  if(el.dim != dim - 1) return ;
  
   /* courant : i.n = 2(V_c - V_a) */
  for(i=0;i<el.fi->np;i++) g1[i] = 2*(f_1[i*NVI+1] - f_1[i*NVI]) ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_i] = - dt*rb[i] ;
  
  /* o2 : w_o2.n = 0.5*V_c */
  for(i=0;i<el.fi->np;i++) g1[i] = 0.5*f_1[i*NVI+1] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_o2] = - dt*rb[i] ;
  
  /* fer : w_fe.n = -V_a */
  for(i=0;i<el.fi->np;i++) g1[i] = - f_1[i*NVI] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fe] = - dt*rb[i] ;
  
  /* metal : 2(Vc - V_a) = 0 */
  for(i=0;i<el.fi->np;i++) g1[i] = 2*(f_1[i*NVI+1] - f_1[i*NVI]) ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_m] = - dt*rb[i] ;

}


int so25(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)

/* Les valeurs exploitees (s) */
{
  double v_a,v_c ;
  double pl,pc,c_o2,psi,psi_m,sl ;
  int    i,j,nso ;
  double h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;
  
  if(el.dim != dim - 1) return(0) ;

  /* Donnees */
  i0_a    = el.mat->pr[pm("i_a")] ;
  i0_c    = el.mat->pr[pm("i_c")] ;
  tafel_a = el.mat->pr[pm("tafel_a")] ;
  tafel_c = el.mat->pr[pm("tafel_c")] ;
  phi     = el.mat->pr[pm("phi")] ;
  c0_o2   = el.mat->pr[pm("c_o2")] ;
  
  /* initialisation */
  nso = 7 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pression */
  pl = param(u,h_s,el.nn,I_p_l) ;
  pc = p_g - pl ;
  /* saturation */
  sl = courbe(pc,el.mat->cb[0]) ;
  /* concentrations */
  c_o2    = param(u,h_s,el.nn,I_c_o2) ;
  /* potentiel electrique */
  psi    = param(u,h_s,el.nn,I_psi) ;
  psi_m  = param(u,h_s,el.nn,I_psi_m) ;  
  /* vitesses des reactions partielles */
  v_a  = i0_a*phi*sl*exp((psi_m - psi)/tafel_a)  ;
  v_c  = i0_c*c_o2/c0_o2*phi*sl*exp(-(psi_m - psi)/tafel_c) ;
  /* quantites exploitees par element */
  strcpy(r[0].text,"pression_liquide") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"molarite_O2") ; r[1].n = 1 ;
  r[1].v[0] = c_o2 ;
  strcpy(r[2].text,"potentiel_electrique du beton") ; r[2].n = 1 ;
  r[2].v[0] = psi ;
  strcpy(r[3].text,"potentiel_electrique du metal") ; r[3].n = 1 ;
  r[3].v[0] = psi_m ;
  strcpy(r[4].text,"Vitesse anodique") ; r[4].n = 1 ;
  r[4].v[0] = v_a ;
  strcpy(r[5].text,"vitesse cathodique") ; r[5].n = 1 ;
  r[5].v[0] = v_c ;
  strcpy(r[6].text,"saturation_surfacique") ; r[6].n = 1 ;
  r[6].v[0] = sl ;  
  return(nso) ;

}
int c25(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
  int    dec ;
  double *h,*c1 ;
  double v_a,v_c ;
  double pl,pc,sl ;
  double c_o2,psi,psi_m ;
  int    i,p ;
  double zero = 0. ;
  double dvasdpl,dvasdpsi,dvasdpsim,dvcsdpl,dvcsdpsi,dvcsdpsim,dvcsdco2,dslsdpc ;

  /* Donnees */
  i0_a    = el.mat->pr[pm("i_a")] ;
  i0_c    = el.mat->pr[pm("i_c")] ;
  tafel_a = el.mat->pr[pm("tafel_a")] ;
  tafel_c = el.mat->pr[pm("tafel_c")] ;
  phi     = el.mat->pr[pm("phi")] ;
  c0_o2   = el.mat->pr[pm("c_o2")] ;

  dec = NEQ*NEQ ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pression */
    pl  = param(u_1,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* concentrations */
    c_o2    = param(u_1,h,el.nn,I_c_o2) ;
    /* potentiel electrique */
    psi    = param(u_1,h,el.nn,I_psi) ;
    psi_m  = param(u_1,h,el.nn,I_psi_m) ;  
    /* vitesses partielles */
    v_a    = i0_a*phi*sl*exp((psi_m - psi)/tafel_a) ;
    v_c    = i0_c*c_o2/c0_o2*phi*sl*exp(-(psi_m - psi)/tafel_c) ;
    /* derivees */
    dslsdpc   = dcourbe(pc,el.mat->cb[0]) ;
    dvasdpl   = -i0_a*phi*dslsdpc*exp((psi_m - psi)/tafel_a) ;
    dvasdpsi  = -v_a/tafel_a ;
    dvasdpsim = v_a/tafel_a ;

    dvcsdpl   = -i0_c*c_o2/c0_o2*phi*dslsdpc*exp(-(psi_m - psi)/tafel_c) ;
    dvcsdpsi  = v_c/tafel_c ;
    dvcsdco2  = v_c/c_o2 ;
    dvcsdpsim = -v_c/tafel_c ;    

    /* courant : i.n = 2(V_c - V_a) */
    c1[E_i*NEQ+I_p_l]    = 2*(dvcsdpl - dvasdpl) ; 
    c1[E_i*NEQ+I_psi]    = 2*(dvcsdpsi - dvasdpsi) ; 
    c1[E_i*NEQ+I_c_o2]   = 2*(dvcsdco2) ;
    c1[E_i*NEQ+I_psi_m]  = 2*(dvcsdpsim - dvasdpsim) ;

    /* o2 : w_o2.n = 0.5V_c */
    c1[E_o2*NEQ+I_p_l]   = 0.5*dvcsdpl ;
    c1[E_o2*NEQ+I_psi]   = 0.5*dvcsdpsi ;
    c1[E_o2*NEQ+I_c_o2]  = 0.5*dvcsdco2 ;
    c1[E_o2*NEQ+I_psi_m] = 0.5*dvcsdpsim ;

    /* fer : w_fer.n = - V_a */
    c1[E_fe*NEQ+I_p_l]   = - dvasdpl ;
    c1[E_fe*NEQ+I_psi]   = - dvasdpsi ;
    c1[E_fe*NEQ+I_psi_m] = - dvasdpsim ;

    /* metal : 2(V_c - V_a) = 0  (i.n = 0) */
    c1[E_m*NEQ+I_p_l]    = 2*(dvcsdpl - dvasdpl) ; 
    c1[E_m*NEQ+I_psi]    = 2*(dvcsdpsi - dvasdpsi) ; 
    c1[E_m*NEQ+I_c_o2]   = 2*(dvcsdco2) ;
    c1[E_m*NEQ+I_psi_m]  = 2*(dvcsdpsim - dvasdpsim) ;

  }
  return(dec) ;
}
