#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "model.h"

#define MODELINDEX  15
#define TITLE "Sols non satures"
#define AUTHORS "Dangla"

#include "OldMethods.h"

/* Macros */
#define NEQ   (1+dim)
#define NVE   (27)
#define NVA   (1)
#define E_liq (dim)
#define E_mec (0)
#define I_p_l (dim)
#define I_u   (0)
/* Fonctions */
static int    pm(char *) ;
static double pie(double,double,crbe_t) ;
static int    c15(double **,double **,double **,double *,double *,double *,elem_t,int,geom_t,double,double *) ;
static int    k15(double **,double **,double **,double *,double *,double *,elem_t,int,geom_t,double,double *) ;
static double rn15(double *,double *,double *,double *,mate_t) ;
static double cr15(double *,double,double *,double *,double *,mate_t) ;
/* Parametres */
static double gravite,kappa,mu,phi0,k_int,mu_l,rho_l,p_l0,p_g,rho_s,lambda,m,sig0,p_co0 ;

int pm(char *s)
{
  if(strcmp(s,"gravite") == 0) return (0) ;
  else if(strcmp(s,"kappa") == 0) return (1) ;
  else if(strcmp(s,"mu") == 0) return (2) ;
  else if(strcmp(s,"phi") == 0) return (3) ;
  else if(strcmp(s,"rho_l") == 0) return (4) ;
  else if(strcmp(s,"p_l0") == 0) return (5) ;
  else if(strcmp(s,"k_int") == 0) return (6) ;
  else if(strcmp(s,"mu_l") == 0) return (7) ;
  else if(strcmp(s,"lambda") == 0) return (8) ;
  else if(strcmp(s,"M") == 0) return (9) ;
  else if(strcmp(s,"p_g") == 0) return (10) ;
  else if(strcmp(s,"rho_s") == 0) return (11) ;
  else if(strcmp(s,"p_co") == 0) return (12) ;
  else if(strcmp(s,"sig0") == 0) return (13) ;
  else if(strcmp(s,"courbes") == 0) return (14) ;
  else return(-1) ;
}

int dm15(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int  i,n_donnees = 15 ;
  
  mat->neq = NEQ ;
  for(i=0;i<dim;i++) sprintf(mat->eqn[E_mec+i],"meca_%d",i+1) ;
  for(i=0;i<dim;i++) sprintf(mat->inc[I_u+i],"u_%d",i+1) ;
  strcpy(mat->eqn[E_liq],"liq") ;
  strcpy(mat->inc[I_p_l],"p_l") ;
  
  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}

int qm15(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est formee de (1 + dim) equations :\n\
\t 1. Conservation de la masse d\'eau (p_l)\n\
\t 2. Equilibre mecanique (u_1,u_2,u_3)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"gravite = 0    # gravite\n") ;
  fprintf(ficd,"rho_s = 2000   # masse volumique du sol sec\n") ;
  fprintf(ficd,"kappa = 0.004  # pente de la courbe de decharge elastique\n") ;
  fprintf(ficd,"lambda = 0.037 # pente de la courbe de consolidation vierge\n") ;
  fprintf(ficd,"mu = 1.e8      # module de cisaillement\n") ;
  fprintf(ficd,"M = 1.2        # pente de la courbe d\'etat critique\n") ;
  fprintf(ficd,"p_co = 18000   # pression de consolidation initiale\n") ;
  fprintf(ficd,"sig0 = -1000   # contrainte moyenne initiale\n") ;
  fprintf(ficd,"phi = 0.25     # porosite\n") ;
  fprintf(ficd,"rho_l = 1000   # masse volumique du fluide\n") ;
  fprintf(ficd,"p_l0 = 0       # pression initiale du fluide\n") ;
  fprintf(ficd,"k_int = 1e-20  # permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 0.001   # viscosite du liquide\n") ;
  fprintf(ficd,"p_g = 0        # pression de gaz\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier p_c S_l k_rl l\n") ;
  
  return(NEQ) ;
}

void tb15(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVE*el->fi->np ;
  el->n_ve = NVA*el->fi->np ;
}

void ch15(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  /* hydraulique */
  if(isdigit(cg.eqn[0]) && (atoi(cg.eqn) - 1) == E_liq) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  } else if(!strcmp(cg.eqn,el.mat->eqn[E_liq])) {
    for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
  }
}

void in15(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double eps[9],sig[9],phi,pl,pc,tre,m_l,w_l[3],k_l,gpl[3],sl,pp,pp_0,ppi ;
  double dfsds[3][3],dgsds[3][3],hm,crit ;
  int    i,p ;
  double *h,*dh ;
  double zero = 0. ;

  if(el.dim < dim) return ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  kappa   = el.mat->pr[pm("kappa")] ;
  mu      = el.mat->pr[pm("mu")] ;
  lambda  = el.mat->pr[pm("lambda")] ;
  m       = el.mat->pr[pm("M")] ;
  p_co0   = el.mat->pr[pm("p_co")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  sig0    = el.mat->pr[pm("sig0")] ;
  pp_0    = pie(p_l0,p_g,el.mat->cb[0]) ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    tre  = eps[0] + eps[4] + eps[8] ;
    /* pressions */
    pl  = param(u,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    pp  = pie(pl,p_g,el.mat->cb[0]) ;
    /* porosite */
    phi  = phi0 + tre ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* contraintes effectives */
    for(i=0;i<9;i++) sig[i] = 2*mu*eps[i] ;
    sig[0] += sig0 - 2*mu*tre/3. + pp_0 ;
    sig[4] += sig0 - 2*mu*tre/3. + pp_0 ;
    sig[8] += sig0 - 2*mu*tre/3. + pp_0 ;
    /* pression de consolidation */
    ppi = p_co0*courbe(pc,el.mat->cb[2]) ;
    /* critere */
    crit = cr15(sig,ppi,dfsds[0],dgsds[0],&hm,*el.mat) ;
    /* contraintes totales */
    sig[0] -= pp ;
    sig[4] -= pp ;
    sig[8] -= pp ;
    /* masse liquide */
    m_l = rho_l*phi*sl ;
    /* coefficient de transfert */
    k_l = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    /* flux */
    grad(x,u,dh,gpl,el.nn,dim,I_p_l) ;
    for(i=0;i<3;i++) w_l[i] = - k_l*gpl[i] ;
    w_l[dim-1] += k_l*rho_l*gravite ;
    /* rangement dans f */
    f[p*NVE] = m_l ;
    for(i=0;i<3;i++) f[p*NVE+1+i]  = w_l[i] ;
    for(i=0;i<9;i++) f[p*NVE+4+i]  = sig[i] ;
    for(i=0;i<3;i++) f[p*NVE+13+i] = zero ;
    f[p*NVE+13+dim-1] = (rho_s+m_l)*gravite ;
    for(i=0;i<9;i++) f[p*NVE+16+i]  = 0. ; /* def. pl. */
    f[p*NVE+25] = p_co0 ;
    f[p*NVE+26] = crit ;
  }
}

int ex15(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
  double pl,pc,k_l ;
  int    p ;
  double *h,*dh ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pressions */
    pl  = param(u,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* coefficient de transfert */
    k_l = rho_l*k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    /* rangement dans va */
    va[p*NVA] = k_l ;
  }
  return(0) ;
}

int ct15(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double pp,pp_n,pl,pl_n,pc ;
  double phi,m_l,w_l[3],k_l,gpl[3],sl,ppi,lpc ;
  double eps[9],sig[9],tre,sigm,p_co,eps_p[9],deps[9],deps_p[9] ;
  double eps_n[9],sig_n[9],trde,sigm_n ;
  double crit ;
  int    i,p ;
  double *h,*dh ;
  double zero = 0. ;
  
  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  kappa   = el.mat->pr[pm("kappa")] ;
  mu      = el.mat->pr[pm("mu")] ;
  lambda  = el.mat->pr[pm("lambda")] ;
  m       = el.mat->pr[pm("M")] ;
  p_co0   = el.mat->pr[pm("p_co")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u_1,h,dh,eps,el.nn,dim,geom,I_u) ;
    def(x,u_n,h,dh,eps_n,el.nn,dim,geom,I_u) ;
    for(i=0;i<9;i++) deps[i] =  eps[i] - eps_n[i] ;
    tre  = eps[0] + eps[4] + eps[8] ;
    trde = deps[0] + deps[4] + deps[8] ;
    /* pressions */
    pl  = param(u_1,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    pp  = pie(pl,p_g,el.mat->cb[0]) ;
    pl_n= param(u_n,h,el.nn,I_p_l) ;
    pp_n= pie(pl_n,p_g,el.mat->cb[0]) ;
    /* porosite */
    phi  = phi0 + tre ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* contraintes effectives */
    for(i=0;i<9;i++) sig_n[i] = f_n[p*NVE+4+i] ;
    sig_n[0] += pp_n ;
    sig_n[4] += pp_n ;
    sig_n[8] += pp_n ;
    sigm_n  = (sig_n[0] + sig_n[4] + sig_n[8])/3. ;
    sigm    = sigm_n*exp(-trde/((1.-phi0)*kappa)) ;
    for(i=0;i<9;i++) sig[i] = sig_n[i] + 2*mu*deps[i] ;
    sig[0] += sigm - sigm_n - 2*mu*trde/3. ;
    sig[4] += sigm - sigm_n - 2*mu*trde/3. ;
    sig[8] += sigm - sigm_n - 2*mu*trde/3. ;
    /* deformations plastiques */
    for(i=0;i<9;i++) eps_p[i] = f_n[p*NVE+16+i] ;
    /* pression de consolidation */
    p_co = f_n[p*NVE+25] ;
    lpc  = courbe(pc,el.mat->cb[2]) ;
    ppi  = p_co*lpc ;
    /* critere */
    crit = rn15(sig,sig_n,&ppi,deps_p,*el.mat) ;
    if(crit > 0.) {
      for(i=0;i<9;i++) eps_p[i] += deps_p[i] ;
      p_co = ppi/lpc ;
    }
    /* contraintes totales */
    sig[0] -= pp ;
    sig[4] -= pp ;
    sig[8] -= pp ;
    /* masse liquide */
    m_l = rho_l*phi*sl ;
    /* coefficient de transfert */
    k_l = va[p*NVA] ;
    /* flux */
    grad(x,u_1,dh,gpl,el.nn,dim,I_p_l) ;
    for(i=0;i<3;i++) w_l[i] = - k_l*gpl[i] ;
    w_l[dim-1] += k_l*rho_l*gravite ;
    /* rangement dans f_1 */
    f_1[p*NVE] = m_l ;
    for(i=0;i<3;i++) f_1[p*NVE+1+i] = w_l[i] ;
    for(i=0;i<9;i++) f_1[p*NVE+4+i] = sig[i] ;
    for(i=0;i<3;i++) f_1[p*NVE+13+i] = zero ;
    f_1[p*NVE+13+dim-1] = (rho_s+m_l)*gravite ;
    for(i=0;i<9;i++) f_1[p*NVE+16+i]  = eps_p[i] ;
    f_1[p*NVE+25] = p_co ;
    f_1[p*NVE+26] = crit ;
    /*
    printf("\ndeformations\n") ;
    for(i=0;i<9;i++) printf("%e ",eps[i]) ;
    printf("\npl = %e pp = %e",pl,pp) ;
    printf("\ncontraintes\n") ;
    for(i=0;i<9;i++) printf("%e ",sig[i]) ;
    */
  }
  return(0) ;
}

int mx15(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])
  int    i,j,n,m,dec ;
  double c[MAX_PGAUSS*100] ;
  double kb[9*MAX_NOEUDS*MAX_NOEUDS] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;

  /*
  ** Matrice de comportement
  */
  dec = c15(x,u_1,u_n,f_1,f_n,va,el,dim,geom,dt,c) ;
  /* 
  ** 1.  Mecanique
  */
  /* 1.2 rigidite du squelette */
  mxrig(kb,x,*el.fi,c,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) for(m=0;m<el.nn;m++) for(j=0;j<dim;j++) {
    K(E_mec+i+n*NEQ,I_u+j+m*NEQ) = kb[(i+n*dim)*el.nn*dim+j+m*dim] ;
  }
  /* 1.3 couplage mecanique */
  mxcpl(kb,x,*el.fi,c+81,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) for(m=0;m<el.nn;m++) {
    K(E_mec+i+n*NEQ,I_p_l+m*NEQ) = kb[(i+n*dim)*el.nn+m] ;
  }
  /*
  ** 2.  Hydraulique
  */
  /* 2.1 couplage diffusif */
  mxcpl(kb,x,*el.fi,c+90,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) for(m=0;m<el.nn;m++) {
    K(E_liq+m*NEQ,I_u+i+n*NEQ) = kb[(i+n*dim)*el.nn+m] ;
  }
  /* 2.2 accumulation */
  mxmass(kb,x,*el.fi,c+99,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(m=0;m<el.nn;m++) {
    K(E_liq+n*NEQ,I_p_l+m*NEQ) = kb[n*el.nn+m] ;
  }

  /*
  ** Matrice de conduction
  */
  dec = k15(x,u_1,u_n,f_1,f_n,va,el,dim,geom,dt,c) ;
  mxcond(kb,x,*el.fi,c,dim,dec,geom) ;
  for(n=0;n<el.nn;n++) for(m=0;m<el.nn;m++) {
    K(E_liq+n*NEQ,I_p_l+m*NEQ) += dt*kb[n*el.nn+m] ;
  }

  return(0) ;
#undef K
}

void rs15(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
#define R(n,i)    (r[(n)*NEQ+(i)])
  int    i,n ;
  double rb[3*MAX_NOEUDS],g1[MAX_PGAUSS] ;
  double zero = 0. ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = zero ;

  if(el.dim < dim) return ;

  /* 1. Mecanique */
  rscont(rb,x,*el.fi,f_1+4,dim,NVE,geom) ;
  for(n=0;n<el.nn;n++) for(i=0;i<dim;i++) R(n,E_mec+i) -= rb[i+n*dim] ;
  rsmass(rb,x,*el.fi,f_1+13+dim-1,dim,NVE,geom) ;
  for(n=0;n<el.nn;n++) R(n,E_mec+dim-1) -= -rb[n] ;
  /* 2. Hydraulique */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVE] - f_n[i*NVE] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_liq) -= rb[i] ;
  rsflux(rb,x,*el.fi,f_1+1,dim,NVE,geom) ;
  for(i=0;i<el.nn;i++) R(i,E_liq) -= -dt*rb[i] ;
#undef R
}

int so15(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double eps[9],sig[9],pl,pc,pp,w_l[3],sl,tre ;
  int    i,j,p,nso ;
  double *h,*dh,h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    donnees
  */
  phi0    = el.mat->pr[pm("phi")] ;
  
  /* initialisation */
  nso = 8 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* quantites exploitees en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;
  /* pressions */
  pl  = param(u,h_s,el.nn,I_p_l) ;
  pc  = p_g - pl ;
  pp  = pie(pl,p_g,el.mat->cb[0]) ;
  /* saturation */
  sl  = courbe(pc,el.mat->cb[0]) ;
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] =  pl ;
  strcpy(r[1].text,"deplacements") ; r[1].n = 3 ;
  for(i=0;i<dim;i++) r[1].v[i] = param(u,h_s,el.nn,I_u+i) ;
  strcpy(r[4].text,"saturation") ; r[4].n = 1 ;
  r[4].v[0] = sl ;
  strcpy(r[5].text,"pression-pi") ; r[5].n = 1 ;
  r[5].v[0] = pp ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* deformations */
    def(x,u,h,dh,eps,el.nn,dim,geom,I_u) ;
    tre  = eps[0] + eps[4] + eps[8] ;
    /* contraintes */
    for(i=0;i<9;i++) sig[i] = f[p*NVE+4+i] ;
    /* flux */
    for(i=0;i<3;i++) w_l[i] = f[p*NVE+1+i] ;
    /* quantites exploitees par element */
    strcpy(r[2].text,"flux-liquide") ; r[2].n = 3 ;
    for(i=0;i<dim;i++) r[2].v[i] += w_l[i]/el.fi->np ;
    strcpy(r[3].text,"contraintes") ; r[3].n = 9 ;
    for(i=0;i<9;i++) r[3].v[i] += sig[i]/el.fi->np ;
    strcpy(r[6].text,"indice-des-vides") ; r[6].n = 1 ;
    r[6].v[0] += tre/(1.-phi0)/el.fi->np ;
    strcpy(r[7].text,"sig_m") ; r[7].n = 1 ;
    r[7].v[0] += (sig[0] + sig[4] + sig[8])/3./el.fi->np ;
  }
  return (nso) ;
}

double pie(double pl,double pg,crbe_t cb)
{
  int    i,n_i = cb.n - 1 ;
  double pc,sl,sg,u ;
  double pc1 = cb.a[0],pc2 = cb.a[1],dpc = (pc2 - pc1)/n_i ;
  double zero = 0.,un = 1. ;

  pc  = pg - pl ;
  /* U */
  if(pc < pc1) {
    u = zero ;
    sl = cb.f[0] ;
  } else {
    if(pc > pc2) pc = pc2 ;
    u  = zero ;
    for(i=0;pc1+(i+1)*dpc<pc;i++) u += cb.f[i] + cb.f[i+1] ;
    u *= dpc*0.5 ;
    u += cb.f[0]*pc1 ;
    sl  = cb.f[i] + (cb.f[i+1] - cb.f[i])/dpc*(pc - pc1 - i*dpc) ;
    u += (cb.f[i] + sl)*0.5*(pc - pc1 - i*dpc) - sl*pc ;
  }
  /* pi */
  sg  = un - sl ;
  return(sl*pl + sg*pg - 2./3.*u) ;
}

int c15(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
#define C1(i,j,k,l)  (c1[(((i)*3+(j))*3+(k))*3+(l)])
#define B1(i,j)      (c1[(i)*3+(j)])
  double pl,pc,sl,dslsdpc,pp,lpc ;
  double sig[9],sigm,ppi,hm,crit ;
  double ppi_n ;
  double lame ;
  double dfsds[3][3],dgsds[3][3],fc[3][3],cg[3][3],fcg ;
  int    dec ;
  int    i,j,k,l,p ;
  double *h,*c1 ;
  double zero = 0. ;
  
  /*
    Donnees
  */
  gravite = el.mat->pr[pm("gravite")] ;
  kappa   = el.mat->pr[pm("kappa")] ;
  mu      = el.mat->pr[pm("mu")] ;
  lambda  = el.mat->pr[pm("lambda")] ;
  m       = el.mat->pr[pm("M")] ;
  p_co0   = el.mat->pr[pm("p_co")] ;
  phi0    = el.mat->pr[pm("phi")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  rho_l   = el.mat->pr[pm("rho_l")] ;
  rho_s   = el.mat->pr[pm("rho_s")] ;
  p_l0    = el.mat->pr[pm("p_l0")] ;
  p_g     = el.mat->pr[pm("p_g")] ;
  
  dec = 100 ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pressions */
    pl  = param(u_1,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    pp  = pie(pl,p_g,el.mat->cb[0]) ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    dslsdpc = dcourbe(pc,el.mat->cb[0]) ;
    /* contraintes effectives */
    for(i=0;i<9;i++) sig[i] = f_1[p*NVE+4+i] ;
    sig[0] += pp ;
    sig[4] += pp ;
    sig[8] += pp ;
    sigm  = (sig[0] + sig[4] + sig[8])/3. ;
    /* pression de consolidation */
    lpc    = courbe(pc,el.mat->cb[2]) ;
    ppi    = f_1[p*NVE+25]*lpc ;
    ppi_n  = f_n[p*NVE+25]*lpc ;

    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* mecanique */
    lame = -sigm/((1.-phi0)*kappa) - 2*mu/3. ;
    for(i=0;i<3;i++) for(j=0;j<3;j++) {
      C1(i,i,j,j) += lame ;
      C1(i,j,i,j) += mu ;
      C1(i,j,j,i) += mu ;
    }
    /* critere */
    crit = f_1[p*NVE+26] ;
    if(crit >= 0.) {
      cr15(sig,ppi,dfsds[0],dgsds[0],&hm,*el.mat) ;
      fcg = hm ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) {
	fc[i][j] = 0. ;
	cg[i][j] = 0. ;
	for(k=0;k<3;k++) for(l=0;l<3;l++) {
	  fc[i][j] += dfsds[k][l]*C1(l,k,i,j) ;
	  cg[i][j] += C1(i,j,k,l)*dgsds[l][k] ;
	  fcg += dfsds[j][i]*C1(i,j,k,l)*dgsds[l][k] ;
	}
      }
      fcg = 1./fcg ;
      for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) for(l=0;l<3;l++) {
	C1(i,j,k,l) -= cg[i][j]*fc[k][l]*fcg ;
      }
    }
    c1 += 81 ;
    for(i=0;i<3;i++) B1(i,i) = -(sl+pc*dslsdpc/3.) ;
    
    c1 += 9 ;
    /* hydraulique */
    for(i=0;i<3;i++) B1(i,i) = rho_l*sl ;
    c1 += 9 ;
    c1[0] = -rho_l*phi0*dslsdpc ;
  }
  return(dec) ;
  
#undef C1
#undef B1
}

int k15(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  int    dec ;
  double *c1 ;
  int    i,p ;
  double zero = 0. ;
  
  dec = 9 ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* permeabilite */
    c1[0] = va[p*NVA] ;
    c1[4] = va[p*NVA] ;
    c1[8] = va[p*NVA] ;
  }
  return(dec) ;
}


double cr15(double *sig,double pc,double *dfsds,double *dgsds,double *hm,mate_t mat)
/* Critere de Cam-Clay */
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit,m2,v,dev ;
  int    i ;
  /*
    Donnees
  */
  kappa   = mat.pr[pm("kappa")] ;
  lambda  = mat.pr[pm("lambda")] ;
  m       = mat.pr[pm("M")] ;
  phi0    = mat.pr[pm("phi")] ;
  m2      = m*m ;
  v       = 1./(lambda-kappa) ;
  /* 
     Le critere
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  /*
    Les gradients
  */
  for(i=0;i<9;i++) {
    dev      = sig[i] - p*id[i] ;
    dfsds[i] = (2*p+pc)*id[i]/3. + 3./m2*dev ;
    dgsds[i] = dfsds[i] ;
  }
  /* Le module d'ecrouissage */
  *hm = v/(1 - phi0)*p*(2*p + pc)*pc ;
  return(crit) ;
}


double rn15(double *sig,double *sig_n,double *p_co,double *eps_p,mate_t mat)
/* Critere de Cam-Clay : return mapping algorithm (sig,p_co,eps_p)
   (d apres Borja & Lee 1990, modifie par Dangla)
*/
{
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double tol = 1.e-8 ;
  double p_t,q_t,p,q,p_n,pc,pc_n,crit,fcrit,m2,v,dev,a ;
  double df,dfsdp,dfsdq,dfsdpc,dl,dlsdp,dqsdp,dpcsdp,dfsds ;
  int    i,nf ;
  /*
    Donnees
  */
  kappa   = mat.pr[pm("kappa")] ;
  mu      = mat.pr[pm("mu")] ;
  lambda  = mat.pr[pm("lambda")] ;
  m       = mat.pr[pm("M")] ;
  phi0    = mat.pr[pm("phi")] ;
  m2      = m*m ;
  v       = 1./(lambda-kappa) ;
  /* 
     Le critere
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*j2(sig)) ;
  pc   = *p_co ;
  crit = q*q/m2 + p*(p + pc) ;
  /*
    Algorithme de projection (closest point projection)
    Une seule boucle iterative pour le calcul de p, racine de
    q*q/m2 + p*(p + pc) = 0
    Les autres variables (pc,q,dl) sont donnees explicitement par p.
  */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  pc_n  = pc ;
  p_n   = (sig_n[0] + sig_n[4] + sig_n[8])/3. ;
  if(crit > 0.) {
    fcrit = crit ;
    nf    = 0 ;
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      dfsdp  = 2*p + pc ;
      dfsdq  = 2*q/m2 ;
      dfsdpc = p ;
      dpcsdp = -v*kappa*pc/p ;
      dlsdp  = ((1.-phi0)*kappa/p - dl*(2+dpcsdp))/dfsdp ;
      dqsdp  = -q*6*mu/(m2+6*mu*dl)*dlsdp ;
      df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      p     -= fcrit/df ;
      /* Les variables (pc,dl,q) sont explicites en p */
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1.-phi0)*kappa*log(p/p_t)/(2*p+pc) ;
      q      = q_t*m2/(m2+6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      if(nf++ > 20) {
	printf("pas de convergence (rn15)") ;
	exit(0) ;
      }
    }
  }
  /*
    Les contraintes et deformations plastiques
  */
  a = 1./(1. + 6*mu/m2*dl) ;
  for(i=0;i<9;i++) {
    dev      = a*(sig[i] - p_t*id[i]) ;
    sig[i]   = p*id[i] + dev ;
    dfsds    = (2*p+pc)*id[i]/3. + 3./m2*dev ;
    eps_p[i] = dl*dfsds ;
  }
  /* La pression de consolidation */
  *p_co = pc ;
  return(crit) ;
}

#undef NEQ
#undef NVE
#undef NVA
#undef E_liq
#undef E_mec
#undef I_p_l
#undef I_u

