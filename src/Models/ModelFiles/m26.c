/*
	Modele (ancien 64) : quatre especes ionique (Cl, Na, K, OH, H)
	Diffusion: Loi de Frick + Migration + reaction chimique (H2O = H+ + OH-)
	Isotherme de fixation: Freundlich
	D est calculee D_oh=D_cl/Do_cl*Do_oh
	Tenir en compte: Activitee de la solution 
	Non stature
	Activite
	pression gaz constante
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"

#define MODELINDEX  26
#define TITLE "Sechage isotherme avec sels"
#define AUTHORS "Nguyen"

#include "OldMethods.h"

/* Macros */
#define NEQ 	(2)
#define NVI     (11)
#define NVE     (12)

#define E_VW	(0)
#define E_Cl	(1)

#define I_Pv    (0)
#define I_XS    (1) /* Salt */

#define P_v(n)      (u[(n)][I_Pv])
#define XS(n)       (u[(n)][I_XS])

#define N_V(n)      (f[(n)])
#define N_W(n)      (f[(n+2)])
#define N_Cl(n)     (f[(n+4)])
#define W_V         (f[6])
#define W_W         (f[7])
#define W_Cl        (f[8])
#define C_Cl(n)     (f[(9+n)])

#define N_Vn(n)      (f_n[(n)])
#define N_Wn(n)      (f_n[(n+2)])
#define N_Cln(n)     (f_n[(n+4)])
#define C_Cln(n)     (f_n[(9+n)])

#define KD_W        (va[(0)])
#define KF_VA       (va[(1)])
#define KF_Cl       (va[(2)])
#define KF_Na       (va[(3)])
#define KD_Cl       (va[(4)])
#define KD_Na       (va[(5)])
#define WA_Cl       (va[(6)])
#define WA_Na       (va[(7)])
#define c_e(n)      (va[(8+n)])
#define lna_es(n)   (va[(10+n)])

/* valence */
#define z_cl    (-1.)
#define z_na    (+1.)

/* volumes molaires */
#define v_nacl  (0.)/*(24.5e-6)*/
#define v_h2o   (1.80e-5)
#define v_na    (1.87e-5)
#define v_cl    (2.52e-6)

/* coefficients de diffusion dans l'eau dilue */
#define do_cl   (2.032e-9)
#define do_na   (1.334e-9)

/* constante d'equilibre */
#define k_s     (6.e3)

/* constante physique */
#define R_g     (8.3143)

/* viscosites */
#define mu_g    (1.8e-5)
#define mu_w    (1.002e-3)

/* Masses molaires */
#define M_v     (1.8e-2)
#define M_a     (2.896e-2)

/* autres */
#define p_atm   (1.01325e5)
#define rho_w   (1.0e3)
#define gravite (9.81)

/* Fonctions */
static int    pm(const char *s) ;
static double activite(double,double,double,double,double *) ;
static double lng_LinLee(double,double,double,double,double,double) ;
static double lna_i(double,double,double,double,double,double) ;
static double lng_TQN(double,double,double,double,double,double,double,double) ;

/* Parametres */
static double phio,r_d,T,k_int,p_vs,n,anpha,m,A,S_max,Pc_min ;
static double s_srf,s_so = 0. ;
static int    noc_phi = 0 ;
static double d_cl,d_va ;

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0) return (0) ;
  else if(strcmp(s,"D_Cl") == 0) return (1);  
  else if(strcmp(s,"r_d") == 0) return (2);  
  else if(strcmp(s,"T") == 0) return (3);
  else if(strcmp(s,"k_int") == 0) return (4);
  else if(strcmp(s,"d_va") == 0) return (5);
  else if(strcmp(s,"B_h") == 0) return (6);
  else if(strcmp(s,"n") == 0) return (7);
  else if(strcmp(s,"anpha") == 0) return (8);
  else if(strcmp(s,"S_max") == 0) return (9);
  else if(strcmp(s,"Pc_min") == 0) return (10);
  else if(strcmp(s,"phi_x") == 0) return (11);
  else if(strcmp(s,"s_srf") == 0) return (12) ;
  else if(strcmp(s,"s_so") == 0) return (13);
  else {
    printf("donnee \"%s\" non connue (pm26)\n",s) ; exit(0) ;
  }
}

int dm26(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    nd = 14,nc = 0;

  if(dim > 1) arret("dm26 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_VW],"liq") ;
  strcpy(mat->eqn[E_Cl],"sel") ;
  strcpy(mat->inc[I_Pv],"p_v") ;
  strcpy(mat->inc[I_XS],"x_s") ;

  lit_mate(mat,ficd,pm,nd,nc) ;

  /*
  mat->n = nd ;
  mat->nc = 0 ;
  for(i=0;i<nd+nc-ns;i++) {
    if(!fgets(line,sizeof(line),ficd)) arret("dmat (1) : erreur ou fin de fichier") ;
    sscanf(line," %[^= ] =",mot) ;
    p = strchr(line,'=') + 1 ;
    if(!strncmp(mot,"courbes",6)) {
  *//*
      fscanf(ficd,"%s",nom) ; 
      fscanf(ficd,"%d",&n_pa) ;	
      for (j=0;j<n_pa;j++) fscanf(ficd,"%lf",para+j);
      if (n_pa!=0) ecrit_courbe(nom,n_pa,para);	
      lit_courbe(mat,nom) ;
    *//*
      ecrit_courbe(p);	
      lit_courbe(mat,line) ;
    } else sscanf(p,"%lf",mat->pr+pm(mot)) ;
  }
      */

  return(mat->n) ;
}

int qm26(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de 2 equations :\n\
\t 1. Conservation de la masse d\'eau (p_v)\n\
\t 2. Conservation de la masse de sel (x_s)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = # Porosite\n") ;
  fprintf(ficd,"D_Cl = # Diffusion effective de Cl\n") ;
  fprintf(ficd,"T = # Temperature\n") ;
  fprintf(ficd,"r_d = # Rapport des tortuosites des anions et des cations\n") ;
  fprintf(ficd,"activite = # model de calcul activite de la solution\n") ;
  fprintf(ficd,"isotherme = # type de l'isotherme\n") ;

  return(NEQ) ;
}

void tb26(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ; /* implicite */
  el->n_ve = NVE ; /* explicite */
}

void ch26(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
#define R(n,i)     (r[(n)*NEQ+(i)])

  double B_h,Mt,xs,ys,s_s,phi;
  int    i,ieq ;
  double zero = 0. ;

 /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0.; 


  phio    = el.mat->pr[pm("porosite")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  k_int	  = el.mat->pr[pm("k_int")] ;
  d_va	  = el.mat->pr[pm("d_va")];
  B_h	  = el.mat->pr[pm("B_h")];

  s_so    = el.mat->pr[pm("s_so")] ;
  s_srf   = el.mat->pr[pm("s_srf")] ;
  noc_phi = floor(el.mat->pr[pm("phi_x")]+0.5) ;  

  xs   = XS(0); ys   = fabs(xs);
  s_s  = 0.5*(xs+ys)*s_srf;	 		

  {
    Function_t *fn = Material_GetFunction(Element_GetMaterial(&el)) ;
    phi = ((noc_phi<1) ? phio : fonction(x[0][0],fn[noc_phi-1])) + (s_so-s_s)*24.5e-6;
  }
  if (phi<0.) phi  = 0.;

  /* on calcule le numero de l'equation */
  if(isdigit(cg.eqn[0])) { /* donne sous forme numerique */
    ieq  = atoi(cg.eqn) - 1 ;
  } else {                 /* donne sous forme alphabetique */
    for(ieq=0;ieq<NEQ;ieq++) if(!strcmp(cg.eqn,el.mat->eqn[ieq])) break ;
    if(ieq == NEQ) arret("ch26 (1) : equation non connue") ;
  }

  if(ieq == E_VW) {
    if (B_h==zero) {
      Mt = fonction(t,*cg.fn) - fonction(t-dt,*cg.fn); 
      R(0,E_VW) =  Mt*champ(x[0],dim,*cg.ch);
    } else {
      R(0,E_VW) = - rho_w*B_h*dt*(P_v(0)-champ(x[0],dim,*cg.ch))/p_vs*pow(phi/phio,3.)*pow((1.-phio)/(1.-phi),2.);
    }

  } else {printf("\nchargement non prevu (ch26)\n") ; exit(0) ;}
 
#undef R
}

void in26(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define lna_cl(n)    (lna[n][0])
#define lna_na(n)    (lna[n][1])

  double s_w,s_g,c_cl,c_na,xs,ys,s_s;
  double grd_pw,grd_pvg, RT,grd_cl,grd_na,grd_si;
  double dx,fa,top,ton,g;
  double p_v,p_l[2],p_c, k_rw,k_re;
  int    i,ipc_sw;
  double lna[2][2],S_Dcz2,phi[2],phif,B_h;
  double zero = 0.,un = 1.,deux = 2. ;
  Function_t *fn = Material_GetFunction(Element_GetMaterial(&el)) ;

  if(el.dim < dim) return ;
  /*
    Donnees
  */
  phio     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ; 
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  k_int	  = el.mat->pr[pm("k_int")] ;
  d_va	  = el.mat->pr[pm("d_va")];

  n	  = el.mat->pr[pm("n")];
  anpha	  = el.mat->pr[pm("anpha")];
  S_max	  = el.mat->pr[pm("S_max")];
  Pc_min  = el.mat->pr[pm("Pc_min")];
  m       = un-un/n;
  A       = rho_w*gravite/anpha;
  RT      = R_g*T;

  s_srf    = el.mat->pr[pm("s_srf")] ;
  noc_phi = floor(el.mat->pr[pm("phi_x")]+0.5) ;  

  B_h	  = el.mat->pr[pm("B_h")];
  if (B_h==zero) {
    ipc_sw = 0;
  } else  {
    ipc_sw = floor(B_h+0.5) ;  
  }
  
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    xs   = XS(i);    ys   = fabs(xs);
    C_Cl(i) = k_s*(0.5*(xs-ys)+1.);
    c_cl = C_Cl(i);
    c_na = c_cl;
    s_s  = 0.5*(xs+ys)*s_srf;	 		

    el.mat->pr[pm("s_so")]        = s_s;
    s_so = el.mat->pr[pm("s_so")] ;

    c_e(i)    = (un - (c_cl*v_cl + c_na*v_na))/v_h2o ;   
    lna_es(i) = activite(c_cl,c_na,c_e(i),T,lna[i]);
    
    phi[i] = ((noc_phi<1) ? phio : fonction(x[i][0],fn[noc_phi-1])) + (s_so-s_s)*v_nacl;
    

    p_v = P_v(i);
    p_l[i] = RT*c_e(i)*(log(p_v/p_vs)-lna_es(i))+p_atm;
    p_c = p_atm-p_l[i];

    if (ipc_sw<1) {
      if (p_c>=0) {
        s_w  = pow(pow(p_c/A,n)+S_max,-m);
      } else {
        s_w  =pow(S_max,-m)-p_c*(un-pow(S_max,-m))/fabs(Pc_min);
      }
    } else {
      s_w  = fonction(p_c,fn[ipc_sw-1]);
    }
    s_g  = un -s_w;
        
    N_V(i)   = phi[i]*M_v/RT*p_v*s_g;  
    N_W(i)   = phi[i]*s_w*c_e(i)*M_v;  

    N_Cl(i)  = phi[i]*s_w*c_cl+s_s;  

    if(s_w<0.||s_w>1.||xs<=-1.) {
      printf("\n x   = %e\n s_w     = %e\n p_c     = %e\n",x[i][0],s_w,p_c) ;
    }

  }

  /* Coefficient de transfert */

    p_v     = (P_v(0)+P_v(1))/deux ;
    p_c     = p_atm-0.5*(p_l[0]+p_l[1]);
    
    if (ipc_sw<1) {
      if (p_c>=0) {
	s_w  = pow(pow(p_c/A,n)+S_max,-m);
      } else {
	s_w  =pow(S_max,-m)-p_c*(un-pow(S_max,-m))/fabs(Pc_min);
      }
    } else {
      s_w  = fonction(p_c,fn[ipc_sw-1]);
    }    
  
    s_g     = un -s_w;    if (s_g<zero) s_g = zero;

    phif    = (phi[0]+phi[1])/deux ;
    k_rw    = pow(s_w,0.5)*pow(un-pow(un-pow(s_w,un/m),m),deux);
    /*  k_re    = un/(un+pow(-10.*log(s_w),2.)); */
    k_re =k_rw;
    fa      = pow(s_g,10./3.)*pow(phif,4./3.);

    /*    k_int = k_int*pow(phif/phio,3.)*pow((un-phio)/(un-phif),2.); */

    KD_W  = 0.5*(c_e(0)+c_e(1))*M_v*k_int/mu_w*k_rw;  
    KF_VA = d_va*fa/RT/p_atm;  

    c_cl     = (C_Cl(0)+C_Cl(1))/deux ;
    c_na     = c_cl;

    KD_Cl = k_int/mu_w*k_rw*c_cl;
    KD_Na = k_int/mu_w*k_rw*c_na;

    g = (0.001+0.07*pow(phif,2.)+1.8*((phif<=0.18) ? 0.: 1.)*pow(phif-0.18,2.));
    g = phif/phio*g/(0.001+0.07*pow(phio,2.)+1.8*((phio<=0.18) ? 0.: 1.)*pow(phio-0.18,2.));

    d_cl = d_cl*g;

    ton = d_cl/do_cl*k_re;
    top = ton/r_d;

    KF_Cl  = ton*do_cl;
    KF_Na  = top*do_na;

 /* Flux */
    dx     = x[1][0] - x[0][0];
    WA_Cl   = - c_cl*KF_Cl*(lna_cl(1)-lna_cl(0))/dx;
    WA_Na   = - c_na*KF_Na*(lna_na(1)-lna_na(0))/dx;

    grd_cl = (C_Cl(1) - C_Cl(0))/dx ;
    grd_na = grd_cl; 
    S_Dcz2 = KF_Cl*pow(z_cl,2.)+KF_Na*pow(z_na,2.);  
    grd_si =-(z_cl*(KF_Cl*grd_cl-WA_Cl)+z_na*(KF_Na*grd_na-WA_Na))/S_Dcz2;
		
    grd_pw = (p_l[1]-p_l[0])/dx;
    grd_pvg= (P_v(1) - P_v(0))/dx;

    W_W    = - KD_W*grd_pw;
    W_V    = - KF_VA*M_v*grd_pvg;

/*    W_Cl   = - KF_Cl*grd_cl  - KD_Cl*grd_pw;*/

    W_Cl   = - KF_Cl*grd_cl - KF_Cl*z_cl*grd_si - KD_Cl*grd_pw;

#undef lna_cl
#undef lna_na
}


int ex26(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define lna_cl(n)    (lna[n][0])
#define lna_na(n)    (lna[n][1])

  double RT,c_cl,c_na,lna[2][2],xs,ys,s_s;
  double s_w,s_g,k_rw,k_re,p_c,p_v,p_l[2];
  double fa,dx,ton,top,phi[2],phif,g,B_h;
  double zero = 0.,un = 1.,deux = 2. ;
  int i,ipc_sw;
  Function_t *fn = Material_GetFunction(Element_GetMaterial(&el)) ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phio     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d    = el.mat->pr[pm("r_d")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  k_int	  = el.mat->pr[pm("k_int")] ;
  d_va	  = el.mat->pr[pm("d_va")] ;

  n	  = el.mat->pr[pm("n")];
  anpha	  = el.mat->pr[pm("anpha")];
  S_max	  = el.mat->pr[pm("S_max")];
  Pc_min  = el.mat->pr[pm("Pc_min")];

  m       = un-un/n;
  A       = rho_w*gravite/anpha;

  RT = R_g*T;

  s_so = el.mat->pr[pm("s_so")] ;
  s_srf    = el.mat->pr[pm("s_srf")] ;

  noc_phi = floor(el.mat->pr[pm("phi_x")]+0.5) ;  

  B_h	  = el.mat->pr[pm("B_h")];
  if (B_h==zero) {
    ipc_sw = 0;
  } else  {
    ipc_sw = floor(B_h+0.5) ;  
  }

  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    xs   = XS(i);    ys   = fabs(xs);
    c_cl = C_Cl(i);
    c_na = c_cl;

    c_e(i)   = (un - (c_cl*v_cl+ c_na*v_na))/v_h2o ;   
    lna_es(i) = activite(c_cl,c_na,c_e(i),T,lna[i]);	

    s_s  = 0.5*(xs+ys)*s_srf;	
    
    phi[i] = ((noc_phi<1) ? phio : fonction(x[i][0],fn[noc_phi-1])) + (s_so-s_s)*v_nacl;

    p_v = P_v(i);
    p_l[i] = RT*c_e(i)*(log(p_v/p_vs)-lna_es(i))+p_atm;
  }
 
 /* Coefficient de transfert */

    p_v     = (P_v(0)+P_v(1))/deux ;
    p_c     = p_atm-0.5*(p_l[0]+p_l[1]);

    if (ipc_sw<1) {
      if (p_c>=0) {
	s_w  = pow(pow(p_c/A,n)+S_max,-m);
      } else {
	s_w  =pow(S_max,-m)-p_c*(un-pow(S_max,-m))/fabs(Pc_min);
      }
    } else {
      s_w  = fonction(p_c,fn[ipc_sw-1]);
    }    

    s_g     = un -s_w;    if (s_g<zero) s_g = zero;

    phif    = (phi[0]+phi[1])/deux ;

    k_rw    = pow(s_w,0.5)*pow(un-pow(un-pow(s_w,un/m),m),deux);
    /*    k_re    = un/(un+pow(-10.*log(s_w),2.)); */
    k_re=k_rw;
    fa      = pow(s_g,10./3.)*pow(phif,4./3.);

    /*    k_int = k_int*pow(phif/phio,3.)*pow((un-phio)/(un-phif),2.); */

    KD_W  = 0.5*(c_e(0)+c_e(1))*M_v*k_int/mu_w*k_rw;  
    KF_VA = d_va*fa/RT/p_atm;  

    c_cl     = (C_Cl(0)+C_Cl(1))/deux ;
    c_na     = c_cl;

    KD_Cl = k_int/mu_w*k_rw*c_cl;
    KD_Na = k_int/mu_w*k_rw*c_na;

    g = (0.001+0.07*pow(phif,2.)+1.8*((phif<=0.18) ? 0.: 1.)*pow(phif-0.18,2.));
    g = phif/phio*g/(0.001+0.07*pow(phio,2.)+1.8*((phio<=0.18) ? 0.: 1.)*pow(phio-0.18,2.));

    d_cl = d_cl*g;

    ton = d_cl/do_cl*k_re;
    top = ton/r_d;

    KF_Cl  = ton*do_cl;
    KF_Na  = top*do_na;

    dx     = x[1][0] - x[0][0];
    WA_Cl   = - c_cl*KF_Cl*(lna_cl(1)-lna_cl(0))/dx;
    WA_Na   = - c_na*KF_Na*(lna_na(1)-lna_na(0))/dx;
  
    return(0) ;

#undef lna_cl
#undef lna_na

}

int ct26(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double s_w,s_g,c_cl,c_na,xs,ys,s_s;
  double grd_pw,grd_pvg, RT,grd_cl,grd_na,grd_si;
  double p_v,p_l[2],p_c,B_h;
  int    i,ipc_sw;
  double dx;
  double S_Dcz2, phi;
  double zero = 0.,un = 1. ;
  Function_t *fn = Material_GetFunction(Element_GetMaterial(&el)) ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */

  phio     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ; 
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  k_int	  = el.mat->pr[pm("k_int")] ;
  d_va	  = el.mat->pr[pm("d_va")] ;

  n	  = el.mat->pr[pm("n")];
  anpha	  = el.mat->pr[pm("anpha")];
  S_max	  = el.mat->pr[pm("S_max")];
  Pc_min  = el.mat->pr[pm("Pc_min")];
  m       = un-un/n;
  A       = rho_w*gravite/anpha;

  RT =R_g*T;

  s_so = el.mat->pr[pm("s_so")] ;
  s_srf    = el.mat->pr[pm("s_srf")] ;
  noc_phi = floor(el.mat->pr[pm("phi_x")]+0.5) ;  

  B_h	  = el.mat->pr[pm("B_h")];
  if (B_h==zero) {
    ipc_sw = 0;
  } else  {
    ipc_sw = floor(B_h+0.5) ;  
  }
   
  /* Contenus molaires */

 for(i=0;i<el.nn;i++) {
    xs   = XS(i);    ys   = fabs(xs);
    C_Cl(i) = k_s*(0.5*(xs-ys)+1.);
    c_cl = C_Cl(i);
    c_na = c_cl;
    s_s  = 0.5*(xs+ys)*s_srf;	 		

    phi = ((noc_phi<1) ? phio : fonction(x[i][0],fn[noc_phi-1])) + (s_so-s_s)*v_nacl;

    p_v = P_v(i);
    p_l[i] = RT*c_e(i)*(log(p_v/p_vs)-lna_es(i))+p_atm;
    p_c = p_atm-p_l[i];

    if (ipc_sw<1) {
      if (p_c>=0) {
	s_w  = pow(pow(p_c/A,n)+S_max,-m);
      } else {
	s_w  =pow(S_max,-m)-p_c*(un-pow(S_max,-m))/fabs(Pc_min);
      }
    } else {
      s_w  = fonction(p_c,fn[ipc_sw-1]);
    }    

    s_g  = un -s_w;

    N_V(i)   = phi*M_v/RT*p_v*s_g;  
    N_W(i)   = phi*s_w*c_e(i)*M_v;  

    N_Cl(i)  = phi*s_w*c_cl+s_s;

    if(s_w<0.||s_w>1.||p_v<=0||phi<=0.) {
      printf("\n\
x       = %e\n\
p_v     = %e\n\
s_w     = %e\n\
p_c     = %e\n\
x_s     = %e\n\
phi     = %e\n",x[i][0],p_v,s_w,p_c,xs,phi) ;
      return(-1) ;
    }
 }

    /* Flux */

    dx     = x[1][0] - x[0][0];

    grd_cl = (C_Cl(1) - C_Cl(0))/dx ;
    grd_na = grd_cl; 
    S_Dcz2 = KF_Cl*pow(z_cl,2.)+KF_Na*pow(z_na,2.);  
    grd_si =-(z_cl*(KF_Cl*grd_cl-WA_Cl)+z_na*(KF_Na*grd_na-WA_Na))/S_Dcz2;
	
    grd_pw = (p_l[1]-p_l[0])/dx;
    grd_pvg= (P_v(1) - P_v(0))/dx;

    W_W    = - KD_W*grd_pw;
    W_V    = - KF_VA*M_v*grd_pvg;

/*    W_Cl   = - KF_Cl*grd_cl - KD_Cl*grd_pw;*/

    W_Cl   = - KF_Cl*grd_cl - KF_Cl*z_cl*grd_si - KD_Cl*grd_pw;

    return(0) ;
}

int mx26(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])

  double tr,RT ;
  double s_w,p_c,p_v,s_g,p_l[2],xs,ys,c_cl,c_na,s_s;
  double dx,xm, dsw_dpc,phi,dccl_xs[2],dys_xs;
  double volume[2],surf,S_Dcz2,B_h;
  int    i,ipc_sw;
  double zero = 0.,un = 1.,deux = 2. ;
  Function_t *fn = Material_GetFunction(Element_GetMaterial(&el)) ;
  
  /* initialisation */
  for(i=0;i<el.nn*el.nn*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  
  phio    = el.mat->pr[pm("porosite")] ; 
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ; 
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  k_int	  = el.mat->pr[pm("k_int")] ;
  d_va	  = el.mat->pr[pm("d_va")] ;

  n	  = el.mat->pr[pm("n")];
  anpha	  = el.mat->pr[pm("anpha")];
  S_max	  = el.mat->pr[pm("S_max")];
  Pc_min  = el.mat->pr[pm("Pc_min")];

  m       = un-un/n;
  A       = rho_w*gravite/anpha;

  RT = R_g*T;
  s_so = el.mat->pr[pm("s_so")] ;
  s_srf    = el.mat->pr[pm("s_srf")] ;
  noc_phi = floor(el.mat->pr[pm("phi_x")]+0.5) ;  

  B_h	  = el.mat->pr[pm("B_h")];
  if (B_h==zero) {
    ipc_sw = 0;
  } else  {
    ipc_sw = floor(B_h+0.5) ;  
  }

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
    xs   = XS(i);    ys   = fabs(xs);
    c_cl = C_Cl(i);
    c_na = c_cl;
    s_s   = 0.5*(xs+ys)*s_srf;	 		
    if (xs>zero) dys_xs = un; else dys_xs = - un; 
    
    phi = ((noc_phi<1) ? phio : fonction(x[i][0],fn[noc_phi-1])) + (s_so-s_s)*v_nacl;

    p_v = P_v(i);
    p_l[i] = RT*c_e(i)*(log(p_v/p_vs)-lna_es(i))+p_atm;
    p_c = p_atm-p_l[i];
 
    dccl_xs[i]= k_s*0.5*(un-dys_xs);

    if (ipc_sw<1) {
      if (p_c>=0) {
        s_w  = pow(pow(p_c/A,n)+S_max,-m);
      } else {
        s_w  =pow(S_max,-m)-p_c*(un-pow(S_max,-m))/fabs(Pc_min);
      }
    } else {
      s_w  = fonction(p_c,fn[ipc_sw-1]);
    }    

    s_g  = un -s_w;
    if (ipc_sw<1){
    if (p_c<zero) {
      dsw_dpc = -(un-pow(S_max,-m))/fabs(Pc_min);
    } else {
      dsw_dpc = -m*pow(pow(p_c/A,n)+un,-m-un)*n*pow(p_c/A,n-un)/A;
    }
    }else {
      dsw_dpc = (fonction(p_c+1.,fn[ipc_sw-1])-s_w)/1.;
    }

    /*
      Conservation de Cl (chlore) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
    */
    K(i*NEQ+E_Cl,i*NEQ+I_XS) += volume[i]*(phi*s_w*dccl_xs[i]+0.5*(un+dys_xs)*s_srf) ;
    K(i*NEQ+E_Cl,i*NEQ+I_Pv) += volume[i]*(phi*c_cl)*(-dsw_dpc*c_e(i)*RT/p_v) ;
    K(i*NEQ+E_Cl,i*NEQ+I_XS) += -volume[i]*(s_w*c_cl*v_nacl*0.5*(un+dys_xs)*s_srf) ;

    /*
      Conservation de Vapeur (V) : (n_V1 - n_Vn) + dt * div(w_V) = 0
    */

    K(i*NEQ+E_VW,i*NEQ+I_Pv) += volume[i]*(phi*M_v/RT)*(s_g+p_v*dsw_dpc*c_e(i)*RT/p_v) ;
    K(i*NEQ+E_VW,i*NEQ+I_XS) += volume[i]*(M_v/RT*p_v*s_g)*(-v_nacl*0.5*(un+dys_xs)*s_srf) ;

   /*
      Conservation de Humidite (W) : (n_W1 - n_Wn) + dt * div(w_W) = 0
    */

      K(i*NEQ+E_VW,i*NEQ+I_Pv) += volume[i]*(phi*M_v*c_e(i))*(-dsw_dpc*c_e(i)*RT/p_v)  ;
      K(i*NEQ+E_VW,i*NEQ+I_XS) += volume[i]*(s_w*M_v*c_e(i))*(-v_nacl*0.5*(un+dys_xs)*s_srf)  ;
  }

  /* termes d'ecoulement */
  tr     = dt*surf/dx ;
  S_Dcz2  = KF_Cl*pow(z_cl,deux)+KF_Na*pow(z_na,deux);   

 /*
    Conservation de Cl (chlore) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
/*
    K(E_Cl,I_XS)             += + tr*KF_Cl*dccl_xs[0];
    K(E_Cl,I_XS+NEQ)         += - tr*KF_Cl*dccl_xs[1];
    K(E_Cl+NEQ,I_XS)         += - tr*KF_Cl*dccl_xs[0];
    K(E_Cl+NEQ,I_XS+NEQ)     += + tr*KF_Cl*dccl_xs[1];
  */
    
    K(E_Cl,I_XS)             += + tr*KF_Cl*dccl_xs[0] - tr*KF_Cl*pow(z_cl,2.)*KF_Cl/S_Dcz2*dccl_xs[0];
    K(E_Cl,I_XS+NEQ)         += - tr*KF_Cl*dccl_xs[1] + tr*KF_Cl*pow(z_cl,2.)*KF_Cl/S_Dcz2*dccl_xs[1];
    K(E_Cl+NEQ,I_XS)         += - tr*KF_Cl*dccl_xs[0] + tr*KF_Cl*pow(z_cl,2.)*KF_Cl/S_Dcz2*dccl_xs[0];
    K(E_Cl+NEQ,I_XS+NEQ)     += + tr*KF_Cl*dccl_xs[1] - tr*KF_Cl*pow(z_cl,2.)*KF_Cl/S_Dcz2*dccl_xs[1];

    K(E_Cl,I_XS)             += - tr*KF_Cl*z_cl*z_na*KF_Na/S_Dcz2*dccl_xs[0];
    K(E_Cl,I_XS+NEQ)         += + tr*KF_Cl*z_cl*z_na*KF_Na/S_Dcz2*dccl_xs[1];
    K(E_Cl+NEQ,I_XS)         += + tr*KF_Cl*z_cl*z_na*KF_Na/S_Dcz2*dccl_xs[0];
    K(E_Cl+NEQ,I_XS+NEQ)     += - tr*KF_Cl*z_cl*z_na*KF_Na/S_Dcz2*dccl_xs[1];


    K(E_Cl,I_Pv)              += + tr*KD_Cl*c_e(0)*RT/P_v(0);
    K(E_Cl,I_Pv+NEQ)          += - tr*KD_Cl*c_e(1)*RT/P_v(1);
    K(E_Cl+NEQ,I_Pv)          += - tr*KD_Cl*c_e(0)*RT/P_v(0);
    K(E_Cl+NEQ,I_Pv+NEQ)      += + tr*KD_Cl*c_e(1)*RT/P_v(1);

   /*
    Conservation de Vapeur (V) : (n_V1 - n_Vn) + dt * div(w_V) = 0
  */


  K(E_VW,I_Pv)             += +tr*M_v*KF_VA;
  K(E_VW,I_Pv+NEQ)         += -tr*M_v*KF_VA;
  K(E_VW+NEQ,I_Pv)         += -tr*M_v*KF_VA;
  K(E_VW+NEQ,I_Pv+NEQ)     += +tr*M_v*KF_VA; 


  /*
    Conservation de Water (W) : (n_W1 - n_Wn) + dt * div(w_W) = 0
  */

  K(E_VW,I_Pv)             += +tr*KD_W*c_e(0)*RT/P_v(0);
  K(E_VW,I_Pv+NEQ)         += -tr*KD_W*c_e(1)*RT/P_v(1);
  K(E_VW+NEQ,I_Pv)         += -tr*KD_W*c_e(0)*RT/P_v(0);
  K(E_VW+NEQ,I_Pv+NEQ)     += +tr*KD_W*c_e(1)*RT/P_v(1); 
 
  return(0) ;
#undef K
}

void rs26(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
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
/*
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*(W_Cl);
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*(W_Cl);
*/

  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*(W_Cl + WA_Cl) ;
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*(W_Cl + WA_Cl) ;

 /*
      Conservation de VW (V) : [(n_V1 - n_Vn)] + dt * div(w_V) = 0
    */
  R(0,E_VW) -= volume[0]*(N_V(0)+N_W(0) - N_Vn(0)-N_Wn(0)) + dt*surf*(W_V+W_W) ;
  R(1,E_VW) -= volume[1]*(N_V(1)+N_W(1) - N_Vn(1)-N_Wn(1)) - dt*surf*(W_V+W_W);
#undef R
}

int so26(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso,ipc_sw ;
  double RT,phi,c_e,lna_es;
  double s_w,s_g, p_v,p_l,p_c,xs,ys,c_cl,c_na,s_s;
  double zero = 0.,un = 1.,deux = 2. ;
  double lna[2],B_h;
  Function_t *fn = Material_GetFunction(Element_GetMaterial(&el)) ;

  /*
    Donnees
  */

  phio    = el.mat->pr[pm("porosite")] ; 
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  r_d	  = el.mat->pr[pm("r_d")] ; 
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  k_int	  = el.mat->pr[pm("k_int")] ;
  d_va	  = el.mat->pr[pm("d_va")] ;

  n	  = el.mat->pr[pm("n")];
  anpha	  = el.mat->pr[pm("anpha")];
  S_max	  = el.mat->pr[pm("S_max")];
  Pc_min  = el.mat->pr[pm("Pc_min")];

  m       = un-un/n;
  A       = rho_w*gravite/anpha;

  RT = R_g*T;
  s_so = el.mat->pr[pm("s_so")] ;
  s_srf    = el.mat->pr[pm("s_srf")] ;
  noc_phi = floor(el.mat->pr[pm("phi_x")]+0.5) ;  

  B_h	  = el.mat->pr[pm("B_h")];
  if (B_h==zero) {
    ipc_sw = 0;
  } else  {
    ipc_sw = floor(B_h+0.5) ;  
  }
 
  /* initialisation */
  nso = 11;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* concentration */
  if(s[0] < (x[0][0] + x[1][0])/deux) {
    p_v  = P_v(0);
    c_cl = C_Cl(0);
    xs   = XS(0);	
    phi = ((noc_phi<1) ? phio : fonction(x[0][0],fn[noc_phi-1])) ;
  } else {
    p_v  = P_v(1);
    c_cl = C_Cl(1);
    xs   = XS(1);
    phi = ((noc_phi<1) ? phio : fonction(x[1][0],fn[noc_phi-1])) ;
  }
  c_na =c_cl;
  c_e   = (un - (c_cl*v_cl+ c_na*v_na))/v_h2o ; 
  lna_es = activite(c_cl,c_na,c_e,T,lna);  
  
  p_l = RT*c_e*(log(p_v/p_vs)-lna_es)+p_atm;
  p_c = p_atm-p_l;

    if (ipc_sw<1) {
      if (p_c>=0) {
	s_w  = pow(pow(p_c/A,n)+S_max,-m);
      } else {
	s_w  =pow(S_max,-m)-p_c*(un-pow(S_max,-m))/fabs(Pc_min);
      }
    } else {
      s_w  = fonction(p_c,fn[ipc_sw-1]);
    }    

  s_g = un - s_w;
  ys   = fabs(xs);
  s_s   = 0.5*(xs+ys)*s_srf;	 		

  phi = phi + (s_so-s_s)*v_nacl;

  /* quantites exploitees */
  strcpy(r[0].text,"h_r") ; r[0].n = 1 ;
  r[0].v[0] = p_v/p_vs;
  strcpy(r[1].text,"tener en eau") ; r[1].n = 1 ;
  r[1].v[0] = s_w*phi;
  strcpy(r[2].text,"Saturation") ; r[2].n = 1 ;
  r[2].v[0] = s_w;
  strcpy(r[3].text,"p_w") ; r[3].n = 1 ;
  r[3].v[0] = p_l;
  strcpy(r[4].text,"flux eau") ; r[4].n = 1 ;
  r[4].v[0] = W_W;
  strcpy(r[5].text,"flux vapeur") ; r[5].n = 1 ;
  r[5].v[0] = W_V;
  strcpy(r[6].text,"Cl libre") ; r[6].n = 1 ;
  r[6].v[0] = c_cl;
  strcpy(r[7].text,"NaCl") ; r[7].n = 1 ;
  r[7].v[0] = s_s;
  strcpy(r[8].text,"Cl total") ; r[8].n = 1 ;
  r[8].v[0] = phi*s_w*c_cl+s_s;
  strcpy(r[9].text,"phi") ; r[9].n = 1 ;
  r[9].v[0] = phi;
  strcpy(r[10].text,"Flux Cl") ; r[10].n = 1 ;
  r[10].v[0] = W_Cl+WA_Cl;
  return(nso) ;
}


double activite(double c_cl,double c_na,double c_w,double T,double *lna)
/* L'activite chimique de l'eau d'une solution de NaCl */
{
#define lna_cl    (lna[0])
#define lna_na    (lna[1])

/* valences */
#define z_cl   (-1.)
#define z_na   (+1.)

/* masse molaire */
#define M_h2o   (18.e-3)

  double m_cl,m_na,m_T ;
  double I,A,lna_w,epsi ;

  double T_0 = 273.15 ;
  double b0 = sqrt(M_h2o),S0 = pow(M_h2o,1.29) ; /* references */
  double b_na = 4.352/b0,b_cl = 1.827/b0 ; /* donnees intrinseques */
  double S_na = 26.448/S0,S_cl = 19.245/S0 ;

  double zero = 0. ;

  /* molarites */
  if(c_cl < zero) c_cl = zero ;
  if(c_na < zero) c_na = zero ;
  
  epsi = 0.0007*(T - T_0)*(T - T_0) - 0.3918*(T - T_0) + 87.663 ;
  A = 1398779.816/pow(epsi*T,1.5)/b0 ;
  
  /* molalites*M_h2o (en moles/mole) */
  m_cl = c_cl/c_w ;
  m_na = c_na/c_w ;

  /* la force ionique */
  I = 0.5*(z_cl*z_cl*m_cl + z_na*z_na*m_na) ;
  
  if (I > zero) {
    m_T =  m_cl + m_na ;

    lna_w = m_cl*lna_i(T,I,z_cl,b_cl,S_cl,A)
          + m_na*lna_i(T,I,z_na,b_na,S_na,A) ;

    /* selon Lin & Lee */
    lna_cl = lng_LinLee(T,I,z_cl,b_cl,S_cl,A) ;
    lna_na = lng_LinLee(T,I,z_na,b_na,S_na,A) ;

    /* selon TQN */
    /*
    lna_cl = lng_TQN(T,I,z_cl,b_cl,S_cl,A,lna_w,m_T) ;
    lna_na = lng_TQN(T,I,z_na,b_na,S_na,A,lna_w,m_T) ;
    */

  }  else {
    lna_cl = 0.;
    lna_na = 0.;
    lna_w  = 0.;
  }

  return(lna_w) ;

#undef lna_cl
#undef lna_na

#undef z_cl
#undef z_na

#undef M_h2o
}



double lng_LinLee(double T,double I,double z,double b,double S,double A)
/* Le log du coefficient d'activite d'un ion d'apres Lin & Lee */ 
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*(II/(1 + b*II) + 2*log(1 + b*II)/b) + S*pow(I,alpha)/T ;
  
  return(lng*z*z) ;
}

double lng_TQN(double T,double I,double z,double b,double S,double A,double lna_w,double m_t)
/* Le log du coefficient d'activite d'un ion (T.Q Nguyen) :
   lng_i = dGamma/dm_i = (dGamma/dm_i)_I - 0.5*z_i*z_i*(lna_w + m_t)/I 
   lna_w = - m_t - sum_i ( m_i*lng_i ) + Gamma */
{
  double alpha = 1.29,II = sqrt(I) ;
  double lng ;
  
  lng = - A*2*log(1 + b*II)/b + S*pow(I,alpha)/(1+alpha)/T - 0.5*(lna_w + m_t)/I ;
  
  return(lng*z*z) ;
}

double lna_i(double T,double I,double z,double b,double S,double A)
/* Contribution de chaque ion au log de l'activite du solvant 
   lna_w = sum_i ( m_i*lna_i ) (T.Q Nguyen) */ 
{
  double alpha = 1.29,a1 = alpha/(1+alpha),II = sqrt(I) ;
  double lna ;
  
  lna = A*II/(1 + b*II) - a1*S*pow(I,alpha)/T ;
  
  return(-1 + lna*z*z) ;
}
