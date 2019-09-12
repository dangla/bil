/*
	Modele (ancien 150): quatre especes ionique (Cl, Na, K, OH, H)
	Diffusion: Loi de Frick + Migration + reaction chimique (H2O = H+ + OH-)
	Isotherme de fixation: Freundlich
	D est calculee D_oh=D_cl/Do_cl*Do_oh
	Tenir en compte: Activitee de la solution 
	Non stature
	Activite
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "CommonModel.h"

#define MODELINDEX  44
#define TITLE "Chlorures dans les betons non satures"
#define AUTHORS "Nguyen"

#include "OldMethods.h"

/* Macros */


#define NEQ     (7)
#define NVI     (24)
#define NVE     (24)

#define E_Cl    (0)
#define E_OH    (1)
#define E_Na    (2)
#define E_K     (3)
#define E_PS    (4)
#define E_WV    (5)
#define E_A     (6)

#define I_Cl    (0)
#define I_OH    (1)
#define I_Na    (2)
#define I_K     (3)
#define I_SI    (4)
#define I_Pv    (5)
#define I_Pg    (6)

#define C_Cl(n)     (u[(n)][I_Cl])
#define C_OH(n)     (u[(n)][I_OH])
#define C_Na(n)     (u[(n)][I_Na])
#define C_K(n)      (u[(n)][I_K])
#define SI(n)       (u[(n)][I_SI])
#define P_v(n)      (u[(n)][I_Pv])
#define P_g(n)      (u[(n)][I_Pg])

#define N_Cl(n)     (f[(n)])
#define N_OH(n)     (f[(n+2)])
#define N_Na(n)     (f[(n+4)])
#define N_K(n)      (f[(n+6)])
#define N_PS(n)     (f[(n+8)])
#define N_V(n)      (f[(n+10)])
#define N_A(n)      (f[(n+12)])
#define N_W(n)      (f[(n+14)])
#define W_Cl        (f[16])
#define W_OH        (f[17])
#define W_Na        (f[18])
#define W_K         (f[19])
#define W_PS        (f[20])
#define W_V         (f[21])
#define W_A         (f[22])
#define W_W         (f[23])

#define N_Cln(n)     (f_n[(n)])
#define N_OHn(n)     (f_n[(n+2)])
#define N_Nan(n)     (f_n[(n+4)])
#define N_Kn(n)      (f_n[(n+6)])
#define N_PSn(n)     (f_n[(n+8)])
#define N_Vn(n)      (f_n[(n+10)])
#define N_An(n)      (f_n[(n+12)])
#define N_Wn(n)      (f_n[(n+14)])

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

#define KD_Cl       (va[(10)])
#define KD_OH       (va[(11)])
#define KD_Na       (va[(12)])
#define KD_K        (va[(13)])
#define KD_H        (va[(14)])

#define KD_W        (va[(15)])
#define KD_V        (va[(16)])
#define KD_A        (va[(17)])
#define KF_VA       (va[(18)])

#define KF_PS       (va[(19)])

#define JA_Cl       (va[(20)])
#define JA_OH       (va[(21)])
#define JA_Na       (va[(22)])
#define JA_K        (va[(23)])

/* valences */
#define z_cl    (-1.)
#define z_oh    (-1.)
#define z_na    (+1.)
#define z_k     (+1.)
#define z_h     (+1.)

/* coefficients de diffusion dans l'eau dilue */
#define do_cl   (2.032e-9)
#define do_oh   (5.273e-9)
#define do_na   (1.334e-9)
#define do_k    (1.957e-9)
#define do_h    (9.310e-9)

/* constante d'equilibre */
#define k_e     (0.) /* (1.e-8) */

/* constantes physiques */
#define R_g	(8.3143)    /* Gaz parfait */
#define F   	(9.64846e4) /* Faraday */
#define epsi_o 	(8.854e-12) /* permeabilite electrique de l'eau */

#define p_atm   (1.01325e5) /* Pression atmospherique */
#define mu_g    (1.8e-5)    /* Viscosite de l'air (Pa.s) */
#define mu_w    (1.002e-3)  /* Viscosite de l'eau (Pa.s) */

/* Volumes molaires (m3/mol) */
#define v_h2o   (18.00e-6)
#define v_na    (22.47e-6)
#define v_k     (54.11e-6)
#define v_h     (-5.50e-6)
#define v_oh    (23.50e-6)
#define v_cl    (-0.35e-6)

/* Masses molaires (kg/mol) */
#define M_Cl    (35.45e-3)
#define M_OH    (17.00e-3)
#define M_Na    (22.95e-3)
#define M_K     (39.10e-3)
#define M_v     (1.8e-2)
#define M_a     (2.896e-2)

/* Densite initiale de radicaux hydroxyles */
#define s_oho   (1000.)

/* Fonctions */
static int    pm(const char *s) ;
static double activite(double,double,double,double,double,int,double*) ;
static double lng_LinLee(double,double,double,double,double,double) ;
static double lna_i(double,double,double,double,double,double) ;
static double lng_TQN(double,double,double,double,double,double,double,double) ;

/* Parametres */
static int    noc_Iso,mode_a ;
static double phi,d_cl,T, kw_int,kg_int,p_vs,d_va,klinkenberg;
static double m,n,a ;

int pm(const char *s)
{
  if(strcmp(s,"porosite") == 0) return (0) ;
  else if(strcmp(s,"D_Cl") == 0) return (1);
  else if(strcmp(s,"T") == 0) return (2);
  else if(strcmp(s,"activite") == 0) return (3);
  else if(strcmp(s,"kw_int") == 0) return (4); 
  else if(strcmp(s,"kg_int") == 0) return (5); 
  else if(strcmp(s,"klinkenberg") == 0) return (6);
  else if(strcmp(s,"d_va") == 0) return (7);
  else if(strcmp(s,"B_h") == 0) return (8);
  else if(strcmp(s,"a") == 0) return (9);
  else if(strcmp(s,"m") == 0) return (10);
  else if(strcmp(s,"n") == 0) return (11);
  else if(strcmp(s,"Isotherme") == 0) return (12);
  else if(strcmp(s,"courbes") == 0) return (13);
  else {
    printf("donnee \"%s\" non connue (pm44)\n",s) ; exit(0) ;
  }
}

int dm44(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int    nd = 13,nc = 1;

  if(dim > 1) arret("dm44 : dimension > 1 non prevue") ;

  mat->neq = NEQ ;
  strcpy(mat->eqn[E_Cl],"E_Cl") ;
  strcpy(mat->eqn[E_OH],"E_OH") ;
  strcpy(mat->eqn[E_Na],"E_Na") ;
  strcpy(mat->eqn[E_K],"E_K") ;
  strcpy(mat->eqn[E_PS],"E_PS") ;
  strcpy(mat->eqn[E_WV],"E_WV") ;
  strcpy(mat->eqn[E_A],"E_A") ;

  strcpy(mat->inc[I_Cl],"c_cl") ;
  strcpy(mat->inc[I_OH],"c_oh") ;
  strcpy(mat->inc[I_Na],"c_na") ;
  strcpy(mat->inc[I_K],"c_k") ;
  strcpy(mat->inc[I_SI],"si") ;
  strcpy(mat->inc[I_Pv],"p_v") ;
  strcpy(mat->inc[I_Pg],"p_g") ;

  lit_mate(mat,ficd,pm,nd,nc) ;

  /*
  mat->n = nd ;
  mat->nc = 0 ;
  for(i=0;i<nd+nc;i++) {
    if(!fgets(line,sizeof(line),ficd)) arret("dmat (1) : erreur ou fin de fichier") ;
    sscanf(line," %[^= ] =",mot) ;
    p = strchr(line,'=') + 1 ;
    if(!strncmp(mot,"courbes",6)) {
  *//*
      fscanf(ficd,"%s",nom) ; 
      fscanf(ficd,"%d",&n_pa) ;	
      for (j=0;j<n_pa;j++) fscanf(ficd,"%lf",para+j);
      if (n_pa!=0) ecrit_courbe(nom,n_pa,para);	
    *//*
      ecrit_courbe(p);	
      lit_courbe(mat,line) ;
    } else sscanf(p,"%lf",mat->pr+pm(mot)) ;
  }
      */


  return(mat->n) ;
}


int qm44(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;
  
  printf("\n\n\
Le systeme est forme de 7 equations :\n\
\t 1. Conservation de la masse de Cl (c_cl)\n\
\t 2. Conservation de la masse de OH (c_oh)\n\
\t 3. Conservation de la masse de Na (c_na)\n\
\t 4. Conservation de la masse de K  (c_k)\n\
\t 5. Conservation de la masse d\'eau (p_v)\n\
\t 6. Conservation de la masse d\'air (p_g)\n\
\t 7. Equilibre chimique             (si)\n") ;

  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"porosite = # porosite\n") ;
  fprintf(ficd,"D_Cl = # coefficient de diffusion effective de Cl\n") ;
  fprintf(ficd,"T = # Temperature T\n") ;
  fprintf(ficd,"kw_int = # \n") ;
  fprintf(ficd,"kg_int = # \n") ;
  fprintf(ficd,"klikenberg = # \n") ;
  fprintf(ficd,"d_va = # \n") ;
  fprintf(ficd,"B_h = # \n") ;
  fprintf(ficd,"a = # \n") ;
  fprintf(ficd,"m = # \n") ;
  fprintf(ficd,"n = # \n") ;
  fprintf(ficd,"activite = # model de calcul activite de la solution\n") ;
  fprintf(ficd,"Isotherme = # type de l'isotherme\n") ;
  fprintf(ficd,"courbes = my_file # Nom du fichier \n") ;

  return(NEQ) ;

}

void tb44(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI ; /* implicite */
  el->n_ve = NVE ; /* explicite */
}

void ch44(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
#define R(n,i)     (r[(n)*NEQ+(i)])

  double B_h;
  int    i,ieq ;

 /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0.; 

  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  kw_int  = el.mat->pr[pm("kw_int")] ;
  d_va	  = el.mat->pr[pm("d_va")];
  B_h	  = el.mat->pr[pm("B_h")];

  /* on calcule le numero de l'equation */
  if(isdigit(cg.eqn[0])) { /* donne sous forme numerique */
    ieq  = atoi(cg.eqn) - 1 ;
  } else {                 /* donne sous forme alphabetique */
    for(ieq=0;ieq<NEQ;ieq++) if(!strcmp(cg.eqn,el.mat->eqn[ieq])) break ;
    if(ieq == NEQ) arret("ch44 (1) : equation non connue") ;
  }

  if(ieq == E_WV) {
    R(0,E_WV) = - B_h*dt*(P_v(0)-champ(x[0],dim,*cg.ch)) ;
   } else {printf("\nchargement non prevu (ch44)\n") ; exit(0) ;}
 
#undef R
}

void in44(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
#define lna_cl(n)    (lna[n][0])
#define lna_oh(n)    (lna[n][1])
#define lna_na(n)    (lna[n][2])
#define lna_k(n)     (lna[n][3])

  double c_cl,c_oh,c_na,c_k,c_h,s_w,s_g;
  double grd_cl,grd_oh,grd_na,grd_k,grd_h,grd_si,grd_pw,grd_pag,grd_pvg,grd_pg,RT;
  double s_cl,s_oh ;
  double dx,ton,top,fa ;
  double lna[2][4],epsi,lna_es[2],c_e[2];
  double p_g,p_v,p_a,p_l[2],p_c, k_rw,k_re,k_rg;
  int    i;
  double J_Cl,J_OH,J_Na,J_K,J_W;
  double cv[2],ca[2],rho_v,rho_a,rho_g,rho_w;
  double zero = 0.,un = 1.,deux = 2. ;

  if(el.dim < dim) return ;
  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  kw_int  = el.mat->pr[pm("kw_int")] ;
  kg_int  = el.mat->pr[pm("kg_int")] ;
  d_va	  = el.mat->pr[pm("d_va")]; /* 2.17e-5*/
  klinkenberg  = el.mat->pr[pm("klinkenberg")];
  a  = el.mat->pr[pm("a")];
  n  = el.mat->pr[pm("n")];
  m  = el.mat->pr[pm("m")];

  noc_Iso = floor(el.mat->pr[pm("Isotherme")]+0.5) ;
  mode_a = floor(el.mat->pr[pm("activite")]+0.5) ;

  RT      = R_g*T;
  epsi = 0.0007*pow(T-273.15,2.)-0.3918*(T-273.15)+87.663;
 
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    c_cl = C_Cl(i);
    c_oh = C_OH(i);
    c_na = C_Na(i) ;
    c_h	 = k_e/c_oh;
    c_k  = C_K(i) ;

    c_e[i]    = (un - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_h*v_h + c_k*v_k))/v_h2o ;   
    lna_es[i] = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna[i]);	

    p_v = P_v(i);
    p_g = P_g(i);
    p_a = p_g-p_v;
    p_l[i] = RT*c_e[i]*(log(p_v/p_vs)-lna_es[i])+p_atm;
    p_c = p_g-p_l[i];
 
    s_w  = pow(pow(p_c/a,n)+1.,-m);
    s_g  = un -s_w;

    s_cl    = (noc_Iso<1) ? zero : courbe(c_cl,el.mat->cb[noc_Iso-1]); 
    s_oh     = s_oho - s_cl;

    N_Cl(i)  = phi*s_w*c_cl + s_cl;
    N_OH(i)  = phi*s_w*(c_oh-c_h) + s_oh;
    N_Na(i)  = phi*s_w*c_na;
    N_K(i)   = phi*s_w*c_k;
    N_PS(i)  = F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*s_w;  

    rho_v = M_v/RT*p_v;
    rho_a = M_a/RT*p_a;
    rho_g = rho_v + rho_a;
    rho_w = c_e[i]*M_v;
    cv[i] = rho_v/rho_g;
    ca[i] = un - cv[i];

    N_V(i)   = phi*s_g*rho_v;  
    N_A(i)   = phi*s_g*rho_a;  
    N_W(i)   = phi*s_w*rho_w;  

    if(c_cl <0. || c_oh <= 0. || c_na <0. || c_k < 0.|| c_h < 0.||s_w<0.||s_w>1.) {
      printf("\n\
x       = %e\n\
c_cl    = %e\n\
c_oh    = %e\n\
c_na    = %e\n\
c_k     = %e\n\
s_w     = %e\n",x[i][0],c_cl,c_oh,c_na,c_k,s_w) ;
    }

  }

  /* Coefficient de transfert */

    c_cl     = (C_Cl(0)+C_Cl(1))/deux ;
    c_oh     = (C_OH(0)+C_OH(1))/deux ;
    c_na     = (C_Na(0)+C_Na(1))/deux ;
    c_k      = (C_K(0)+C_K(1))/deux ;
    c_h      = k_e/c_oh ;
    
    p_g     = (P_g(0)+P_g(1))/deux ;
    p_v     = (P_v(0)+P_v(1))/deux ;
    p_a     = p_g - p_v;
    p_c     = p_g-0.5*(p_l[0]+p_l[1]);
    s_w     = pow(pow(p_c/a,n)+1.,-m);
    s_g     = un -s_w;    if (s_g<zero) s_g = zero;

    k_rw    = pow(s_w,0.5)*pow(1-pow(1-pow(s_w,1/m),m),2.);
    k_re    = pow(s_w,6.);  /*fittage*/
    k_rg    = pow(1-s_w,5.5)*pow(1-pow(s_w,1/m),2.*m);  /*van Genuchten*/
    fa      = pow(s_g,4.2)*pow(phi,2.74);

    KD_Cl = kw_int/mu_w*k_rw*c_cl;
    KD_OH = kw_int/mu_w*k_rw*c_oh;
    KD_Na = kw_int/mu_w*k_rw*c_na;
    KD_K  = kw_int/mu_w*k_rw*c_k;
    KD_H  = kw_int/mu_w*k_rw*c_h;

    rho_v = M_v/RT*p_v;
    rho_a = M_a/RT*p_a;
    rho_g = rho_v + rho_a;
    rho_w = 0.5*(c_e[0]+c_e[1])*M_v;

    KD_W  = rho_w*kw_int/mu_w*k_rw;  
    KD_V  = rho_v*kg_int*(1.+klinkenberg/p_g)/mu_g*k_rg;  
    KD_A  = rho_a*kg_int*(1.+klinkenberg/p_g)/mu_g*k_rg;  
    KF_VA = rho_g*fa*d_va*p_atm/p_g*pow(T/273.,1.88);  

    ton = d_cl/do_cl*k_re;
    top = ton;

    KF_Cl  = ton*do_cl;
    KF_OH  = ton*do_oh ;
    KF_Na  = top*do_na;
    KF_K   = top*do_k;
    KF_H   = top*do_h;   

    dx     = x[1][0] - x[0][0];
    JA_Cl   = - c_cl*KF_Cl*(lna_cl(1)-lna_cl(0))/dx;
    JA_OH   = - c_oh*KF_OH*(lna_oh(1)-lna_oh(0))/dx;
    JA_Na   = - c_na*KF_Na*(lna_na(1)-lna_na(0))/dx;
    JA_K    = - c_k*KF_K*(lna_k(1)-lna_k(0))/dx;
   
    KF_PS  = -ton*s_w*epsi*epsi_o;

    KM_Cl  = z_cl*KF_Cl*c_cl*F/RT;	
    KM_OH  = z_oh*KF_OH*c_oh*F/RT;
    KM_Na  = z_na*KF_Na*c_na*F/RT;
    KM_K   = z_k*KF_K*c_k*F/RT;	
    KM_H   = z_h*KF_H*c_h*F/RT;

 /* Flux */
    
    dx     = x[1][0] - x[0][0];
		
    grd_cl = (C_Cl(1) - C_Cl(0))/dx ;
    grd_oh = (C_OH(1) - C_OH(0))/dx ;
    grd_na = (C_Na(1) - C_Na(0))/dx ;
    grd_k  = (C_K(1) - C_K(0))/dx ;
    grd_h  =  k_e*(un/C_OH(1) - un/C_OH(0))/dx ;

    grd_si = (SI(1)-SI(0))/dx;

    grd_pg = (P_g(1) - P_g(0))/dx ;
    grd_pw = (p_l[1]-p_l[0])/dx;
    grd_pvg= (cv[1] - cv[0])/dx;
    grd_pag= -grd_pvg;
    
    J_Cl =  - KF_Cl*grd_cl - KM_Cl*grd_si;
    J_OH =  - KF_OH*grd_oh - KM_OH*grd_si;
    J_Na =  - KF_Na*grd_na - KM_Na*grd_si;
    J_K  =  - KF_K*grd_k - KM_K*grd_si;
  
    J_W = -(M_Cl*(J_Cl+JA_Cl)+M_OH*(J_OH+JA_OH)+M_Na*(J_Na+JA_Na)+M_K*(J_K+JA_K));

    W_Cl   = - KF_Cl*grd_cl - KM_Cl*grd_si - KD_Cl*grd_pw;
    W_OH   = - KF_OH*grd_oh + KF_H*grd_h   - (KM_OH-KM_H)*grd_si - (KD_OH-KD_H)*grd_pw;
    W_Na   = - KF_Na*grd_na - KM_Na*grd_si - KD_Na*grd_pw;
    W_K    = - KF_K*grd_k   - KM_K*grd_si - KD_K*grd_pw;
    W_W    = - KD_W*grd_pw + J_W;
    W_V    = - KD_V*grd_pg-KF_VA*grd_pvg;
    W_A    = - KD_A*grd_pg-KF_VA*grd_pag;

    W_PS   = - KF_PS*grd_si;

#undef lna_cl
#undef lna_oh
#undef lna_na
#undef lna_k
}


int ex44(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t) 
/* Termes explicites (va)  */
{
#define lna_cl(n)    (lna[n][0])
#define lna_oh(n)    (lna[n][1])
#define lna_na(n)    (lna[n][2])
#define lna_k(n)     (lna[n][3])

  double c_cl,c_oh,c_na,c_k,c_h;
  double lna[2][4],RT,epsi,c_e[2],lna_es[2];
  double s_w,s_g,k_rw,k_re,p_c,p_v,p_a,p_g,p_l[2],k_rg;
  double ton,top,fa,dx;
  double zero = 0.,un = 1.,deux = 2. ;
  int i;
  double rho_v,rho_a,rho_g,rho_w;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  kw_int  = el.mat->pr[pm("kw_int")] ;
  kg_int  = el.mat->pr[pm("kg_int")] ;
  d_va	  = el.mat->pr[pm("d_va")]; /* 2.17e-3*/
  klinkenberg  = el.mat->pr[pm("klinkenberg")];
  a  = el.mat->pr[pm("a")];
  n  = el.mat->pr[pm("n")];
  m  = el.mat->pr[pm("m")];

  noc_Iso = floor(el.mat->pr[pm("Isotherme")]+0.5) ;
  mode_a = floor(el.mat->pr[pm("activite")]+0.5) ;

  RT      = R_g*T;
  epsi = 0.0007*pow(T-273.15,2.)-0.3918*(T-273.15)+87.663;

 
  /* Contenus molaires */
  for(i=0;i<el.nn;i++) {
    c_cl = C_Cl(i);
    c_oh = C_OH(i);
    c_na = C_Na(i) ;
    c_h	 = k_e/c_oh;
    c_k  = C_K(i) ;

    c_e[i]   = (un - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_h*v_h + c_k*v_k))/v_h2o ;   
    lna_es[i] = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna[i]);	

    p_v = P_v(i);
    p_g = P_g(i);
    p_a = p_g-p_v;
    p_l[i] = RT*c_e[i]*(log(p_v/p_vs)-lna_es[i])+p_atm;

  }
 
 /* Coefficient de transfert */

    c_cl     = (C_Cl(0)+C_Cl(1))/deux ;
    c_oh     = (C_OH(0)+C_OH(1))/deux ;
    c_na     = (C_Na(0)+C_Na(1))/deux ;
    c_k      = (C_K(0)+C_K(1))/deux ;
    c_h      = k_e/c_oh ;	

    p_g     = (P_g(0)+P_g(1))/deux ;
    p_v     = (P_v(0)+P_v(1))/deux ;
    p_a     = p_g -p_v;
    p_c     = p_g-0.5*(p_l[0]+p_l[1]);
    s_w     = pow(pow(p_c/a,n)+1.,-m);
    s_g     = un -s_w;    if (s_g<zero) s_g = zero;

    k_rw    = pow(s_w,0.5)*pow(1-pow(1-pow(s_w,1/m),m),2.);
    k_re    = pow(s_w,6.);  /*fittage*/
    k_rg    = pow(1-s_w,5.5)*pow(1-pow(s_w,1/m),2.*m);  /*van Genuchten*/
    fa      = pow(s_g,4.2)*pow(phi,2.74);

    KD_Cl = kw_int/mu_w*k_rw*c_cl;
    KD_OH = kw_int/mu_w*k_rw*c_oh;
    KD_Na = kw_int/mu_w*k_rw*c_na;
    KD_K  = kw_int/mu_w*k_rw*c_k;
    KD_H  = kw_int/mu_w*k_rw*c_h;

    rho_v = M_v/RT*p_v;
    rho_a = M_a/RT*p_a;
    rho_g = rho_v + rho_a;
    rho_w = 0.5*(c_e[0]+c_e[1])*M_v;

    KD_W  = rho_w*kw_int/mu_w*k_rw;  
    KD_V  = rho_v*kg_int*(1.+klinkenberg/p_g)/mu_g*k_rg;  
    KD_A  = rho_a*kg_int*(1.+klinkenberg/p_g)/mu_g*k_rg;  
    KF_VA = rho_g*fa*d_va*p_atm/p_g*pow(T/273.,1.88);  

    ton = d_cl/do_cl*k_re;
    top = ton;

    KF_Cl  = ton*do_cl;
    KF_OH  = ton*do_oh ;
    KF_Na  = top*do_na;
    KF_K   = top*do_k;
    KF_H   = top*do_h;   

    dx     = x[1][0] - x[0][0];
    JA_Cl   = - c_cl*KF_Cl*(lna_cl(1)-lna_cl(0))/dx;
    JA_OH   = - c_oh*KF_OH*(lna_oh(1)-lna_oh(0))/dx;
    JA_Na   = - c_na*KF_Na*(lna_na(1)-lna_na(0))/dx;
    JA_K    = - c_k*KF_K*(lna_k(1)-lna_k(0))/dx;
    
    KF_PS  = -ton*s_w*epsi*epsi_o;

    KM_Cl  = z_cl*KF_Cl*c_cl*F/RT;	
    KM_OH  = z_oh*KF_OH*c_oh*F/RT;
    KM_Na  = z_na*KF_Na*c_na*F/RT;
    KM_K   = z_k*KF_K*c_k*F/RT;	
    KM_H   = z_h*KF_H*c_h*F/RT;
    
    return(0) ;

#undef lna_cl
#undef lna_oh
#undef lna_na
#undef lna_k

}

int ct44(double **x,double **u,double **u_n,double *f,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double c_cl,c_oh,c_na,c_k,c_h,s_w,s_g;
  double grd_cl,grd_oh,grd_na,grd_k,grd_h,grd_si,grd_pw,grd_pag,grd_pvg,grd_pg,RT;
  double s_cl,s_oh ;
  double dx,lna[2][4];
  double epsi;
  double p_g,p_v,p_a,p_l[2],p_c,c_e[2],lna_es[2];
  int    i;
  double J_Cl,J_OH,J_Na,J_K,J_W;
  double cv[2],ca[2],rho_v,rho_a,rho_g,rho_w;
  double zero = 0.,un = 1. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  kw_int  = el.mat->pr[pm("kw_int")] ;
  kg_int  = el.mat->pr[pm("kg_int")] ;
  d_va	  = el.mat->pr[pm("d_va")]; /* 2.17e-3*/
  klinkenberg  = el.mat->pr[pm("klinkenberg")];
  a  = el.mat->pr[pm("a")];
  n  = el.mat->pr[pm("n")];
  m  = el.mat->pr[pm("m")];

  noc_Iso = floor(el.mat->pr[pm("Isotherme")]+0.5) ;
  mode_a = floor(el.mat->pr[pm("activite")]+0.5) ;

  RT      = R_g*T;
  epsi = 0.0007*pow(T-273.15,2.)-0.3918*(T-273.15)+87.663;

   
  /* Contenus molaires */

 for(i=0;i<el.nn;i++) {
    c_cl = C_Cl(i);
    c_oh = C_OH(i);
    c_na = C_Na(i) ;
    c_h	 = k_e/c_oh;
    c_k  = C_K(i) ;

    p_v = P_v(i);
    p_g = P_g(i);
    p_a = p_g-p_v;

    c_e[i]   = (un - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_h*v_h + c_k*v_k))/v_h2o ;   
    lna_es[i] = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna[i]);	

    p_l[i] = RT*c_e[i]*(log(p_v/p_vs)-lna_es[i])+p_atm;
    p_c = p_g-p_l[i];

    s_w  = pow(pow(p_c/a,n)+1.,-m);
    s_g  = un -s_w;

    s_cl    = (noc_Iso<1) ? zero : courbe(c_cl,el.mat->cb[noc_Iso-1]);
    s_oh     = s_oho - s_cl;

    N_Cl(i)  = phi*s_w*c_cl + s_cl ;
    N_OH(i)  = phi*s_w*(c_oh-c_h) + s_oh;
    N_Na(i)  = phi*s_w*c_na;
    N_K(i)   = phi*s_w*c_k;
    N_PS(i)  = F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*s_w;  

    rho_v = M_v/RT*p_v;
    rho_a = M_a/RT*p_a;
    rho_g = rho_v + rho_a;
    rho_w = c_e[i]*M_v;
    cv[i] = rho_v/rho_g;
    ca[i] = un - cv[i];

    N_V(i)   = phi*s_g*rho_v;  
    N_A(i)   = phi*s_g*rho_a;  
    N_W(i)   = phi*s_w*rho_w;  

    if(s_w<0.||s_w>1.||p_v<=0||p_c<0||p_g<0) {
      printf("\n x       = %e\n c_cl    = %e\n c_oh    = %e\n c_na    = %e\n c_k     = %e\n p_v     = %e\n c_e     = %e\n lna_es  = %e\n s_w     = %e\n p_c  = %e\n",x[i][0],c_cl,c_oh,c_na,c_k,p_v,c_e[i],lna_es[i],s_w,p_c) ;
      return(-1);
    }

 }

    /* Flux */

    dx     = x[1][0] - x[0][0];
		
    grd_cl = (C_Cl(1) - C_Cl(0))/dx ;
    grd_oh = (C_OH(1) - C_OH(0))/dx ;
    grd_na = (C_Na(1) - C_Na(0))/dx ;
    grd_k  = (C_K(1) - C_K(0))/dx ;
    grd_h  =  k_e*(un/C_OH(1) - un/C_OH(0))/dx ;

    grd_si = (SI(1)-SI(0))/dx;

    grd_pg = (P_g(1) - P_g(0))/dx ;
    grd_pw = (p_l[1]-p_l[0])/dx;
    grd_pvg= (cv[1] - cv[0])/dx;
    grd_pag= -grd_pvg;

    J_Cl =  - KF_Cl*grd_cl - KM_Cl*grd_si;
    J_OH =  - KF_OH*grd_oh - KM_OH*grd_si;
    J_Na =  - KF_Na*grd_na - KM_Na*grd_si;
    J_K  =  - KF_K*grd_k - KM_K*grd_si;
  
    J_W = -(M_Cl*(J_Cl+JA_Cl)+M_OH*(J_OH+JA_OH)+M_Na*(J_Na+JA_Na)+M_K*(J_K+JA_K));

    W_Cl   = - KF_Cl*grd_cl - KM_Cl*grd_si - KD_Cl*grd_pw;
    W_OH   = - KF_OH*grd_oh + KF_H*grd_h   - (KM_OH-KM_H)*grd_si - (KD_OH-KD_H)*grd_pw;
    W_Na   = - KF_Na*grd_na - KM_Na*grd_si - KD_Na*grd_pw;
    W_K    = - KF_K*grd_k   - KM_K*grd_si - KD_K*grd_pw;
    W_W    = - KD_W*grd_pw + J_W;
    W_V    = - KD_V*grd_pg-KF_VA*grd_pvg;
    W_A    = - KD_A*grd_pg-KF_VA*grd_pag;
    W_PS   = - KF_PS*grd_si;

    return(0) ;
}

int mx44(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define K(i,j)    (k[(i)*el.nn*NEQ+(j)])

  double c_cl,c_oh,c_na,c_k,c_h,s_cl,dfsdc,dc,tr,RT,epsi;
  double s_w,p_c,p_v,s_g, p_a,p_g, p_l[2],c_e[2],lna_es[2];
  double dlna_ccl[2],dlna_coh[2],dlna_cna[2],dlna_ck[2];
  double dx,xm,dsw_dpc,lna[2][4];
  double volume[2],surf;
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
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  kw_int  = el.mat->pr[pm("kw_int")] ;
  kg_int  = el.mat->pr[pm("kg_int")] ;
  d_va	  = el.mat->pr[pm("d_va")]; /* 2.17e-3*/
  klinkenberg  = el.mat->pr[pm("klinkenberg")];
  a  = el.mat->pr[pm("a")];
  n  = el.mat->pr[pm("n")];
  m  = el.mat->pr[pm("m")];

  noc_Iso = floor(el.mat->pr[pm("Isotherme")]+0.5) ;
  mode_a = floor(el.mat->pr[pm("activite")]+0.5) ;

  RT      = R_g*T;
  epsi = 0.0007*pow(T-273.15,2.)-0.3918*(T-273.15)+87.663;
 
  dc = (noc_Iso<1) ? un : (el.mat->cb[noc_Iso-1].a[1]-el.mat->cb[noc_Iso-1].a[0])/(el.mat->cb[noc_Iso-1].n-1)/10.;
  

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
    c_cl = C_Cl(i);
    c_oh = C_OH(i);
    c_na = C_Na(i) ;
    c_h	 = k_e/c_oh;
    c_k  = C_K(i) ;

    p_v = P_v(i);
    p_g = P_g(i);
    p_a = p_g-p_v;

    c_e[i]   = (un - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_h*v_h + c_k*v_k))/v_h2o ;   
    lna_es[i] = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna[i]);	

    p_l[i] = RT*c_e[i]*(log(p_v/p_vs)-lna_es[i])+p_atm;
    p_c = p_g-p_l[i];

    s_w  = pow(pow(p_c/a,n)+1.,-m); 
    s_g  = un - s_w;
  
    s_cl     = (noc_Iso<1) ? zero : courbe(c_cl,el.mat->cb[noc_Iso-1]); 
    
    dfsdc    = ((noc_Iso<1) ? zero : courbe(c_cl+dc,el.mat->cb[noc_Iso-1]) - s_cl)/dc ;

    dsw_dpc = -m*n/a*pow(p_c/a,n-1.)*pow(pow(p_c/a,n)+1.,-m-1.);

    dlna_ccl[i] = (activite(c_cl+dc,c_oh,c_na,c_k,T,mode_a,lna[i])-lna_es[i])/dc;
    dlna_coh[i] = (activite(c_cl,c_oh+dc,c_na,c_k,T,mode_a,lna[i])-lna_es[i])/dc;
    dlna_cna[i] = (activite(c_cl,c_oh,c_na+dc,c_k,T,mode_a,lna[i])-lna_es[i])/dc;
    dlna_ck[i]  = (activite(c_cl,c_oh,c_na,c_k+dc,T,mode_a,lna[i])-lna_es[i])/dc;
  
    /*
      Conservation de Cl (chlore) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
    */
      K(i*NEQ+E_Cl,i*NEQ+I_Cl) += volume[i]*(phi*s_w + dfsdc) ;
      K(i*NEQ+E_Cl,i*NEQ+I_Pv) += volume[i]*(phi*c_cl)*(-dsw_dpc*c_e[i]*RT/p_v) ;
      K(i*NEQ+E_Cl,i*NEQ+I_Pg) += volume[i]*(phi*c_cl)*(dsw_dpc) ;

      K(i*NEQ+E_Cl,i*NEQ+I_Cl) += volume[i]*(phi*c_cl)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_Cl,i*NEQ+I_OH) += volume[i]*(phi*c_cl)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_Cl,i*NEQ+I_Na) += volume[i]*(phi*c_cl)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_Cl,i*NEQ+I_K)  += volume[i]*(phi*c_cl)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);

    /*
      Conservation de OH (oxy) : [(n_O1 - n_On)] + dt * div(w_O) = 0
    */
      K(i*NEQ+E_OH,i*NEQ+I_Cl) +=  volume[i]*(-dfsdc);
      K(i*NEQ+E_OH,i*NEQ+I_OH) +=  volume[i]*(phi*s_w*(1.+k_e/pow(c_oh,2.))) ;
      K(i*NEQ+E_OH,i*NEQ+I_Pv) +=  volume[i]*(phi*(c_oh-c_h))*(-dsw_dpc*c_e[i]*RT/p_v) ;
      K(i*NEQ+E_OH,i*NEQ+I_Pg) +=  volume[i]*(phi*(c_oh-c_h))*(dsw_dpc) ;  

      K(i*NEQ+E_OH,i*NEQ+I_Cl) += volume[i]*(phi*(c_oh-c_h))*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_OH,i*NEQ+I_OH) += volume[i]*(phi*(c_oh-c_h))*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_OH,i*NEQ+I_Na) += volume[i]*(phi*(c_oh-c_h))*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_OH,i*NEQ+I_K)  += volume[i]*(phi*(c_oh-c_h))*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);

     /*
      Conservation de Na (natri) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
    */
      K(i*NEQ+E_Na,i*NEQ+I_Na) +=  volume[i]*(phi*s_w) ;
      K(i*NEQ+E_Na,i*NEQ+I_Pv) +=  volume[i]*(phi*c_na)*(-dsw_dpc*c_e[i]*RT/p_v) ;
      K(i*NEQ+E_Na,i*NEQ+I_Pg) +=  volume[i]*(phi*c_na)*(dsw_dpc) ;

      K(i*NEQ+E_Na,i*NEQ+I_Cl) += volume[i]*(phi*c_na)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_Na,i*NEQ+I_OH) += volume[i]*(phi*c_na)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_Na,i*NEQ+I_Na) += volume[i]*(phi*c_na)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_Na,i*NEQ+I_K)  += volume[i]*(phi*c_na)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);

    /*
      Conservation de K (kali) : (n_K1 - n_Kn) + dt * div(w_K) = 0
    */
      K(i*NEQ+E_K,i*NEQ+I_K) +=  volume[i]*(phi*s_w) ;
      K(i*NEQ+E_K,i*NEQ+I_Pv) += volume[i]*(phi*c_k)*(-dsw_dpc*c_e[i]*RT/p_v) ;
      K(i*NEQ+E_K,i*NEQ+I_Pg) += volume[i]*(phi*c_k)*(dsw_dpc) ;

      K(i*NEQ+E_K,i*NEQ+I_Cl) += volume[i]*(phi*c_k)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_K,i*NEQ+I_OH) += volume[i]*(phi*c_k)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_K,i*NEQ+I_Na) += volume[i]*(phi*c_k)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_K,i*NEQ+I_K)  += volume[i]*(phi*c_k)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);

    /*
      Conservation de Vapeur (V) : (n_V1 - n_Vn) + dt * div(w_V) = 0
    */
 
      K(i*NEQ+E_WV,i*NEQ+I_Pv) += volume[i]*(phi*M_v/RT)*(s_g+p_v*dsw_dpc*c_e[i]*RT/p_v) ;
      K(i*NEQ+E_WV,i*NEQ+I_Pg) += volume[i]*(phi*M_v/RT*p_v)*(-dsw_dpc) ;

      K(i*NEQ+E_WV,i*NEQ+I_Cl) += volume[i]*(phi*M_v)*(-p_v*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_WV,i*NEQ+I_OH) += volume[i]*(phi*M_v)*(-p_v*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_WV,i*NEQ+I_Na) += volume[i]*(phi*M_v)*(-p_v*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_WV,i*NEQ+I_K)  += volume[i]*(phi*M_v)*(-p_v*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);


   /*
      Conservation de Vapeur (A) : (n_A1 - n_An) + dt * div(w_A) = 0
    */

      K(i*NEQ+E_A,i*NEQ+I_Pg) += volume[i]*(phi*M_a/RT)*(s_g-p_a*dsw_dpc) ;
      K(i*NEQ+E_A,i*NEQ+I_Pv) += volume[i]*(phi*M_a/RT)*(-s_g+p_a*dsw_dpc*c_e[i]*RT/p_v) ;

      K(i*NEQ+E_A,i*NEQ+I_Cl) += volume[i]*(phi*M_a)*(-p_a*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_A,i*NEQ+I_OH) += volume[i]*(phi*M_a)*(-p_a*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_A,i*NEQ+I_Na) += volume[i]*(phi*M_a)*(-p_a*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_A,i*NEQ+I_K)  += volume[i]*(phi*M_a)*(-p_a*dsw_dpc)*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);

   /*
      Conservation de Humidite (W) : (n_W1 - n_Wn) + dt * div(w_W) = 0
    */

      K(i*NEQ+E_WV,i*NEQ+I_Pv) += volume[i]*(phi*M_v*c_e[i])*(-dsw_dpc*c_e[i]*RT/p_v)  ;
      K(i*NEQ+E_WV,i*NEQ+I_Pg) += volume[i]*(phi*M_v*c_e[i])*(dsw_dpc)  ;

      K(i*NEQ+E_WV,i*NEQ+I_Cl) += volume[i]*(phi*c_e[i]*M_v)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_WV,i*NEQ+I_OH) += volume[i]*(phi*c_e[i]*M_v)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]);
      K(i*NEQ+E_WV,i*NEQ+I_Na) += volume[i]*(phi*c_e[i]*M_v)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]);
      K(i*NEQ+E_WV,i*NEQ+I_K)  += volume[i]*(phi*c_e[i]*M_v)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]);

      K(i*NEQ+E_WV,i*NEQ+I_Cl) += volume[i]*(phi*s_w*M_v)*(-v_cl/v_h2o);
      K(i*NEQ+E_WV,i*NEQ+I_OH) += volume[i]*(phi*s_w*M_v)*(-v_oh/v_h2o);
      K(i*NEQ+E_WV,i*NEQ+I_Na) += volume[i]*(phi*s_w*M_v)*(-v_na/v_h2o);
      K(i*NEQ+E_WV,i*NEQ+I_K)  += volume[i]*(phi*s_w*M_v)*(-v_k/v_h2o);
 
    /*
      Conservation de Charge (S) : (n_Sw1 - n_Swn) + dt * div(w_Sw) = 0
    */

      K(i*NEQ+E_PS,i*NEQ+I_Cl) +=  volume[i]*F*s_w*z_cl;
      K(i*NEQ+E_PS,i*NEQ+I_OH) +=  volume[i]*F*s_w*(z_oh-z_h*k_e/pow(c_oh,2.));
      K(i*NEQ+E_PS,i*NEQ+I_Na) +=  volume[i]*F*s_w*z_na;
      K(i*NEQ+E_PS,i*NEQ+I_K)  +=  volume[i]*F*s_w*z_k;      
      K(i*NEQ+E_PS,i*NEQ+I_Pv) +=  volume[i]*F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*(-dsw_dpc*c_e[i]*RT/p_v) ;
      K(i*NEQ+E_PS,i*NEQ+I_Pg) +=  volume[i]*F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*(dsw_dpc) ;

      K(i*NEQ+E_PS,i*NEQ+I_Cl) += volume[i]*F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_cl/v_h2o+ c_e[i]*dlna_ccl[i]);
      K(i*NEQ+E_PS,i*NEQ+I_OH) += volume[i]*F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_oh/v_h2o+ c_e[i]*dlna_coh[i]); 
      K(i*NEQ+E_PS,i*NEQ+I_Na) += volume[i]*F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_na/v_h2o+ c_e[i]*dlna_cna[i]); 
      K(i*NEQ+E_PS,i*NEQ+I_K)  += volume[i]*F*(c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*(dsw_dpc)*RT*((log(p_v/p_vs)-lna_es[i])*v_k/v_h2o+ c_e[i]*dlna_ck[i]); 

  }

  /* termes d'ecoulement */
  tr     = dt*surf/dx ;
   
 /*
    Conservation de Cl (chlore) : (n_Cl1 - n_Cln) + dt * div(w_Cl) = 0
  */
      
  K(E_Cl,I_Cl)             += + tr*KF_Cl ;
  K(E_Cl,I_Cl+NEQ)         += - tr*KF_Cl ;
  K(E_Cl+NEQ,I_Cl)         += - tr*KF_Cl ;
  K(E_Cl+NEQ,I_Cl+NEQ)     += + tr*KF_Cl ;

  K(E_Cl,I_SI)              += + tr*KM_Cl;
  K(E_Cl,I_SI+NEQ)          += - tr*KM_Cl;
  K(E_Cl+NEQ,I_SI)          += - tr*KM_Cl;
  K(E_Cl+NEQ,I_SI+NEQ)      += + tr*KM_Cl;

  K(E_Cl,I_Pv)              += + tr*KD_Cl*c_e[0]*RT/P_v(0);
  K(E_Cl,I_Pv+NEQ)          += - tr*KD_Cl*c_e[1]*RT/P_v(1);
  K(E_Cl+NEQ,I_Pv)          += - tr*KD_Cl*c_e[0]*RT/P_v(0);
  K(E_Cl+NEQ,I_Pv+NEQ)      += + tr*KD_Cl*c_e[1]*RT/P_v(1);

  K(E_Cl,I_Cl)              += - tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_Cl,I_Cl+NEQ)          += + tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);
  K(E_Cl+NEQ,I_Cl)          += + tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_Cl+NEQ,I_Cl+NEQ)      += - tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);

  K(E_Cl,I_OH)              += - tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_Cl,I_OH+NEQ)          += + tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);
  K(E_Cl+NEQ,I_OH)          += + tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_Cl+NEQ,I_OH+NEQ)      += - tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);

  K(E_Cl,I_Na)              += - tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_Cl,I_Na+NEQ)          += + tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);
  K(E_Cl+NEQ,I_Na)          += + tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_Cl+NEQ,I_Na+NEQ)      += - tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);

  K(E_Cl,I_K)              += - tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_Cl,I_K+NEQ)          += + tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);
  K(E_Cl+NEQ,I_K)          += + tr*KD_Cl*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_Cl+NEQ,I_K+NEQ)      += - tr*KD_Cl*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);


     /*
      Conservation de OH (Oxy) : [(n_O1 - n_On)] + dt * div(w_O) = 0
    */

  K(E_OH,I_OH)         += + tr*(KF_OH+KF_H*k_e/pow(C_OH(0),2));
  K(E_OH,I_OH+NEQ)     += - tr*(KF_OH+KF_H*k_e/pow(C_OH(1),2));
  K(E_OH+NEQ,I_OH)     += - tr*(KF_OH+KF_H*k_e/pow(C_OH(0),2));
  K(E_OH+NEQ,I_OH+NEQ) += + tr*(KF_OH+KF_H*k_e/pow(C_OH(1),2));

  K(E_OH,I_SI)              += + tr*(KM_OH-KM_H);
  K(E_OH,I_SI+NEQ)          += - tr*(KM_OH-KM_H);
  K(E_OH+NEQ,I_SI)          += - tr*(KM_OH-KM_H);
  K(E_OH+NEQ,I_SI+NEQ)      += + tr*(KM_OH-KM_H);

  K(E_OH,I_Pv)             += + tr*(KD_OH-KD_H)*c_e[0]*RT/P_v(0);
  K(E_OH,I_Pv+NEQ)         += - tr*(KD_OH-KD_H)*c_e[1]*RT/P_v(1);
  K(E_OH+NEQ,I_Pv)         += - tr*(KD_OH-KD_H)*c_e[0]*RT/P_v(0);
  K(E_OH+NEQ,I_Pv+NEQ)     += + tr*(KD_OH-KD_H)*c_e[1]*RT/P_v(1);

  K(E_OH,I_Cl)              += - tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_OH,I_Cl+NEQ)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);
  K(E_OH+NEQ,I_Cl)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_OH+NEQ,I_Cl+NEQ)      += - tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);

  K(E_OH,I_OH)              += - tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_OH,I_OH+NEQ)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);
  K(E_OH+NEQ,I_OH)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_OH+NEQ,I_OH+NEQ)      += - tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);

  K(E_OH,I_Na)              += - tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_OH,I_Na+NEQ)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);
  K(E_OH+NEQ,I_Na)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_OH+NEQ,I_Na+NEQ)      += - tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);

  K(E_OH,I_K)              += - tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_OH,I_K+NEQ)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);
  K(E_OH+NEQ,I_K)          += + tr*(KD_OH-KD_H)*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_OH+NEQ,I_K+NEQ)      += - tr*(KD_OH-KD_H)*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);

   /*
    Conservation de Na (Natri) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  K(E_Na,I_Na)             += + tr*KF_Na;
  K(E_Na,I_Na+NEQ)         += - tr*KF_Na;
  K(E_Na+NEQ,I_Na)         += - tr*KF_Na;
  K(E_Na+NEQ,I_Na+NEQ)     += + tr*KF_Na;

  K(E_Na,I_SI)              += + tr*KM_Na;
  K(E_Na,I_SI+NEQ)          += - tr*KM_Na;
  K(E_Na+NEQ,I_SI)          += - tr*KM_Na;
  K(E_Na+NEQ,I_SI+NEQ)      += + tr*KM_Na;

  K(E_Na,I_Pv)             += +tr*KD_Na*c_e[0]*RT/P_v(0);
  K(E_Na,I_Pv+NEQ)         += -tr*KD_Na*c_e[1]*RT/P_v(1);
  K(E_Na+NEQ,I_Pv)         += -tr*KD_Na*c_e[0]*RT/P_v(0);
  K(E_Na+NEQ,I_Pv+NEQ)     += +tr*KD_Na*c_e[1]*RT/P_v(1); 

  K(E_Na,I_Cl)              += - tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_Na,I_Cl+NEQ)          += + tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);
  K(E_Na+NEQ,I_Cl)          += + tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_Na+NEQ,I_Cl+NEQ)      += - tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);

  K(E_Na,I_OH)              += - tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_Na,I_OH+NEQ)          += + tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);
  K(E_Na+NEQ,I_OH)          += + tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_Na+NEQ,I_OH+NEQ)      += - tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);

  K(E_Na,I_Na)              += - tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_Na,I_Na+NEQ)          += + tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);
  K(E_Na+NEQ,I_Na)          += + tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_Na+NEQ,I_Na+NEQ)      += - tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);

  K(E_Na,I_K)              += - tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_Na,I_K+NEQ)          += + tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);
  K(E_Na+NEQ,I_K)          += + tr*KD_Na*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_Na+NEQ,I_K+NEQ)      += - tr*KD_Na*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);

   /*
    Conservation de K (kali) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */

  K(E_K,I_K)              += + tr*KF_K;
  K(E_K,I_K+NEQ)          += - tr*KF_K;
  K(E_K+NEQ,I_K)          += - tr*KF_K;
  K(E_K+NEQ,I_K+NEQ)      += + tr*KF_K;

  K(E_K,I_SI)              += + tr*KM_K;
  K(E_K,I_SI+NEQ)          += - tr*KM_K;
  K(E_K+NEQ,I_SI)          += - tr*KM_K;
  K(E_K+NEQ,I_SI+NEQ)      += + tr*KM_K; 

  K(E_K,I_Pv)             += +tr*KD_K*c_e[0]*RT/P_v(0);
  K(E_K,I_Pv+NEQ)         += -tr*KD_K*c_e[1]*RT/P_v(1);
  K(E_K+NEQ,I_Pv)         += -tr*KD_K*c_e[0]*RT/P_v(0);
  K(E_K+NEQ,I_Pv+NEQ)     += +tr*KD_K*c_e[1]*RT/P_v(1); 

  K(E_K,I_Cl)              += - tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_K,I_Cl+NEQ)          += + tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);
  K(E_K+NEQ,I_Cl)          += + tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_K+NEQ,I_Cl+NEQ)      += - tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);

  K(E_K,I_OH)              += - tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_K,I_OH+NEQ)          += + tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);
  K(E_K+NEQ,I_OH)          += + tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_K+NEQ,I_OH+NEQ)      += - tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);

  K(E_K,I_Na)              += - tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_K,I_Na+NEQ)          += + tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);
  K(E_K+NEQ,I_Na)          += + tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_K+NEQ,I_Na+NEQ)      += - tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);

  K(E_K,I_K)              += - tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_K,I_K+NEQ)          += + tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);
  K(E_K+NEQ,I_K)          += + tr*KD_K*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_K+NEQ,I_K+NEQ)      += - tr*KD_K*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);

  /*
     Poisson  : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */

  K(E_PS,I_SI)             += +surf/dx*KF_PS;
  K(E_PS,I_SI+NEQ)         += -surf/dx*KF_PS;
  K(E_PS+NEQ,I_SI)         += -surf/dx*KF_PS;
  K(E_PS+NEQ,I_SI+NEQ)     += +surf/dx*KF_PS; 

   /*
    Conservation de Vapeur (V) : (n_V1 - n_Vn) + dt * div(w_V) = 0
  */

  K(E_WV,I_Pg)             += +tr*(KD_V + KF_VA*(-M_v*M_a*P_v(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.)));
  K(E_WV,I_Pg+NEQ)         += -tr*(KD_V + KF_VA*(-M_v*M_a*P_v(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.)));
  K(E_WV+NEQ,I_Pg)         += -tr*(KD_V + KF_VA*(-M_v*M_a*P_v(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.)));
  K(E_WV+NEQ,I_Pg+NEQ)     += +tr*(KD_V + KF_VA*(-M_v*M_a*P_v(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.))); 

  K(E_WV,I_Pv)             += +tr*KF_VA*(M_v*M_a*P_g(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.));
  K(E_WV,I_Pv+NEQ)         += -tr*KF_VA*(M_v*M_a*P_g(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.));
  K(E_WV+NEQ,I_Pv)         += -tr*KF_VA*(M_v*M_a*P_g(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.));
  K(E_WV+NEQ,I_Pv+NEQ)     += +tr*KF_VA*(M_v*M_a*P_g(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.)); 

  /*
    Conservation de Aire (A) : (n_A1 - n_An) + dt * div(w_A) = 0
  */

  K(E_A,I_Pg)             += +tr*(KD_A + KF_VA*(M_v*M_a*P_v(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.)));
  K(E_A,I_Pg+NEQ)         += -tr*(KD_A + KF_VA*(M_v*M_a*P_v(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.)));
  K(E_A+NEQ,I_Pg)         += -tr*(KD_A + KF_VA*(M_v*M_a*P_v(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.)));
  K(E_A+NEQ,I_Pg+NEQ)     += +tr*(KD_A + KF_VA*(M_v*M_a*P_v(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.))); 

  K(E_A,I_Pv)             += +tr*KF_VA*(-M_v*M_a*P_g(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.));
  K(E_A,I_Pv+NEQ)         += -tr*KF_VA*(-M_v*M_a*P_g(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.));
  K(E_A+NEQ,I_Pv)         += -tr*KF_VA*(-M_v*M_a*P_g(0)/pow(P_g(0)*M_a+(M_v-M_a)*P_v(0),2.));
  K(E_A+NEQ,I_Pv+NEQ)     += +tr*KF_VA*(-M_v*M_a*P_g(1)/pow(P_g(1)*M_a+(M_v-M_a)*P_v(1),2.)); 

  /*
    Conservation de Water (W) : (n_W1 - n_Wn) + dt * div(w_W) = 0
  */

  K(E_WV,I_Pv)             += +tr*KD_W*c_e[0]*RT/P_v(0);
  K(E_WV,I_Pv+NEQ)         += -tr*KD_W*c_e[1]*RT/P_v(1);
  K(E_WV+NEQ,I_Pv)         += -tr*KD_W*c_e[0]*RT/P_v(0);
  K(E_WV+NEQ,I_Pv+NEQ)     += +tr*KD_W*c_e[1]*RT/P_v(1);

  K(E_WV,I_Cl)              += - tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_WV,I_Cl+NEQ)          += + tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);
  K(E_WV+NEQ,I_Cl)          += + tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_cl/v_h2o+ c_e[0]*dlna_ccl[0]);
  K(E_WV+NEQ,I_Cl+NEQ)      += - tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_cl/v_h2o+ c_e[1]*dlna_ccl[1]);

  K(E_WV,I_OH)              += - tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_WV,I_OH+NEQ)          += + tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);
  K(E_WV+NEQ,I_OH)          += + tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_oh/v_h2o+ c_e[0]*dlna_coh[0]);
  K(E_WV+NEQ,I_OH+NEQ)      += - tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_oh/v_h2o+ c_e[1]*dlna_coh[1]);

  K(E_WV,I_Na)              += - tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_WV,I_Na+NEQ)          += + tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);
  K(E_WV+NEQ,I_Na)          += + tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_na/v_h2o+ c_e[0]*dlna_cna[0]);
  K(E_WV+NEQ,I_Na+NEQ)      += - tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_na/v_h2o+ c_e[1]*dlna_cna[1]);

  K(E_WV,I_K)              += - tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_WV,I_K+NEQ)          += + tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);
  K(E_WV+NEQ,I_K)          += + tr*KD_W*RT*((log(P_v(0)/p_vs)-lna_es[0])*v_k/v_h2o+ c_e[0]*dlna_ck[0]);
  K(E_WV+NEQ,I_K+NEQ)      += - tr*KD_W*RT*((log(P_v(1)/p_vs)-lna_es[1])*v_k/v_h2o+ c_e[1]*dlna_ck[1]);
  
  K(E_WV,I_Cl)             += - M_Cl*tr*KF_Cl ;
  K(E_WV,I_Cl+NEQ)         += + M_Cl*tr*KF_Cl ;
  K(E_WV+NEQ,I_Cl)         += + M_Cl*tr*KF_Cl ;
  K(E_WV+NEQ,I_Cl+NEQ)     += - M_Cl*tr*KF_Cl ;

  K(E_WV,I_SI)              += - M_Cl*tr*KM_Cl;
  K(E_WV,I_SI+NEQ)          += + M_Cl*tr*KM_Cl;
  K(E_WV+NEQ,I_SI)          += + M_Cl*tr*KM_Cl;
  K(E_WV+NEQ,I_SI+NEQ)      += - M_Cl*tr*KM_Cl;

  K(E_WV,I_OH)         += - M_OH*tr*(KF_OH+KF_H*k_e/pow(C_OH(0),2));
  K(E_WV,I_OH+NEQ)     += + M_OH*tr*(KF_OH+KF_H*k_e/pow(C_OH(1),2));
  K(E_WV+NEQ,I_OH)     += + M_OH*tr*(KF_OH+KF_H*k_e/pow(C_OH(0),2));
  K(E_WV+NEQ,I_OH+NEQ) += - M_OH*tr*(KF_OH+KF_H*k_e/pow(C_OH(1),2));

  K(E_WV,I_SI)              += - M_OH*tr*(KM_OH-KM_H);
  K(E_WV,I_SI+NEQ)          += + M_OH*tr*(KM_OH-KM_H);
  K(E_WV+NEQ,I_SI)          += + M_OH*tr*(KM_OH-KM_H);
  K(E_WV+NEQ,I_SI+NEQ)      += - M_OH*tr*(KM_OH-KM_H);

  K(E_WV,I_Na)             += - M_Na*tr*KF_Na;
  K(E_WV,I_Na+NEQ)         += + M_Na*tr*KF_Na;
  K(E_WV+NEQ,I_Na)         += + M_Na*tr*KF_Na;
  K(E_WV+NEQ,I_Na+NEQ)     += - M_Na*tr*KF_Na;

  K(E_WV,I_SI)              += - M_Na*tr*KM_Na;
  K(E_WV,I_SI+NEQ)          += + M_Na*tr*KM_Na;
  K(E_WV+NEQ,I_SI)          += + M_Na*tr*KM_Na;
  K(E_WV+NEQ,I_SI+NEQ)      += - M_Na*tr*KM_Na;

  K(E_WV,I_K)              += - M_K*tr*KF_K;
  K(E_WV,I_K+NEQ)          += + M_K*tr*KF_K;
  K(E_WV+NEQ,I_K)          += + M_K*tr*KF_K;
  K(E_WV+NEQ,I_K+NEQ)      += - M_K*tr*KF_K;

  K(E_WV,I_SI)              += - M_K*tr*KM_K;
  K(E_WV,I_SI+NEQ)          += + M_K*tr*KM_K;
  K(E_WV+NEQ,I_SI)          += + M_K*tr*KM_K;
  K(E_WV+NEQ,I_SI+NEQ)      += - M_K*tr*KM_K; 
  
  return(0) ;

#undef K
}

void rs44(double **x,double **u,double **u_n,double *f,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
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
  R(0,E_Cl) -= volume[0]*(N_Cl(0) - N_Cln(0)) + dt*surf*(W_Cl +JA_Cl);
  R(1,E_Cl) -= volume[1]*(N_Cl(1) - N_Cln(1)) - dt*surf*(W_Cl +JA_Cl);

  /*
      Conservation de O (oxy) : [(n_O1 - n_On)] + dt * div(w_O) = 0
    */
  R(0,E_OH) -= volume[0]*(N_OH(0) - N_OHn(0)) + dt*surf*(W_OH+JA_OH) ;
  R(1,E_OH) -= volume[1]*(N_OH(1) - N_OHn(1)) - dt*surf*(W_OH+JA_OH) ;

  /*
    Conservation de Na (natrie) : (n_Na1 - n_Nan) + dt * div(w_Na) = 0
  */
  R(0,E_Na) -= volume[0]*(N_Na(0) - N_Nan(0)) + dt*surf*(W_Na+JA_Na) ;
  R(1,E_Na) -= volume[1]*(N_Na(1) - N_Nan(1)) - dt*surf*(W_Na+JA_Na);
 /*
    Conservation de K (kali) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_K) -= volume[0]*(N_K(0) - N_Kn(0)) + dt*surf*(W_K+JA_K) ;
  R(1,E_K) -= volume[1]*(N_K(1) - N_Kn(1)) - dt*surf*(W_K+JA_K);

 /*
    Conservation de charge (kali) : (n_K1 - n_Kn) + dt * div(w_K) = 0
  */
  R(0,E_PS) -= volume[0]*N_PS(0)+surf*W_PS;
  R(1,E_PS) -= volume[1]*N_PS(1)-surf*W_PS;
 
  /*
      Conservation de A (A) : [(n_A1 - n_An)] + dt * div(w_A) = 0
    */
  R(0,E_A) -= volume[0]*(N_A(0) - N_An(0)) + dt*surf*W_A;
  R(1,E_A) -= volume[1]*(N_A(1) - N_An(1)) - dt*surf*W_A;

  /*
      Conservation de W (W) : [(n_W1 - n_Wn)] + dt * div(w_W) = 0
    */
  R(0,E_WV) -= volume[0]*(N_W(0)+N_V(0) - N_Wn(0)-N_Vn(0)) + dt*surf*(W_W +W_V) ;
  R(1,E_WV) -= volume[1]*(N_W(1)+N_V(1) - N_Wn(1)-N_Vn(1)) - dt*surf*(W_W +W_V);

#undef R
}

int so44(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  int    i,j,nso ;
  double c_cl,s_cl,c_oh,c_na,c_k,c_h,si,lna_es,RT;
  double s_w,s_g,p_g,p_v,p_l, p_c,p_w[2],c_e;
  double zero = 0.,un = 1.,deux = 2. ;
  double lna[4],r_Cl;
  double grd_cl,grd_oh,grd_na,grd_k,grd_h,grd_si,grd_pw;
  double dx,epsi;
  double J_Cl,J_OH,J_Na,J_K,J_W;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("porosite")] ;
  d_cl    = el.mat->pr[pm("D_Cl")] ;
  T	  = el.mat->pr[pm("T")] ;
  p_vs	  = 609.14*pow(10.,7.45*(T-273.)/(T-38.));
  kw_int  = el.mat->pr[pm("kw_int")] ;
  kg_int  = el.mat->pr[pm("kg_int")] ;
  d_va	  = el.mat->pr[pm("d_va")]; /* 2.17e-3*/
  klinkenberg  = el.mat->pr[pm("klinkenberg")];
  a  = el.mat->pr[pm("a")];
  n  = el.mat->pr[pm("n")];
  m  = el.mat->pr[pm("m")];

  noc_Iso = floor(el.mat->pr[pm("Isotherme")]+0.5) ;
  mode_a = floor(el.mat->pr[pm("activite")]+0.5) ;

  RT      = R_g*T;
  epsi = 0.0007*pow(T-273.15,2.)-0.3918*(T-273.15)+87.663;

 
  /* initialisation */
  nso = 18;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* concentration */

  for(i=0;i<el.nn;i++) {
    c_cl = C_Cl(i);
    c_oh = C_OH(i);
    c_na = C_Na(i) ;
    c_h	 = k_e/c_oh;
    c_k  = C_K(i) ;
    c_e   = (un - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_h*v_h + c_k*v_k))/v_h2o ;   
    lna_es = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna);
    p_v = P_v(i);
    p_g = P_g(i);
    p_w[i] = RT*c_e*(log(p_v/p_vs)-lna_es)+p_atm;
  }

  if(s[0] < (x[0][0] + x[1][0])/deux) 
    {
      c_cl = C_Cl(0) ;
      c_oh = C_OH(0) ;
      c_na = C_Na(0) ;
      c_k  = C_K(0) ;
      si  = SI(0);
      p_v  = P_v(0);
      p_g  = P_g(0);
      lna_es = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna);     
   }
  else 
    {
      c_cl = C_Cl(1) ;
      c_oh = C_OH(1) ;
      c_na = C_Na(1) ;
      c_k  = C_K(1) ;
      si   = SI(1);
      p_v  = P_v(1);
      p_g  = P_g(1);
      lna_es = activite(c_cl,c_oh,c_na,c_k,T,mode_a,lna);     
   }

  c_h   = k_e/c_oh;
  c_e   = (un - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_h*v_h + c_k*v_k))/v_h2o ;   
  p_l = RT*c_e*(log(p_v/p_vs)-lna_es)+p_atm;
  p_c = p_g-p_l;

  s_w  = pow(pow(p_c/a,n)+1.,-m); 

  s_cl   = (noc_Iso<1) ? zero : courbe(c_cl,el.mat->cb[noc_Iso-1]);
  s_g = un - s_w;

    dx     = x[1][0] - x[0][0];
		
    grd_cl = (C_Cl(1) - C_Cl(0))/dx ;
    grd_oh = (C_OH(1) - C_OH(0))/dx ;
    grd_na = (C_Na(1) - C_Na(0))/dx ;
    grd_k  = (C_K(1) - C_K(0))/dx ;
    grd_h  =  k_e*(un/C_OH(1) - un/C_OH(0))/dx ;

    grd_si = (SI(1)-SI(0))/dx;
  

   grd_pw = (p_w[1]-p_w[0])/dx;

    J_Cl =  - KF_Cl*grd_cl - KM_Cl*grd_si;
    J_OH =  - KF_OH*grd_oh - KM_OH*grd_si;
    J_Na =  - KF_Na*grd_na - KM_Na*grd_si;
    J_K  =  - KF_K*grd_k - KM_K*grd_si;
  
    J_W = -(M_Cl*(J_Cl+JA_Cl)+M_OH*(J_OH+JA_OH)+M_Na*(J_Na+JA_Na)+M_K*(J_K+JA_K));

    if (J_Cl+JA_Cl ==0.) r_Cl =0. ; else  r_Cl = -KD_Cl*grd_pw/(JA_Cl+J_Cl);
      
  /* quantites exploitees */
  strcpy(r[0].text,"Cl libre") ; r[0].n = 1 ;
  r[0].v[0] = c_cl ;
  strcpy(r[1].text,"Cl totaux") ; r[1].n = 1 ;
  r[1].v[0] = s_cl+s_w*phi*c_cl; 
  strcpy(r[2].text,"c_OH") ; r[2].n = 1 ;
  r[2].v[0] = c_oh ;
  strcpy(r[3].text,"c_Na") ; r[3].n = 1 ;
  r[3].v[0] = c_na ;
  strcpy(r[4].text,"c_K") ; r[4].n = 1 ;
  r[4].v[0] = c_k ;
  strcpy(r[5].text,"c_h2o") ; r[5].n = 1 ;
  r[5].v[0] = c_e; 
  strcpy(r[6].text,"charge d'elctrique") ; r[6].n = 1 ;
  r[6].v[0] = (c_cl*z_cl+c_oh*z_oh+c_na*z_na+c_k*z_k+c_h*z_h)*F;  
  strcpy(r[7].text,"Saturation") ; r[7].n = 1 ;
  r[7].v[0] = s_w;
  strcpy(r[8].text,"SI") ; r[8].n = 1 ;
  r[8].v[0] = si;
  strcpy(r[9].text,"I") ; r[9].n = 1 ;
  r[9].v[0] = ((W_Cl+JA_Cl)*z_cl+(W_Na+JA_Na)*z_na+(W_K+JA_K)*z_k+(W_OH+JA_OH)*z_oh)*F;
  strcpy(r[10].text,"h_r") ; r[10].n = 1 ;
  r[10].v[0] = p_v/p_vs;
  strcpy(r[11].text,"p_g") ; r[11].n = 1 ;
  r[11].v[0] = p_g/p_atm;
  strcpy(r[12].text,"p_w") ; r[12].n = 1 ;
  r[12].v[0] = p_l;
  strcpy(r[13].text,"Flux vapeur") ; r[13].n = 1 ;
  r[13].v[0] = W_V;
  strcpy(r[14].text,"Flux liquide") ; r[14].n = 1 ;
  r[14].v[0] = W_W;
  strcpy(r[15].text,"H2O diffusion") ; r[15].n = 1 ;
  r[15].v[0] = J_W;
  strcpy(r[16].text,"Flux toltal") ; r[16].n = 1 ;
  r[16].v[0] = W_W+W_A+W_V;
  strcpy(r[17].text,"Peclet") ; r[17].n = 1 ;
  r[17].v[0] = fabs(r_Cl);

  return(nso) ;
}


double activite(double c_cl,double c_oh,double c_na,double c_k,double T,int model,double *lna)
/* L'activite chimique de l'eau d'une solution electrolytique 
   composee de NaCl, KCl et NaOH (cas du beton)
   (les concentrations c_i doivent etre en moles/m3) */
{
#define lna_cl    (lna[0])
#define lna_oh    (lna[1])
#define lna_na    (lna[2])
#define lna_k     (lna[3])

/* valences */
#define z_cl    (-1.)
#define z_oh    (-1.)
#define z_na    (+1.)
#define z_k     (+1.)

/* volumes molaires (m3/mol) */
#define v_h2o   (18.00e-6)
#define v_na    (22.47e-6)
#define v_k     (54.11e-6)
#define v_oh    (23.50e-6)
#define v_cl    (-0.35e-6)

/* masse molaire */
#define M_h2o   (18.e-3)

  double c_w ;
  double m_cl,m_oh,m_na,m_k,m_T ;
  double I,A,epsi,lna_w ;
  double b_cl,b_oh,b_na,b_k,S_cl,S_oh,S_na,S_k ;
  double y_cl,y_oh,y_na,y_k ;

  double T_0 = 273.15 ;
  double b0 = sqrt(M_h2o),S0 = pow(M_h2o,1.29) ; /* references */
  double b_na_nacl = 4.352/b0,b_cl_nacl = 1.827/b0 ;
  double b_na_naoh = 0.971/b0,b_oh_naoh = 6.052/b0 ;
  double b_k_kcl   = 1.243/b0,b_cl_kcl  = 3.235/b0 ;
  double b_k_koh   = 0.002/b0,b_oh_koh  = 22.347/b0 ;
  double S_na_nacl = 26.448/S0,S_cl_nacl = 19.245/S0 ;
  double S_na_naoh = 59.306/S0,S_oh_naoh = 15.685/S0 ;
  double S_k_kcl   = 13.296/S0,S_cl_kcl  = 11.158/S0 ;
  double S_k_koh   = 36.479/S0,S_oh_koh  = 159.038/S0 ;

  double zero = 0. ;

  /* molarites */
  if(c_cl < zero) c_cl = zero ;
  if(c_oh < zero) c_oh = zero ;
  if(c_na < zero) c_na = zero ;
  if(c_k  < zero) c_k  = zero ;

  epsi = 0.0007*(T - T_0)*(T - T_0) - 0.3918*(T - T_0) + 87.663 ;
  A = 1398779.816/pow(epsi*T,1.5)/b0 ;

  /* concentration molaire du solvant (les c_i doivent etre en mole/m3) */
  c_w = (1. - (c_cl*v_cl + c_oh*v_oh + c_na*v_na + c_k*v_k))/v_h2o ;
  
  /* molalites*M_h2o (moles/mole) */
  m_cl = c_cl/c_w ;
  m_oh = c_oh/c_w ;
  m_na = c_na/c_w ;
  m_k  = c_k/c_w ;
  
  /* force ionique (divisee par M_h2o) */
  I = 0.5*(z_cl*z_cl*m_cl + z_oh*z_oh*m_oh + z_na*z_na*m_na + z_k*z_k*m_k) ;
  
  if (I > zero && model > 0) {
    m_T = m_cl + m_oh + m_na + m_k ;
    
    y_cl = ((c_na + c_k)  > 0) ? c_cl/(c_na + c_k)  : 0. ;
    y_oh = ((c_na + c_k)  > 0) ? c_oh/(c_na + c_k)  : 0. ;
    y_na = ((c_cl + c_oh) > 0) ? c_na/(c_cl + c_oh) : 0. ;
    y_k  = ((c_cl + c_oh) > 0) ? c_k/(c_cl + c_oh)  : 0. ;
    
    b_cl = b_cl_nacl*y_na + b_cl_kcl*y_k ;
    b_oh = b_oh_naoh*y_na + b_oh_koh*y_k ;
    b_na = b_na_naoh*y_oh + b_na_nacl*y_cl ;
    b_k  = b_k_koh*y_oh   + b_k_kcl*y_cl ;
    
    S_cl = S_cl_nacl*y_na + S_cl_kcl*y_k ;
    S_oh = S_oh_naoh*y_na + S_oh_koh*y_k ;
    S_na = S_na_naoh*y_oh + S_na_nacl*y_cl ;
    S_k  = S_k_koh*y_oh   + S_k_kcl*y_cl ;

    lna_w = m_cl*lna_i(T,I,z_cl,b_cl,S_cl,A) 
          + m_oh*lna_i(T,I,z_oh,b_oh,S_oh,A)
          + m_na*lna_i(T,I,z_na,b_na,S_na,A)
          + m_k*lna_i(T,I,z_k,b_k,S_k,A) ;

    if(model == 1) {
      /* selon Lin & Lee */
      lna_cl = lng_LinLee(T,I,z_cl,b_cl,S_cl,A) ;
      lna_oh = lng_LinLee(T,I,z_oh,b_oh,S_oh,A) ;
      lna_na = lng_LinLee(T,I,z_na,b_na,S_na,A) ;
      lna_k  = lng_LinLee(T,I,z_k,b_k,S_k,A) ;
      
      /* selon TQN */
    } else if(model == 2) {
      lna_cl = lng_TQN(T,I,z_cl,b_cl,S_cl,A,lna_w,m_T) ;
      lna_na = lng_TQN(T,I,z_na,b_na,S_na,A,lna_w,m_T) ;
      lna_oh = lng_TQN(T,I,z_oh,b_oh,S_oh,A,lna_w,m_T) ;
      lna_k  = lng_TQN(T,I,z_k,b_k,S_k,A,lna_w,m_T) ;
    } else {
      arret("activite : model non connu") ;
    }
  } else {
    lna_cl = 0.;
    lna_oh = 0.;
    lna_na = 0.;
    lna_k  = 0.;
    lna_w  = 0.;
  }
  
  return(lna_w) ;

#undef lna_cl
#undef lna_oh
#undef lna_na
#undef lna_k 

#undef z_cl
#undef z_oh
#undef z_na
#undef z_k

#undef v_h2o
#undef v_na 
#undef v_k  
#undef v_oh 
#undef v_cl 

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
