#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CommonModel.h"

#define MODELINDEX  24
#define TITLE "Transport ionique dans le beton (3D)"
#define AUTHORS "Dridi"

#include "OldMethods.h"

/* Macros */
#define NEQ       (6)

#define NVI       (6+6*dim)
#define NVE       (31)

#define E_h       (0)
#define E_fe      (1)
#define E_i       (2)
#define E_o2      (3)
#define E_cat     (4)
#define E_ani     (5)

#define I_p_l     (0)
#define I_c_fe    (1)
#define I_psi     (2)
#define I_c_o2    (3)
#define I_c_cat   (4)
#define I_c_ani   (5)

/* Fonctions */
static int    pm(const char *s) ;
static int    c24(double **,double **,double **,double *,double *,double *,elem_t,geom_t,double,double *) ;
static int    k24(double **,double **,double **,double *,double *,double *,elem_t,geom_t,double,double *) ;
static double tortuosite_l(double,mate_t*) ;
static double tortuosite_g(double,mate_t*) ;
static double concentration_oh(double,double,double,double) ;
/* Parametres */
static double phi,c_h2o,k_int,mu_l,RT_0,k_eau,k_feoh2,farad,k_hen,d0_oh,d0_h,d0_o2_l,d0_o2_g,d0_fe,d0_feoh2,d0_cat,d0_ani,p_g=0 ;
static const double M_h = 1,M_o = 16,M_fe = 55.85,M_cat = 40,M_ani = 40 ;
static const double M_h2o = 18,M_o2 = 32,M_oh = 17,M_feoh2 = 89.85 ;

int pm(const char *s)
{
if(strcmp(s,"phi") == 0) return (0) ;
else if(strcmp(s,"c_h2o") == 0) return (1) ;
else if(strcmp(s,"k_int") == 0) return (2) ;
else if(strcmp(s,"mu_l") == 0) return (3) ;
else if(strcmp(s,"RT_0") == 0) return (4) ;
else if(strcmp(s,"K_eau") == 0) return (5) ;
else if(strcmp(s,"K_Fe(OH)2") == 0) return (6) ;
else if(strcmp(s,"K_Far") == 0) return (7) ;
else if(strcmp(s,"K_Hen") == 0) return (8) ;
else if(strcmp(s,"D_oh") == 0) return (9) ;
else if(strcmp(s,"D_h") == 0) return (10) ;
else if(strcmp(s,"D_o2_l") == 0) return (11) ;
else if(strcmp(s,"D_fe") == 0) return (12) ;
else if(strcmp(s,"D_feoh2") == 0) return (13) ;
else if(strcmp(s,"D_cat") == 0) return (14) ;
else if(strcmp(s,"D_o2_g") == 0) return (15) ;
else if(strcmp(s,"D_ani") == 0) return (16) ;
else if(strcmp(s,"courbes") == 0) return (17) ;
else
{ printf("donnee \"%s\" non connue (pm24)\n",s) ; exit(0) ; }
}


int dm24(int dim,mate_t *mat,FILE *ficd)
/* Lecture des donnees materiaux dans le fichier ficd */
{
  int n_donnees = 18 ;

  mat->neq = NEQ ;

  strcpy(mat->eqn[E_h]    ,"E_h") ;
  strcpy(mat->eqn[E_fe]   ,"E_fe") ;
  strcpy(mat->eqn[E_i]    ,"E_i") ;
  strcpy(mat->eqn[E_o2]   ,"E_o2") ;
  strcpy(mat->eqn[E_cat]  ,"E_cat") ;
  strcpy(mat->eqn[E_ani]  ,"E_ani") ;

  strcpy(mat->inc[I_p_l]  ,"p_l") ;
  strcpy(mat->inc[I_c_fe] ,"c_fe") ;
  strcpy(mat->inc[I_psi]  ,"psi") ;
  strcpy(mat->inc[I_c_o2] ,"c_o2") ;
  strcpy(mat->inc[I_c_cat],"c_cat") ;
  strcpy(mat->inc[I_c_ani],"c_ani") ;

  dmat(mat,ficd,pm,n_donnees) ;

  return(mat->n) ;
}


int qm24(int dim,FILE *ficd)
/* Saisie des donnees materiaux */
{
  
  printf(TITLE) ;
  
  if(!ficd) return(NEQ) ;

  printf("\n\n\
Le systeme est forme de :\n\
\t 1. Conservation de la masse de O2 (c_o2) \n\
\t 2. Conservation de la masse de H  (p_l)  \n\
\t 3. Conservation de la masse de Fe (c_fe) \n\
\t 4. Conservation de la charge Q    (psi) \n\
\t 5. Conservation de la masse de K+ (c_cat) \n\
\t 6. Conservation de la masse de A- (c_ani)\n") ;
  
  printf("\n\
Exemple de donnees\n\n") ;

  fprintf(ficd,"phi = 0.28            # La porosite\n") ;
  fprintf(ficd,"c_h2o = 5.55e4        # Concentration molaire de l\'eau\n") ;
  fprintf(ficd,"k_int = 3.e-21        # Permeabilite intrinseque\n") ;
  fprintf(ficd,"mu_l = 1.e-3          # Viscosite du liquide\n") ;
  fprintf(ficd,"RT_0 = 2479           # Constante des gaz parfaits fois la temperature\n") ;
  fprintf(ficd,"K_eau = 1.e-8         # Constante d\'equilibre pour l\'eau\n") ;
  fprintf(ficd,"K_Fe(OH)2 = 3.1429e-2 # Constante d\'equilibre pour Fe(OH)2\n") ;
  fprintf(ficd,"S_Fe(OH)2 = 1.995e-5  # Solubilite de Fe(OH)2 solide\n") ;
  fprintf(ficd,"K_Far = 96485         # Constante de Faraday\n") ;
  fprintf(ficd,"K_Hen = 1.25e-5       # Constante de Henry pour O2\n") ;
  fprintf(ficd,"D_oh = 5.24e-9        # Diffusion moleculaire de OH-\n") ;
  fprintf(ficd,"D_h = 9.32e-9         # Diffusion moleculaire de H+\n") ;
  fprintf(ficd,"D_o2_l = 2.51e-9      # Diffusion moleculaire de O2 liquide\n") ;
  fprintf(ficd,"D_o2_g = 1.e-7        # Diffusion moleculaire de O2 gazeux\n") ;
  fprintf(ficd,"D_fe = 0.72e-9        # Diffusion moleculaire de Fe2+\n") ;
  fprintf(ficd,"D_feoh2 = 0.72e-9     # Diffusion moleculaire de Fe(OH)2\n") ;
  fprintf(ficd,"D_cat = 1.33e-9       # Diffusion moleculaire des cations\n") ;
  fprintf(ficd,"D_cat = 1.18e-9       # Diffusion moleculaire des anions\n") ;
  fprintf(ficd,"courbes = my_file     # Nom du fichier : p_c S_l k_rl tau_l/tau_l_sat tau_g/tau_g_sat\n") ;

return(NEQ) ;
}


void tb24(elem_t *el,inte_t *fi,unsigned int *n_fi,int dim)
{
  el->n_vi = NVI*el->fi->np ;
  el->n_ve = NVE*el->fi->np ;
}



void ch24(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t,char_t cg)
/* Residu du aux chargements (r) */
{
  int    i ;

  chsurf(x,r,dim,geom,dt,t,cg,el,el.fi) ;
  for(i=0;i<NEQ*el.nn;i++) r[i] = -r[i] ;
}


void in24(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom)
/* Initialise les variables du systeme (f,va) */ 
{
  double un = 1. ;
  double pl,sl,pc,sg ;
  double c_oh,c_h,c_o2_l,c_cat,c_ani,c_fe,c_feoh2,c_o2_g ;
  double dc_hsdc_fe,dc_hsdc_cat,dc_hsdc_ani ;
  double dc_ohsdc_fe,dc_ohsdc_cat,dc_ohsdc_ani ;
  double dc_feoh2sdc_fe,dc_feoh2sdc_cat,dc_feoh2sdc_ani ;
  double n_h,n_o,n_fe,n_o2,n_cat,n_ani,n_h2o ;
  double grd_p_l[3],grd_oh[3],grd_o2[3],grd_cat[3],grd_ani[3],grd_psi[3],grd_h[3],grd_fe[3],grd_feoh2[3] ;
  double w_h[3],w_fe[3],w_q[3],w_o2[3],w_cat[3],w_ani[3],w_h2o[3] ;
  double k_l,kd_h2o,kd_h,kd_oh,kd_fe,kd_feoh2,kd_o2,kd_cat,kd_ani ;
  double d_h,d_oh,d_fe,d_feoh2,d_o2,d_o2_l,d_cat,d_ani,d_o2_g,tau_l,tau_g ;
  double ke_h,ke_oh,ke_fe,ke_cat,ke_ani,ke_h2o ;
  int    i,p ;
  double *h,*dh ;

  if(el.dim < dim) return ;

  /* Donnees */
  phi      = el.mat->pr[pm("phi")] ;
  c_h2o    = el.mat->pr[pm("c_h2o")] ;
  k_int    = el.mat->pr[pm("k_int")] ;
  mu_l     = el.mat->pr[pm("mu_l")] ;
  k_eau    = el.mat->pr[pm("K_eau")] ;
  k_feoh2  = el.mat->pr[pm("K_Fe(OH)2")] ;
  farad    = el.mat->pr[pm("K_Far")] ;
  d0_oh    = el.mat->pr[pm("D_oh")] ;
  d0_h     = el.mat->pr[pm("D_h")] ;
  d0_o2_l  = el.mat->pr[pm("D_o2_l")] ;
  d0_fe    = el.mat->pr[pm("D_fe")] ;
  d0_feoh2 = el.mat->pr[pm("D_feoh2")] ;
  d0_cat   = el.mat->pr[pm("D_cat")] ;
  d0_ani   = el.mat->pr[pm("D_ani")] ;
  d0_o2_g  = el.mat->pr[pm("D_o2_g")] ;
  RT_0     = el.mat->pr[pm("RT_0")] ;
  k_hen    = el.mat->pr[pm("K_Hen")]*RT_0 ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    pl  = param(u,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    sg  = un - sl ; 
    /* concentrations */
    c_fe    = param(u,h,el.nn,I_c_fe) ;
    c_o2_l  = param(u,h,el.nn,I_c_o2) ;
    c_cat   = param(u,h,el.nn,I_c_cat) ;
    c_ani   = param(u,h,el.nn,I_c_ani) ;
    c_o2_g  = c_o2_l/k_hen ;
    c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
    c_h     = k_eau/c_oh ;
    c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;
    /* derivees */
    dc_ohsdc_fe     = 2*c_oh/(c_h + c_oh) ;
    dc_ohsdc_cat    =  c_oh/(c_h + c_oh) ;
    dc_ohsdc_ani    = -c_oh/(c_h + c_oh) ;
    
    dc_hsdc_fe      = -k_eau/(c_oh*c_oh)*dc_ohsdc_fe ;
    dc_hsdc_cat     = -k_eau/(c_oh*c_oh)*dc_ohsdc_cat ;
    dc_hsdc_ani     = -k_eau/(c_oh*c_oh)*dc_ohsdc_ani ;
    
    dc_feoh2sdc_fe  = (c_oh*c_oh + 2*c_fe*c_oh*dc_ohsdc_fe)/k_feoh2 ;
    dc_feoh2sdc_cat = 2*c_fe*c_oh*dc_ohsdc_cat/k_feoh2 ;
    dc_feoh2sdc_ani = 2*c_fe*c_oh*dc_ohsdc_ani/k_feoh2 ;
    /* contenus molaires */
    n_h2o = phi*sl*c_h2o ;
    n_h   = 2*n_h2o + phi*sl*(c_h + c_oh + 2*c_feoh2) ;
    n_o2  = phi*(sl*c_o2_l + sg*c_o2_g) ;
    n_o   = phi*sl*(0.5*c_oh - 0.5*c_h + c_feoh2) ;
    n_fe  = phi*sl*(c_fe + c_feoh2) ;
    n_cat = phi*sl*c_cat ;
    n_ani = phi*sl*c_ani ;
    /* tortuosites */
    tau_l = tortuosite_l(pc,el.mat) ;
    tau_g = tortuosite_g(pc,el.mat) ;
    /* coefficient de transfert */
    k_l      = k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    kd_h2o   = c_h2o*k_l ;
    kd_h     = c_h*k_l ;
    kd_oh    = c_oh*k_l ;
    kd_fe    = c_fe*k_l ;
    kd_feoh2 = c_feoh2*k_l ;
    kd_o2    = c_o2_l*k_l ;
    kd_cat   = c_cat*k_l ;
    kd_ani   = c_ani*k_l ;
    d_h      = phi*sl*tau_l*d0_h ;
    d_oh     = phi*sl*tau_l*d0_oh ;
    d_o2_l   = phi*sl*tau_l*d0_o2_l ;
    d_fe     = phi*sl*tau_l*d0_fe ;
    d_feoh2  = phi*sl*tau_l*d0_feoh2 ;
    d_cat    = phi*sl*tau_l*d0_cat ;
    d_ani    = phi*sl*tau_l*d0_ani ;
    d_o2_g   = phi*sg*tau_g*d0_o2_g ;
    d_o2     = d_o2_l + d_o2_g/k_hen ;
    ke_h     = c_h*d_h*farad/RT_0 ;
    ke_oh    = - c_oh*d_oh*farad/RT_0 ;
    ke_fe    = 2.*c_fe*d_fe*farad/RT_0 ;
    ke_cat   = c_cat*d_cat*farad/RT_0 ;
    ke_ani   = - c_ani*d_ani*farad/RT_0 ;
    ke_h2o   = (M_oh*ke_oh + M_h*ke_h + M_fe*ke_fe + M_cat*ke_cat + M_ani*ke_ani)/M_h2o ;
    /* flux */
    grad(x,u,dh,grd_p_l,el.nn,dim,I_p_l) ;
    grad(x,u,dh,grd_fe,el.nn,dim,I_c_fe) ;
    grad(x,u,dh,grd_o2,el.nn,dim,I_c_o2) ;
    grad(x,u,dh,grd_cat,el.nn,dim,I_c_cat) ;
    grad(x,u,dh,grd_ani,el.nn,dim,I_c_ani) ;
    grad(x,u,dh,grd_psi,el.nn,dim,I_psi) ;
    for(i=0;i<dim;i++){
      grd_h[i]     = dc_hsdc_fe*grd_fe[i]
	           + dc_hsdc_cat*grd_cat[i] 
	           + dc_hsdc_ani*grd_ani[i] ;
      grd_oh[i]    = dc_ohsdc_fe*grd_fe[i] 
	           + dc_ohsdc_cat*grd_cat[i] 
	           + dc_ohsdc_ani*grd_ani[i] ;
      grd_feoh2[i] = dc_feoh2sdc_fe*grd_fe[i] 
	           + dc_feoh2sdc_cat*grd_cat[i] 
                   + dc_feoh2sdc_ani*grd_ani[i] ;
    }
    /* H2O */
    for(i=0;i<3;i++) {
      w_h2o[i]  = - kd_h2o*grd_p_l[i]
                + M_o2/M_h2o*d_o2_l*grd_o2[i]
                + M_h/M_h2o*d_h*grd_h[i]
                + M_oh/M_h2o*d_oh*grd_oh[i]
	        + M_fe/M_h2o*d_fe*grd_fe[i]
                + M_feoh2/M_h2o*d_feoh2*grd_feoh2[i]
                + M_cat/M_h2o*d_cat*grd_cat[i]
	        + M_ani/M_h2o*d_ani*grd_ani[i]
	        + ke_h2o*grd_psi[i] ;
    }
    /* H = 2*H2O + H+ + OH- + 2*Fe(OH)2 */
    for(i=0;i<3;i++) {
      w_h[i]   = 2*w_h2o[i]
	       - (kd_h + kd_oh + 2*kd_feoh2)*grd_p_l[i]
	       - d_h*grd_h[i] - d_oh*grd_oh[i] - 2.*d_feoh2*grd_feoh2[i]
               - (ke_h + ke_oh)*grd_psi[i] ;
    }
    /* O2 */
    for(i=0;i<3;i++) {
      w_o2[i]  = - kd_o2*grd_p_l[i] - d_o2*grd_o2[i] ;
    }
    /* Fe = Fe2+ + Fe(OH)2 */
    for(i=0;i<3;i++) {
      w_fe[i]  = - (kd_fe + kd_feoh2)*grd_p_l[i]
               - d_fe*grd_fe[i] - d_feoh2*grd_feoh2[i]
	       - ke_fe*grd_psi[i] ;
    }
    /* Charge Q = H+ - OH- + 2Fe2+ + K+ - A- */
    for(i=0;i<3;i++) {
      w_q[i]   = - d_h*grd_h[i] + d_oh*grd_oh[i] - 2*d_fe*grd_fe[i] - d_cat*grd_cat[i] + d_ani*grd_ani[i]
               - (ke_h - ke_oh + 2*ke_fe + ke_cat - ke_ani)*grd_psi[i] ;
    }
    /* Cations */
    for(i=0;i<3;i++) {
      w_cat[i] = - kd_cat*grd_p_l[i] - d_cat*grd_cat[i] - ke_cat*grd_psi[i] ;
    }
    /* Anions */
    for(i=0;i<3;i++) {
      w_ani[i] = - kd_ani*grd_p_l[i] - d_ani*grd_ani[i] - ke_ani*grd_psi[i] ;
    }

    /* rangement dans f */
    f[p*NVI+E_h]   = n_h ;
    f[p*NVI+E_fe]  = n_fe ;
    f[p*NVI+E_i]   = 0. ;
    f[p*NVI+E_o2]  = n_o2 ;
    f[p*NVI+E_cat] = n_cat ;
    f[p*NVI+E_ani] = n_ani ;
    for(i=0;i<dim;i++) {
      f[p*NVI+NEQ+E_h*dim+i]    = w_h[i] ;
      f[p*NVI+NEQ+E_fe*dim+i]   = w_fe[i] ;
      f[p*NVI+NEQ+E_i*dim+i]    = w_q[i] ;
      f[p*NVI+NEQ+E_o2*dim+i]   = w_o2[i] ;
      f[p*NVI+NEQ+E_cat*dim+i]  = w_cat[i] ;
      f[p*NVI+NEQ+E_ani*dim+i]  = w_ani[i] ;
    }
    /* rangement dans va */
    va[p*NVE]    = kd_h2o ;
    va[p*NVE+1]  = kd_h ;
    va[p*NVE+2]  = kd_oh ;
    va[p*NVE+3]  = kd_fe ;
    va[p*NVE+4]  = kd_feoh2 ;
    va[p*NVE+5]  = kd_o2 ;
    va[p*NVE+6]  = kd_cat ;
    va[p*NVE+7]  = kd_ani ;
    va[p*NVE+8]  = d_oh ;
    va[p*NVE+9]  = d_h ;
    va[p*NVE+10] = d_fe ;
    va[p*NVE+11] = d_feoh2 ;
    va[p*NVE+12] = d_o2_l ;
    va[p*NVE+13] = d_o2_g ;
    va[p*NVE+14] = d_cat ;
    va[p*NVE+15] = d_ani ;
    va[p*NVE+16] = ke_h ;
    va[p*NVE+17] = ke_oh ;
    va[p*NVE+18] = ke_fe ;
    va[p*NVE+19] = ke_cat ;
    va[p*NVE+20] = ke_ani ;
    va[p*NVE+21] = dc_hsdc_fe ;
    va[p*NVE+22] = dc_hsdc_cat ;
    va[p*NVE+23] = dc_hsdc_ani ;
    va[p*NVE+24] = dc_ohsdc_fe ;
    va[p*NVE+25] = dc_ohsdc_cat ;
    va[p*NVE+26] = dc_ohsdc_ani ;
    va[p*NVE+27] = dc_feoh2sdc_fe ;
    va[p*NVE+28] = dc_feoh2sdc_cat ;
    va[p*NVE+29] = dc_feoh2sdc_ani ;
    va[p*NVE+30] = ke_h2o ;
  }
}

int ex24(double **x,double **u,double *f,double *va,elem_t el,int dim,geom_t geom,double t)
/* Termes explicites (va)  */
{
  int    p ;
  double un = 1. ;
  double pl,pc,sl,sg ;
  double c_oh,c_h,c_o2_l,c_cat,c_ani,c_fe,c_feoh2 ;
  double dc_hsdc_fe,dc_hsdc_cat,dc_hsdc_ani ;
  double dc_ohsdc_fe,dc_ohsdc_cat,dc_ohsdc_ani ;
  double dc_feoh2sdc_fe,dc_feoh2sdc_cat,dc_feoh2sdc_ani ;
  double k_l,kd_h2o,kd_h,kd_oh,kd_fe,kd_feoh2,kd_o2,kd_cat,kd_ani ;
  double d_h,d_oh,d_fe,d_feoh2,d_o2_l,d_o2_g,d_cat,d_ani,tau_l,tau_g ;
  double ke_h,ke_oh,ke_fe,ke_cat,ke_ani,ke_h2o ;
  double *h ;

  if(el.dim < dim) return(0) ;

  /* Donnees */
  phi      = el.mat->pr[pm("phi")] ;
  c_h2o    = el.mat->pr[pm("c_h2o")] ;
  k_int    = el.mat->pr[pm("k_int")] ;
  mu_l     = el.mat->pr[pm("mu_l")] ;
  k_eau    = el.mat->pr[pm("K_eau")] ;
  k_feoh2  = el.mat->pr[pm("K_Fe(OH)2")] ;
  farad    = el.mat->pr[pm("K_Far")] ;
  d0_oh    = el.mat->pr[pm("D_oh")] ;
  d0_h     = el.mat->pr[pm("D_h")] ;
  d0_o2_l  = el.mat->pr[pm("D_o2_l")] ;
  d0_fe    = el.mat->pr[pm("D_fe")] ;
  d0_feoh2 = el.mat->pr[pm("D_feoh2")] ;
  d0_cat   = el.mat->pr[pm("D_cat")] ;
  d0_ani   = el.mat->pr[pm("D_ani")] ;
  d0_o2_g  = el.mat->pr[pm("D_o2_g")] ;
  RT_0     = el.mat->pr[pm("RT_0")] ;
  k_hen    = el.mat->pr[pm("K_Hen")]*RT_0 ;
  
  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* pression */
    pl  = param(u,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    sg  = un - sl ;
    /* concentrations */
    c_fe    = param(u,h,el.nn,I_c_fe) ;
    c_o2_l  = param(u,h,el.nn,I_c_o2) ;
    c_cat   = param(u,h,el.nn,I_c_cat) ;
    c_ani   = param(u,h,el.nn,I_c_ani) ;
    c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
    c_h     = k_eau/c_oh ;
    c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;

    if(c_oh < 0. || c_fe < 0. || c_cat < 0. || c_ani < 0. || c_o2_l < 0.) {
      printf("\n\
      c_oh    = %e\n\
      c_fe    = %e\n\
      c_cat   = %e\n\
      c_ani   = %e\n\
      c_o2    = %e\n",c_oh,c_fe,c_cat,c_ani,c_o2_l) ;
      return(1) ;
    }
    /* derivees */
    dc_ohsdc_fe     = 2*c_oh/(c_h + c_oh) ;
    dc_ohsdc_cat    =  c_oh/(c_h + c_oh) ;
    dc_ohsdc_ani    = -c_oh/(c_h + c_oh) ;
    
    dc_hsdc_fe      = -k_eau/(c_oh*c_oh)*dc_ohsdc_fe ;
    dc_hsdc_cat     = -k_eau/(c_oh*c_oh)*dc_ohsdc_cat ;
    dc_hsdc_ani     = -k_eau/(c_oh*c_oh)*dc_ohsdc_ani ;
    
    dc_feoh2sdc_fe  = (c_oh*c_oh + 2*c_fe*c_oh*dc_ohsdc_fe)/k_feoh2 ;
    dc_feoh2sdc_cat = 2*c_fe*c_oh*dc_ohsdc_cat/k_feoh2 ;
    dc_feoh2sdc_ani = 2*c_fe*c_oh*dc_ohsdc_ani/k_feoh2 ;
    /* tortuosites */
    tau_l = tortuosite_l(pc,el.mat) ;
    tau_g = tortuosite_g(pc,el.mat) ;
    /* coefficients de transfert */
    k_l      = k_int/mu_l*courbe(pc,el.mat->cb[1]) ;
    kd_h2o   = c_h2o*k_l ;
    kd_h     = c_h*k_l ;
    kd_oh    = c_oh*k_l ;
    kd_fe    = c_fe*k_l ;
    kd_feoh2 = c_feoh2*k_l ;
    kd_o2    = c_o2_l*k_l ;
    kd_cat   = c_cat*k_l ;
    kd_ani   = c_ani*k_l ;
    d_h      = phi*sl*tau_l*d0_h ;
    d_oh     = phi*sl*tau_l*d0_oh ;
    d_o2_l   = phi*sl*tau_l*d0_o2_l ;
    d_fe     = phi*sl*tau_l*d0_fe ;
    d_feoh2  = phi*sl*tau_l*d0_feoh2 ;
    d_cat    = phi*sl*tau_l*d0_cat ;
    d_ani    = phi*sl*tau_l*d0_ani ;
    d_o2_g   = phi*sg*tau_g*d0_o2_g ;
    ke_h     = c_h*d_h*farad/RT_0 ;
    ke_oh    = - c_oh*d_oh*farad/RT_0 ;
    ke_fe    = 2*c_fe*d_fe*farad/RT_0 ;
    ke_cat   = c_cat*d_cat*farad/RT_0 ;
    ke_ani   = - c_ani*d_ani*farad/RT_0 ;
    ke_h2o   = (M_oh*ke_oh + M_h*ke_h + M_fe*ke_fe + M_cat*ke_cat + M_ani*ke_ani)/M_h2o ;
    /* rangement dans va */
    va[p*NVE]    = kd_h2o ;
    va[p*NVE+1]  = kd_h ;
    va[p*NVE+2]  = kd_oh ;
    va[p*NVE+3]  = kd_fe ;
    va[p*NVE+4]  = kd_feoh2 ;
    va[p*NVE+5]  = kd_o2 ;
    va[p*NVE+6]  = kd_cat ;
    va[p*NVE+7]  = kd_ani ;
    va[p*NVE+8]  = d_oh ;
    va[p*NVE+9]  = d_h ;
    va[p*NVE+10] = d_fe ;
    va[p*NVE+11] = d_feoh2 ;
    va[p*NVE+12] = d_o2_l ;
    va[p*NVE+13] = d_o2_g ;
    va[p*NVE+14] = d_cat ;
    va[p*NVE+15] = d_ani ;
    va[p*NVE+16] = ke_h ;
    va[p*NVE+17] = ke_oh ;
    va[p*NVE+18] = ke_fe ;
    va[p*NVE+19] = ke_cat ;
    va[p*NVE+20] = ke_ani ;
    va[p*NVE+21] = dc_hsdc_fe ;
    va[p*NVE+22] = dc_hsdc_cat ;
    va[p*NVE+23] = dc_hsdc_ani ;
    va[p*NVE+24] = dc_ohsdc_fe ;
    va[p*NVE+25] = dc_ohsdc_cat ;
    va[p*NVE+26] = dc_ohsdc_ani ;
    va[p*NVE+27] = dc_feoh2sdc_fe ;
    va[p*NVE+28] = dc_feoh2sdc_cat ;
    va[p*NVE+29] = dc_feoh2sdc_ani ;
    va[p*NVE+30] = ke_h2o ;
  }
  return(0) ;
}

int ct24(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,int dim,geom_t geom,double dt,double t)
/* Les variables donnees par la loi de comportement (f_1) */
{
  double un = 1. ;
  double pl,sl,pc,sg ;
  double c_oh,c_h,c_o2_l,c_cat,c_ani,c_fe,c_feoh2,c_o2_g ;
  double dc_hsdc_fe,dc_hsdc_cat,dc_hsdc_ani ;
  double dc_ohsdc_fe,dc_ohsdc_cat,dc_ohsdc_ani ;
  double dc_feoh2sdc_fe,dc_feoh2sdc_cat,dc_feoh2sdc_ani ;
  double n_h,n_o,n_fe,n_o2,n_cat,n_ani,n_h2o ;
  double grd_p_l[3],grd_oh[3],grd_o2[3],grd_cat[3],grd_ani[3],grd_psi[3],grd_h[3],grd_fe[3],grd_feoh2[3] ;
  double w_h[3],w_fe[3],w_q[3],w_o2[3],w_cat[3],w_ani[3],w_h2o[3] ;
  double kd_h2o,kd_h,kd_oh,kd_fe,kd_feoh2,kd_o2,kd_cat,kd_ani ;
  double d_h,d_oh,d_fe,d_feoh2,d_o2,d_o2_l,d_o2_g,d_cat,d_ani ;
  double ke_h,ke_oh,ke_fe,ke_cat,ke_ani,ke_h2o ;
  int    i,p ;
  double *h,*dh ;

  if(el.dim < dim) return (0) ;

  /* Donnees */
  phi      = el.mat->pr[pm("phi")] ;
  c_h2o    = el.mat->pr[pm("c_h2o")] ;
  RT_0     = el.mat->pr[pm("RT_0")] ;
  k_eau    = el.mat->pr[pm("K_eau")] ;
  k_feoh2  = el.mat->pr[pm("K_Fe(OH)2")] ;
  k_hen    = el.mat->pr[pm("K_Hen")]*RT_0 ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    pl  = param(u_1,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    sg  = un - sl ;
    /* concentrations */
    c_fe    = param(u_1,h,el.nn,I_c_fe) ;
    c_o2_l  = param(u_1,h,el.nn,I_c_o2) ;
    c_cat   = param(u_1,h,el.nn,I_c_cat) ;
    c_ani   = param(u_1,h,el.nn,I_c_ani) ;
    c_o2_g  = c_o2_l/k_hen ;
    c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
    c_h     = k_eau/c_oh ;
    c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;
    /* contenus molaires */
    n_h2o = phi*sl*c_h2o ;
    n_h   = 2*n_h2o + phi*sl*(c_h + c_oh + 2*c_feoh2) ;
    n_o2  = phi*(sl*c_o2_l + sg*c_o2_g) ;
    n_o   = phi*sl*(0.5*c_oh - 0.5*c_h + c_feoh2) ;
    n_fe  = phi*sl*(c_fe + c_feoh2) ;
    n_cat = phi*sl*c_cat ;
    n_ani = phi*sl*c_ani ;
    /* coefficient de transfert */
    kd_h2o   = va[p*NVE] ;
    kd_h     = va[p*NVE+1] ;
    kd_oh    = va[p*NVE+2] ;
    kd_fe    = va[p*NVE+3] ;
    kd_feoh2 = va[p*NVE+4] ;
    kd_o2    = va[p*NVE+5] ;
    kd_cat   = va[p*NVE+6] ;
    kd_ani   = va[p*NVE+7] ;
    d_oh     = va[p*NVE+8] ;
    d_h      = va[p*NVE+9] ;
    d_fe     = va[p*NVE+10] ;
    d_feoh2  = va[p*NVE+11] ;
    d_o2_l   = va[p*NVE+12] ;
    d_o2_g   = va[p*NVE+13] ;
    d_cat    = va[p*NVE+14] ;
    d_ani    = va[p*NVE+15] ;
    ke_h     = va[p*NVE+16] ;
    ke_oh    = va[p*NVE+17] ;
    ke_fe    = va[p*NVE+18] ;
    ke_cat   = va[p*NVE+19] ;
    ke_ani   = va[p*NVE+20] ;
    ke_h2o   = va[p*NVE+30] ;
    d_o2     = d_o2_l + d_o2_g/k_hen ;
    /* derivees */
    dc_hsdc_fe      = va[p*NVE+21] ;
    dc_hsdc_cat     = va[p*NVE+22] ;
    dc_hsdc_ani     = va[p*NVE+23] ;
    dc_ohsdc_fe     = va[p*NVE+24] ;
    dc_ohsdc_cat    = va[p*NVE+25] ;
    dc_ohsdc_ani    = va[p*NVE+26] ;
    dc_feoh2sdc_fe  = va[p*NVE+27] ;
    dc_feoh2sdc_cat = va[p*NVE+28] ;
    dc_feoh2sdc_ani = va[p*NVE+29] ;
    /* flux */
    grad(x,u_1,dh,grd_p_l,el.nn,dim,I_p_l) ;
    grad(x,u_1,dh,grd_fe,el.nn,dim,I_c_fe) ;
    grad(x,u_1,dh,grd_o2,el.nn,dim,I_c_o2) ;
    grad(x,u_1,dh,grd_cat,el.nn,dim,I_c_cat) ;
    grad(x,u_1,dh,grd_ani,el.nn,dim,I_c_ani) ;
    grad(x,u_1,dh,grd_psi,el.nn,dim,I_psi) ;
    for(i=0;i<dim;i++){
      grd_h[i]     = dc_hsdc_fe*grd_fe[i]
	           + dc_hsdc_cat*grd_cat[i] 
	           + dc_hsdc_ani*grd_ani[i] ;
      grd_oh[i]    = dc_ohsdc_fe*grd_fe[i] 
	           + dc_ohsdc_cat*grd_cat[i] 
	           + dc_ohsdc_ani*grd_ani[i] ;
      grd_feoh2[i] = dc_feoh2sdc_fe*grd_fe[i] 
	           + dc_feoh2sdc_cat*grd_cat[i] 
                   + dc_feoh2sdc_ani*grd_ani[i] ;
    }
    /* H2O */
    for(i=0;i<3;i++) {
      w_h2o[i]  = - kd_h2o*grd_p_l[i]
                + M_o2/M_h2o*d_o2_l*grd_o2[i]
                + M_h/M_h2o*d_h*grd_h[i]
                + M_oh/M_h2o*d_oh*grd_oh[i]
	        + M_fe/M_h2o*d_fe*grd_fe[i]
                + M_feoh2/M_h2o*d_feoh2*grd_feoh2[i]
                + M_cat/M_h2o*d_cat*grd_cat[i]
	        + M_ani/M_h2o*d_ani*grd_ani[i]
	        + ke_h2o*grd_psi[i] ;
    }
    /* H = 2*H2O + H+ + OH- + 2*Fe(OH)2 */
    for(i=0;i<3;i++) {
      w_h[i]   = 2*w_h2o[i]
	       - (kd_h + kd_oh + 2*kd_feoh2)*grd_p_l[i]
	       - d_h*grd_h[i] - d_oh*grd_oh[i] - 2*d_feoh2*grd_feoh2[i]
               - (ke_h + ke_oh)*grd_psi[i] ;
    }
    /* O2 */
    for(i=0;i<3;i++) {
      w_o2[i]  = - kd_o2*grd_p_l[i] - d_o2*grd_o2[i] ;
    }
    /* Fe = Fe2+ + Fe(OH)2 */
    for(i=0;i<3;i++) {
      w_fe[i]  = - (kd_fe + kd_feoh2)*grd_p_l[i]
               - d_fe*grd_fe[i] - d_feoh2*grd_feoh2[i]
	       - ke_fe*grd_psi[i] ;
    }
    /* Charge Q = H+ - OH- + 2Fe2+ + K+ - A- */
    for(i=0;i<3;i++) {
      w_q[i]   = - d_h*grd_h[i] + d_oh*grd_oh[i] - 2*d_fe*grd_fe[i] - d_cat*grd_cat[i] + d_ani*grd_ani[i]
               - (ke_h - ke_oh + 2*ke_fe + ke_cat - ke_ani)*grd_psi[i] ;
    }
    /* Cations */
    for(i=0;i<3;i++) {
      w_cat[i] = - kd_cat*grd_p_l[i] - d_cat*grd_cat[i] - ke_cat*grd_psi[i] ;
    }
    /* Anions */
    for(i=0;i<3;i++) {
      w_ani[i] = - kd_ani*grd_p_l[i] - d_ani*grd_ani[i] - ke_ani*grd_psi[i] ;
    }

    /* rangement dans f_1 */
    f_1[p*NVI+E_h]   = n_h ;
    f_1[p*NVI+E_fe]  = n_fe ;
    f_1[p*NVI+E_i]   = 0. ;
    f_1[p*NVI+E_o2]  = n_o2 ;
    f_1[p*NVI+E_cat] = n_cat ;
    f_1[p*NVI+E_ani] = n_ani ;
    for(i=0;i<dim;i++) {
      f_1[p*NVI+NEQ+E_h*dim+i]   = w_h[i] ;
      f_1[p*NVI+NEQ+E_fe*dim+i]  = w_fe[i] ;
      f_1[p*NVI+NEQ+E_i*dim+i]   = w_q[i] ;
      f_1[p*NVI+NEQ+E_o2*dim+i]  = w_o2[i] ;
      f_1[p*NVI+NEQ+E_cat*dim+i] = w_cat[i] ;
      f_1[p*NVI+NEQ+E_ani*dim+i] = w_ani[i] ;
    }
  }
  return(0) ;
}

int mx24(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *k,elem_t el,int dim,geom_t geom,double dt,double t)
/* Matrice (k) */
{
#define NN        (el.nn)
  int    i,dec ;
  double c[MAX_PGAUSS*9*NEQ*NEQ] ;
  double kb[NEQ*NEQ*MAX_NOEUDS*MAX_NOEUDS] ;
  double zero = 0. ;
  
  /* initialisation */
  for(i=0;i<NN*NN*NEQ*NEQ;i++) k[i] = zero ;

  if(el.dim < dim) return(0) ;

  /*
  ** Matrice de comportement
  */
  dec = c24(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxcmss(k,x,*el.fi,c,dim,dec,geom,NEQ) ;
  /*
  ** Matrice de conduction
  */
  dec = k24(x,u_1,u_n,f_1,f_n,va,el,geom,dt,c) ;
  mxccnd(kb,x,*el.fi,c,dim,dec,geom,NEQ) ;
  for(i=0;i<NN*NN*NEQ*NEQ;i++) k[i] += dt*kb[i] ;
  return(0) ;
#undef NN
}


void rs24(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,double *r,elem_t el,int dim,geom_t geom,double dt,double t)
/* Residu (r) */
{
  int    i ;
  double rb[MAX_NOEUDS],g1[MAX_PGAUSS],*h1 ;

  /* initialisation */
  for(i=0;i<el.nn*NEQ;i++) r[i] = 0. ;

  if(el.dim < dim) return ;
  
  /* hydrogene */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_h] - f_n[i*NVI+E_h] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_h] = - rb[i] ;
  h1 = f_1 + NEQ + E_h*dim ;
  rsflux(rb,x,*el.fi,h1,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_h] += dt*rb[i] ;
  
  /* fer */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_fe] - f_n[i*NVI+E_fe] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fe] = - rb[i] ;
  h1 = f_1 + NEQ + E_fe*dim ;
  rsflux(rb,x,*el.fi,h1,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_fe] += dt*rb[i] ;
  
  /* courant : div i = 0 */
  h1 = f_1 + NEQ + E_i*dim ;
  rsflux(rb,x,*el.fi,h1,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_i] += dt*rb[i] ;
  
  /* o2 */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_o2] - f_n[i*NVI+E_o2] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_o2] = - rb[i] ;
  h1 = f_1 + NEQ + E_o2*dim ;
  rsflux(rb,x,*el.fi,h1,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_o2] += dt*rb[i] ;
  
  /* cations */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_cat] - f_n[i*NVI+E_cat] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_cat] = - rb[i] ;
  h1 = f_1 + NEQ + E_cat*dim ;
  rsflux(rb,x,*el.fi,h1,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_cat] += dt*rb[i] ;
  
  /* anions */
  for(i=0;i<el.fi->np;i++) g1[i] = f_1[i*NVI+E_ani] - f_n[i*NVI+E_ani] ;
  rsmass(rb,x,*el.fi,g1,dim,1,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_ani] = - rb[i] ;
  h1 = f_1 + NEQ + E_ani*dim ;
  rsflux(rb,x,*el.fi,h1,dim,NVI,geom) ;
  for(i=0;i<el.nn;i++) r[i*NEQ+E_ani] += dt*rb[i] ;
}


int so24(double **x,double **u,double *f,double *va,double *s,resu_t *r,elem_t el,int dim,geom_t geom,double t)
/* Les valeurs exploitees (s) */
{
  double pl,sl,pc ;
  double c_oh,c_h,c_o2_l,c_cat,c_ani,c_fe,c_feoh2,psi ;
  double dc_hsdc_fe,dc_hsdc_cat,dc_hsdc_ani ;
  double dc_ohsdc_fe,dc_ohsdc_cat,dc_ohsdc_ani ;
  double dc_feoh2sdc_fe,dc_feoh2sdc_cat,dc_feoh2sdc_ani ;
  double grd_p_l[3],grd_oh[3],grd_o2[3],grd_cat[3],grd_ani[3],grd_psi[3],grd_h[3],grd_fe[3],grd_feoh2[3] ;
  double w_h2o[3],w_h_p[3],w_oh[3],w_fe_2p[3],w_feoh2[3],w_o2[3],w_cat[3],w_ani[3],i_ionic[3],i_ohmic[3] ;
  double kd_h2o,kd_h,kd_oh,kd_fe,kd_feoh2,kd_o2,kd_cat,kd_ani ;
  double d_h,d_oh,d_fe,d_feoh2,d_o2,d_o2_l,d_o2_g,d_cat,d_ani ;
  double ke_h,ke_oh,ke_fe,ke_cat,ke_ani,ke_h2o,r_elec ;
  int    i,j,p,nso ;
  double *h,*dh,h_s[MAX_NOEUDS],dh_s[3*MAX_NOEUDS] ;
  double zero = 0. ;

  if(el.dim < dim) return(0) ;

  /*
    Donnees
  */
  phi     = el.mat->pr[pm("phi")] ;
  c_h2o   = el.mat->pr[pm("c_h2o")] ;
  k_int   = el.mat->pr[pm("k_int")] ;
  mu_l    = el.mat->pr[pm("mu_l")] ;
  k_eau   = el.mat->pr[pm("K_eau")] ;
  k_feoh2 = el.mat->pr[pm("K_Fe(OH)2")] ;
  RT_0    = el.mat->pr[pm("RT_0")] ;
  farad   = el.mat->pr[pm("K_Far")] ;
  k_hen   = el.mat->pr[pm("K_Hen")]*RT_0 ;
  
  /* initialisation */
  nso = 21 ;
  for(i=0;i<nso;i++) for(j=0;j<9;j++) r[i].v[j] = zero ;

  /* fonctions d'interpolation en s */
  fint_abs(dim,el.nn,x,s,el.dim,h_s,dh_s) ;

  /* pression */
  pl = param(u,h_s,el.nn,I_p_l) ;
  pc = p_g - pl ;
  /* saturation */
  sl  = courbe(pc,el.mat->cb[0]) ;
  /* concentrations */
  c_fe    = param(u,h_s,el.nn,I_c_fe) ;
  c_o2_l  = param(u,h_s,el.nn,I_c_o2) ;
  c_cat   = param(u,h_s,el.nn,I_c_cat) ;
  c_ani   = param(u,h_s,el.nn,I_c_ani) ;
  c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
  c_h     = k_eau/c_oh ;
  c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;
  /* potentiel electrique */
  psi     = param(u,h_s,el.nn,I_psi) ;

  /* quantites exploitees */
  strcpy(r[0].text,"pression-liquide") ; r[0].n = 1 ;
  r[0].v[0] = pl ;
  strcpy(r[1].text,"saturation") ; r[1].n = 1 ;
  r[1].v[0] = sl ;
  strcpy(r[2].text,"molarite_O2") ; r[2].n = 1 ;
  r[2].v[0] = c_o2_l ;
  strcpy(r[3].text,"molarite_H+") ; r[3].n = 1 ;
  r[3].v[0] = c_h ;
  strcpy(r[4].text,"molarite_OH-") ; r[4].n = 1 ;
  r[4].v[0] = c_oh ;
  strcpy(r[5].text,"molarite_Fe2+") ; r[5].n = 1 ;
  r[5].v[0] = c_fe ;
  strcpy(r[6].text,"molarite_FeOH2") ; r[6].n = 1 ;
  r[6].v[0] = c_feoh2 ;
  strcpy(r[7].text,"molarite_cations") ; r[7].n = 1 ;
  r[7].v[0] = c_cat ;
  strcpy(r[8].text,"molarite_anions") ; r[8].n = 1 ;
  r[8].v[0] = c_ani ;
  strcpy(r[9].text,"potentiel_electrique") ; r[9].n = 1 ;
  r[9].v[0] = psi ;

  /* boucle sur les points d'integration */
  for(p=0;p<el.fi->np;p++) {
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    dh = el.fi->dh + p*dim*el.nn ;
    /* pression */
    pl  = param(u,h,el.nn,I_p_l) ;
    pc  = p_g - pl ;
    /* saturation */
    sl  = courbe(pc,el.mat->cb[0]) ;
    /* concentrations */
    c_fe    = param(u,h,el.nn,I_c_fe) ;
    c_o2_l  = param(u,h,el.nn,I_c_o2) ;
    c_cat   = param(u,h,el.nn,I_c_cat) ;
    c_ani   = param(u,h,el.nn,I_c_ani) ;
    c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
    c_h     = k_eau/c_oh ;
    c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;
    /* coefficient de transfert */
    kd_h2o   = va[p*NVE] ;
    kd_h     = va[p*NVE+1] ;
    kd_oh    = va[p*NVE+2] ;
    kd_fe    = va[p*NVE+3] ;
    kd_feoh2 = va[p*NVE+4] ;
    kd_o2    = va[p*NVE+5] ;
    kd_cat   = va[p*NVE+6] ;
    kd_ani   = va[p*NVE+7] ;
    d_oh     = va[p*NVE+8] ;
    d_h      = va[p*NVE+9] ;
    d_fe     = va[p*NVE+10] ;
    d_feoh2  = va[p*NVE+11] ;
    d_o2_l   = va[p*NVE+12] ;
    d_o2_g   = va[p*NVE+13] ;
    d_cat    = va[p*NVE+14] ;
    d_ani    = va[p*NVE+15] ;
    ke_h     = va[p*NVE+16] ;
    ke_oh    = va[p*NVE+17] ;
    ke_fe    = va[p*NVE+18] ;
    ke_cat   = va[p*NVE+19] ;
    ke_ani   = va[p*NVE+20] ;
    ke_h2o   = va[p*NVE+30] ;
    d_o2     = d_o2_l + d_o2_g/k_hen ;
    r_elec   = 1./farad/(ke_h + 4.* ke_fe + ke_cat - ke_oh - ke_ani) ;
    /* derivees */
    dc_hsdc_fe      = va[p*NVE+21] ;
    dc_hsdc_cat     = va[p*NVE+22] ;
    dc_hsdc_ani     = va[p*NVE+23] ;
    dc_ohsdc_fe     = va[p*NVE+24] ;
    dc_ohsdc_cat    = va[p*NVE+25] ;
    dc_ohsdc_ani    = va[p*NVE+26] ;
    dc_feoh2sdc_fe  = va[p*NVE+27] ;
    dc_feoh2sdc_cat = va[p*NVE+28] ;
    dc_feoh2sdc_ani = va[p*NVE+29] ;
    /* flux */
    grad(x,u,dh,grd_p_l,el.nn,dim,I_p_l) ;
    grad(x,u,dh,grd_fe,el.nn,dim,I_c_fe) ;
    grad(x,u,dh,grd_o2,el.nn,dim,I_c_o2) ;
    grad(x,u,dh,grd_cat,el.nn,dim,I_c_cat) ;
    grad(x,u,dh,grd_ani,el.nn,dim,I_c_ani) ;
    grad(x,u,dh,grd_psi,el.nn,dim,I_psi) ;
    for(i=0;i<dim;i++){
      grd_h[i]     = dc_hsdc_fe*grd_fe[i]
	           + dc_hsdc_cat*grd_cat[i] 
	           + dc_hsdc_ani*grd_ani[i] ;
      grd_oh[i]    = dc_ohsdc_fe*grd_fe[i] 
	           + dc_ohsdc_cat*grd_cat[i] 
	           + dc_ohsdc_ani*grd_ani[i] ;
      grd_feoh2[i] = dc_feoh2sdc_fe*grd_fe[i] 
	           + dc_feoh2sdc_cat*grd_cat[i] 
                   + dc_feoh2sdc_ani*grd_ani[i] ;
    }
    for(i=0;i<3;i++) {
      w_h2o[i]   = - kd_h2o*grd_p_l[i]
                 + M_o2/M_h2o*d_o2_l*grd_o2[i]
                 + M_h/M_h2o*d_h*grd_h[i]
                 + M_oh/M_h2o*d_oh*grd_oh[i]
	         + M_fe/M_h2o*d_fe*grd_fe[i]
                 + M_feoh2/M_h2o*d_feoh2*grd_feoh2[i]
                 + M_cat/M_h2o*d_cat*grd_cat[i]
	         + M_ani/M_h2o*d_ani*grd_ani[i]
	         + ke_h2o*grd_psi[i] ;
      w_h_p[i]   = - kd_h*grd_p_l[i] - d_h*grd_h[i] - ke_h*grd_psi[i] ;
      w_oh[i]    = - kd_oh*grd_p_l[i] - d_oh*grd_oh[i] - ke_oh*grd_psi[i] ;
      w_fe_2p[i] = - kd_fe*grd_p_l[i] - d_fe*grd_fe[i] - ke_fe*grd_psi[i] ;
      w_feoh2[i] = - kd_feoh2*grd_p_l[i] - d_feoh2*grd_feoh2[i] ;
      w_o2[i]    = - kd_o2*grd_p_l[i] - d_o2*grd_o2[i] ;
      w_cat[i]   = - kd_cat*grd_p_l[i] - d_cat*grd_cat[i] - ke_cat*grd_psi[i] ;
      w_ani[i]   = - kd_ani*grd_p_l[i] - d_ani*grd_ani[i] - ke_ani*grd_psi[i] ;
      i_ionic[i] = farad*(2*w_fe_2p[i] + w_h_p[i] + w_cat[i] - w_oh[i] - w_ani[i]) ;
      i_ohmic[i] = - grd_psi[i]/r_elec ;
    }
    /* quantites exploitees par element */
    strcpy(r[10].text,"flux_H2O") ; r[10].n = 3 ;
    for(i=0;i<dim;i++) r[10].v[i] += w_h2o[i]/el.fi->np ;
    strcpy(r[11].text,"flux_H+") ; r[11].n = 3 ;
    for(i=0;i<dim;i++) r[11].v[i] += w_h_p[i]/el.fi->np ;
    strcpy(r[12].text,"flux_OH-") ; r[12].n = 3 ;
    for(i=0;i<dim;i++) r[12].v[i] += w_oh[i]/el.fi->np ;
    strcpy(r[13].text,"flux_Fe2+") ; r[13].n = 3 ;
    for(i=0;i<dim;i++) r[13].v[i] += w_fe_2p[i]/el.fi->np ;
    strcpy(r[14].text,"flux_FeOH2") ; r[14].n = 3 ;
    for(i=0;i<dim;i++) r[14].v[i] += w_feoh2[i]/el.fi->np ;
    strcpy(r[15].text,"flux_O2") ; r[15].n = 3 ;
    for(i=0;i<dim;i++) r[15].v[i] += w_o2[i]/el.fi->np ;
    strcpy(r[16].text,"flux_cations") ; r[16].n = 3 ;
    for(i=0;i<dim;i++) r[16].v[i] += w_cat[i]/el.fi->np ;
    strcpy(r[17].text,"flux_anions") ; r[17].n = 3 ;
    for(i=0;i<dim;i++) r[17].v[i] += w_ani[i]/el.fi->np ;
    strcpy(r[18].text,"courant_ionique") ; r[18].n = 3 ;
    for(i=0;i<dim;i++) r[18].v[i] += i_ionic[i]/el.fi->np ;
    strcpy(r[19].text,"courant_ohmique") ; r[19].n = 3 ;
    for(i=0;i<dim;i++) r[19].v[i] += i_ohmic[i]/el.fi->np ;
    strcpy(r[20].text,"resistivite_electrique") ; r[20].n = 1 ;
    r[20].v[0] = r_elec ;
  }
  return (nso) ;
}


int c24(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de comportement (c) et decalage (dec)
*/
{
  double pl,pc,sl,sg,dslsdpc ;
  double c_oh,c_h,c_o2_l,c_cat,c_ani,c_fe,c_feoh2,c_o2_g ;
  double dc_hsdc_fe,dc_hsdc_cat,dc_hsdc_ani ;
  double dc_ohsdc_fe,dc_ohsdc_cat,dc_ohsdc_ani ;
  double dc_feoh2sdc_fe,dc_feoh2sdc_cat,dc_feoh2sdc_ani ;
  int    i,dec,p ;
  double *h,*c1 ;
  double zero = 0.,un = 1. ;
 
 /* Donnees */
  phi      = el.mat->pr[pm("phi")] ;
  c_h2o    = el.mat->pr[pm("c_h2o")] ;
  k_eau    = el.mat->pr[pm("K_eau")] ;
  k_feoh2  = el.mat->pr[pm("K_Fe(OH)2")] ;
  RT_0     = el.mat->pr[pm("RT_0")] ;
  k_hen    = el.mat->pr[pm("K_Hen")]*RT_0 ;
  
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
    sg  = un - sl ;
    dslsdpc = dcourbe(pc,el.mat->cb[0]) ;
    /* concentrations */
    c_fe    = param(u_1,h,el.nn,I_c_fe) ;
    c_o2_l  = param(u_1,h,el.nn,I_c_o2) ;
    c_cat   = param(u_1,h,el.nn,I_c_cat) ;
    c_ani   = param(u_1,h,el.nn,I_c_ani) ;
    c_o2_g  = c_o2_l/k_hen ;
    c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
    c_h     = k_eau/c_oh ;
    c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;
    /* derivees */
    dc_ohsdc_fe     = 2*c_oh/(c_h + c_oh) ;
    dc_ohsdc_cat    =  c_oh/(c_h + c_oh) ;
    dc_ohsdc_ani    = -c_oh/(c_h + c_oh) ;
    
    dc_hsdc_fe      = -k_eau/(c_oh*c_oh)*dc_ohsdc_fe ;
    dc_hsdc_cat     = -k_eau/(c_oh*c_oh)*dc_ohsdc_cat ;
    dc_hsdc_ani     = -k_eau/(c_oh*c_oh)*dc_ohsdc_ani ;
    
    dc_feoh2sdc_fe  = (c_oh*c_oh + 2*c_fe*c_oh*dc_ohsdc_fe)/k_feoh2 ;
    dc_feoh2sdc_cat = 2*c_fe*c_oh*dc_ohsdc_cat/k_feoh2 ;
    dc_feoh2sdc_ani = 2*c_fe*c_oh*dc_ohsdc_ani/k_feoh2 ;

    /* hydrogene */
    c1[E_h*NEQ+I_p_l]   = phi*(-dslsdpc)*(2*c_h2o + c_h + c_oh + 2*c_feoh2) ;
    c1[E_h*NEQ+I_c_fe]  = phi*sl*(dc_hsdc_fe + dc_ohsdc_fe + 2*dc_feoh2sdc_fe) ;
    c1[E_h*NEQ+I_c_cat] = phi*sl*(dc_hsdc_cat + dc_ohsdc_cat + 2*dc_feoh2sdc_cat) ;
    c1[E_h*NEQ+I_c_ani] = phi*sl*(dc_hsdc_ani + dc_ohsdc_ani + 2*dc_feoh2sdc_ani) ;

    /* fe */
    c1[E_fe*NEQ+I_p_l]   = phi*(-dslsdpc)*(c_fe + c_feoh2) ;
    c1[E_fe*NEQ+I_c_fe]  = phi*sl*(1. + dc_feoh2sdc_fe) ;
    c1[E_fe*NEQ+I_c_cat] = phi*sl*(dc_feoh2sdc_cat) ;
    c1[E_fe*NEQ+I_c_ani] = phi*sl*(dc_feoh2sdc_ani) ;

    /* charge Q = 0 */

    /* o2 */
    c1[E_o2*NEQ+I_p_l]   = phi*(-dslsdpc)*(c_o2_l - c_o2_g) ;
    c1[E_o2*NEQ+I_c_o2]  = phi*(sl + sg/k_hen) ;

    /* cations */
    c1[E_cat*NEQ+I_p_l]   = phi*(-dslsdpc)*c_cat ;
    c1[E_cat*NEQ+I_c_cat] = phi*sl ;

    /* anions */
    c1[E_ani*NEQ+I_p_l]   = phi*(-dslsdpc)*c_ani ;
    c1[E_ani*NEQ+I_c_ani] = phi*sl ;
  }
  return(dec) ;
}


int k24(double **x,double **u_1,double **u_n,double *f_1,double *f_n,double *va,elem_t el,geom_t geom,double dt,double *c)
/*
**  Matrice de conduction (c) et decalage (dec)
*/
{
  int    dec ;
  double *h,*c1,c2 ;
  double c_oh,c_cat,c_ani,c_h,c_fe,c_feoh2 ;
  double dc_hsdc_fe,dc_hsdc_cat,dc_hsdc_ani ;
  double dc_ohsdc_fe,dc_ohsdc_cat,dc_ohsdc_ani ;
  double dc_feoh2sdc_fe,dc_feoh2sdc_cat,dc_feoh2sdc_ani ;
  double kd_h2o,kd_h,kd_oh,kd_fe,kd_feoh2,kd_o2,kd_cat,kd_ani ;
  double d_h,d_oh,d_fe,d_feoh2,d_o2,d_o2_l,d_o2_g,d_cat,d_ani ;
  double ke_h,ke_oh,ke_fe,ke_cat,ke_ani,ke_h2o ;
  int    i,p ;
  double zero = 0. ;

 /* Donnees */
  k_eau    = el.mat->pr[pm("K_eau")] ;
  k_feoh2  = el.mat->pr[pm("K_Fe(OH)2")] ;
  RT_0     = el.mat->pr[pm("RT_0")] ;
  k_hen    = el.mat->pr[pm("K_Hen")]*RT_0 ;  

  dec = 9*NEQ*NEQ ;
  for(p=0;p<el.fi->np;p++) {
    c1 = c + p*dec ;
    /* initialisation */
    for(i=0;i<dec;i++) c1[i] = zero ;
    /* fonctions d'interpolation */
    h  = el.fi->h  + p*el.nn ;
    /* concentrations */
    c_fe    = param(u_1,h,el.nn,I_c_fe) ;
    c_cat   = param(u_1,h,el.nn,I_c_cat) ;
    c_ani   = param(u_1,h,el.nn,I_c_ani) ;
    c_oh    = concentration_oh(c_fe,c_cat,c_ani,k_eau) ;
    c_h     = k_eau/c_oh ;
    c_feoh2 = c_fe*c_oh*c_oh/k_feoh2 ;
    /* coefficient de transfert */
    kd_h2o   = va[p*NVE] ;
    kd_h     = va[p*NVE+1] ;
    kd_oh    = va[p*NVE+2] ;
    kd_fe    = va[p*NVE+3] ;
    kd_feoh2 = va[p*NVE+4] ;
    kd_o2    = va[p*NVE+5] ;
    kd_cat   = va[p*NVE+6] ;
    kd_ani   = va[p*NVE+7] ;
    d_oh     = va[p*NVE+8] ;
    d_h      = va[p*NVE+9] ;
    d_fe     = va[p*NVE+10] ;
    d_feoh2  = va[p*NVE+11] ;
    d_o2_l   = va[p*NVE+12] ;
    d_o2_g   = va[p*NVE+13] ;
    d_cat    = va[p*NVE+14] ;
    d_ani    = va[p*NVE+15] ;
    ke_h     = va[p*NVE+16] ;
    ke_oh    = va[p*NVE+17] ;
    ke_fe    = va[p*NVE+18] ;
    ke_cat   = va[p*NVE+19] ;
    ke_ani   = va[p*NVE+20] ;
    ke_h2o   = va[p*NVE+30] ;
    d_o2     = d_o2_l + d_o2_g/k_hen ;
    /* derivees */
    dc_hsdc_fe      = va[p*NVE+21] ;
    dc_hsdc_cat     = va[p*NVE+22] ;
    dc_hsdc_ani     = va[p*NVE+23] ;
    dc_ohsdc_fe     = va[p*NVE+24] ;
    dc_ohsdc_cat    = va[p*NVE+25] ;
    dc_ohsdc_ani    = va[p*NVE+26] ;
    dc_feoh2sdc_fe  = va[p*NVE+27] ;
    dc_feoh2sdc_cat = va[p*NVE+28] ;
    dc_feoh2sdc_ani = va[p*NVE+29] ;

    /* hydrogene */
    c2 = 2*kd_h2o + kd_h + kd_oh + 2*kd_feoh2 ;
    c1[(E_h*NEQ+I_p_l)*9+0]    = c2 ;
    c1[(E_h*NEQ+I_p_l)*9+4]    = c2 ;
    c1[(E_h*NEQ+I_p_l)*9+8]    = c2 ;

    c2 = d_h*dc_hsdc_fe + d_oh*dc_ohsdc_fe + 2*d_feoh2*dc_feoh2sdc_fe 
       - 2*M_oh/M_h2o*d_oh*dc_ohsdc_fe 
       - 2*M_h/M_h2o*d_h*dc_hsdc_fe
       - 2*M_feoh2/M_h2o*d_feoh2*dc_feoh2sdc_fe
       - 2*M_fe/M_h2o*d_fe ;
    c1[(E_h*NEQ+I_c_fe)*9+0]   = c2 ;
    c1[(E_h*NEQ+I_c_fe)*9+4]   = c2 ;
    c1[(E_h*NEQ+I_c_fe)*9+8]   = c2 ;

    c2 = - 2*M_o2/M_h2o*d_o2_l ;
    c1[(E_h*NEQ+I_c_o2)*9+0]   = c2 ;
    c1[(E_h*NEQ+I_c_o2)*9+4]   = c2 ;
    c1[(E_h*NEQ+I_c_o2)*9+8]   = c2 ;

    c2 = d_h*dc_hsdc_cat + d_oh*dc_ohsdc_cat + 2*d_feoh2*dc_feoh2sdc_cat
       - 2*M_oh/M_h2o*d_oh*dc_ohsdc_cat 
       - 2*M_h/M_h2o*d_h*dc_hsdc_cat
       - 2*M_cat/M_h2o*d_cat
       - 2*M_feoh2/M_h2o*d_feoh2*dc_feoh2sdc_cat ;
    c1[(E_h*NEQ+I_c_cat)*9+0]  = c2 ;
    c1[(E_h*NEQ+I_c_cat)*9+4]  = c2 ;
    c1[(E_h*NEQ+I_c_cat)*9+8]  = c2 ;

    c2 = d_h*dc_hsdc_ani + d_oh*dc_ohsdc_ani + 2*d_feoh2*dc_feoh2sdc_ani
       - 2*M_oh/M_h2o*d_oh*dc_ohsdc_ani 
       - 2*M_h/M_h2o*d_h*dc_hsdc_ani
       - 2*M_ani/M_h2o*d_ani
       - 2*M_feoh2/M_h2o*d_feoh2*dc_feoh2sdc_ani ;
    c1[(E_h*NEQ+I_c_ani)*9+0]  = c2 ;
    c1[(E_h*NEQ+I_c_ani)*9+4]  = c2 ;
    c1[(E_h*NEQ+I_c_ani)*9+8]  = c2 ;

    c2 = ke_h + ke_oh - 2*ke_h2o ;
    c1[(E_h*NEQ+I_psi)*9+0]    = c2 ;
    c1[(E_h*NEQ+I_psi)*9+4]    = c2 ;
    c1[(E_h*NEQ+I_psi)*9+8]    = c2 ;

    /* fe */
    c2 = kd_fe + kd_feoh2 ;
    c1[(E_fe*NEQ+I_p_l)*9+0]    = c2 ;
    c1[(E_fe*NEQ+I_p_l)*9+4]    = c2 ;
    c1[(E_fe*NEQ+I_p_l)*9+8]    = c2 ;

    c2 = d_fe + d_feoh2*dc_feoh2sdc_fe ;
    c1[(E_fe*NEQ+I_c_fe)*9+0]   = c2 ;
    c1[(E_fe*NEQ+I_c_fe)*9+4]   = c2 ;
    c1[(E_fe*NEQ+I_c_fe)*9+8]   = c2 ;

    c2 = d_feoh2*dc_feoh2sdc_cat ;
    c1[(E_fe*NEQ+I_c_cat)*9+0]  = c2 ;
    c1[(E_fe*NEQ+I_c_cat)*9+4]  = c2 ;
    c1[(E_fe*NEQ+I_c_cat)*9+8]  = c2 ;

    c2 = d_feoh2*dc_feoh2sdc_ani ;
    c1[(E_fe*NEQ+I_c_ani)*9+0]  = c2 ;
    c1[(E_fe*NEQ+I_c_ani)*9+4]  = c2 ;
    c1[(E_fe*NEQ+I_c_ani)*9+8]  = c2 ;

    c2 = ke_fe ;
    c1[(E_fe*NEQ+I_psi)*9+0]    = c2 ;
    c1[(E_fe*NEQ+I_psi)*9+4]    = c2 ;
    c1[(E_fe*NEQ+I_psi)*9+8]    = c2 ;

    /* courant i = j_h - j_oh + 2j_fe + j_cat - j_ani */
    c2 = d_h*dc_hsdc_fe - d_oh*dc_ohsdc_fe + 2*d_fe ;
    c1[(E_i*NEQ+I_c_fe)*9+0]   = c2 ;
    c1[(E_i*NEQ+I_c_fe)*9+4]   = c2 ;
    c1[(E_i*NEQ+I_c_fe)*9+8]   = c2 ;

    c2 = d_h*dc_hsdc_cat - d_oh*dc_ohsdc_cat + d_cat ;
    c1[(E_i*NEQ+I_c_cat)*9+0]  = c2 ;
    c1[(E_i*NEQ+I_c_cat)*9+4]  = c2 ;
    c1[(E_i*NEQ+I_c_cat)*9+8]  = c2 ;

    c2 = d_h*dc_hsdc_ani - d_oh*dc_ohsdc_ani - d_ani ;
    c1[(E_i*NEQ+I_c_ani)*9+0]  = c2 ;
    c1[(E_i*NEQ+I_c_ani)*9+4]  = c2 ;
    c1[(E_i*NEQ+I_c_ani)*9+8]  = c2 ;

    c2 = ke_h - ke_oh + 2*ke_fe + ke_cat - ke_ani ;
    c1[(E_i*NEQ+I_psi)*9+0]    = c2 ;
    c1[(E_i*NEQ+I_psi)*9+4]    = c2 ;
    c1[(E_i*NEQ+I_psi)*9+8]    = c2 ;

    /* o2 */
    c1[(E_o2*NEQ+I_p_l)*9+0]   = kd_o2 ;
    c1[(E_o2*NEQ+I_p_l)*9+4]   = kd_o2 ;
    c1[(E_o2*NEQ+I_p_l)*9+8]   = kd_o2 ;

    c1[(E_o2*NEQ+I_c_o2)*9+0]  = d_o2 ;
    c1[(E_o2*NEQ+I_c_o2)*9+4]  = d_o2 ;
    c1[(E_o2*NEQ+I_c_o2)*9+8]  = d_o2 ;

    /* cations */
    c1[(E_cat*NEQ+I_p_l)*9+0]  = kd_cat ;
    c1[(E_cat*NEQ+I_p_l)*9+4]  = kd_cat ;
    c1[(E_cat*NEQ+I_p_l)*9+8]  = kd_cat ;

    c1[(E_cat*NEQ+I_c_cat)*9+0]= d_cat ;
    c1[(E_cat*NEQ+I_c_cat)*9+4]= d_cat ;
    c1[(E_cat*NEQ+I_c_cat)*9+8]= d_cat ;

    c1[(E_cat*NEQ+I_psi)*9+0]  = ke_cat ;
    c1[(E_cat*NEQ+I_psi)*9+4]  = ke_cat ;
    c1[(E_cat*NEQ+I_psi)*9+8]  = ke_cat ;

    /* anions */
    c1[(E_ani*NEQ+I_p_l)*9+0]  = kd_ani ;
    c1[(E_ani*NEQ+I_p_l)*9+4]  = kd_ani ;
    c1[(E_ani*NEQ+I_p_l)*9+8]  = kd_ani ;

    c1[(E_ani*NEQ+I_c_ani)*9+0]= d_ani ;
    c1[(E_ani*NEQ+I_c_ani)*9+4]= d_ani ;
    c1[(E_ani*NEQ+I_c_ani)*9+8]= d_ani ;

    c1[(E_ani*NEQ+I_psi)*9+0]  = ke_ani ;
    c1[(E_ani*NEQ+I_psi)*9+4]  = ke_ani ;
    c1[(E_ani*NEQ+I_psi)*9+8]  = ke_ani ;
  }
  return(dec) ;
}

double tortuosite_l(double p_c,mate_t *mat)
{
  double phi = mat->pr[pm("phi")] ;
  double tau_l_sat = 0.296e-3*exp(9.95*phi)/phi ;
  return(tau_l_sat*courbe(p_c,mat->cb[2])) ;
}

double tortuosite_g(double p_c,mate_t *mat)
{
  double phi = mat->pr[pm("phi")] ;
  double tau_g_sat = pow(phi,1.74) ;
  return(tau_g_sat*courbe(p_c,mat->cb[3])) ;
}

double concentration_oh(double c_fe, double c_cat, double c_ani, double k_eau)
{
  double d = c_ani - c_cat - 2*c_fe ;
  return (0.5*(-d + sqrt(d*d + 4*k_eau))) ;
}
