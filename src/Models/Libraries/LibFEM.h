#ifndef LIBFEM_H
#define LIBFEM_H

#include "Elements.h"
#include "IntFcts.h"
#include "Loads.h"
#include "Geometry.h"

/* Fonctions */

extern void   def(double **,double **,double *,double *,double *,int,int,Symmetry_t,int) ;
extern void   grad(double **,double **,double *,double *,int,int,int) ;
extern void   mxcpl(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   mxmass(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   mxcond(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   mxcmss(double *,double **,IntFct_t,double *,int,int,Symmetry_t,int) ;
extern void   mxccnd(double *,double **,IntFct_t,double *,int,int,Symmetry_t,int) ;
extern void   mxrig(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   mxbiot(double *,double **,IntFct_t,double *,int,int,Symmetry_t,int,int,int,int,int) ;
extern double param(double **,double *,int,int) ;
extern void   rscont(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   rsflux(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   rsmass(double *,double **,IntFct_t,double *,int,int,Symmetry_t) ;
extern void   chsurf(double **,double *,int,Symmetry_t,double,double,Load_t,Element_t,IntFct_t*) ;



#endif
