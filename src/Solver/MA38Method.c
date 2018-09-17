#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include "MA38Method.h"
#include "Solver.h"
#include "Message.h"
#include "Matrix.h"


#include "CoordinateFormat.h"




/* 
 * HSL - MA38 method
 * -----------------
 */
int   MA38Method_Solve(Solver_t* solver)
/** Solve a sparse unsymmetric system a.x = b 
 *  by a multifrontal approach through ma38 from HSL 
 */
{
  Matrix_t* a = Solver_GetMatrix(solver) ;
  CoordinateFormat_t* ac = (CoordinateFormat_t*) Matrix_GetStorage(a) ;
  double* b = Solver_GetRHS(solver) ;
  double* x = Solver_GetSolution(solver) ;
  int n = Solver_GetNbOfColumns(solver) ;
  int ne = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
  double* w = (double*) Matrix_GetWorkSpace(a) ;
  
  int lvalue    = CoordinateFormat_GetLengthOfArrayValue(ac) ;
  int lindex    = CoordinateFormat_GetLengthOfArrayIndex(ac) ;
  double* value = CoordinateFormat_GetNonZeroValue(ac) ;
  int* index    = CoordinateFormat_GetIndex(ac) ;
  
  int    keep[20] ;
  int    icntl[20] ;
  double cntl[10] ;
  int    info[40] ;
  double rinfo[20] ;
  
  /* Initialize controls */
  ma38id_(keep,cntl,icntl) ;
  
  /* Printing controls:
   * 1: Print only error messages, 
   * 2: + warnings, 
   * 3: + terse diagnostics 
   */
  {
    Options_t* options = CoordinateFormat_GetOptions(ac) ;
    
    icntl[2] = 1 ;
    
    if(Options_IsToPrintOutAtEachIteration(options)) {
      icntl[2] = 2 ;
    }
  }
  
  /* Factorize the matrix */
  {
    int job = 0 ;
    bool transa = false ;
  
    ma38ad_(&n,&ne,&job,&transa,&lvalue,&lindex,value,index,keep,cntl,icntl,info,rinfo) ;
    
    if(info[0] < 0) {
      return(-1) ;
    }
  }
  
  /* Solve a * x = b */
  {
    int job = 0 ;
    bool transc = false ;
    
    ma38cd_(&n,&job,&transc,&lvalue,&lindex,value,index,keep,b,x,w,cntl,icntl,info,rinfo) ;
    
    if(info[0] < 0) {
      return(-1) ;
    }
  }
  
  return(0) ;
}







/* Not used */
#ifdef NOTDEFINED

int ludcmp1(LDUSKLFormat_t* a,int n)
/*
Remplace une matrice par sa decomposition LU (d apres NR).
"a" est une matrice ligne de ciel donnee sous la forme "LDU":

  D = a->d ,  U = a->u , L = a->l

  D    pointe sur le premier terme de la diagonale.
  U[0] pointe sur le premier terme non nul de la matrice triangulaire sup.
  L[0] pointe sur le premier terme non nul de la matrice triangulaire inf.
  U[j] - U[j-1] donne la hauteur de la colonne j de la matrice sup.
  L[i] - L[i-1] donne la longueur de la ligne i de la matrice inf.

diagonale   : a(i,i) = D[i]
matrice sup : a(i,j) = *(U[j]-j+i) si i1 <= i < j , a(i,j) = 0 sinon
matrice inf : a(i,j) = *(L[i]-i+j) si j1 <= j < i , a(i,j) = 0 sinon
avec          i1 = (j - U[j] + U[j-1])  et  j1 = (i - L[i] + L[i-1]).
*/
{
#define PU(i,j)     (LDUSKLFormat_GetPointerToUpperColumn(a)[j] - j + i)
#define PL(i,j)     (LDUSKLFormat_GetPointerToLowerRow(a)[i] - i + j)
  double** u = LDUSKLFormat_GetPointerToUpperColumn(a) ;
  double** l = LDUSKLFormat_GetPointerToLowerRow(a);
  double*  d = LDUSKLFormat_GetDiagonal(a) ;
  int    i,j,k,i1 = 0,j1,k1 ;
  double dum = 0.,*p,*pl,*pu ;
  double zero = 0. ;
  
  /* boucle sur les colonnes (algorithme de Crout) */
  for(j=0;j<(int) n;j++) {
    /* 1. U(i,j) = a(i,j) - sum_1^{i-1}L(i,k)*U(k,j) pour i<=j */
    if(j>0) { /* U(0,0) inchange */
      /* 1.a pour i<j */
      i1 = j - (u[j] - u[j-1]) ;   /* 1ere ligne non nulle de la colonne j */
      for(i=i1+1;i<j;i++) {          /* a(i1,j) inchange */
	j1 = i - (l[i] - l[i-1]) ; /* 1ere colonne non nulle de la ligne i */
	k1 = (i1 > j1) ? i1 : j1 ;
	p  = PU(i,j) ;  /* pointe sur a(i,j) i<j */
	pl = PL(i,k1) ; /* pointe sur a(i,k1) i>k1 */
	pu = PU(k1,j) ; /* pointe sur a(k1,j) k1<j */
	for(k=k1;k<i;k++) *p -= (*pl++)*(*pu++) ;
      }
      /* 1.b pour i=j */
      j1 = j - (l[j] - l[j-1]) ; /* 1ere colonne non nulle de la ligne j */
      k1 = (i1 > j1) ? i1 : j1 ;
      p  = d + j ;    /* pointe sur a(j,j) */
      pl = PL(i,k1) ; /* pointe sur a(i,k1) i>k1 */
      pu = PU(k1,j) ; /* pointe sur a(k1,j) k1<j */
      for(k=k1;k<j;k++) *p -= (*pl++)*(*pu++) ;
    }
    if(d[j] == zero) {
      /* arret("Le pivot est nul (ludcmp)") ; */
      return(-1) ;
    } else dum = 1./d[j] ;
    /* 2. L(i,j) = (a(i,j) - sum_1^{j-1}L(i,k)*U(k,j))/U(j,j) pour i>j */
    if(j==0) i1 = 0 ;
    for(i=j+1;i<(int) n;i++) {
      j1 = i - (l[i] - l[i-1]) ; /* 1ere colonne non nulle de la ligne i */
      if(j1 > j) continue ; /* a(i,j) i>j est nul */
      k1 = (i1 > j1) ? i1 : j1 ;
      p  = PL(i,j) ;  /* pointe sur a(i,j)  i>j */
      pl = PL(i,k1) ; /* pointe sur a(i,k1) i>k1 */
      pu = PU(k1,j) ; /* pointe sur a(k1,j) k1<j */
      for(k=k1;k<j;k++) *p -= (*pl++)*(*pu++) ;
      *p *= dum ;
    }
  }
  return(0) ;
#undef PU
#undef PL
}


void lubksb1(LDUSKLFormat_t* a,double* x,double* b,int n)
/* Resout a*x=b. a est donnee sous la forme LU (voir ludcmp) */
{
#define U(i,j)     (*(LDUSKLFormat_GetPointerToUpperColumn(a)[j] - j + i))
#define L(i,j)     (*(LDUSKLFormat_GetPointerToLowerRow(a)[i] - i + j))
  double** u = LDUSKLFormat_GetPointerToUpperColumn(a) ;
  double** l = LDUSKLFormat_GetPointerToLowerRow(a);
  double*  d = LDUSKLFormat_GetDiagonal(a) ;
  int    i,j,i1,j1 ;
  double* p ;
  
  /* descente : y(i) = b(i) - sum_1^{i-1}L(i,j)*b(j) */
  for(i=1;i<(int) n;i++) { /* b[0] inchangee */
    j1 = i - (l[i] - l[i-1]) ; /* 1ere colonne non nulle de la ligne i */
    p  = b+i ;
    for(j=j1;j<i;j++) *p -= L(i,j)*b[j] ;
  }
  /* remontee : x(i) = (y(i) - sum_{i+1}^nU(i,j)*x(j))/U(i,i) */
  for(i=(int) n-1;i>=0;i--) {
    p = b+i ;
    for(j=i+1;j<(int) n;j++) {
      i1 = j - (u[j] - u[j-1]) ; /* 1ere ligne non nulle de la colonne j */
      if(i1 > i) continue ; /* a(i,j) i<j est nul */
      *p -= U(i,j)*b[j] ;
    }
    *p /= d[i] ;
  }
#undef U
#undef L
}

#endif
