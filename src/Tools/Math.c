#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>

#include "Message.h"
#include "Exception.h"
#include "Math.h"


/* 
   Extern Functions 
*/



double Math_ComputeSecondDeviatoricStressInvariant(const double* sig)
/** Second invariant du deviateur des contraintes */
{
#define SIG(i,j) (sig[3*(i)+(j)])
  double j2a = (SIG(0,0) - SIG(1,1))*(SIG(0,0) - SIG(1,1))
             + (SIG(1,1) - SIG(2,2))*(SIG(1,1) - SIG(2,2))
             + (SIG(2,2) - SIG(0,0))*(SIG(2,2) - SIG(0,0)) ;
  double j2b = SIG(0,1)*SIG(1,0) + SIG(1,2)*SIG(2,1) + SIG(2,0)*SIG(0,2) ;
  return(j2a/6. + j2b) ;
#undef SIG
}



double* Math_SolveByGaussElimination(double* a,double* b,int n)
/** Solve the system a.x = b by Gaussian elimination with 
 *  backsubstitution with partial pivoting (interchange of rows).
 *  The rhs vector is replaced by the solution vector. */
{
#define  A(i,j)  (a[(i)*n+(j)])
#define  SWAP(a,b) {double tmp=(a);(a)=(b);(b)=tmp;}
  int    i ;
  int    irow = 0 ;

  for(i = 0 ; i < n ; i++) {
    double big = 0. ; /* Initialize for the search for largest pivot */
    int j ;
    
    /* Search for largest pivot */
    for(j = i ; j < n ; j++) {
      if(fabs(A(j,i)) >= big) {
	      big  = fabs(A(j,i)) ;
	      irow = j ;
      }
    }
    
    /* Do we need to interchange rows? */
    if(irow != i) {
      for(j = i ; j < n ; j++) SWAP(A(irow,j),A(i,j)) ;
      SWAP(b[irow],b[i]) ;
    }
    
    if(A(i,i) == 0.) {
      arret("Math_SolveByGaussElimination: pivot nul") ;
    } else {
      A(i,i) = 1./A(i,i) ;
    }
    
    for(j = i + 1 ; j < n ; j++) {
      int k ;
      A(j,i) *= A(i,i) ; /* Divide by the pivot */
      for(k = i + 1 ; k < n ; k++) A(j,k) -= A(j,i) * A(i,k) ;
      b[j] -= A(j,i) * b[i] ; /* Forward substitution */
    }
  }
  
  /* Backsubstitution */
  for(i = n - 1 ; i >= 0 ; i--) {
    int j ;
    b[i] *= A(i,i) ;
    for(j = i + 1 ; j < n ; j++) b[i] -= A(i,j) * b[j] ;
  }
  
  return(b) ;
#undef A
#undef SWAP
}



#ifdef NOTDEFINED  /* Not used */
extern double* (Math_SolveByGaussJordanElimination)(double*,double*,int,int) ;

double* Math_SolveByGaussJordanElimination(double* a,double* b,int n,int m)
/** Return matrix inverse and replace rhs vectors by solution vectors.
 *  We use Gauss-Jordan elimination with full pivoting 
 *  (after W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery. 
 *  Numerical Recipes, Cambridge University Press, 2007)
 *  The input matrix is a[n*n]. b[n*m] is input containing the m right-hand 
 *  side vectors. On output, a is replaced by its matrix inverse, and b is 
 *  replaced by the corresponding set of solution vectors.
 *  Used with m = 0, it only replaces a by its inverse.
 */
{
#define  A(i,j)  (a[(i)*n+(j)])
#define  B(i,j)  (b[(i)*m+(j)])
#define  SWAP(a,b) {double tmp=(a);(a)=(b);(b)=tmp;}
#define  MaxNbOfRows (10)
  int    indxc[MaxNbOfRows],indxr[MaxNbOfRows],ipiv[MaxNbOfRows] ;
  /* These integer arrays are used for bookkeeping on the pivoting. */
  int    icol = 0,irow = 0 ;
  int    i ;

  if(n > MaxNbOfRows) arret("gaussj (1)") ;

  for(i = 0 ; i < n ; i++) ipiv[i] = 0 ;
  
  for(i = 0 ; i < n ; i++) { 
    /* This is the main loop over the columns to be reduced */
    double big = 0. ;
    int j,l,ll ;
    
    for(j = 0 ; j < n ; j++) { 
      /* This is the outer loop of the search for a pivot element */
      if(ipiv[j] != 1) {
        int k ;
        
        for(k = 0 ; k < n ; k++) {
          if(ipiv[k] == 0) {
            if(fabs(A(j,k)) >= big) {
              big  = fabs(A(j,k)) ;
              irow = j ;
              icol = k ;
            }
          } else if(ipiv[k] > 1) {
            arret("gaussj (2) : matrice singuliere") ;
          }
        }
      }
    }
    
    ++(ipiv[icol]) ;
    
    /* We now have the pivot element, so we interchange rows, if needed, to put the pivot element on the diagonal. The columns are not physically interchanged, only relabeled: indxc[i], the column of the (i+1)th pivot element, is the (i+1)th column that is reduced, while indxr[i] is the row in which that pivot element was originally located. If indxr[i] != indxc[i], there is an implied column interchange. With this form of bookkeeping, the solution b’s will end up in the correct order, and the inverse matrix will be scrambled by columns. */
    
    if(irow != icol) {
      for(l = 0 ; l < n ; l++) SWAP(A(irow,l),A(icol,l)) ;
      for(l = 0 ; l < m ; l++) SWAP(B(irow,l),B(icol,l)) ;
    }
    
    indxr[i] = irow ; 
    indxc[i] = icol ;
    
    /* We are now ready to divide the pivot row by the pivot element, 
     * located at irow and icol. */
    
    if(A(icol,icol) == 0.) arret("gaussj (3) : matrice singuliere") ;
    
    {
      double pivinv = 1./A(icol,icol) ;
      
      A(icol,icol) = 1. ;
    
      for(l = 0 ; l < n ; l++) A(icol,l) *= pivinv ;
      for(l = 0 ; l < m ; l++) B(icol,l) *= pivinv ;
    }
    
    for(ll = 0 ; ll < n ; ll++) { /* Next, we reduce the rows... */
      if(ll != icol) {         /* ...except for the pivot one, of course. */
        double dum = A(ll,icol) ;
        
        A(ll,icol) = 0. ;
        
        for(l = 0 ; l < n ; l++) A(ll,l) -= A(icol,l)*dum ;
        for(l = 0 ; l < m ; l++) B(ll,l) -= B(icol,l)*dum ;
      }
    }
    
    /* This is the end of the main loop over columns of the reduction. It only remains to unscramble the solution in view of the column interchanges. We do this by interchanging pairs of columns in the reverse order that the permutation was built up. */

    for(l = n-1 ; l >= 0 ; l--) {
      if(indxr[l] != indxc[l]) {
        int k ;
        for(k = 0 ; k < n ; k++) SWAP(A(k,indxr[l]),A(k,indxc[l])) ;
      }
    }
  }
  
  return(a) ;
#undef A
#undef B
#undef SWAP
}
#endif



/*
 * Solve polynomial equations of degrees 1 to 4
 */
static int Math_ComputePolynomialEquationRootsOfDegree4(double*) ;
static int Math_ComputePolynomialEquationRootsOfDegree3(double*) ;
static int Math_ComputePolynomialEquationRootsOfDegree2(double*) ;

int Math_ComputePolynomialEquationRoots(double* x,int n)
/** Real Roots of Polynomial Equation of Degree n
 *  Return the nb of real solutions, sorted 
 *  from max to min */
{
#define MAX_DEGREE  (4)
  //double y[MAX_DEGREE + 1] ;
  int m = 0 ;
  
  /*
  {
    int i ;
    
    for(i = 0 ; i < n + 1 ; i++)  y[i] = x[i] ;
  }
  */
  
  if(x[0] != 0) {
    if(n == 0) {
      arret("Math_ComputePolynomialEquationRoots: degree 0") ;
    } else if(n == 1) {
      x[0] = - x[1]/x[0] ;
      return(1) ;
    } else if(n == 2) {
      m = Math_ComputePolynomialEquationRootsOfDegree2(x) ;
      return(m) ;
    } else if(n == 3) {
      m = Math_ComputePolynomialEquationRootsOfDegree3(x) ;
    } else if(n == 4) {
      m = Math_ComputePolynomialEquationRootsOfDegree4(x) ;
    } else {
      arret("Math_ComputePolynomialEquationRoots: degree too big") ;
    }
  } else {
    int i ;
    
    m = Math_ComputePolynomialEquationRoots(x + 1,n - 1) ;
    for(i = 0 ; i < m ; i++) x[i] = x[i + 1] ;
    return(m) ;
  }
  
  /* Polish the solutions */
  /*
  {
    double tol = 1.e-10 ;
    int i ;

    for(i = 0 ; i < m ; i++) {
      int k = Math_PolishPolynomialEquationRoot(y,n,x+i,fabs(z)*tol,10) ;
    }
  }
  */
  
  return(m) ;
#undef MAX_DEGREE
}



int Math_PolishPolynomialEquationRoot(double* x,int n,double* proot,double tol,int iterations)
{
  double root = proot[0] ;
  int it ;
  
  for(it = 0 ; it < iterations ; it++) {
    int i ;
    double error = x[0] ;
    double derivative = 0 ;
    double droot ;
      
    for(i = 0 ; i < n ; i++) {
      double a = error ;
      double b = x[i + 1] ;
      
      error = a*root + b ;
    }
      
    for(i = 0 ; i < n ; i++) {
      double a = derivative ;
      double b = (n - i)*x[i] ;
      
      derivative = a*root + b ;
    }
	    
    if(derivative == 0) {
      proot[0] = root ;
      return(0) ;
    }
	    
    droot = - error / derivative ;
    root += droot ;
      
    if(fabs(droot) < tol) {
      proot[0] = root ;
      return(root) ;
    }
  }
  
  /* Raise an interrupt signal instead of exit */
  Message_Warning("Math_PolishPolynomialEquationRoot: no convergence") ;
  {
    int i ;
    
    Message_Direct("\nthe %dth order polynomial equation:\n",n) ;
    
    for(i = 0 ; i <= n ; i++) {
      Message_Direct("%d order coeffcient: %lf\n",n-i,x[i]) ;
    }
    
    Message_Direct("\nhas not converged for the root %lf\n",root) ;
  }
  //Exception_Interrupt ;
  /*
  arret("Math_PolishPolynomialEquationRoot: not converged ") ;
  */
  return(-1) ;
}



/*
 * Evaluate an expression in a string
 * From snippets.
 * C snippets from "The snippets Collection"
 * http://www.brokersys.com/snippets/
 */

#include "Libraries/Snippets/eval.c"
#include "Libraries/Snippets/rmallws.c"
#include "Libraries/Snippets/strupr.c"

#undef strdup

//typedef enum {R_ERROR = -2 /* range */, ERROR /* syntax */, SUCCESS} STATUS;

double (Math_EvaluateExpression)(char *line)
{
  double val ;
  int i = evaluate(line,&val) ;
  
  if(i < 0) {
    /*
    if(i == ERROR) {
      Message_RuntimeError("Math_EvaluateExpression: syntax error\n%s",line) ;
    } else if(i == R_ERROR) {
      Message_RuntimeError("Math_EvaluateExpression: range error\n%s",line) ;
    }
    */
    Message_RuntimeError("Math_EvaluateExpression: syntax or range error\n%s",line) ;
  }
  
  return(val) ;
}



/*
 * Evaluate expressions in a string
 * C code created by a parser generator from AnaGram
 * http://www.parsifalsoft.com/
 */
 
#include "Libraries/evaluateExpression/evalkern.c"
#include "Libraries/evaluateExpression/evalwrap.c"

double (Math_EvaluateExpressions)(char* variablename,char* expressionstrings)
{
  double val ;
  int errorFlag = evaluateExpression(expressionstrings) ;
  
  if(errorFlag) {
    Message_RuntimeError("Math_EvaluateExpression1: syntax error\n\
                          %s at line %d, column %d\n",\
                        errorRecord.message,\
                        errorRecord.line,\
                        errorRecord.column) ;
  }
  
  {
    int i = 0 ;
    
    while(strncmp(variablename,variable[i].name,strlen(variablename))) i++ ;
    
    val = variable[i].value ;
  }
  
  return(val) ;
}



/*
 * Intern functions
 * ================
 * 
 */



/*
   Polynomial roots
*/
int Math_ComputePolynomialEquationRootsOfDegree2(double* x)
{
  double b = x[1]/x[0] ;
  double c = x[2]/x[0] ;
  double delta = b*b - 4*c ;
  
  if(c == 0) {
    if(b == 0) {
      x[0] = 0 ;
      return(1) ;
    } else {
      x[0] = MAX(-b,0) ;
      x[1] = MIN(-b,0) ;
      return(2) ;
    }
  } else if(delta == 0) {
    double temp1 = - b*0.5 ;
    x[0] = temp1 ;
    return(1) ;
  } else if(delta > 0) {
    double temp1 = - b*0.5;
    double temp2 = sqrt(delta)*0.5 ;
    if(b < 0) {
      x[0] = temp1 + temp2 ;
      x[1] = c/x[0] ;
    } else {
      x[1] = temp1 - temp2 ;
      x[0] = c/x[1] ;
    }
    return(2) ;
  }
  
  return(0) ;
}



int Math_ComputePolynomialEquationRootsOfDegree3(double* x)
{
  double b = x[1]/x[0] ;
  double c = x[2]/x[0] ;
  double d = x[3]/x[0] ;
  /* Depressed cubic equation : x3 + px + q = 0 */
  double p  = c - b*b/3 ;
  double q  = d + 2*b*b*b/27 - b*c/3 ;
  double p3 = p*p*p ;
  double q2 = q*q ;
  double r2 = p3/27 + q2/4 ;
  int n = 0 ;
  
  if(q == 0) {
    if(p >= 0) {
      x[0] = 0 ;
      n = 1 ;
    } else {
      double r = sqrt(-p) ;
      x[0] = r ;
      x[1] = 0 ;
      x[2] = - r ;
      n = 3 ;
    }
    
  } else if(r2 == 0) { /* q != 0 and p3/27 = -q2/4 (< 0) */
    double x1 = 3*q/p ;
    double x2 = - 0.5*x1 ; /* double root */
    x[0] = MAX(x1,x2) ;
    x[1] = MIN(x1,x2) ;
    n = 2 ;
    
  } else if(r2 > 0) {
    double r  = sqrt(r2) ;
    double r1 = (q < 0) ? r : -r ;
    double u3 = - 0.5*q + r1 ; /* so u3 is != 0 even for p = 0 */
    double u  = (u3 > 0) ? pow(u3,1/3.) : - pow(-u3,1/3.) ;
    x[0] = u - p/(3*u) ;
    n = 1 ;
    
  /* Trigonometric solution */
  } else { /* q != 0 and p3/27 < -q2/4 (< 0) */
    double srp3 = sqrt(-p/3) ;
    double k    = acos(1.5*q/(p*srp3))/3 ; /* 0 < k < pi/3 */
    double ck   = srp3*cos(k) ;
    double sk   = srp3*sin(k) ;
    double sr3  = sqrt(3.) ;
    double x1   = 2*ck ;
    double x2   = - ck + sr3*sk ;
    double x3   = - ck - sr3*sk ;
    x[0] = x1 ;
    x[1] = x2 ;
    x[2] = x3 ;
    n = 3 ;
  }
  
  {
    int i ;
    for(i = 0 ; i < n ; i++) x[i] -= b/3 ;
  }
  
  return(n) ;
}



int Math_ComputePolynomialEquationRootsOfDegree4(double* x)
{
  double b = x[1]/x[0] ;
  double c = x[2]/x[0] ;
  double d = x[3]/x[0] ;
  double e = x[4]/x[0] ;
  /* Depressed quartic equation: x4 + px2 + qx + r = 0 */
  double b2 = b*b ;
  double b3 = b2*b ;
  double b4 = b2*b2 ;
  double p  = - 0.375*b2 + c ;
  double q  =   0.125*b3 - 0.5*b*c + d ;
  double r  = - 0.01171875*b4 + 0.0625*b2*c - 0.25*b*d + e ;
  int n = 0 ;
  
  if(e == 0) {
    int i ;
    n = Math_ComputePolynomialEquationRootsOfDegree3(x) ;
    x[n] = 0 ;
    for(i = n ; i > 0 ; i--) {
      if(x[i - 1] >= 0) break ;
      x[i] = x[i - 1] ;
      x[i - 1] = 0 ;
    }
    return(n + 1) ;
    
  /* Biquadractic equations */
  } else if(q == 0) { 
    double y[3] = {1,0,0} ;
    int i,m ;
    y[1] = p ;
    y[2] = r ;
    m = Math_ComputePolynomialEquationRootsOfDegree2(y) ;
    /* First keep only positive roots */
    for(i = 0 ; i < m ; i++) {
      if(y[i] >= 0) {
        x[i] = sqrt(y[i]) ;
        n++ ;
      }
    }
    /* Second put the negative roots behind */
    for(i = 0 ; i < n ; i++) {
        x[2*n - 1 - i] =  - x[i] ;
    }
    n *= 2 ;
    
  /* Ferrari's solution */
  } else { 
    /* The nested depressed cubic equation: x3 + p1x + q1 */
    double p1 = - p*p/12 - r ;
    double q1 = - p*p*p/108 + p*r/3 - 0.125*q*q ;
    double y[4] = {1,0,0,0} ;
    y[2] = p1 ;
    y[3] = q1 ;
    /* Solving for the nested depressed cubic equation */
    Math_ComputePolynomialEquationRootsOfDegree3(y) ;
    {
      double s = y[0] - 5./6*p ;
      double w2 = p + 2*s ;
      if(w2 <= 0) { /* Theoretically it's never met */
        arret("Math_ComputePolynomialEquationRootsOfDegree4: never met?") ;
        return(0) ;
      } else {
        double w = sqrt(w2) ;
        double h1 = - (3*p + 2*s + 2*q/w) ;
        double h2 = - (3*p + 2*s - 2*q/w) ;
        if(h1 == 0) {
          x[0] = 0.5*w ;
          n = 1 ;
        } else if(h1 > 0) {
          x[0] = 0.5*(w + sqrt(h1)) ;
          x[1] = 0.5*(w - sqrt(h1)) ;
          n = 2 ;
        }
        if(h2 == 0) {
          x[n++] = - 0.5*w ;
        } else if(h2 > 0) {
          x[n++] = 0.5*(- w + sqrt(h2)) ;
          x[n++] = 0.5*(- w - sqrt(h2)) ;
        }
      }
    }
  }
  
  {
    int i ;
    for(i = 0 ; i < n ; i++) x[i] -= 0.25*b ;
  }
  
  return(n) ;
}
