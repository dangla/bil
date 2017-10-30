#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Solver.h"
#include "Message.h"
#include "BilLib.h"

#ifdef SLU_DIR
  #define val(x) #x
  #define xval(x) val(x)
  #include xval(SLU_DIR/SRC/dsp_defs.h)
  #undef val
  #undef xval
#endif


static int    solveCROUT(Solver_t*) ;

#ifdef SLU_DIR
static void   resolSLU(Matrix_t*,double*,double*,int) ;
static int    solveSLU(Solver_t*) ;
#endif

static int    ludcmp(LDUSKLformat_t*,int) ;

static void   lubksb(LDUSKLformat_t*,double*,double*,int) ;


/*
  Extern functions
*/

Solver_t*  Solver_Create(Mesh_t* mesh,Options_t* options)
{
  Solver_t* solver = (Solver_t*) malloc(sizeof(Solver_t)) ;
  
  if(!solver) arret("Solver_Create") ;
  
  
  /*  Method */
  {
    char* method = Options_GetResolutionMethod(options) ;
    
    if(!strcmp(method,"crout")) {
    
      Solver_GetResolutionMethod(solver) = CROUT ;
      Solver_GetSolve(solver) = solveCROUT ;
    
    #ifdef SLU_DIR
    } else if(!strcmp(method,"slu")) {
    
      Solver_GetResolutionMethod(solver) = SLU ;
      Solver_GetSolve(solver) = solveSLU ;
    #endif

    } else {
      arret("Solver_Create(1) : methode non connue") ;
    }
  }
  
  
  /* Nb of rows/columns */
  {
    int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
    
    Solver_GetNbOfColumns(solver) = n_col ;
  }
  
  
  /* Allocation of space for the matrix */
  {
    Solver_GetMatrix(solver) = Matrix_Create(mesh,options) ;
  }
  
  
  /* Allocation of space for the right hand side */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    double* rhs = (double*) malloc(n_col*sizeof(double)) ;
    
    if(!rhs) arret("Solver_Create") ;
    
    Solver_GetRHS(solver) = rhs ;
  }
  
  
  /* Allocation of space for the solution */
  {
    int n_col = Solver_GetNbOfColumns(solver) ;
    double* sol = (double*) malloc(n_col*sizeof(double)) ;
    
    if(!sol) arret("Solver_Create") ;
    
    Solver_GetSolution(solver) = sol ;
  }
  
  return(solver) ;
}



void Solver_Print(Solver_t* solver,char* keyword)
{
  static int i_debug=0 ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"debug(%d)\n",i_debug++) ;
  fprintf(stdout,"-----\n") ;
            
  if(!strcmp(keyword,"residu")) {
    double*  rhs = Solver_GetRHS(solver) ;
    int n_col = Solver_GetNbOfColumns(solver) ;
    int i ;
    
    fprintf(stdout,"\n") ;
    fprintf(stdout,"residu:\n") ;
    fprintf(stdout,"n = %d\n",n_col) ;
            
    for(i = 0 ; i < n_col ; i++) {
      fprintf(stdout,"res %d: % e\n",i,rhs[i]) ;
    }
  }
  
  if(!strncmp(keyword,"matrix",6)) {
    Matrix_t*  a = Solver_GetMatrix(solver) ;
    
    Matrix_PrintMatrix(a,keyword) ;
  }

}


/* Intern functions */

int ludcmp(LDUSKLformat_t* a,int n)
/**
 * Replace the given matrix a by the LDU decomposition (after NR).
 * "a" is arranged as a skyline matrix given as "LDU".

 * D    pointer to the diagonal.
 * U    pointer to the strictly upper triangular matrix.
 * L    pointer to the strictly lower triangular matrix.
 * 
 * U[j] pointer to the first non zero value of the column j + 1.
 * L[i] pointer to the first non zero value of the row i + 1.
 * U[j] - U[j - 1] is the height of the column j of the upper matrix.
 * L[i] - L[i - 1] is the height of the row i of the lower matrix.
 *
 * diagonal     : a(i,i) = D[i]
 * upper matrix : a(i,j) = (U[j] - j)[i] if i1 <= i < j , a(i,j) = 0 otherwise
 * lower matrix : a(i,j) = (L[i] - i)[j] if j1 <= j < i , a(i,j) = 0 otherwise
 * with           i1 = (j - U[j] + U[j-1])
 * and            j1 = (i - L[i] + L[i-1]).
 */
{
#define UpperColumn(j)      (LDUSKLformat_GetUpperColumn(a,j))
#define LowerRow(i)         (LDUSKLformat_GetLowerRow(a,i))
  double* diag = LDUSKLformat_GetDiagonal(a) ;
  double dum = 0. ;
  double zero = 0. ;
  int    j ;
  
  
  /* Loop on columns (Crout's algorithm) */
  for(j = 0 ; j < n ; j++) {
    double* colj = UpperColumn(j) ;
    /* i1 is the row index which starts the column j */
    int i1 = (j > 0) ? LDUSKLformat_RowIndexStartingColumn(a,j) : 0 ;
    int i ;
    
    /* 1. U_ij = a_ij - sum_1^{i-1} L_ik * U_kj for i <= j */
    if(j > 0) { /* U(0,0) unchanged */
    
      /* 1.a For i < j */
      for(i = i1 + 1 ; i < j ; i++) { /* a(i1,j) unchanged */
        /* j1 is the column index which starts the row i */
        int j1 = LDUSKLformat_ColumnIndexStartingRow(a,i) ;
        int k1 = (i1 > j1) ? i1 : j1 ;
        double* rowi = LowerRow(i) ;
        int k ;
        
        for(k = k1 ; k < i ; k++) {
          colj[i] -= rowi[k]*colj[k] ;
        }
      }
      
      /* 1.b For i = j */
      {
        /* Column index which starts the row j */
        int j1 = LDUSKLformat_ColumnIndexStartingRow(a,j) ;
        int k1 = (i1 > j1) ? i1 : j1 ;
        double* rowj = LowerRow(j) ;
        int k ;
        
        for(k = k1 ; k < j ; k++) {
          diag[j] -= rowj[k]*colj[k] ;
        }
      }
    }
    
    if(diag[j] == zero) {
      Message_Direct("\nludcmp: diagonal term (pivot) is zero at row/col %d",j) ;
      return(-1) ;

    } else dum = 1./diag[j] ;
    
    /* 2. L(i,j) = (a(i,j) - sum_1^{j-1}L(i,k)*U(k,j))/U(j,j) pour i > j */
    for(i = j + 1 ; i < n ; i++) {
      /* Column index which starts the row i */
      int j1 = LDUSKLformat_ColumnIndexStartingRow(a,i) ;
      
      if(j1 <= j) { /* a(i,j) i > j is zero */
        int k1 = (i1 > j1) ? i1 : j1 ;
        double* rowi = LowerRow(i) ;
        int k ;
      
        for(k = k1 ; k < j ; k++) {
          rowi[j] -= rowi[k]*colj[k] ;
        }
      
        rowi[j] *= dum ;
      }
    }
  }
  
  return(0) ;
#undef UpperColumn
#undef LowerRow
}


void lubksb(LDUSKLformat_t* a,double* x,double* b,int n)
/* Solve a*x = b. 
 * Here a is input as its LU decomposition determined by ludcmp.
 */
{
#define UpperColumn(j)      (LDUSKLformat_GetUpperColumn(a,j))
#define LowerRow(i)         (LDUSKLformat_GetLowerRow(a,i))
  double*  diag = LDUSKLformat_GetDiagonal(a) ;
  int    i ;
  
  if(x != b) {
    for(i = 0 ; i < n ; i++) x[i] = b[i] ;
  }
  
  /* Forward substitution : y_i = b_i - sum_0^{i-1} L_ij * y_j */
  for(i = 1 ; i < n ; i++) { /* x[0] unchanged */
    double* rowi = LowerRow(i) ;
    /* j1 is the column index which starts the row i */
    int j1 = LDUSKLformat_ColumnIndexStartingRow(a,i) ;
    int j ;
    
    for(j = j1 ; j < i ; j++) {
      x[i] -= rowi[j]*x[j] ;
    }
  }
  
  /* Backsubstitution : x_i = (y_i - sum_{i+1}^n U_ij * x_j)/D_ii */
  for(i = n - 1 ; i >= 0 ; i--) {
    int j ;
    
    for(j = i + 1 ; j < n ; j++) {
      /* i1 is the row index which starts the column j */
      int i1 = LDUSKLformat_RowIndexStartingColumn(a,j) ;
      
      /* a(i,j) = 0 for i < i1 */
      if(i >= i1) {
        double* colj = UpperColumn(j) ;
      
        x[i] -= colj[i]*x[j] ;
      }
    }
    
    x[i] /= diag[i] ;
  }
#undef UpperColumn
#undef LowerRow
}



int   solveCROUT(Solver_t* solver)
/** Resolution of a.x = b by Crout's method */
{
  Matrix_t* a = Solver_GetMatrix(solver) ;
  double* b = Solver_GetRHS(solver) ;
  double* x = Solver_GetSolution(solver) ;
  int n = Solver_GetNbOfColumns(solver) ;
  LDUSKLformat_t* askl = (LDUSKLformat_t*) Matrix_GetStorage(a) ;

  int i = ludcmp(askl,n) ;
  if(i < 0) return(i) ;
  lubksb(askl,x,b,n) ;
  return(0) ;
}



#ifdef SLU_DIR
int   solveSLU(Solver_t* solver)
/** Resolution of a.x = b by SuperLU's method */
{
  Matrix_t* a = Solver_GetMatrix(solver) ;
  double* b = Solver_GetRHS(solver) ;
  double* x = Solver_GetSolution(solver) ;
  int n = Solver_GetNbOfColumns(solver) ;
  
  resolSLU(a,x,b,n) ;
  return(0) ;
}


void   resolSLU(Matrix_t* a,double* x,double* b,int n)
/* Resolution of a.x = b by SuperLU's method */
{
  /* pour dgssv */
  SuperMatrix_t* Aslu = Matrix_GetStorage(a) ;
  SuperMatrix_t L,U,B ;
  DNformat    Bstore ;
  superlu_options_t options ;
  SuperLUStat_t stat ;
  int   info ;
  /* int   *perm_r = a.work,*perm_c = (int *) a.work + n ; */
  static int*   perm_r ;
  static int*   perm_c ;
  /* pour dgssvx */
  char          equed[1] ;
  double        rpg ;
  double        rcond ;
  mem_usage_t   mem_usage ;
  static int    iresol = 0 ;
  static int*   etree ;
  SuperMatrix_t X ;
  DNformat    Xstore ;
  static double* R ;
  static double* C ;
  static void*   work ;
  static int    lwork = 0 ;
  static double ferr[1] ;
  static double berr[1] ;
  
  /* allocation de memoire pour les variables statiques */
  if(!iresol) {
    perm_r = (int*) malloc(n*sizeof(int)) ;
    if(!perm_r) arret("resolSLU(1) : allocation impossible") ;
    perm_c = (int*) malloc(n*sizeof(int)) ;
    if(!perm_c) arret("resolSLU(2) : allocation impossible") ;
    etree = (int*) malloc(n*sizeof(int)) ;
    if(!etree) arret("resolSLU(3) : allocation impossible") ;
    R = (double*) malloc(n*sizeof(double)) ;
    if(!R) arret("resolSLU(5) : allocation impossible") ;
    C = (double*) malloc(n*sizeof(double)) ;
    if(!C) arret("resolSLU(6) : allocation impossible") ;
    if(lwork > 0) {
      work = (double*) malloc(lwork*sizeof(double)) ;
      if(!work) arret("resolSLU(7) : allocation impossible") ;
    }
  }

  /* le second membre */
  B.Stype = SLU_DN ;
  B.Dtype = SLU_D ;
  B.Mtype = SLU_GE ;
  B.nrow  = n ;
  B.ncol  = 1 ;
  B.Store = (DNformat *) &Bstore ;
  Bstore.lda = n ;
  Bstore.nzval = (double*) b ;
  /* dCreate_Dense_Matrix(&B,n,1,b,n,SLU_DN,SLU_D,SLU_GE) ; */

  /* la solution */
  X.Stype = SLU_DN ;
  X.Dtype = SLU_D ;
  X.Mtype = SLU_GE ;
  X.nrow  = n ;
  X.ncol  = 1 ;
  X.Store = (DNformat *) &Xstore ;
  Xstore.lda = n ;
  Xstore.nzval = (double*) x ;
  
  /* les options */
  /* par defaut 
     options->Fact = DOFACT;
     options->Equil = YES;
     options->ColPerm = COLAMD;
     options->DiagPivotThresh = 1.0;
     options->Trans = NOTRANS;
     options->IterRefine = NOREFINE;
     options->SymmetricMode = NO;
     options->PivotGrowth = NO;
     options->ConditionNumber = NO;
     options->PrintStat = YES;
  */
  set_default_options(&options) ;
  if(iresol++) options.Fact = SamePattern ;
  options.ColPerm = NATURAL ;
  options.PrintStat = NO ;
  /* options.DiagPivotThresh = 0. ; */
  /* options.IterRefine = EXTRA ; */

  /* les statistiques */
  StatInit(&stat) ;

  /* dgssv(&options,Aslu,perm_c,perm_r,&L,&U,&B,&stat,&info) ; */

  dgssvx(&options,Aslu,perm_c,perm_r,etree,equed,R,C,
	 &L,&U,work,lwork,&B,&X,&rpg,&rcond,ferr,berr,
	 &mem_usage,&stat,&info);

  /* if(options.PrintStat) StatPrint(&stat) ; */
}
#endif







/* Not used */
#ifdef NOTDEFINED

int ludcmp1(LDUSKLformat_t* a,int n)
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
#define PU(i,j)     (LDUSKLformat_GetPointerToUpperColumn(a)[j] - j + i)
#define PL(i,j)     (LDUSKLformat_GetPointerToLowerRow(a)[i] - i + j)
  double** u = LDUSKLformat_GetPointerToUpperColumn(a) ;
  double** l = LDUSKLformat_GetPointerToLowerRow(a);
  double*  d = LDUSKLformat_GetDiagonal(a) ;
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


void lubksb1(LDUSKLformat_t* a,double* x,double* b,int n)
/* Resout a*x=b. a est donnee sous la forme LU (voir ludcmp) */
{
#define U(i,j)     (*(LDUSKLformat_GetPointerToUpperColumn(a)[j] - j + i))
#define L(i,j)     (*(LDUSKLformat_GetPointerToLowerRow(a)[i] - i + j))
  double** u = LDUSKLformat_GetPointerToUpperColumn(a) ;
  double** l = LDUSKLformat_GetPointerToLowerRow(a);
  double*  d = LDUSKLformat_GetDiagonal(a) ;
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
