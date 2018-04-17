#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "BilLib.h"
#include "Matrix.h"


#ifdef SLU_DIR
  #define val(x) #x
  #define xval(x) val(x)
  #include xval(SLU_DIR/SRC/dsp_defs.h)
  #undef val
  #undef xval
#endif


/* Extern functions */

Matrix_t*   Matrix_Create(Mesh_t* mesh,Options_t* options)
/* Alloue la memoire de la matrice */
{
  Matrix_t* a = (Matrix_t*) malloc(sizeof(Matrix_t)) ;
  
  if(!a) assert(a) ;


  /*  Nb of rows and columns */
  {
    Matrix_GetNbOfRows(a) = Nodes_GetNbOfMatrixRows(Mesh_GetNodes(mesh)) ;
    Matrix_GetNbOfColumns(a) = Nodes_GetNbOfMatrixColumns(Mesh_GetNodes(mesh)) ;
  }
  
  
  /*  Matrix format */
  {
    char* method = Options_GetResolutionMethod(options) ;
    
    if(!strcmp(method,"crout")) {
      Matrix_GetMatrixFormat(a) = LDUSKL ;
    
    } else if(!strcmp(method,"slu")) {
      Matrix_GetMatrixFormat(a) = HBSLU ;
    
    } else {
      arret("Matrix_Create(1): unknown method") ;
    }
  }


  /* Allocation of space for the matrix */
  {
    MatrixFormat_t fmt = Matrix_GetMatrixFormat(a) ;
    
    /* Skyline format */
    if(fmt == LDUSKL) {
      LDUSKLformat_t* askl = LDUSKLformat_Create(mesh) ;
    
      Matrix_GetStorage(a) = (void*) askl ;
      Matrix_GetNbOfNonZeroValues(a) = LDUSKLformat_GetNbOfNonZeroValues(askl) ;
      Matrix_GetNonZeroValue(a) = LDUSKLformat_GetNonZeroValue(askl) ;
    
    /* SuperMatrix format */
    #ifdef SLU_DIR
    } else if(Matrix_GetMatrixFormat(a) == HBSLU) {
      SuperMatrix_t* aslu = SuperMatrix_Create(mesh) ;
      NCformat_t*  asluNC = SuperMatrix_GetStorage(aslu) ;
    
      Matrix_GetStorage(a) = (void*) aslu ;
      Matrix_GetNbOfNonZeroValues(a) = NCformat_GetNbOfNonZeroValues(asluNC) ;
      Matrix_GetNonZeroValue(a) = NCformat_GetNonZeroValue(asluNC) ;

      /*  Work space for NCformat_AssembleElementMatrix */
      {
        int n_col = Matrix_GetNbOfColumns(a) ;
        void* work = (void*) malloc(n_col*sizeof(int)) ;
      
        if(!work) {
          arret("Matrix_Create(4): allocation impossible") ;
        }
      
        Matrix_GetWorkSpace(a) = work ;
      }
    #endif

    } else {
      arret("Matrix_Create(2): unknown format") ;
    }
  }


  /* State of the matrix */
  {
    Matrix_GetState(a) = 0 ;
  }
  
  return(a) ;
}



void Matrix_Delete(void* self)
{
  Matrix_t** a = (Matrix_t**) self ;
  MatrixFormat_t fmt = Matrix_GetMatrixFormat(*a) ;
  void* storage = Matrix_GetStorage(*a) ;


  /* Skyline format */
  if(fmt == LDUSKL) {
    LDUSKLformat_t* askl = (LDUSKLformat_t*) storage ;
      
    LDUSKLformat_Delete(&askl) ;
    
  /* SuperMatrix format */
  #ifdef SLU_DIR
  } else if(fmt == HBSLU) {
    SuperMatrix_t* aslu = (SuperMatrix_t*) storage ;
      
    SuperMatrix_Delete(&aslu) ;

    free(Matrix_GetWorkSpace(*a)) ;
  #endif

  } else {
    arret("Matrix_Delete(2): unknown format") ;
  }
}



void Matrix_AssembleElementMatrix_(Matrix_t* a,double* ke,int* cole,int* lige,int n)
/* Assemblage de la matrice elementaire ke dans la matrice globale a */
{
  /* Skyline format */
  if(Matrix_GetMatrixFormat(a) == LDUSKL) {
    LDUSKLformat_t* askl = (LDUSKLformat_t*) Matrix_GetStorage(a) ;
    
    LDUSKLformat_AssembleElementMatrix(askl,ke,cole,lige,n) ;
    return ;
    
#ifdef SLU_DIR
  /* CCS format (or Harwell-Boeing format) used in SuperLU */
  } else if(Matrix_GetMatrixFormat(a) == HBSLU) {
    SuperMatrix_t* aslu   = (SuperMatrix_t*) Matrix_GetStorage(a) ;
    NCformat_t*    asluNC = (NCformat_t*) SuperMatrix_GetStorage(aslu) ;
    int*        rowind = (int*) Matrix_GetWorkSpace(a) ;
    int         nrow = SuperMatrix_GetNbOfRows(aslu) ;
    
    NCformat_AssembleElementMatrix(asluNC,ke,cole,lige,n,rowind,nrow) ;
    return ;
#endif
  }
  
  arret("Matrix_AssembleElementMatrix : format non connu") ;
}



void Matrix_PrintMatrix(Matrix_t* a,const char* keyword)
{
  /* format Sky Line */
  if(Matrix_GetMatrixFormat(a) == LDUSKL) {
    LDUSKLformat_t* askl = (LDUSKLformat_t*) Matrix_GetStorage(a) ;
    int nrows = Matrix_GetNbOfRows(a) ;
    
    LDUSKLformat_PrintMatrix(askl,nrows,keyword) ;
    
#ifdef SLU_DIR
  /* format Harwell-Boeing de SuperLU */
  } else if(Matrix_GetMatrixFormat(a) == HBSLU) {
    SuperMatrix_t* aslu   = (SuperMatrix_t*) Matrix_GetStorage(a) ;
    NCformat_t*    asluNC = (NCformat_t*) SuperMatrix_GetStorage(aslu) ;
    int nrows = Matrix_GetNbOfRows(a) ;
    
    NCformat_PrintMatrix(asluNC,nrows,keyword) ;
#endif

  } else {
      arret("Matrix_PrintMatrix : format non connu") ;
  }

  fflush(stdout) ;
}
