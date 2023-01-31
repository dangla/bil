#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "BilExtraLibs.h"
#include "Mry.h"
#include "Matrix.h"
#include "MatrixStorageFormat.h"
#include "LDUSKLFormat.h"
#include "NCFormat.h"
#include "CoordinateFormat.h"

#if defined (SUPERLULIB) || defined (SUPERLUMTLIB) || defined (SUPERLUDISTLIB)
  #define SUPERLU
  #include "SuperLUFormat.h"
#else
  #undef SUPERLU
#endif


/* Extern functions */

Matrix_t*   (Matrix_Create)(Mesh_t* mesh,Options_t* options,const int imatrix)
{
  Matrix_t* a = (Matrix_t*) Mry_New(Matrix_t) ;
  
  if(imatrix >= Mesh_GetNbOfMatrices(mesh)) {
    arret("Matrix_Create") ;
  }
  
  Matrix_GetMatrixIndex(a) = imatrix ;


  /*  Nb of rows and columns */
  {
    Matrix_GetNbOfRows(a)    = Mesh_GetNbOfMatrixRows(mesh)[imatrix] ;
    Matrix_GetNbOfColumns(a) = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;
  }
  
  
  /* Allocate space for the permutation of rows */
  {
    int nrows = Matrix_GetNbOfRows(a) ;
    int* rowperm = (int*) Mry_New(int[nrows]) ;
    
    Matrix_GetRowPermutation(a) = rowperm ;
  }
  
  
  /* Allocate space for the permutation of columns */
  {
    int ncols = Matrix_GetNbOfColumns(a) ;
    int* colperm = (int*) Mry_New(int[ncols]) ;
    
    Matrix_GetColumnPermutation(a) = colperm ;
  }
  
  
  /* Matrix storage format */
  {
    Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Create(options) ;
  }


  /* Allocation of space for the matrix */
  {
    /* Skyline format */
    if(Matrix_StorageFormatIs(a,LDUSKL)) {
      LDUSKLFormat_t* askl = LDUSKLFormat_Create(mesh,imatrix) ;
    
      Matrix_GetStorage(a) = (void*) askl ;
      Matrix_GetNbOfNonZeroValues(a) = LDUSKLFormat_GetNbOfNonZeroValues(askl) ;
      Matrix_GetNonZeroValue(a) = LDUSKLFormat_GetNonZeroValue(askl) ;
    
    /* SuperLU format */
    #ifdef SUPERLU
    } else if(Matrix_StorageFormatIs(a,SuperLU)) {
      SuperLUFormat_t* aslu = SuperLUFormat_Create(mesh,imatrix) ;
      NCFormat_t*  asluNC = (NCFormat_t*) SuperLUFormat_GetStorage(aslu) ;
    
      Matrix_GetStorage(a) = (void*) aslu ;
      Matrix_GetNbOfNonZeroValues(a) = NCFormat_GetNbOfNonZeroValues(asluNC) ;
      Matrix_GetNonZeroValue(a) = (double*) NCFormat_GetNonZeroValue(asluNC) ;

      /*  Work space for NCFormat_AssembleElementMatrix */
      {
        int n_col = Matrix_GetNbOfColumns(a) ;
        void* work = Mry_New(int[n_col]) ;
        GenericData_t* gw = GenericData_Create(n_col,work,int,"rowptr") ;
      
        Matrix_AppendGenericWorkSpace(a,gw) ;
      }
    #endif
    
    } else if(Matrix_StorageFormatIs(a,Coordinate)) {
      CoordinateFormat_t* ac = CoordinateFormat_Create(mesh,options,imatrix) ;
    
      Matrix_GetStorage(a) = (void*) ac ;
      Matrix_GetNbOfNonZeroValues(a) = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
      Matrix_GetNonZeroValue(a) = CoordinateFormat_GetNonZeroValue(ac) ;

    } else {
      arret("Matrix_Create: unknown format") ;
    }
  }


  /* State of the matrix */
  {
    Matrix_SetToUnfactorizedState(a) ;
  }


  /* Pattern state of the matrix */
  {
    Matrix_UnsetSameSparsityPattern(a) ;
  }
  
  return(a) ;
}



void (Matrix_Delete)(void* self)
{
  Matrix_t* a = (Matrix_t*) self ;
  
  {
    int* rowperm = Matrix_GetRowPermutation(a) ;
    
    if(rowperm) {
      free(rowperm) ;
    }
  }
  
  {
    int* colperm = Matrix_GetColumnPermutation(a) ;
    
    if(colperm) {
      free(colperm) ;
    }
  }


  {
    void* storage = Matrix_GetStorage(a) ;
  
    /* Skyline format */
    if(Matrix_StorageFormatIs(a,LDUSKL)) {
      LDUSKLFormat_t* askl = (LDUSKLFormat_t*) storage ;
    
      if(askl) {
        LDUSKLFormat_Delete(askl) ;
        free(askl) ;
      }
    
    /* SuperLU format */
    #ifdef SUPERLU
    } else if(Matrix_StorageFormatIs(a,SuperLU)) {
      SuperLUFormat_t* aslu = (SuperLUFormat_t*) storage ;
    
      if(aslu) {
        SuperLUFormat_Delete(aslu) ;
        free(aslu) ;
      }
    #endif
  
    } else if(Matrix_StorageFormatIs(a,Coordinate)) {
      CoordinateFormat_t* ac = (CoordinateFormat_t*) storage ;
    
      if(ac) {
        CoordinateFormat_Delete(ac) ;
        free(ac) ;
      }

    } else {
      arret("Matrix_Delete: unknown format") ;
    }
  }


  {
    MatrixStorageFormat_t* msf = Matrix_GetMatrixStorageFormat(a) ;

    if(msf) {
      MatrixStorageFormat_Delete(msf) ;
      free(msf) ;
    }
  }
  
  {
    GenericData_t* genericwork = Matrix_GetGenericWorkSpace(a) ;
    
    if(genericwork) {
      GenericData_Delete(genericwork) ;
      free(genericwork) ;
      Matrix_GetGenericWorkSpace(a) = NULL ;
    }
  }
}


void Matrix_AssembleElementMatrix(Matrix_t* a,Element_t* el,double* ke)
/** Assemble the element matrix ke in the global matrix a */
{
  int imatrix = Matrix_GetMatrixIndex(a) ;
  int  ndof = Element_GetNbOfDOF(el) ;
  int* row = Element_ComputeSelectedMatrixRowAndColumnIndices(el,imatrix) ;
  int* col = row + ndof ;

  /* Skyline format */
  if(Matrix_StorageFormatIs(a,LDUSKL)) {
    LDUSKLFormat_t* askl = (LDUSKLFormat_t*) Matrix_GetStorage(a) ;
    
    LDUSKLFormat_AssembleElementMatrix(askl,ke,col,row,ndof) ;
    return ;
    
#ifdef SUPERLU
  /* CCS format (or Harwell-Boeing format) used in SuperLU */
  } else if(Matrix_StorageFormatIs(a,SuperLU)) {
    SuperLUFormat_t* aslu   = (SuperLUFormat_t*) Matrix_GetStorage(a) ;
    NCFormat_t*    asluNC = (NCFormat_t*) SuperLUFormat_GetStorage(aslu) ;
    GenericData_t* gw = Matrix_GetGenericWorkSpace(a) ;
    int*        rowptr = GenericData_FindData(gw,int,"rowptr") ;
    int         nrow = SuperLUFormat_GetNbOfRows(aslu) ;
    
    NCFormat_AssembleElementMatrix(asluNC,ke,col,row,ndof,rowptr,nrow) ;
    return ;
#endif
  
  } else if(Matrix_StorageFormatIs(a,Coordinate)) {
    CoordinateFormat_t* ac = (CoordinateFormat_t*) Matrix_GetStorage(a) ;
    int len = Matrix_GetNbOfEntries(a) ;
    
    len = CoordinateFormat_AssembleElementMatrix(ac,ke,col,row,ndof,len) ;
    
    Matrix_GetNbOfEntries(a) = len ;
    return ;
  }
  
  arret("Matrix_AssembleElementMatrix: unknown format") ;
}




void Matrix_PrintMatrix(Matrix_t* a,const char* keyword)
{
  /* format Sky Line */
  if(Matrix_StorageFormatIs(a,LDUSKL)) {
    LDUSKLFormat_t* askl = (LDUSKLFormat_t*) Matrix_GetStorage(a) ;
    int nrows = Matrix_GetNbOfRows(a) ;
    
    LDUSKLFormat_PrintMatrix(askl,nrows,keyword) ;
    
#ifdef SUPERLU
  /* format Harwell-Boeing de SuperLU */
  } else if(Matrix_StorageFormatIs(a,SuperLU)) {
    SuperLUFormat_t* aslu   = (SuperLUFormat_t*) Matrix_GetStorage(a) ;
    NCFormat_t*    asluNC = (NCFormat_t*) SuperLUFormat_GetStorage(aslu) ;
    int nrows = Matrix_GetNbOfRows(a) ;
    
    NCFormat_PrintMatrix(asluNC,nrows,keyword) ;
#endif
  
  } else if(Matrix_StorageFormatIs(a,Coordinate)) {
    CoordinateFormat_t* ac = (CoordinateFormat_t*) Matrix_GetStorage(a) ;
    int nrows = Matrix_GetNbOfRows(a) ;
    
    CoordinateFormat_PrintMatrix(ac,nrows,keyword) ;
    return ;

  } else {
      arret("Matrix_PrintMatrix: unknown format") ;
  }

  fflush(stdout) ;
}
