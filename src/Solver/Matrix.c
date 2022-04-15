#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "BilLib.h"
#include "Mry.h"
#include "Matrix.h"
#include "MatrixStorageFormat.h"
#include "LDUSKLFormat.h"
#include "NCFormat.h"
#include "CoordinateFormat.h"

#ifdef SUPERLULIB
  #include "SuperLUFormat.h"
#endif


/* Extern functions */

Matrix_t*   (Matrix_Create)(Mesh_t* mesh,Options_t* options,const int imatrix)
/* Alloue la memoire de la matrice */
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
  
  
  /*  Matrix format */
  {
    char* method = Options_GetResolutionMethod(options) ;
    
    if(!strncmp(method,"crout",5)) {
      Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Type(LDUSKL) ;
    
    } else if(!strcmp(method,"slu")) {
      Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Type(SuperLU) ;

    } else if(!strcmp(method,"ma38")) {
      Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Type(Coordinate) ;
    
    } else {
      arret("Matrix_Create: unknown method") ;
    }
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
    #ifdef SUPERLULIB
    } else if(Matrix_StorageFormatIs(a,SuperLU)) {
      SuperLUFormat_t* aslu = SuperLUFormat_Create(mesh,imatrix) ;
      NCFormat_t*  asluNC = (NCFormat_t*) SuperLUFormat_GetStorage(aslu) ;
    
      Matrix_GetStorage(a) = (void*) aslu ;
      Matrix_GetNbOfNonZeroValues(a) = NCFormat_GetNbOfNonZeroValues(asluNC) ;
      Matrix_GetNonZeroValue(a) = (double*) NCFormat_GetNonZeroValue(asluNC) ;

      /*  Work space for NCFormat_AssembleElementMatrix */
      {
        int n_col = Matrix_GetNbOfColumns(a) ;
        void* work = (void*) Mry_New(int[n_col]) ;
      
        Matrix_GetWorkSpace(a) = work ;
      }
    #endif
    
    } else if(Matrix_StorageFormatIs(a,Coordinate)) {
      CoordinateFormat_t* ac = CoordinateFormat_Create(mesh,options,imatrix) ;
    
      Matrix_GetStorage(a) = (void*) ac ;
      Matrix_GetNbOfNonZeroValues(a) = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
      Matrix_GetNonZeroValue(a) = CoordinateFormat_GetNonZeroValue(ac) ;

      /*  Work space for ma38cd */
      {
        int n_col = Matrix_GetNbOfColumns(a) ;
        int lwork = 4 * n_col ;
        void* work = (void*) Mry_New(double[lwork]) ;
      
        Matrix_GetWorkSpace(a) = work ;
      }

    } else {
      arret("Matrix_Create: unknown format") ;
    }
  }


  /* State of the matrix */
  {
    Matrix_GetState(a) = 0 ;
  }
  
  return(a) ;
}


#if 0
Matrix_t*   (Matrix_Create)(Mesh_t* mesh,Options_t* options)
/* Alloue la memoire de la matrice */
{
  Matrix_t* a = (Matrix_t*) Mry_New(Matrix_t) ;
  
  Matrix_GetMatrixIndex(a) = 0 ;


  /*  Nb of rows and columns */
  {
    Matrix_GetNbOfRows(a)    = Mesh_GetNbOfMatrixRows(mesh)[0] ;
    Matrix_GetNbOfColumns(a) = Mesh_GetNbOfMatrixColumns(mesh)[0] ;
  }
  
  
  /*  Matrix format */
  {
    char* method = Options_GetResolutionMethod(options) ;
    
    if(!strcmp(method,"crout")) {
      Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Type(LDUSKL) ;
    
    } else if(!strcmp(method,"slu")) {
      Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Type(SuperLU) ;

    } else if(!strcmp(method,"ma38")) {
      Matrix_GetMatrixStorageFormat(a) = MatrixStorageFormat_Type(Coordinate) ;
    
    } else {
      arret("Matrix_Create(1): unknown method") ;
    }
  }


  /* Allocation of space for the matrix */
  {
    /* Skyline format */
    if(Matrix_StorageFormatIs(a,LDUSKL)) {
      LDUSKLFormat_t* askl = LDUSKLFormat_Create(mesh) ;
    
      Matrix_GetStorage(a) = (void*) askl ;
      Matrix_GetNbOfNonZeroValues(a) = LDUSKLFormat_GetNbOfNonZeroValues(askl) ;
      Matrix_GetNonZeroValue(a) = LDUSKLFormat_GetNonZeroValue(askl) ;
    
    /* SuperLU format */
    #ifdef SUPERLULIB
    } else if(Matrix_StorageFormatIs(a,SuperLU)) {
      SuperLUFormat_t* aslu = SuperLUFormat_Create(mesh) ;
      NCFormat_t*  asluNC = (NCFormat_t*) SuperLUFormat_GetStorage(aslu) ;
    
      Matrix_GetStorage(a) = (void*) aslu ;
      Matrix_GetNbOfNonZeroValues(a) = NCFormat_GetNbOfNonZeroValues(asluNC) ;
      Matrix_GetNonZeroValue(a) = (double*) NCFormat_GetNonZeroValue(asluNC) ;

      /*  Work space for NCFormat_AssembleElementMatrix */
      {
        int n_col = Matrix_GetNbOfColumns(a) ;
        void* work = (void*) Mry_New(int[n_col]) ;
      
        Matrix_GetWorkSpace(a) = work ;
      }
    #endif
    
    } else if(Matrix_StorageFormatIs(a,Coordinate)) {
      CoordinateFormat_t* ac = CoordinateFormat_Create(mesh,options) ;
    
      Matrix_GetStorage(a) = (void*) ac ;
      Matrix_GetNbOfNonZeroValues(a) = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
      Matrix_GetNonZeroValue(a) = CoordinateFormat_GetNonZeroValue(ac) ;

      /*  Work space for ma38cd */
      {
        int n_col = Matrix_GetNbOfColumns(a) ;
        int lwork = 4 * n_col ;
        void* work = (void*) Mry_New(double[lwork]) ;
      
        Matrix_GetWorkSpace(a) = work ;
      }

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
#endif



void (Matrix_Delete)(void* self)
{
  Matrix_t* a = (Matrix_t*) self ;
  void* storage = Matrix_GetStorage(a) ;


  /* Skyline format */
  if(Matrix_StorageFormatIs(a,LDUSKL)) {
    LDUSKLFormat_t* askl = (LDUSKLFormat_t*) storage ;
    
    if(askl) {
      LDUSKLFormat_Delete(askl) ;
      free(askl) ;
    }
    
  /* SuperLU format */
  #ifdef SUPERLULIB
  } else if(Matrix_StorageFormatIs(a,SuperLU)) {
    SuperLUFormat_t* aslu = (SuperLUFormat_t*) storage ;
    
    if(aslu) {
      SuperLUFormat_Delete(aslu) ;
      free(aslu) ;
    }
    
    {
      void* work = Matrix_GetWorkSpace(a) ;

      free(work) ;
    }
  #endif
  
  } else if(Matrix_StorageFormatIs(a,Coordinate)) {
    CoordinateFormat_t* ac = (CoordinateFormat_t*) storage ;
    
    if(ac) {
      CoordinateFormat_Delete(ac) ;
      free(ac) ;
    }
    
    {
      void* work = Matrix_GetWorkSpace(a) ;

      free(work) ;
    }

  } else {
    arret("Matrix_Delete(2): unknown format") ;
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
    
#ifdef SUPERLULIB
  /* CCS format (or Harwell-Boeing format) used in SuperLU */
  } else if(Matrix_StorageFormatIs(a,SuperLU)) {
    SuperLUFormat_t* aslu   = (SuperLUFormat_t*) Matrix_GetStorage(a) ;
    NCFormat_t*    asluNC = (NCFormat_t*) SuperLUFormat_GetStorage(aslu) ;
    int*        rowptr = (int*) Matrix_GetWorkSpace(a) ;
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
    
#ifdef SUPERLULIB
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
