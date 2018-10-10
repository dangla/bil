#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "BilLib.h"
#include "CoordinateFormat.h"



CoordinateFormat_t* (CoordinateFormat_Create)(Mesh_t* mesh,Options_t* options)
/** Create a matrix in CoordinateFormat format with duplicate entries */
{
  CoordinateFormat_t* ac = (CoordinateFormat_t*) malloc(sizeof(CoordinateFormat_t)) ;

  assert(ac) ;
  
  {
    /* Nb of entries */
    {
      int nnz = Mesh_ComputeNbOfMatrixEntries(mesh) ;
      
      CoordinateFormat_GetNbOfNonZeroValues(ac) = nnz ;
    }
    
    /* Allocate memory space for the values */
    {
      int nnz = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
      /* The length required by ma38 must not be lower than 2*nnz */
      int ff = Options_GetFillFactor(options) ;
      int lv  = ff*(2*nnz) ;
      size_t sz = lv * sizeof(double) ;
      double* v = (double*) malloc(sz) ;
      
      assert(v) ;
      
      CoordinateFormat_GetLengthOfArrayValue(ac) = lv ;
      CoordinateFormat_GetNonZeroValue(ac) = v ;
    }
    
    /* Allocate memory space for the indices, namely:
     * the row and column indices;
     * some other informations provided by ma38 */
    {
      int nnz = CoordinateFormat_GetNbOfNonZeroValues(ac) ;
      int n = Mesh_GetNbOfMatrixColumns(mesh) ;
      /* The length required by ma38 must not be lower than 3*nnz+2*n+1 */
      int ff = Options_GetFillFactor(options) ;
      int lindex = ff*(3*nnz + 2*n + 1) ;
      size_t sz = lindex * sizeof(int) ;
      int* index = (int*) malloc(sz) ;
      
      assert(index) ;
      
      CoordinateFormat_GetLengthOfArrayIndex(ac) = lindex ;
      CoordinateFormat_GetIndex(ac)              = index ;
    }
  }
  
  CoordinateFormat_GetOptions(ac) = options ;
  
  return(ac) ;
}



void (CoordinateFormat_Delete)(void* self)
{
  CoordinateFormat_t** pac = (CoordinateFormat_t**) self ;
  CoordinateFormat_t*   ac = *pac ;
  
  free(CoordinateFormat_GetNonZeroValue(ac)) ;
  free(CoordinateFormat_GetIndex(ac)) ;
  free(ac) ;
  *pac = NULL ;
}







int CoordinateFormat_AssembleElementMatrix(CoordinateFormat_t* a,Element_t* el,double* ke, int nnz)
/** Assemble the local matrix ke into the global matrix a 
 *  Return the updated nb of entries */
{
#define KE(i,j) (ke[(i)*ndof+(j)])
  int   ndof   = Element_GetNbOfDOF(el) ;
  double*  val = CoordinateFormat_GetNonZeroValue(a) ;
  int*  colind = CoordinateFormat_GetColumnIndexOfValue(a) ;
  int*  rowind = CoordinateFormat_GetRowIndexOfValue(a) ;
  int*  row    = Element_ComputeMatrixRowAndColumnIndices(el) ;
  int*  col    = row + ndof ;
  
  
  {
    int   je ;

    for(je = 0 ; je < ndof ; je++) {
      int jcol = col[je] ;
      int ie ;
    
      if(jcol < 0) continue ;

      for(ie = 0 ; ie < ndof ; ie++) {
        int irow = row[ie] ;
      
        if(irow < 0) continue ;
      
        /* Row (and column) indices outside the range 1,nnz are ignored in ma38 */
        if(KE(ie,je) == 0.) irow = -1 ;
        
        val[nnz]    = KE(ie,je) ;
        colind[nnz] = jcol + 1 ;
        rowind[nnz] = irow + 1 ;
        
        nnz++ ;
      }
    }
  }
  
  return(nnz) ;

#undef KE
}




void CoordinateFormat_PrintMatrix(CoordinateFormat_t* a,const int n_col,const char* keyword)
{
  double* val  = CoordinateFormat_GetNonZeroValue(a) ;
  int*  colind = CoordinateFormat_GetColumnIndexOfValue(a) ;
  int*  rowind = CoordinateFormat_GetRowIndexOfValue(a) ;
  int   nnz    = CoordinateFormat_GetNbOfNonZeroValues(a) ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"Matrix in coordinate format with duplicate entries:\n") ;
  fprintf(stdout,"order = %u ; nb of entries = %d\n",n_col,nnz) ;
  fprintf(stdout,"\n") ;
  
  {
    int k ;
    
    for(k = 0 ; k < nnz ; k++) {
    
      fprintf(stdout,"%d: K(%d,%d) = %e",k,rowind[k],colind[k],val[k]) ;
    
      fprintf(stdout,"\n") ;
    }
  }
}
