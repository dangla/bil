#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "LDUSKLFormat.h"



/* Extern functions */


LDUSKLFormat_t* LDUSKLFormat_Create(Mesh_t* mesh)
/** Create a matrix in LDU Skyline format */
{
  LDUSKLFormat_t* a = (LDUSKLFormat_t*) malloc(sizeof(LDUSKLFormat_t)) ;
  
  assert(a) ;


  /* Allocation of space */
  {
    int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
    /*  les hauteurs de colonne (hc) */
    int*  hc = (int*) malloc(n_col*sizeof(int)) ;
    
    assert(hc) ;


    {
      int i ;
      
      for(i = 0 ; i < n_col ; i++) hc[i] = 0 ;
    }
    

    {
      int n_ddl ;
      
      {
        int n_no = Mesh_GetNbOfNodes(mesh) ;
        Node_t* no = Mesh_GetNode(mesh) ;
        int i ;
      
        n_ddl = 0 ;
        for(i = 0 ; i < n_no ; i++) n_ddl += Node_GetNbOfEquations(no + i) ;
      }
      
      /* The upper column heights: hc[i] = height of the upper column i */
      {
        int n_el = Mesh_GetNbOfElements(mesh) ;
        Element_t* el = Mesh_GetElement(mesh) ;
        int ie ;
        
        for(ie = 0 ; ie < n_el ; ie++) {
          if(Element_GetMaterial(el + ie)) {
            int   nn  = Element_GetNbOfNodes(el + ie) ;
            int   neq = Element_GetNbOfEquations(el + ie) ;
            int   k0 = n_ddl ;
            int i ;
    
            for(i = 0 ; i < nn ; i++) {
              Node_t* node_i = Element_GetNode(el + ie,i) ;
              int j ;
      
              for(j = 0 ; j < neq ; j++) {
                int ij = i*neq + j ;
                int ii = Element_GetUnknownPosition(el + ie)[ij] ;
        
                if(ii >= 0) {
                  int k  = Node_GetMatrixColumnIndex(node_i)[ii] ;
                  if(k >= 0 && k < k0) k0 = k ;
                }
              }
            }
    
            for(i = 0 ; i < nn ; i++) {
              Node_t* node_i = Element_GetNode(el + ie,i) ;
              int j ;
      
              for(j = 0 ; j < neq ; j++) {
                int ij = i*neq + j ;
                int ii = Element_GetUnknownPosition(el + ie)[ij] ;
        
                if(ii >= 0) {
                  int k  = Node_GetMatrixColumnIndex(node_i)[ii] ;
                  if(k >= 0 && k - k0 > hc[k]) hc[k] = k - k0 ;
                }
              }
            }
          }
        }
      }
    }
  
  
    {
      int   nnz_l ;
      
      /* Number of non zero values in the upper triangular matrix */
      {
        int i ;
      
        nnz_l = 0 ;
        for(i = 0 ; i < n_col ; i++) nnz_l += hc[i] ;
      }


      /* Allocation of space for the non zeros */
      {
        int nnz = 2*nnz_l + n_col ;
        double* z = (double*) malloc(nnz*sizeof(double)) ;
      
        assert(z) ;
    
        LDUSKLFormat_GetNbOfNonZeroValues(a) = nnz ;
        LDUSKLFormat_GetNonZeroValue(a)  = z ;
      }
  
  
      /* les tableaux de pointeurs de ligne et colonne */
      {
        double* z = LDUSKLFormat_GetNonZeroValue(a) ;
        double** p = (double**) malloc(2*n_col*sizeof(double*)) ;
        int i ;
      
        assert(p) ;
    
        LDUSKLFormat_GetPointerToLowerRow(a) = p ;
        LDUSKLFormat_GetPointerToUpperColumn(a) = p + n_col ;
    
        LDUSKLFormat_GetDiagonal(a) = z ;
    
        z += n_col ;
        LDUSKLFormat_GetPointerToLowerRow(a)[0] = z ;
        LDUSKLFormat_GetPointerToUpperColumn(a)[0] = z + nnz_l ;
      
        for(i = 1 ; i < n_col ; i++) {
          z = LDUSKLFormat_GetPointerToLowerRow(a)[i - 1] ;
          LDUSKLFormat_GetPointerToLowerRow(a)[i] = z + hc[i] ;
          LDUSKLFormat_GetPointerToUpperColumn(a)[i] = z + nnz_l + hc[i] ;
        }
  
      }
    }

    free(hc) ;
  }

  return(a) ;
}




void LDUSKLFormat_Delete(void* self)
{
  LDUSKLFormat_t** a = (LDUSKLFormat_t**) self ;
  
  free(LDUSKLFormat_GetNonZeroValue(*a)) ;
  free(LDUSKLFormat_GetPointerToLowerRow(*a)) ;
  free(*a) ;
  *a = NULL ;
}




void LDUSKLFormat_AssembleElementMatrix(LDUSKLFormat_t* a,double* ke,int* cole,int* lige,int n)
/* Assemblage de la matrice elementaire ke dans la matrice globale k */
{
#define KE(i,j) (ke[(i)*n+(j)])
#define D(i)    (LDUSKLFormat_GetDiagonal(a)[i])
#define U(i,j)  (LDUSKLFormat_GetUpperColumn(a,j)[i])
#define L(i,j)  (LDUSKLFormat_GetLowerRow(a,i)[j])
/*
#define U(i,j)  (*(LDUSKLFormat_GetPointerToUpperColumn(a)[j] - j + i))
#define L(i,j)  (*(LDUSKLFormat_GetPointerToLowerRow(a)[i] - i + j))
*/
  {
    int    ie ;
  
    for(ie = 0 ; ie < n ; ie++) { /* les lignes */
      int i = lige[ie] ;
      int je ;
    
      if(i < 0) continue ;
    
      for(je = 0 ; je < n ; je++) { /* les colonnes */
        int j = cole[je] ;
      
        if(j < 0) continue ;
      
        if(i == j)     D(i)   += KE(ie,je) ;
        else if(i < j) U(i,j) += KE(ie,je) ;
        else if(i > j) L(i,j) += KE(ie,je) ;
      }
    }
  }

#undef KE
#undef D
#undef U
#undef L
}




void LDUSKLFormat_PrintMatrix(LDUSKLFormat_t* a,unsigned int n,const char* keyword)
{
  double*  d = LDUSKLFormat_GetDiagonal(a) ;
  double** u = LDUSKLFormat_GetPointerToUpperColumn(a) ;
  double** l = LDUSKLFormat_GetPointerToLowerRow(a) ;
  int nnz = LDUSKLFormat_GetNbOfNonZeroValues(a) ;
  int    irow,jcol ;

  fprintf(stdout,"\n") ;
  
  fprintf(stdout,"LDU matrix:\n") ;
  fprintf(stdout,"n_col = %u nnz = %d\n",n,nnz) ;

  fprintf(stdout,"\n") ;
  
  fprintf(stdout,"diagonal \"diag\" diag: val\n") ;
  
  for(irow = 0 ; irow < (int) n ; irow++) {
    fprintf(stdout,"diag %d: % e\n",irow,d[irow]) ;
  }
  
  fprintf(stdout,"\n") ;
  
  if(!strcmp(keyword,"matrixdiag")) return ;
  
  if(n < 2) return ;

  fprintf(stdout,"sup matrix \"col\" col: (row)val ...\n") ;
  
  for(jcol = 1 ; jcol < (int) n ; jcol++) {
    double* p ;
    
    fprintf(stdout,"col %d:",jcol) ;
    
    for(p = u[jcol - 1] ; p < u[jcol] ; p++) {
      irow = jcol - (u[jcol] - p) ;
      fprintf(stdout," (%d)% e",irow,*p) ;
    }
    
    fprintf(stdout,"\n") ;
  }

  fprintf(stdout,"\n") ;
  
  fprintf(stdout,"inf matrix \"row\" row: (col)val ...\n") ;
  
  for(irow = 1 ; irow < (int) n ; irow++) {
    double* p ;
    
    fprintf(stdout,"row %d:",irow) ;
    
    for(p = l[irow - 1] ; p < l[irow] ; p++) {
      jcol = irow - (l[irow] - p) ;
      fprintf(stdout," (%d)% e",jcol,*p) ;
    }
    
    fprintf(stdout,"\n") ;
  }
}
