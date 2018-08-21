#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "BilLib.h"
#include "NCformat.h"



/* Extern functions */


#ifdef SLU_DIR

NCformat_t* NCformat_Create(Mesh_t* mesh)
/** Create a matrix in NCformat */
{
  int n_el = Mesh_GetNbOfElements(mesh) ;
  int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  int*   colptr0 ;
  int    ie ;
  int    i,j ;
  int    nnz_max ;
  NCformat_t*    asluNC = (NCformat_t*) malloc(sizeof(NCformat_t)) ;

  assert(asluNC) ;


  /* tableau de travail */
  colptr0 = (int*) malloc((n_col + 1)*sizeof(int)) ;
  
  assert(colptr0) ;

  for(i = 0 ; i < n_col + 1 ; i++) colptr0[i] = 0 ;

  /* nb max de termes par colonne */
  for(ie = 0 ; ie < n_el ; ie++) {
    int neq = Element_GetNbOfEquations(el + ie) ;
    
    for(i = 0 ; i < Element_GetNbOfNodes(el + ie) ; i++) {
      Node_t* node_i = Element_GetNode(el + ie,i) ;
      int ieq ;
      
      for(ieq = 0 ; ieq < neq ; ieq++) {
        int ij = i*neq + ieq ;
        int jj = Element_GetUnknownPosition(el + ie)[ij] ;
        int jcol ;
        
        if(jj < 0) continue ;
        
        jcol = Node_GetMatrixColumnIndex(node_i)[jj] ;
        
        if(jcol < 0) continue ;
        colptr0[jcol+1] += Element_GetNbOfNodes(el + ie)*neq ;
      }
    }
  }


  colptr0[0] = 0 ; /* deja nul ! */
  /* on calcul ou commence chaque colonne */
  for(i = 0 ; i < n_col ; i++) colptr0[i+1] += colptr0[i] ;
  nnz_max = colptr0[n_col] ;


  /* les tableaux colptr et rowind */
  {
    int* colptr = (int*) malloc((n_col+1)*sizeof(int)) ;
    
    assert(colptr) ;
    
    NCformat_GetFirstNonZeroValueIndexOfColumn(asluNC) = colptr ;
  }
    
  {
    int* rowind = (int*) malloc(nnz_max*sizeof(int)) ;
    
    assert(rowind) ;
    
    NCformat_GetRowIndexOfTheNonZeroValue(asluNC) = rowind ;
  }


  /* initialisation */
  {
    int* colptr = NCformat_GetFirstNonZeroValueIndexOfColumn(asluNC) ;
    int* rowind = NCformat_GetRowIndexOfTheNonZeroValue(asluNC) ;
    
    for(i = 0 ; i < n_col + 1 ; i++) colptr[i] = 0 ;

    for(ie = 0 ; ie < n_el ; ie++) {
      int neq = Element_GetNbOfEquations(el + ie) ;
    
      for(i = 0 ; i < Element_GetNbOfNodes(el + ie) ; i++) {
        Node_t* node_i = Element_GetNode(el + ie,i) ;
        int ieq ;
      
        for(ieq = 0 ; ieq < neq ; ieq++) {
          int iieq = i*neq + ieq ;
          int ii = Element_GetUnknownPosition(el + ie)[iieq] ;
          int irow ;
        
          if(ii < 0) continue ;
        
          irow = Node_GetMatrixColumnIndex(node_i)[ii] ;
          if(irow < 0) continue ;
        
          for(j = 0 ; j < Element_GetNbOfNodes(el + ie) ; j++) {
            Node_t* node_j = Element_GetNode(el + ie,j) ;
            int jeq ;
          
            for(jeq = 0 ; jeq < neq ; jeq++) {
              int jjeq = j*neq + jeq ;
              int jj = Element_GetUnknownPosition(el + ie)[jjeq] ;
              int jcol ;
              int k ;
            
              if(jj < 0) continue ;
            
              jcol = Node_GetMatrixColumnIndex(node_j)[jj] ;
              if(jcol < 0) continue ;
            
              /* on verifie que irow n'est pas deja enregistre */
              for(k = 0 ; k < colptr[jcol + 1] ; k++) {
                if(irow == rowind[colptr0[jcol] + k]) break ;
              }
            
              if(k < colptr[jcol + 1]) continue ;
            
              rowind[colptr0[jcol] + colptr[jcol + 1]] = irow ;
              colptr[jcol + 1] += 1 ;
            
              if(irow == jcol) continue ;
            
              rowind[colptr0[irow] + colptr[irow + 1]] = jcol ;
              colptr[irow + 1] += 1 ;
            }
          }
        }
      }
    }
  }


  /* compression de rowind */
  {
    int* colptr = NCformat_GetFirstNonZeroValueIndexOfColumn(asluNC) ;
    int* rowind = NCformat_GetRowIndexOfTheNonZeroValue(asluNC) ;
    
    colptr[0] = 0 ;
    
    for(j = 0 , i = 0 ; i < n_col ; i++) {
      int k ;
      
      for(k = 0 ; k < colptr[i + 1] ; k++) {
        rowind[j++] = rowind[colptr0[i] + k] ;
      }
      
      colptr[i + 1] += colptr[i] ;
    }
  }

  free(colptr0) ;


  {
    int* colptr = NCformat_GetFirstNonZeroValueIndexOfColumn(asluNC) ;
    
    if(nnz_max < colptr[n_col]) arret("NCformat_Create: memory issue") ;
  
    NCformat_GetNbOfNonZeroValues(asluNC) = colptr[n_col] ;
  }


  /* reallocation de la memoire */
  {
    int* rowind = NCformat_GetRowIndexOfTheNonZeroValue(asluNC) ;
    int nnz = NCformat_GetNbOfNonZeroValues(asluNC) ;
    int* rowind1 = (int*) realloc(rowind,nnz*sizeof(int)) ;
    
    if(!rowind1) arret("NCformat_Create: reallocation issue") ;
    
    if(rowind1 != rowind) {
      NCformat_GetRowIndexOfTheNonZeroValue(asluNC) = rowind1 ;
      Message_Warning("NCformat_Create: new memory allocation") ;
    }
  }
  

  /*  1. allocation de l'espace memoire pour la matice */
  {
    int nnz = NCformat_GetNbOfNonZeroValues(asluNC) ;
    double* nzval = (double*) malloc(nnz*sizeof(double)) ;
    
    assert(nzval) ;
    
    NCformat_GetNonZeroValue(asluNC) = nzval ;
  }
  
  return(asluNC) ;
}




void NCformat_Delete(void* self)
{
  NCformat_t** a = (NCformat_t**) self ;
  
  free(NCformat_GetFirstNonZeroValueIndexOfColumn(*a)) ;
  free(NCformat_GetRowIndexOfTheNonZeroValue(*a)) ;
  free(NCformat_GetNonZeroValue(*a)) ;
  free(*a) ;
  *a = NULL ;
}




void NCformat_AssembleElementMatrix(NCformat_t* a,double* ke,int* cole,int* lige,int n,int* rowptr,int n_row)
/* Assemblage de la matrice elementaire ke dans la matrice globale a */
{
#define KE(i,j) (ke[(i)*n+(j)])
  double* nzval  = NCformat_GetNonZeroValue(a) ;
  int*    colptr = NCformat_GetFirstNonZeroValueIndexOfColumn(a) ;
  int*    rowind = NCformat_GetRowIndexOfTheNonZeroValue(a) ;
  int    je,i ;

  for(i = 0 ; i < n_row ; i++) rowptr[i] = -1 ;
  
  for(je = 0 ; je < n ; je++) {
    int jcol = cole[je] ;
    int ie ;
    
    if(jcol < 0) continue ;

    for(i = colptr[jcol] ; i < colptr[jcol+1] ; i++) rowptr[rowind[i]] = i ;

    for(ie = 0 ; ie < n ; ie++) {
      int irow = lige[ie] ;
      
      if(irow < 0) continue ;
      
      if(rowptr[irow] < 0) {
        arret("NCformat_AssembleElementMatrix: assembling not possible") ;
      }
      
      nzval[rowptr[irow]] += KE(ie,je) ;
    }

    for(i = colptr[jcol] ; i < colptr[jcol+1] ; i++) rowptr[rowind[i]] = -1 ;
  }

#undef KE
}




void NCformat_PrintMatrix(NCformat_t* a,unsigned int n_col,const char* keyword)
{
  double* nzval  = NCformat_GetNonZeroValue(a) ;
  int*    colptr = NCformat_GetFirstNonZeroValueIndexOfColumn(a) ;
  int*    rowind = NCformat_GetRowIndexOfTheNonZeroValue(a) ;
  int    nnz = NCformat_GetNbOfNonZeroValues(a) ;
  int    jcol ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"matrice par colonne compressee :\n") ;
  fprintf(stdout,"n_col = %u nnz = %d\n",n_col,nnz) ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"\"col\" col: (lig)val ...\n") ;
  
  for(jcol = 0 ; jcol < (int) n_col ; jcol++) {
    int i ;
    
    fprintf(stdout,"col %d:",jcol) ;
    
    for(i = colptr[jcol] ; i < colptr[jcol+1] ; i++) {
      int irow = rowind[i] ;
      
      fprintf(stdout," (%d)% e",irow,nzval[i]) ;
    }
    
    fprintf(stdout,"\n") ;
  }
}

#endif
