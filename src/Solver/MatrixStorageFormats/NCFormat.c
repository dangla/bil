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
#include "NCFormat.h"



/* Extern functions */



NCFormat_t* NCFormat_CreateSelectedMatrix(Mesh_t* mesh,const int imatrix)
/** Create a matrix in NCFormat */
{
  NCFormat_t* ncformat = (NCFormat_t*) Mry_New(NCFormat_t) ;
  int n_col = Mesh_GetNbOfMatrixColumns(mesh)[imatrix] ;
  /* Working memory */
  int*   colptr0 = (int*) Mry_New(int[n_col + 1]) ;
  
  if(imatrix >= Mesh_GetNbOfMatrices(mesh)) {
    arret("NCFormat_CreateSelectedMatrix") ;
  }


  /* Initialize colptr0 */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int ie ;
    int i ;
    
    for(i = 0 ; i < n_col + 1 ; i++) colptr0[i] = 0 ;

    /* nb max de termes par colonne */
    for(ie = 0 ; ie < n_el ; ie++) {
      int neq = Element_GetNbOfEquations(el + ie) ;
    
      for(i = 0 ; i < Element_GetNbOfNodes(el + ie) ; i++) {
        Node_t* node_i = Element_GetNode(el + ie,i) ;
        int ieq ;
      
        for(ieq = 0 ; ieq < neq ; ieq++) {
          int ij = i*neq + ieq ;
          int ii = Element_GetUnknownPosition(el + ie)[ij] ;
          int jcol ;
        
          if(ii < 0) continue ;
        
          //jcol = Node_GetMatrixColumnIndex(node_i)[ii] ;
          jcol = Node_GetSelectedMatrixColumnIndexOf(node_i,ii,imatrix) ;
        
          if(jcol < 0) continue ;
          colptr0[jcol+1] += Element_GetNbOfNodes(el + ie)*neq ;
        }
      }
    }

    colptr0[0] = 0 ; /* deja nul ! */
    
    /* on calcul ou commence chaque colonne */
    for(i = 0 ; i < n_col ; i++) colptr0[i+1] += colptr0[i] ;
  }


  /* Allocation of space for colptr */
  {
    int* colptr = (int*) Mry_New(int[n_col+1]) ;
    
    NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) = colptr ;
  }
    
  /* Allocation of space for rowind */
  {
    int nnz_max = colptr0[n_col] ;
    int* rowind = (int*) Mry_New(int[nnz_max]) ;
    
    NCFormat_GetRowIndexOfNonZeroValue(ncformat) = rowind ;
  }


  /* Initialize colptr and rowind */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int* colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) ;
    int* rowind = NCFormat_GetRowIndexOfNonZeroValue(ncformat) ;
    int ie ;
    int i ;
    
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
          int j ;
        
          if(ii < 0) continue ;
        
          //irow = Node_GetMatrixColumnIndex(node_i)[ii] ;
          irow = Node_GetSelectedMatrixColumnIndexOf(node_i,ii,imatrix) ;
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
            
              //jcol = Node_GetMatrixColumnIndex(node_j)[jj] ;
              jcol = Node_GetSelectedMatrixColumnIndexOf(node_j,ii,imatrix) ;
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
    int* colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) ;
    int* rowind = NCFormat_GetRowIndexOfNonZeroValue(ncformat) ;
    int i ;
    int j ;
    
    colptr[0] = 0 ;
    
    for(j = 0 , i = 0 ; i < n_col ; i++) {
      int k ;
      
      for(k = 0 ; k < colptr[i + 1] ; k++) {
        rowind[j++] = rowind[colptr0[i] + k] ;
      }
      
      colptr[i + 1] += colptr[i] ;
    }

    if(colptr0[n_col] < colptr[n_col]) {
      arret("NCFormat_Create: memory issue") ;
    }
  }

  free(colptr0) ;


  {
    int* colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) ;
  
    NCFormat_GetNbOfNonZeroValues(ncformat) = colptr[n_col] ;
  }


  /* reallocation de la memoire */
  {
    int* rowind = NCFormat_GetRowIndexOfNonZeroValue(ncformat) ;
    int nnz = NCFormat_GetNbOfNonZeroValues(ncformat) ;
    int* rowind1 = (int*) Mry_Realloc(rowind,nnz*sizeof(int)) ;
    
    if(rowind1 != rowind) {
      NCFormat_GetRowIndexOfNonZeroValue(ncformat) = rowind1 ;
      Message_Warning("NCFormat_Create: new memory allocation") ;
    }
  }
  

  /*  1. allocation de l'espace memoire pour la matice */
  {
    int nnz = NCFormat_GetNbOfNonZeroValues(ncformat) ;
    double* nzval = (double*) Mry_New(double[nnz]) ;
    
    NCFormat_GetNonZeroValue(ncformat) = nzval ;
  }
  
  return(ncformat) ;
}



NCFormat_t* NCFormat_Create(Mesh_t* mesh)
/** Create a matrix in NCFormat */
{
  NCFormat_t* ncformat = (NCFormat_t*) Mry_New(NCFormat_t) ;
  int n_col = Mesh_GetNbOfMatrixColumns(mesh)[0] ;
  /* Working memory */
  int*   colptr0 = (int*) Mry_New(int[n_col + 1]) ;


  /* Initialize colptr0 */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int ie ;
    int i ;
    
    for(i = 0 ; i < n_col + 1 ; i++) colptr0[i] = 0 ;

    /* nb max de termes par colonne */
    for(ie = 0 ; ie < n_el ; ie++) {
      int neq = Element_GetNbOfEquations(el + ie) ;
    
      for(i = 0 ; i < Element_GetNbOfNodes(el + ie) ; i++) {
        Node_t* node_i = Element_GetNode(el + ie,i) ;
        int ieq ;
      
        for(ieq = 0 ; ieq < neq ; ieq++) {
          int ij = i*neq + ieq ;
          int ii = Element_GetUnknownPosition(el + ie)[ij] ;
          int jcol ;
        
          if(ii < 0) continue ;
        
          jcol = Node_GetMatrixColumnIndex(node_i)[ii] ;
        
          if(jcol < 0) continue ;
          colptr0[jcol+1] += Element_GetNbOfNodes(el + ie)*neq ;
        }
      }
    }

    colptr0[0] = 0 ; /* deja nul ! */
    
    /* on calcul ou commence chaque colonne */
    for(i = 0 ; i < n_col ; i++) colptr0[i+1] += colptr0[i] ;
  }


  /* Allocation of space for colptr */
  {
    int* colptr = (int*) Mry_New(int[n_col+1]) ;
    
    NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) = colptr ;
  }
    
  /* Allocation of space for rowind */
  {
    int nnz_max = colptr0[n_col] ;
    int* rowind = (int*) Mry_New(int[nnz_max]) ;
    
    NCFormat_GetRowIndexOfNonZeroValue(ncformat) = rowind ;
  }


  /* Initialize colptr and rowind */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int* colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) ;
    int* rowind = NCFormat_GetRowIndexOfNonZeroValue(ncformat) ;
    int ie ;
    int i ;
    
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
          int j ;
        
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
    int* colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) ;
    int* rowind = NCFormat_GetRowIndexOfNonZeroValue(ncformat) ;
    int i ;
    int j ;
    
    colptr[0] = 0 ;
    
    for(j = 0 , i = 0 ; i < n_col ; i++) {
      int k ;
      
      for(k = 0 ; k < colptr[i + 1] ; k++) {
        rowind[j++] = rowind[colptr0[i] + k] ;
      }
      
      colptr[i + 1] += colptr[i] ;
    }

    if(colptr0[n_col] < colptr[n_col]) {
      arret("NCFormat_Create: memory issue") ;
    }
  }

  free(colptr0) ;


  {
    int* colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(ncformat) ;
  
    NCFormat_GetNbOfNonZeroValues(ncformat) = colptr[n_col] ;
  }


  /* reallocation de la memoire */
  {
    int* rowind = NCFormat_GetRowIndexOfNonZeroValue(ncformat) ;
    int nnz = NCFormat_GetNbOfNonZeroValues(ncformat) ;
    int* rowind1 = (int*) Mry_Realloc(rowind,nnz*sizeof(int)) ;
    
    if(rowind1 != rowind) {
      NCFormat_GetRowIndexOfNonZeroValue(ncformat) = rowind1 ;
      Message_Warning("NCFormat_Create: new memory allocation") ;
    }
  }
  

  /*  1. allocation de l'espace memoire pour la matice */
  {
    int nnz = NCFormat_GetNbOfNonZeroValues(ncformat) ;
    double* nzval = (double*) Mry_New(double[nnz]) ;
    
    NCFormat_GetNonZeroValue(ncformat) = nzval ;
  }
  
  return(ncformat) ;
}




void NCFormat_Delete(void* self)
{
  NCFormat_t** a = (NCFormat_t**) self ;
  
  free(NCFormat_GetFirstNonZeroValueIndexOfColumn(*a)) ;
  free(NCFormat_GetRowIndexOfNonZeroValue(*a)) ;
  free(NCFormat_GetNonZeroValue(*a)) ;
  free(*a) ;
  *a = NULL ;
}




void NCFormat_AssembleElementMatrix(NCFormat_t* a,double* ke,int* cole,int* lige,int n,int* rowptr,int n_row)
/* Assemblage de la matrice elementaire ke dans la matrice globale a */
{
#define KE(i,j) (ke[(i)*n+(j)])
  double* nzval  = (double*) NCFormat_GetNonZeroValue(a) ;
  int*    colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(a) ;
  int*    rowind = NCFormat_GetRowIndexOfNonZeroValue(a) ;
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
        arret("NCFormat_AssembleElementMatrix: assembling not possible") ;
      }
      
      nzval[rowptr[irow]] += KE(ie,je) ;
    }

    for(i = colptr[jcol] ; i < colptr[jcol+1] ; i++) rowptr[rowind[i]] = -1 ;
  }

#undef KE
}




void NCFormat_PrintMatrix(NCFormat_t* a,unsigned int n_col,const char* keyword)
{
  double* nzval  = (double*) NCFormat_GetNonZeroValue(a) ;
  int*    colptr = NCFormat_GetFirstNonZeroValueIndexOfColumn(a) ;
  int*    rowind = NCFormat_GetRowIndexOfNonZeroValue(a) ;
  int    nnz = NCFormat_GetNbOfNonZeroValues(a) ;
  int    jcol ;

  fprintf(stdout,"\n") ;
  fprintf(stdout,"Matrix in compressed column format:\n") ;
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

