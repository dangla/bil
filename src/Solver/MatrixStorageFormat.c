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
#include "MatrixStorageFormat.h"



/* Extern functions */


LDUSKLformat_t* LDUSKLformat_Create(Mesh_t* mesh)
/** Create a matrix in LDU Skyline format */
{
  LDUSKLformat_t* a = (LDUSKLformat_t*) malloc(sizeof(LDUSKLformat_t)) ;
  
  if(!a) {
    assert(a) ;
  }


  /* Allocation of space */
  {
    int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
    /*  les hauteurs de colonne (hc) */
    int*  hc = (int*) malloc(n_col*sizeof(int)) ;
    
    if(!hc) {
      arret("LDUSKLformat_Create (1) : impossible d\'allouer la memoire") ;
    }


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
                int jj = Element_GetUnknownPosition(el + ie)[ij] ;
        
                if(jj >= 0) {
                  int k  = Node_GetMatrixColumnIndex(node_i)[jj] ;
                  if(k >= 0 && k < k0) k0 = k ;
                }
              }
            }
    
            for(i = 0 ; i < nn ; i++) {
              Node_t* node_i = Element_GetNode(el + ie,i) ;
              int j ;
      
              for(j = 0 ; j < neq ; j++) {
                int ij = i*neq + j ;
                int jj = Element_GetUnknownPosition(el + ie)[ij] ;
        
                if(jj >= 0) {
                  int k  = Node_GetMatrixColumnIndex(node_i)[jj] ;
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
      
        if(!z) arret("LDUSKLformat_Create (3) : impossible d\'allouer la memoire") ;
    
        LDUSKLformat_GetNbOfNonZeroValues(a) = nnz ;
        LDUSKLformat_GetNonZeroValue(a)  = z ;
      }
  
  
      /* les tableaux de pointeurs de ligne et colonne */
      {
        double* z = LDUSKLformat_GetNonZeroValue(a) ;
        double** p = (double**) malloc(2*n_col*sizeof(double*)) ;
        int i ;
      
        if(!p) arret("LDUSKLformat_Create (4) : impossible d\'allouer la memoire") ;
    
        LDUSKLformat_GetPointerToLowerRow(a) = p ;
        LDUSKLformat_GetPointerToUpperColumn(a) = p + n_col ;
    
        LDUSKLformat_GetDiagonal(a) = z ;
    
        z += n_col ;
        LDUSKLformat_GetPointerToLowerRow(a)[0] = z ;
        LDUSKLformat_GetPointerToUpperColumn(a)[0] = z + nnz_l ;
      
        for(i = 1 ; i < n_col ; i++) {
          z = LDUSKLformat_GetPointerToLowerRow(a)[i - 1] ;
          LDUSKLformat_GetPointerToLowerRow(a)[i] = z + hc[i] ;
          LDUSKLformat_GetPointerToUpperColumn(a)[i] = z + nnz_l + hc[i] ;
        }
  
      }
    }

    free(hc) ;
  }

  return(a) ;
}



void LDUSKLformat_Delete(void* self)
{
  LDUSKLformat_t** a = (LDUSKLformat_t**) self ;
  free(LDUSKLformat_GetNonZeroValue(*a)) ;
  free(LDUSKLformat_GetPointerToLowerRow(*a)) ;
  free(*a) ;
  *a = NULL ;
}



#ifdef SLU_DIR
SuperMatrix_t* SuperMatrix_Create(Mesh_t* mesh)
/* Create a matrix in SuperMatrix format */
{
  int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
  NCformat_t*    asluNC = NCformat_Create(mesh) ;
  SuperMatrix_t* aslu = (SuperMatrix_t*) malloc(sizeof(SuperMatrix_t)) ;

  if(!aslu) assert(aslu) ;
  
  SuperMatrix_GetStorageType(aslu) = SLU_NC ;
  SuperMatrix_GetDataType(aslu) = SLU_D ;
  SuperMatrix_GetMatrixType(aslu) = SLU_GE ;
  SuperMatrix_GetNbOfColumns(aslu) = n_col ;
  SuperMatrix_GetNbOfRows(aslu) = n_col ;
  SuperMatrix_GetStorage(aslu) = (void*) asluNC ;
  
  return(aslu) ;
}



void SuperMatrix_Delete(void* self)
{
  SuperMatrix_t** aslu = (SuperMatrix_t**) self ;
  NCformat_t* asluNC = (NCformat_t*) SuperMatrix_GetStorage(*aslu) ;
  
  NCformat_Delete(&asluNC) ;
  
  free(*aslu) ;
  *aslu = NULL ;
}



NCformat_t* NCformat_Create(Mesh_t* mesh)
/** Create a matrix in SuperMatrix format */
{
  int n_el = Mesh_GetNbOfElements(mesh) ;
  int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  int*   colptr0 ;
  int*   colptr ;
  int*   rowind ;
  int    ie ;
  int    i,j ;
  int    nnz_max ;
  NCformat_t*    asluNC = (NCformat_t*) malloc(sizeof(NCformat_t)) ;

  if(!asluNC) assert(asluNC) ;


  /* tableau de travail */
  colptr0 = (int*) malloc((n_col + 1)*sizeof(int)) ;
  if(!colptr0) arret("NCformat_Create (1) : impossible d\'allouer la memoire") ;

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
    colptr = (int*) malloc((n_col+1)*sizeof(int)) ;
    
    if(!colptr) arret("NCformat_Create (2) : impossible d\'allouer la memoire") ;
    
    NCformat_GetFirstNonZeroValueIndexOfColumn(asluNC) = colptr ;
    
    rowind = (int*) malloc(nnz_max*sizeof(int)) ;
    if(!rowind) arret("NCformat_Create (3) : impossible d\'allouer la memoire") ;
    
    NCformat_GetRowIndexOfTheNonZeroValue(asluNC) = rowind ;
  }


  /* initialisation */
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


  /* compression de rowind */
  colptr[0] = 0 ;
  for(j = 0 , i = 0 ; i < n_col ; i++) {
    int k ;
    for(k = 0 ; k < colptr[i + 1] ; k++) {
      rowind[j++] = rowind[colptr0[i] + k] ;
    }
    colptr[i + 1] += colptr[i] ;
  }

  free(colptr0) ;


  if(nnz_max < colptr[n_col]) arret("NCformat_Create (4) : pb de memoire") ;
  NCformat_GetNbOfNonZeroValues(asluNC) = colptr[n_col] ;


  /* reallocation de la memoire */
  {
    int nnz = NCformat_GetNbOfNonZeroValues(asluNC) ;
    int* rowind1 = (int*) realloc(rowind,nnz*sizeof(int)) ;
    if(rowind1 != rowind) arret("NCformat_Create (5) : impossible d\'allouer la memoire") ;
  }
  

  /*  1. allocation de l'espace memoire pour la matice */
  {
    int nnz = NCformat_GetNbOfNonZeroValues(asluNC) ;
    double* nzval = (double*) malloc(nnz*sizeof(double)) ;
    
    if(!nzval) arret("SuperMatrix_Create (3) : impossible d\'allouer la memoire") ;
    
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
#endif



void LDUSKLformat_AssembleElementMatrix(LDUSKLformat_t* a,double* ke,int* cole,int* lige,int n)
/* Assemblage de la matrice elementaire ke dans la matrice globale k */
{
#define KE(i,j) (ke[(i)*n+(j)])
#define D(i)    (LDUSKLformat_GetDiagonal(a)[i])
#define U(i,j)  (LDUSKLformat_GetUpperColumn(a,j)[i])
#define L(i,j)  (LDUSKLformat_GetLowerRow(a,i)[j])
/*
#define U(i,j)  (*(LDUSKLformat_GetPointerToUpperColumn(a)[j] - j + i))
#define L(i,j)  (*(LDUSKLformat_GetPointerToLowerRow(a)[i] - i + j))
*/
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

#undef KE
#undef D
#undef U
#undef L
}


#ifdef SLU_DIR
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
      
      if(rowptr[irow] < 0) arret("(assemHBSLU : assemblage impossible") ;
      nzval[rowptr[irow]] += KE(ie,je) ;
    }

    for(i = colptr[jcol] ; i < colptr[jcol+1] ; i++) rowptr[rowind[i]] = -1 ;
  }

#undef KE
}
#endif


void LDUSKLformat_PrintMatrix(LDUSKLformat_t* a,unsigned int n,const char* keyword)
{
  double*  d = LDUSKLformat_GetDiagonal(a) ;
  double** u = LDUSKLformat_GetPointerToUpperColumn(a) ;
  double** l = LDUSKLformat_GetPointerToLowerRow(a) ;
  int nnz = LDUSKLformat_GetNbOfNonZeroValues(a) ;
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


#ifdef SLU_DIR
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
