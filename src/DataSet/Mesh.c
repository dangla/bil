#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include "Elements.h"
#include "Nodes.h"
#include "Geometry.h"
#include "Message.h"
#include "DataFile.h"
#include "Materials.h"
#include "Models.h"
#include "Periodicities.h"
#include "Graph.h"
#include "Solutions.h"
#include "Context.h"
#include "TextFile.h"
#include "Mesh.h"


#if defined(__cplusplus)
  extern "C" {
#endif

extern void   mc40ad_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*) ;
extern void   mc43ad_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*,int*) ;

#if defined(__cplusplus)
  }
#endif

static int*   (Mesh_ReadInversePermutationOfNodes)(DataFile_t*,int) ;
static int*   (Mesh_ComputeInversePermutationOfNodes)(Mesh_t*,const char*) ;
static int*   (Mesh_ComputeInversePermutationOfElements)(Mesh_t*,const char*) ;
static Graph_t*  (Mesh_CreateGraph)(Mesh_t*) ;


static void   mail0d(Mesh_t*) ;
static void   mail1dold(Mesh_t*,FILE*) ;
static void   mail1d(Mesh_t*,char*) ;
static void   maillage(double*,int*,double,int,Node_t*) ;
static int    mesh1d(double*,double*,int,Node_t**) ;
static void   lit_mail_m1d(Mesh_t*,const char*) ;
static void   lit_mail_gmsh(Mesh_t*,const char*) ;
static void   lit_mail_gmsh_1(Mesh_t*,const char*) ;
static void   lit_mail_gmsh_2(Mesh_t*,const char*) ;
static void   lit_mail_cesar(Mesh_t*,const char*) ;
static void   ecrit_mail_msh_1(Mesh_t*,const char*) ;
static void   ecrit_mail_msh_2(Mesh_t*,const char*) ;
static int    gmsh_ElmType(int,int) ;
static int    gmsh_NbNodes(int) ;
static int    gmsh_DimElement(int) ;





/* Extern functions */



Mesh_t*  Mesh_Create(DataFile_t* datafile,Materials_t* materials,Geometry_t* geometry)
{
  Mesh_t* mesh = Mesh_New() ;
  
  Mesh_GetGeometry(mesh) = geometry ;
  
  Mesh_GetElements(mesh) = Elements_New() ;
  
  Mesh_GetNodes(mesh) = Nodes_New() ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"MAIL,MESH,Mesh",",",1) ;
   
  Message_Direct("Enter in %s","Mesh") ;
  Message_Direct("\n") ;


  {
    char* line = DataFile_GetCurrentPositionInFileContent(datafile) ;
    //char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;

    /* 1. Allocation memory for the mesh i.e. 
     *    node, coordinates, element and node numbering
     * ------------------------------------------------*/
  
    if(!Mesh_Scan(mesh,line)) {
      Message_FatalError("Mesh_Create: No such file name") ;
    }
  }
  
  {
    char* filename = DataFile_GetFileName(datafile) ;
    
    ecrit_mail_msh_2(mesh,filename) ;
  }
  
  DataFile_CloseFile(datafile) ;

  /* 2. Allocation memory for unknown and equation positions
   * -------------------------------------------------------*/
  Elements_CreateMore(Mesh_GetElements(mesh),materials) ;

  /* 3. Allocation memory for names of equations and unknowns
   * --------------------------------------------------------*/
  Mesh_CreateEquationContinuity(mesh) ;
  
  /* 4. Allocation memory for matrix row and column indexes
   * ------------------------------------------------------*/
  Nodes_CreateMore(Mesh_GetNodes(mesh)) ;
  
  //Elements_DefineProperties(Mesh_GetElements(mesh)) ;


  return(mesh) ;
}



char*  Mesh_Scan(Mesh_t* mesh,char* line)
{
  char  nom_mail[Mesh_MaxLengthOfFileName] ;

  /* Mesh file ? */
  sscanf(line,"%s",nom_mail) ;
  
  /* Treatment after filename extension */
  if(strstr(nom_mail,".msh")) {
    
    lit_mail_gmsh(mesh,nom_mail) ;
    
  } else if(strstr(nom_mail,".m1d")) {
    
    lit_mail_m1d(mesh,nom_mail) ;
    
  } else if(strstr(nom_mail,".ces")) {
    
    lit_mail_cesar(mesh,nom_mail) ;
    
  } else {
    int dim = Mesh_GetDimension(mesh) ;
    
    if(dim == 0) {
      
      mail0d(mesh) ;
      
    } else if(dim == 1) {
      
      /* Read directly in the data file */
      mail1d(mesh,line) ;
      
    } else {
      
      return(NULL) ;
      
    }
  }
  
  return(line) ;
}



Graph_t* Mesh_CreateGraph(Mesh_t* mesh)
/** Create the graph matrix of node indexes.
 *  Return a pointer to Graph_t.
 **/
{
  Graph_t* graph ;
  
  
  /* Create a graph */
  {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    /* Nb of connections per node (useful to size graph) */
    int*   nnz_no = (int*) calloc(n_no,sizeof(int)) ;
  
    if(!nnz_no) {
      arret("Mesh_CreateGraph(1): impossible d\'allouer la memoire") ;
    }
  
  
    /* An overestimation of nnz_no */
    {
      int n_el = Mesh_GetNbOfElements(mesh) ;
      Element_t* el = Mesh_GetElement(mesh) ;
      int ie ;
    
      for(ie = 0 ; ie < n_el ; ie++) {
        int nn = Element_GetNbOfNodes(el + ie) ;
        int i ;
    
        for(i = 0 ; i < nn ; i++) {
          Node_t* node = Element_GetNode(el + ie,i) ;
          int k = Node_GetNodeIndex(node) ;
      
          nnz_no[k] += nn - 1 ;
        }
      }
    }
    
    graph = Graph_Create(n_no,nnz_no) ;

    free(nnz_no) ;
  }


  /* Compute the graph */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int ie ;

    for(ie = 0 ; ie < n_el ; ie++) {
      int nn = Element_GetNbOfNodes(el + ie) ;
      int i ;
    
      for(i = 0 ; i < nn ; i++) {
        Node_t* node_i = Element_GetNode(el + ie,i) ;
        int  in = Node_GetNodeIndex(node_i) ;
        int* listin = Graph_GetNeighborOfVertex(graph,in) ;
        int j ;
      
        for(j = 0 ; j < nn ; j++) {
          Node_t* node_j = Element_GetNode(el + ie,j) ;
          int  jn = Node_GetNodeIndex(node_j) ;
          int* listjn = Graph_GetNeighborOfVertex(graph,jn) ;
          int  degrjn = Graph_GetDegreeOfVertex(graph,jn) ;
          int  degrin = Graph_GetDegreeOfVertex(graph,in) ;
        
	        if(in == jn) continue ;
        
          /* Has jn been already met? */
          {
            int met = 0 ;
            int k ;
            
	          for(k = 0 ; k < degrin ; k++) {
	            if(jn == listin[k]) {met = 1 ; break ;}
	          }
        
	          if(met) continue ;
          }
        
          /* Not already met. So we increment with jn and in */
          if(listin[degrin] < 0) {
            Graph_GetDegreeOfVertex(graph,in) += 1 ;
	          listin[degrin] = jn ;
          } else {
            arret("Mesh_CreateGraph(3): not enough space") ;
          }
          
          if(listjn[degrjn] < 0) {
            Graph_GetDegreeOfVertex(graph,jn) += 1 ;
	          listjn[degrjn] = in ;
          } else {
            arret("Mesh_CreateGraph(4): not enough space") ;
          }
          
        }
      }
    }
  }
  
    
  /* Update the graph for periodic mesh/BC */
  {
    Periodicities_t* periodicities = Mesh_GetPeriodicities(mesh) ;
  
    if(periodicities) {
      Periodicities_UpdateGraph(mesh,graph) ;
    }
  }
  
  
  /* Nb of edges */
  Graph_UpdateTheNbOfEdges(graph) ;

  
  return(graph) ;
}



void  Mesh_SetMatrixPermutationNumbering(Mesh_t* mesh,BConds_t* bconds,DataFile_t* datafile)
/** Set matrix row (column) index which node equation (unknown) points to
 *  by using the inverse permutation vector 
 **/
{

  /* Reset to 0 */
  Mesh_ResetMatrixNumbering(mesh) ;

  /* Accounting for BC (set indexes to -1) */
  BConds_ResetMatrixNumbering(bconds,mesh) ;

  /* Accounting for periodic mesh/BC (set indexes of slave nodes to -1) */
  Periodicities_ResetMatrixPermutationNumbering(mesh) ;

  /* We set up the system */
  {
    int   n_no = Mesh_GetNbOfNodes(mesh) ;
    Node_t* no = Mesh_GetNode(mesh) ;
    int*  perm = Mesh_ReadInversePermutationOfNodes(datafile,n_no) ;
    int NbOfMatrixRows = 0 ;
    int NbOfMatrixColumns = 0 ;
    int    i ;
    
    
    /* In the order of the permuted nodes */
    for(i = 0 ; i < n_no ; i++) {
      Node_t* node_i = no + perm[i] ;
      int   j ;

      for(j = 0 ; j < Node_GetNbOfUnknowns(node_i) ; j++) {
        int icol = Node_GetMatrixColumnIndex(node_i)[j] ;
        
        if(icol == 0) {
          Node_GetMatrixColumnIndex(node_i)[j] = NbOfMatrixColumns ;
          NbOfMatrixColumns += 1 ;
        }
      }

      for(j = 0 ; j < Node_GetNbOfEquations(node_i) ; j++) {
        int irow = Node_GetMatrixRowIndex(node_i)[j] ;
        
        if(irow == 0) {
          Node_GetMatrixRowIndex(node_i)[j] = NbOfMatrixRows ;
          NbOfMatrixRows += 1 ;
        }
      }
    }
  
    free(perm) ;
  
  
    if(NbOfMatrixColumns != NbOfMatrixRows) {
      arret("Mesh_SetMatrixPermutationNumbering(3)") ;
    }
  
    Nodes_GetNbOfMatrixRows(Mesh_GetNodes(mesh)) = NbOfMatrixRows ;
    Nodes_GetNbOfMatrixColumns(Mesh_GetNodes(mesh)) = NbOfMatrixColumns ;
  }
  
  
  /* Update indexes of slave nodes for periodic mesh */
  Periodicities_UpdateMatrixPermutationNumbering(mesh) ;
}



void  Mesh_ResetMatrixNumbering(Mesh_t* mesh)
/** Initialization to 0 the Matrix Row/Column Indexes */
{

  /* Initialization to arbitrarily negative value (-1) 
   * so as to eliminate dof of isolated nodes or
   * dof of negative position (see below).
   */
  {
    Node_t* no = Mesh_GetNode(mesh) ;
    int n_no = Mesh_GetNbOfNodes(mesh) ;
    int    i ;
    
    for(i = 0 ; i < n_no ; i++) {
      Node_t* node_i = no + i ;
      int   j ;
    
      for(j = 0 ; j < Node_GetNbOfUnknowns(node_i) ; j++)  {
        Node_GetMatrixColumnIndex(node_i)[j] = -1 ;
      }
    
      for(j = 0 ; j < Node_GetNbOfEquations(node_i) ; j++) {
        Node_GetMatrixRowIndex(node_i)[j] = -1 ;
      }
    }
  }


  /* Initialization to 0 for nodes belonging to elements */
  {
    Element_t* el = Mesh_GetElement(mesh) ;
    int n_el = Mesh_GetNbOfElements(mesh) ;
    int    ie ;
  
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t* elt_i = el + ie ;
      Material_t* mat = Element_GetMaterial(elt_i) ;
    
      if(mat) {
        int    nn  = Element_GetNbOfNodes(elt_i) ;
        int    neq = Element_GetNbOfEquations(elt_i) ;
        int i ;

        for(i = 0 ; i < nn ; i++) {
          Node_t* node_i = Element_GetNode(elt_i,i) ;
          int   j ;

          for(j = 0 ; j < neq ; j++) {
            int ij = i*neq + j ;
            int jj ;
          
            /*  columns (unknowns) */
            if((jj = Element_GetUnknownPosition(elt_i)[ij]) >= 0) {
              Node_GetMatrixColumnIndex(node_i)[jj] = 0 ;
            }
          
            /*  rows (equations) */
            if((jj = Element_GetEquationPosition(elt_i)[ij]) >= 0) {
              Node_GetMatrixRowIndex(node_i)[jj] = 0 ;
            }
          }
        }
      }
    }
  }

}


void Mesh_WriteGraph(Mesh_t* mesh,const char* nom,const char* format)
/** Create the graph file "nom.graph" in the format "format"
 */
{
  Graph_t*  graph = Mesh_CreateGraph(mesh) ;
  FILE*  fic_graph ;
  
  
  /* Open the file "nom.graph" */
  {
    char   nom_graph[Mesh_MaxLengthOfFileName] ;
    
    if(strlen(nom) + 6 > Mesh_MaxLengthOfFileName) {
      arret("Mesh_WriteGraph(1): too long name") ;
    }
  
    sprintf(nom_graph,"%s.graph",nom) ;
  
    printf("Creation of the graph file %s\n",nom_graph) ;
  
    fic_graph = fopen(nom_graph,"w") ;
  
    if(!fic_graph) {
      arret("Mesh_WriteGraph(2): cannot open %s\n",nom_graph) ;
    }
  }


  /* Write the first line */
  {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    int    nnz  = Graph_GetNbOfEdges(graph) ;
    
    fprintf(fic_graph,"%d %d\n",n_no,nnz) ;
  }


	/* Format HSL_MC40 */
  if(!strcmp(format,"hsl") || !strcmp(format,"hsl_mc40")) {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    int in ;
    
    for(in = 0 ; in < n_no ; in++) {
      int  degree = Graph_GetDegreeOfVertex(graph,in) ;
      int* list   = Graph_GetNeighborOfVertex(graph,in) ;
      int i ;
      
      for(i = 0 ; i < degree ; i++) {
        int jn = list[i] ;
        
	      if(jn < in) fprintf(fic_graph,"%d %d\n",in+1,jn+1) ;
      }
    }
    
	/* Format METIS */
  } else if(!strcmp(format,"metis")) {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    int in ;
    
    for(in = 0 ; in < n_no ; in++) {
      int  degree = Graph_GetDegreeOfVertex(graph,in) ;
      int* list   = Graph_GetNeighborOfVertex(graph,in) ;
      int i ;
      
      for(i = 0 ; i < degree ; i++) {
        int jn = list[i] ;
        
        fprintf(fic_graph,"%d ",jn+1) ;
      }
      fprintf(fic_graph,"\n") ;
    }
    
  } else {
    arret("Mesh_WriteGraph: non prevu") ;
  }

  /* fermeture du fichier */
  fclose(fic_graph) ;
  
  Graph_Delete(&graph) ;
}


void   Mesh_WriteInversePermutation(Mesh_t* mesh,const char* nom,const char* format)
/** Create the inverse permutation file "nom.graph.iperm" */
{
  char   nom_iperm[Mesh_MaxLengthOfFileName] ;

  if(strlen(nom) + 12 > Mesh_MaxLengthOfFileName) {
    arret("Mesh_WriteInversePermutation(6): too long name") ;
  }
    
  sprintf(nom_iperm,"%s.graph.iperm",nom) ;
  printf("Creation of the Graph.iperm file = %s\n",nom_iperm) ;
    
  {
    FILE*  fic_iperm = fopen(nom_iperm,"w") ;
  
    if(!fic_iperm) {
      arret("Mesh_WriteInversePermutation(7): can't open the file\n") ;
    }
      
    if(!strcmp(format,"hsl") || !strcmp(format,"hsl_mc40")) {
      int    n_no = Mesh_GetNbOfNodes(mesh) ;
      int*   iperm = Mesh_ComputeInversePermutationOfNodes(mesh,format) ;
      int i ;

      for(i = 0 ; i < n_no ; i++) {
        fprintf(fic_iperm,"%d\n",iperm[i] - 1) ;
      }
        
      free(iperm) ;
    } else if(!strcmp(format,"hsl_mc43")) {
      int    nelt = Mesh_GetNbOfNodes(mesh) ;
      int*   norder = Mesh_ComputeInversePermutationOfElements(mesh,format) ;
      int i ;

      for(i = 0 ; i < nelt ; i++) {
        fprintf(fic_iperm,"%d\n",norder[i] - 1) ;
      }
      
      free(norder) ;
    }
  
    fclose(fic_iperm) ;
  }
}


int*   Mesh_ComputeInversePermutationOfNodes(Mesh_t* mesh,const char* format)
{
  int    n_no = Mesh_GetNbOfNodes(mesh) ;
  int*   iperm = (int*) malloc(n_no*sizeof(int)) ;
  
  if(!iperm) {
    arret("Mesh_ComputeInversePermutationOfNodes(1): not enough memory") ;
  }
  
  
  if(!strcmp(format,"hsl") || !strcmp(format,"hsl_mc40")) {
    Graph_t*  graph = Mesh_CreateGraph(mesh) ;
    int    nnz  = Graph_GetNbOfEdges(graph) ;
    int*   irn = (int*) malloc(2*nnz*sizeof(int)) ;
    int*   jcn = (int*) malloc(nnz*sizeof(int)) ;
  
    if(!irn || !jcn) {
      arret("Mesh_ComputeInversePermutationOfNodes(2): not enough memory") ;
    }

    /* Compute the row and column indexes: irn and jcn */
    {
      int  in ;
      int  k ;
  
      for(k = 0 , in = 0 ; in < n_no ; in++) {
        int  degrin = Graph_GetDegreeOfVertex(graph,in) ;
        int* listin = Graph_GetNeighborOfVertex(graph,in) ;
        int i ;
    
        for(i = 0 ; i < degrin ; i++) {
          int jn = listin[i] ;
      
          if(jn < in) {
            irn[k] = in + 1 ;
	          jcn[k] = jn + 1 ;
	          k++;
          }
        }
      }
    }
  
    Graph_Delete(&graph) ;
    
    /* From the package HSL_MC40 */
    {
      int   iflag ;
      int   itype = 1 ;
      int   iprof[2] ;
      int*  icptr = (int*) malloc((n_no + 1)*sizeof(int)) ;
      int*  iw = (int*) malloc((3*n_no + 2)*sizeof(int)) ;
    
      if(!icptr || !iw) {
        arret("Mesh_ComputeInversePermutationOfNodes(4): not enough memory") ;
      }
      
      mc40ad_(&itype,(int*) &n_no,&nnz,irn,jcn,icptr,iperm,iw,iprof,&iflag) ;
      
      if(iflag < 0) {
        Message_FatalError("Mesh_ComputeInversePermutationOfNodes: something wrong happened in mc40ad\n") ;
        Message_FatalError("iflag = %d\n",iflag) ;
      }
      
      Message_Direct("Profile of the matrix:\n") ;
      Message_Direct("    original ordering  %d\n",iprof[0]) ;
      Message_Direct("    new ordering       %d\n",iprof[1]) ;
      /*
      Message_Direct("The profile of the matrix is the number of coefficients\n") ;
      Message_Direct("of the lower triangle when any zero ahead of the first \n") ;
      Message_Direct("entry in its row is excluded.\n") ;
      */
    
      free(icptr) ;
      free(iw) ;
    }
  
    free(irn) ;
    free(jcn) ;
  } else {
    arret("Mesh_ComputeInversePermutationOfNodes(5): format %s unknown",format) ;
  }
  

  return(iperm) ;
}


int*   Mesh_ComputeInversePermutationOfElements(Mesh_t* mesh,const char* format)
{
  int    nelt = Mesh_GetNbOfNodes(mesh) ;
  int*   norder = (int*) malloc(nelt*sizeof(int)) ;
  
  if(!norder) {
    arret("Mesh_ComputeInversePermutationOfElements(1): not enough memory") ;
  }
  
  
  if(!strcmp(format,"hsl_mc43")) {
    
    /* From the package HSL_MC43 */
    {
      /* If icntl = 0 the direct element reordering algorithm is employed; 
       * If icntl = 1 the indirect element reordering algorithm is employed
       */
      int  icntl = 1 ;
      int  nno = Mesh_GetNbOfNodes(mesh) ;
      int  n = nno ;
      Element_t* elt = Mesh_GetElement(mesh) ;
      int  nz ;
      int* eltptr ;
      int* eltvar ;
      int  liw ;
      int* iw ;
      int  mxwave[2] ;
      int  iflag ;
      
      eltptr = (int*) malloc((nelt + 1)*sizeof(int)) ;
      
      if(!eltptr) {
        arret("Mesh_ComputeInversePermutationOfElements(4): not enough memory") ;
      }
      
      {
        int i ;
        
        eltptr[0] = 1 ;
        for(i = 0 ; i < nelt ; i++) {
          int nn = Element_GetNbOfNodes(elt + i) ;
        
          eltptr[i + 1] = eltptr[i] + nn ;
        }
      }
      
      nz = eltptr[nelt] - 1 ;
      
      eltvar = (int*) malloc(nz*sizeof(int)) ;
      
      if(!eltvar) {
        arret("Mesh_ComputeInversePermutationOfElements(4): not enough memory") ;
      }
      
      {
        int i ;
        int k = 0 ;
        
        for(i = 0 ; i < nelt ; i++) {
          int nn = Element_GetNbOfNodes(elt + i) ;
          int j ;
        
          for(j = 0 ; j < nn ; j++) {
            Node_t* no = Element_GetNode(elt + i,j) ;
            int in = Node_GetNodeIndex(no) ;
            
            eltvar[k] = in + 1 ;
            k++ ;
          }
        }
      }
      
      {
        if(icntl == 1) {
          int nnmax = 0 ;
          int i ;
        
          for(i = 0 ; i < nelt ; i++) {
            int nn = Element_GetNbOfNodes(elt + i) ;
        
            if(nn > nnmax) nnmax = nn ;
          }
          
          liw = 3*(n + nz + 1) + nelt*(nnmax*nnmax + 1) ;
        
          if(liw < nelt*64) liw = nelt*64 ;
        } else if(icntl == 0) {
          int maxel = 10 ; /* Need to be calculated */
          
          liw = 3*nz + 2*nelt*(maxel + 3) + n + 3 ;
          
          if(liw < 2*n) liw = 2*n ;
        }
        
        iw = (int*) malloc(liw*sizeof(int)) ;
        
        if(!iw) {
          arret("Mesh_ComputeInversePermutationOfElements(4): not enough memory") ;
        }
      }
      
      
      mc43ad_(&icntl,&n,&nelt,&nz,eltvar,eltptr,norder,&liw,iw,mxwave,&iflag) ;
      
      if(iflag < 0) {
        Message_FatalError("Mesh_ComputeInversePermutationOfElements: something wrong happened in mc43ad\n") ;
        Message_FatalError("iflag = %d\n",iflag) ;
      }
      
      Message_Direct("Maximum wavefront of the matrix:\n") ;
      Message_Direct("    original ordering  %d\n",mxwave[0]) ;
      Message_Direct("    new ordering       %d\n",mxwave[1]) ;
      
      free(eltptr) ;
      free(eltvar) ;
      free(iw) ;
    }
  } else {
    arret("Mesh_ComputeInversePermutationOfElements(5): format %s unknown",format) ;
  }
  

  return(norder) ;
}




void Mesh_CreateEquationContinuity(Mesh_t* mesh)
/** Compute some informations needed to describe the continuity of 
 *  equations and unknowns at nodes, i.e.:
 *  - the number and names of equations at nodes
 *  - the position of equations and unknowns of elements at nodes
 */
{
  Element_t* el = Mesh_GetElement(mesh) ;
  Node_t* no = Mesh_GetNode(mesh) ;
  int n_no = Mesh_GetNbOfNodes(mesh) ;
  int n_el = Mesh_GetNbOfElements(mesh) ;
  int    ie,in ;

  /* Allocate memory space for names of equations and unknwons */
  {
    int    np = 2*n_no*Model_MaxNbOfEquations ;
    char** p  = (char**) malloc(np*sizeof(char*)) ;
    
    if(!p) {
      arret("Mesh_SetEquationContinuity(1): memory allocation impossible") ;
    }
  
    for(in = 0 ; in < n_no ; in++) {
      Node_GetNameOfUnknown(no + in)  = p + in*Model_MaxNbOfEquations ;
      Node_GetNameOfEquation(no + in) = p + (in + n_no)*Model_MaxNbOfEquations ;
    }
  }

  /** Compute
   *  - the number of equations per node: Node_GetNbOfEquations(no)
   *  - the number of unknowns per node:  Node_GetNbOfUnknowns(no)
   *  - the position of equations of each element at nodes: Element_GetEquationPosition(el)
   *  - the position of unknowns of each element at nodes:  Element_GetUnknownPosition(el)
   */
  for(in = 0 ; in < n_no ; in++) {
    Node_GetNbOfUnknowns(no + in) = 0 ;
    Node_GetNbOfEquations(no + in) = 0 ;
  }

  for(ie = 0 ; ie < n_el ; ie++) {
    Material_t* mat = Element_GetMaterial(el + ie) ;
    
    if(mat) {
      int nn = Element_GetNbOfNodes(el + ie) ;
      int neq = Element_GetNbOfEquations(el + ie) ;
      char** mat_name_unk = Material_GetNameOfUnknown(mat) ;
      char** mat_name_eqn = Material_GetNameOfEquation(mat) ;
      Node_t** pnode = Element_GetPointerToNode(el + ie) ;
      short int*  unk_pos = Element_GetUnknownPosition(el + ie) ;
      short int*  eqn_pos = Element_GetEquationPosition(el + ie) ;
      int i ;
      
      for(i = 0 ; i < nn ; i++) {
        Node_t* node_i = pnode[i] ;
        char** node_name_unk = Node_GetNameOfUnknown(node_i) ;
        char** node_name_eqn = Node_GetNameOfEquation(node_i) ;
        int ieq ;
        
        for(ieq = 0 ; ieq < neq ; ieq++) {
          int ij = i*neq + ieq ;
          int jun ;
          int jeq ;

          /* Continuity of unknowns: 
           * - number of unknowns: Node_GetNbOfUnknowns(node_i), 
           * - name of unknowns: Node_GetNameOfUnknown(node_i) 
           * - position of unknowns at nodes: Element_GetUnknownPosition(el + ie)
           */
          for(jun = 0 ; jun < Node_GetNbOfUnknowns(node_i) ; jun++) {
            if(!strcmp(node_name_unk[jun],mat_name_unk[ieq])) break ;
          }
          
          if(unk_pos[ij] >= 0) { /* if < 0 no unknown at this node */
            unk_pos[ij] = jun ;
            if(jun == Node_GetNbOfUnknowns(node_i)) {
              Node_GetNbOfUnknowns(node_i) += 1 ;
              node_name_unk[jun] = mat_name_unk[ieq] ;
            }
          }

          /* Continuity of equations:
           * - number of equations: Node_GetNbOfEquations(node_i)
           * - name of equations: Node_GetNameOfEquation(node_i)
           * - position of equations: Element_GetEquationPosition(el + ie)
           */
          for(jeq = 0 ; jeq < Node_GetNbOfEquations(node_i) ; jeq++) {
            if(!strcmp(node_name_eqn[jeq],mat_name_eqn[ieq])) break ;
          }
          
          if(eqn_pos[ij] >= 0) { /* if < 0 no equation at this node */
            eqn_pos[ij] = jeq ;
            if(jeq == Node_GetNbOfEquations(node_i)) {
              Node_GetNbOfEquations(node_i) += 1 ;
              node_name_eqn[jeq] = mat_name_eqn[ieq] ;
            }
          }
        }
      }
    }
  }


  /* Checking */
  for(in = 0 ; in < n_no ; in++) {
    if(Node_GetNbOfUnknowns(no + in) != Node_GetNbOfEquations(no + in)) {
      arret("Mesh_SetEquationContinuity(2)") ;
    }
  }


  /* Compression of the memory space */
  {
    int    n_u = 0 ;
    char** p = Node_GetNameOfUnknown(no) ;
    
    /* Compress the memory space for Node_GetNameOfUnknown(no + in) */
    for(in = 0 ; in < n_no ; in++) {
      char** node_name_unk = Node_GetNameOfUnknown(no + in) ;
      int j ;
    
      for(j = 0 ; j < Node_GetNbOfUnknowns(no + in) ; j++) {
        p[j] = node_name_unk[j] ;
      }
    
      Node_GetNameOfUnknown(no + in) = p ;
      p   += Node_GetNbOfUnknowns(no + in) ;
      n_u += Node_GetNbOfUnknowns(no + in) ;
    }
  
    /* Compress the memory space for Node_GetNameOfEquation(no + in) */
    for(in = 0 ; in < n_no ; in++) {
      char** node_name_eqn = Node_GetNameOfEquation(no + in) ;
      int j ;
    
      for(j = 0 ; j < Node_GetNbOfEquations(no + in) ; j++) {
        p[j] = node_name_eqn[j] ;
      }
    
      Node_GetNameOfEquation(no + in) = p ;
      p   += Node_GetNbOfEquations(no + in) ;
      n_u += Node_GetNbOfEquations(no + in) ;
    }
  
  
    /* To avoid that realloc crashes when performing tests ! */
    if(n_u == 0) n_u = 1 ; 
  
  
    /* We shrink the allocated memory */
    p = (char**) realloc(Node_GetNameOfUnknown(no),n_u*sizeof(char*)) ;
    
    if(!p) {
      arret("Mesh_SetEquationContinuity(3): memory allocation impossible") ;
    }
    
    if(Node_GetNameOfUnknown(no) != p) {
      /* arret("Mesh_SetEquationContinuity(4): memory allocation impossible") ; */
      Message_Warning("Mesh_SetEquationContinuity(4): new memory allocation") ;
      /* Compress the memory space for Node_GetNameOfUnknown(no + in) */
      for(in = 0 ; in < n_no ; in++) {
        Node_GetNameOfUnknown(no + in) = p ;
        p   += Node_GetNbOfUnknowns(no + in) ;
      }
  
      /* Compress the memory space for Node_GetNameOfEquation(no + in) */
      for(in = 0 ; in < n_no ; in++) {
        Node_GetNameOfEquation(no + in) = p ;
        p   += Node_GetNbOfEquations(no + in) ;
      }
    }
  }
}



void (Mesh_InitializeSolutionPointers)(Mesh_t* mesh,Solutions_t* sols)
/** Initialize the pointers of nodes and elements to
 ** the pointers of the solution sol. */
{
  Solution_t* sol = Solutions_GetSolution(sols) ;
  int n_no = Mesh_GetNbOfNodes(mesh) ;
  Node_t* no = Mesh_GetNode(mesh) ;
  int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  int    i ;

  for(i = 0 ; i < n_el ; i++) {
    Element_GetElementSol(el + i) = Solution_GetElementSol(sol) + i ;
  }
  
  for(i = 0 ; i < n_no ; i++) {
    Node_GetNodeSol(no + i) = Solution_GetNodeSol(sol) + i ;
  }

}



int (Mesh_LoadCurrentSolution)(Mesh_t* mesh,DataFile_t* datafile,double* t)
/** Load the solution from a continuous file (suffix "cont" or "conti"). 
 ** Return either i > 0 if a continuous file was found and a solution
 ** was loaded from it or 0 if no continuous file was found. */
{
  int ires = 0 ;
  FILE* fic_cont ;

  /* ires = 0 no continuous file is found
   * ires > 0 a continuous file is found */
  {
    char* nom = DataFile_GetFileName(datafile) ;
    char   nom_cont[Mesh_MaxLengthOfFileName] ;

    if(strlen(nom) + 6 > Mesh_MaxLengthOfFileName) {
      arret("Mesh_LoadCurrentSolution(1): too long name") ;
    }
    
    /* From the highest to the lowest priority */
    if(ires == 0) {
      sprintf(nom_cont,"%s.cont",nom) ;
      {
        fic_cont = fopen(nom_cont,"rb") ;
    
        if(fic_cont) {
          ires = 2 ;
          /* Set initialization context to no initialization */
          DataFile_ContextSetToNoInitialization(datafile) ;
        }
      }
    }
    
    if(ires == 0) {
      sprintf(nom_cont,"%s.conti",nom) ;
      {
        fic_cont = fopen(nom_cont,"rb") ;
    
        if(fic_cont) {
          ires = 1 ;
          /* Set initialization context to partial initialization */
          DataFile_ContextSetToPartialInitialization(datafile) ;
        }
      }
    }
  }
  

  if(ires == 0) return(ires) ;

  /* Time */
  fscanf(fic_cont,"%lf",t) ;
  
  /* Unknowns */
  {
    int n_no = Mesh_GetNbOfNodes(mesh) ;
    Node_t* no = Mesh_GetNode(mesh) ;
    int    i ;
    
    for(i = 0 ; i < n_no ; i++) {
      double* u = Node_GetCurrentUnknown(no + i) ;
      int    nb_unk = Node_GetNbOfUnknowns(no + i) ;
      int j ;
    
      for(j = 0 ; j < nb_unk ; j++) {
        fscanf(fic_cont,"%lf",u + j) ;
      }
    }
  }
  
  
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int    i ;
    
    /* Implicit terms */
    for(i = 0 ; i < n_el ; i++) {
      double* vi = Element_GetCurrentImplicitTerm(el + i) ;
      int   n_vi = Element_GetNbOfImplicitTerms(el + i) ;
      int j ;
      
      for(j = 0 ; j < n_vi ; j++) {
        fscanf(fic_cont,"%lf",vi + j) ;
      }
    }
  
    /* Explicit terms */
    for(i = 0 ; i < n_el ; i++) {
      double* ve = Element_GetCurrentExplicitTerm(el + i) ;
      int   n_ve = Element_GetNbOfExplicitTerms(el + i) ;
      int j ;
      
      for(j = 0 ; j < n_ve ; j++) {
        fscanf(fic_cont,"%lf",ve + j) ;
      }
    }
  
    /* Constant terms */
    for(i = 0 ; i < n_el ; i++) {
      double* v0 = Element_GetConstantTerm(el + i) ;
      int   n_v0 = Element_GetNbOfConstantTerms(el + i) ;
      int j ;
      
      for(j = 0 ; j < n_v0 ; j++) {
        fscanf(fic_cont,"%lf",v0 + j) ;
      }
    }
  }

  fclose(fic_cont) ;
  
  return(ires) ;
}



int (Mesh_StoreCurrentSolution)(Mesh_t* mesh,DataFile_t* datafile,double t)
/** Store the solution to a storage file (suffix "sto"). */
{
  FILE* fic_sto ;


  {
    char* nom = DataFile_GetFileName(datafile) ;
    char   nom_sto[Mesh_MaxLengthOfFileName] ;

    if(strlen(nom) + 4 > Mesh_MaxLengthOfFileName) {
      arret("Mesh_StoreCurrentSolution(1): too long name") ;
    }
    
    sprintf(nom_sto,"%s.sto",nom) ;
    
    fic_sto = fopen(nom_sto,"rb") ;
    if(fic_sto) {
      Message_Direct("%s has been replaced\n",nom_sto) ;
    }
    
    fic_sto = fopen(nom_sto,"wb") ;
    if(!fic_sto) {
      Message_FatalError("error while opening %s\n",nom_sto) ;
    }
  }

  /* Time */
  fprintf(fic_sto,"%e ",t) ;
  fprintf(fic_sto,"\n") ;
  
  
  /* Unknowns */
  {
    int n_no = Mesh_GetNbOfNodes(mesh) ;
    Node_t* no = Mesh_GetNode(mesh) ;
    int    i ;
    
    for(i = 0 ; i < n_no ; i++) {
      double* u = Node_GetCurrentUnknown(no + i) ;
      int    nb_unk = Node_GetNbOfUnknowns(no + i) ;
      int j ;
      
      for(j = 0 ; j < nb_unk ; j++) {
        fprintf(fic_sto,"%.12e ",u[j]) ;
      }
      
      if(nb_unk) fprintf(fic_sto,"\n") ;
    }
  }
  
  
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
    int    i ;
    
    /* Implicit terms */
    for(i = 0 ; i < n_el ; i++) {
      double* vi = Element_GetCurrentImplicitTerm(el + i) ;
      int   n_vi = Element_GetNbOfImplicitTerms(el + i) ;
      int j ;
      
      for(j = 0 ; j < n_vi ; j++) {
        fprintf(fic_sto,"%.12e ",vi[j]) ;
      }
      
      if(n_vi) fprintf(fic_sto,"\n") ;
    }
    
    /* Explicit terms */
    for(i = 0 ; i < n_el ; i++) {
      double* ve = Element_GetCurrentExplicitTerm(el + i) ;
      int   n_ve = Element_GetNbOfExplicitTerms(el + i) ;
      int j ;
      
      for(j = 0 ; j < n_ve ; j++) {
        fprintf(fic_sto,"%.12e ",ve[j]) ;
      }
      
      if(n_ve) fprintf(fic_sto,"\n") ;
    }
    
    /* Constant terms */
    for(i = 0 ; i < n_el ; i++) {
      double* v0 = Element_GetConstantTerm(el + i) ;
      int   n_v0 = Element_GetNbOfConstantTerms(el + i) ;
      int j ;
      
      for(j = 0 ; j < n_v0 ; j++) {
        fprintf(fic_sto,"%.12e ",v0[j]) ;
      }
      
      if(n_v0) fprintf(fic_sto,"\n") ;
    }
  }

  fclose(fic_sto) ;
  
  return(0) ;
}



void (Mesh_SetCurrentUnknownsWithBoundaryConditions)(Mesh_t* mesh,BConds_t* bconds,double t)
/* Set the current values.. */
{
  /*
	  1. .. with previous ones
  */
  {
    unsigned int   nb_nodes = Mesh_GetNbOfNodes(mesh) ;
    Node_t*        node = Mesh_GetNode(mesh) ;
    unsigned int   i ;
  
    for(i = 0 ; i < nb_nodes ; i++) {
      int  nb_unk = Node_GetNbOfUnknowns(node + i) ;
      double* u_n = Node_GetPreviousUnknown(node + i) ;
      double* u_1 = Node_GetCurrentUnknown(node + i) ;
      int j ;
        
      for(j = 0 ; j < nb_unk ; j++) {
        u_1[j] = u_n[j] ;
      }
    }
  }
      
  /*
    2. .. and with the boundary conditions
  */
  BConds_AssignBoundaryConditions(bconds,mesh,t) ;
}



void (Mesh_UpdateCurrentUnknowns)(Mesh_t* mesh,Solver_t* solver)
{
  double* x = Solver_GetSolution(solver) ;
  Nodes_t* nodes = Mesh_GetNodes(mesh) ;
  Node_t* node = Mesh_GetNode(mesh) ;
  ObVals_t* obvals = Nodes_GetObjectiveValues(nodes) ;
  ObVal_t* obval = ObVals_GetObVal(obvals) ;
  unsigned int nb_nodes = Mesh_GetNbOfNodes(mesh) ;
  unsigned int i ;
          
  for(i = 0 ; i < nb_nodes ; i++) {
    int nin = Node_GetNbOfUnknowns(node + i) ;
    double* u_1 = Node_GetCurrentUnknown(node + i) ;
    int j ;
            
    for(j = 0 ; j < nin ; j++) {
      int   k = Node_GetMatrixColumnIndex(node + i)[j] ;
      int  iobval = Node_GetObValIndex(node + i)[j] ;
      double rfac = ObVal_GetRelaxationFactor(obval + iobval) ;
              
      if(k >= 0) u_1[j] += rfac * x[k] ;
    }
  }
}


/* Intern functions */



#define MESH          (mesh)
#define GEOMETRY      Mesh_GetGeometry(MESH) 
#define DIM           Mesh_GetDimension(MESH)
#define SYMMETRY      Mesh_GetSymmetry(MESH) 
#define COORSYS       Mesh_GetCoordinateSystem(MESH) 
#define N_NO          Mesh_GetNbOfNodes(MESH)
#define NO            Mesh_GetNode(MESH)
#define N_EL          Mesh_GetNbOfElements(MESH)
#define EL            Mesh_GetElement(MESH)



void mail0d(Mesh_t* mesh)
/* MAILLAGE 0D */
{
  /* nombre d'elements */
  N_EL = 1 ;

  /* nombre de noeuds */
  N_NO = 1 ;

  /* les noeuds */
  NO = (Node_t*) malloc(N_NO*sizeof(Node_t)) ;
  if(NO == NULL) arret("mail0d(1) : impossible d\'allouer la memoire") ;

  {
    double* x = (double*) malloc(N_NO*3*sizeof(double)) ;
    if(x == NULL) arret("mail0d(2) : impossible d\'allouer la memoire") ;
  
    Node_GetCoordinate(NO) = x ;
    Node_GetNodeIndex(NO)  = 0 ;
  }

  /* les elements */
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(EL == NULL) {
    fprintf(stderr,"impossible d\'allouer la memoire") ;
    exit(EXIT_FAILURE) ;
  }

  /* Indexes */
  Element_GetElementIndex(EL) = 0 ;

  Element_GetMaterialIndex(EL) = 0 ;
  Element_GetRegionIndex(EL) = 1 ;

  Element_GetNbOfNodes(EL) = 1 ;    /* nb de noeuds */
  Element_GetDimension(EL) = 0 ;

  
  {
    Node_t** no = (Node_t**) malloc(sizeof(Node_t*)) ; /* pointeurs */
    if(!no) arret("mail0d : impossible d\'allouer la memoire") ;
    
    Element_GetPointerToNode(EL) = no ;
  }

  /* numerotation */
  Element_GetNode(EL,0) = NO ;
}



void mail1dold(Mesh_t* mesh,FILE *ficd)
/* MAILLAGE 1D */
{
  int    npt ;
  double* pt ;
  int    *ne,imat ;
  double dx_ini ;
  double* x ;
  int    i,j,k ;

  /* nombre de points */
  fscanf(ficd,"%d",&npt) ;

  pt = (double*) malloc(npt*sizeof(double)) ;
  if(pt == NULL) arret("mail1d (1) : impossible d\'allouer la memoire") ;

  ne = (int*) malloc(npt*sizeof(int)) ;
  if(ne == NULL) arret("mail1d (2) : impossible d\'allouer la memoire") ;

  /* les points */
  for(i = 0 ; i < npt ; i++) fscanf(ficd,"%le",pt + i) ;
  
  /* longueur du premier element */
  fscanf(ficd,"%lf",&dx_ini) ;
  
  /* nombre d'elements de volume */
  for(i = 0 ; i < npt - 1 ; i++) fscanf(ficd,"%d",ne + i) ;

  N_EL = 0 ;
  for(i = 0 ; i < npt - 1 ; i++) N_EL += ne[i] ;

  /* nombre de noeuds */
  if(N_EL > 0) N_NO = N_EL + 1 ; else N_NO = 0 ;

  /* les noeuds */
  NO = (Node_t*) malloc(N_NO*sizeof(Node_t)) ;
  if(NO == NULL) arret("mail1dold(1) : impossible d\'allouer la memoire") ;

  x = (double*) malloc(N_NO*DIM*sizeof(double)) ;
  if(x == NULL) arret("mail1dold(2) : impossible d\'allouer la memoire") ;
  
  for(i = 0 ; i < (int) N_NO ; i++) {
    Node_GetCoordinate(NO + i) = x + i*DIM ;
    Node_GetNodeIndex(NO + i) = i ;
  }

  maillage(pt,ne,dx_ini,npt,NO) ;

  /* les elements */
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(EL == NULL) {
    fprintf(stderr,"impossible d\'allouer la memoire") ;
    exit(EXIT_FAILURE) ;
  }

  /* Indexes */
  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetElementIndex(EL + i) = i ;
  }

  for(k = 0 , i = 0 ; i < npt - 1 ; i++)  {
    fscanf(ficd,"%d",&imat) ; /* les groupes d'elements */
    for(j = k ; j < k + ne[i] ; j++) {
      Element_GetMaterialIndex(EL + j) = imat - 1 ;
      Element_GetRegionIndex(EL + j) = i + 1 ;
    }
    k += ne[i] ;
  }

  free(pt) ;
  free(ne) ;

  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetNbOfNodes(EL + i) = 2 ;    /* nb de noeuds */
    Element_GetDimension(EL + i) = 1 ;
  }

  for(j = 0 , i = 0 ; i < (int) N_EL ; i++) j += Element_GetNbOfNodes(EL + i) ;
  
  {
    Node_t** no = (Node_t**) malloc(j*sizeof(Node_t*)) ; /* pointeurs */
    if(!no) arret("mail1d : impossible d\'allouer la memoire") ;
    
    j = 0 ;
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_GetPointerToNode(EL + i) = no + j ;
      j += Element_GetNbOfNodes(EL + i) ;
    }
  }

  /* numerotation */
  for(i = 0 ; i < (int) N_EL ; i++) {
    for(j = 0 ; j < 2 ; j++) {
      Element_GetNode(EL + i,j) = NO + i + j ;
    }
  }

  /* cas d'elements de surface i=0 et i=n_el-1 */
  if(N_EL > 0) {
    if(Node_GetCoordinate(NO)[0] == Node_GetCoordinate(NO + 1)[0]) {
      Element_GetNbOfNodes(EL) = 1 ;
      Element_GetDimension(EL) = 0 ;
      Element_GetNode(EL,0) = NO + 1 ;
    }
    if(N_EL > 1) {
      if(Node_GetCoordinate(NO + N_EL - 1)[0] == Node_GetCoordinate(NO + N_EL)[0]) {
        Element_GetNbOfNodes(EL + N_EL - 1) = 1 ;
        Element_GetDimension(EL + N_EL - 1) = 0 ;
        Element_GetNode(EL + N_EL - 1,0) = NO + N_EL - 1 ;
      }
    }
  }
}



void mail1d(Mesh_t* mesh,char* str)
/* MAILLAGE 1D */
{
  int    npt ;
  double* pt ;
  int    *ne,imat ;
  double dx_ini ;
  double* x ;
  int    i,j,k ;

  /* nombre de points */
  str += String_Scan(str,"%d",&npt) ;

  pt = (double*) malloc(npt*sizeof(double)) ;
  if(pt == NULL) arret("mail1d (1) : impossible d\'allouer la memoire") ;

  ne = (int*) malloc(npt*sizeof(int)) ;
  if(ne == NULL) arret("mail1d (2) : impossible d\'allouer la memoire") ;

  /* les points */
  for(i = 0 ; i < npt ; i++) {
    str += String_Scan(str,"%le",pt + i) ;
  }
  
  /* longueur du premier element */
  {
    str += String_Scan(str,"%lf",&dx_ini) ;
  }
  
  /* nombre d'elements de volume */
  for(i = 0 ; i < npt - 1 ; i++) {
    str += String_Scan(str,"%d",ne + i) ;
  }

  N_EL = 0 ;
  for(i = 0 ; i < npt - 1 ; i++) N_EL += ne[i] ;

  /* nombre de noeuds */
  if(N_EL > 0) N_NO = N_EL + 1 ; else N_NO = 0 ;

  /* les noeuds */
  NO = (Node_t*) malloc(N_NO*sizeof(Node_t)) ;
  if(NO == NULL) arret("mail1d(1) : impossible d\'allouer la memoire") ;

  x = (double*) malloc(N_NO*DIM*sizeof(double)) ;
  if(x == NULL) arret("mail1d(2) : impossible d\'allouer la memoire") ;
  
  for(i = 0 ; i < (int) N_NO ; i++) {
    Node_GetCoordinate(NO + i) = x + i*DIM ;
    Node_GetNodeIndex(NO + i) = i ;
  }

  maillage(pt,ne,dx_ini,npt,NO) ;

  /* les elements */
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(EL == NULL) {
    fprintf(stderr,"impossible d\'allouer la memoire") ;
    exit(EXIT_FAILURE) ;
  }

  /* Indexes */
  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetElementIndex(EL + i) = i ;
  }

  for(k = 0 , i = 0 ; i < npt - 1 ; i++)  {
    str += String_Scan(str,"%d",&imat) ;
    
    for(j = k ; j < k + ne[i] ; j++) {
      Element_GetMaterialIndex(EL + j) = imat - 1 ;
      Element_GetRegionIndex(EL + j) = i + 1 ;
    }
    k += ne[i] ;
  }

  free(pt) ;
  free(ne) ;

  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetNbOfNodes(EL + i) = 2 ;    /* nb de noeuds */
    Element_GetDimension(EL + i) = 1 ;
  }

  for(j = 0 , i = 0 ; i < (int) N_EL ; i++) j += Element_GetNbOfNodes(EL + i) ;
  
  {
    Node_t** no = (Node_t**) malloc(j*sizeof(Node_t*)) ; /* pointeurs */
    if(!no) arret("mail1d : impossible d\'allouer la memoire") ;
    
    j = 0 ;
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_GetPointerToNode(EL + i) = no + j ;
      j += Element_GetNbOfNodes(EL + i) ;
    }
  }

  /* numerotation */
  for(i = 0 ; i < (int) N_EL ; i++) {
    for(j = 0 ; j < 2 ; j++) {
      Element_GetNode(EL + i,j) = NO + i + j ;
    }
  }

  /* cas d'elements de surface i=0 et i=n_el-1 */
  if(N_EL > 0) {
    if(Node_GetCoordinate(NO)[0] == Node_GetCoordinate(NO + 1)[0]) {
      Element_GetNbOfNodes(EL) = 1 ;
      Element_GetDimension(EL) = 0 ;
      Element_GetNode(EL,0) = NO + 1 ;
    }
    if(N_EL > 1) {
      if(Node_GetCoordinate(NO + N_EL - 1)[0] == Node_GetCoordinate(NO + N_EL)[0]) {
        Element_GetNbOfNodes(EL + N_EL - 1) = 1 ;
        Element_GetDimension(EL + N_EL - 1) = 0 ;
        Element_GetNode(EL + N_EL - 1,0) = NO + N_EL - 1 ;
      }
    }
  }
}



void lit_mail_m1d(Mesh_t* mesh,const char* nom_m1d)
/* MAILLAGE 1D */
{
  int    npt,nreg ;
  double* pt,*lc ;
  int    ireg ;
  int    i,j,k ;
  FILE*  fic_m1d ;

  /* ouverture du fichier */
  fic_m1d = fopen(nom_m1d,"r") ;
  if(fic_m1d == NULL) arret("lit_mail_m1d (1) :impossible d\'ouvrir le fichier") ;

  /* nombre de points */
  fscanf(fic_m1d,"%d",&npt) ;
  nreg = npt - 1 ;

  pt = (double*) malloc(npt*sizeof(double)) ;
  if(!pt) arret("lit_mail_m1d (2) : impossible d\'allouer la memoire") ;

  lc = (double*) malloc(npt*sizeof(double)) ;
  if(!lc) arret("lit_mail_m1d (3) : impossible d\'allouer la memoire") ;

  /* les points et longueurs caracteristiques */
  for(i=0;i<npt;i++) fscanf(fic_m1d,"%le %le",pt+i,lc+i) ;

  N_NO = mesh1d(pt,lc,nreg,&NO) + 1 ;
  
  for(i = 0; i < (int) N_NO; i++) {
    Node_GetNodeIndex(NO + i) = i ;
  }

  /* les elements */
  N_EL = N_NO - 1 ;
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(!EL) arret("lit_mail_m1d (4) : impossible d\'allouer la memoire") ;

  /* les materiaux et les regions */
  for(k = 0 , ireg = 0 ; ireg < nreg ; ireg++)  {
    int imat ;
    fscanf(fic_m1d,"%d",&imat) ; /* materiau */
    while(k < (int) N_EL) {
      Element_GetMaterialIndex(EL + k) = imat - 1 ;
      Element_GetRegionIndex(EL + k) = ireg + 1 ;
      k++ ;
      if(fabs(pt[ireg + 1] - Node_GetCoordinate(NO + k)[0]) < lc[ireg + 1]*0.1) break ;
    }
  }

  free(pt) ;

  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetNbOfNodes(EL + i) = 2 ;    /* nb de noeuds */
    Element_GetDimension(EL + i) = 1 ;
  }

  for(j = 0 , i = 0 ; i < (int) N_EL ; i++) j += Element_GetNbOfNodes(EL + i) ;
  
  {
    Node_t** no = (Node_t**) malloc(j*sizeof(Node_t*)) ;/* pointeurs */
    if(!no) arret("lit_mail_m1d (5) : impossible d\'allouer la memoire") ;
    
    j = 0 ;
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_GetPointerToNode(EL + i) = no + j ;
      j += Element_GetNbOfNodes(EL + i) ;
    }
  }

  /* numerotation */
  for(i = 0 ; i < (int) N_EL ; i++) {
    for(j = 0 ; j < 2 ; j++) {
      Element_GetNode(EL + i,j) = NO + i + j ;
    }
  }
  
  
  /* cas d'elements de surface i=0 et i=n_el-1 */
  if(N_EL > 0) {
    if(Node_GetCoordinate(NO)[0] == Node_GetCoordinate(NO + 1)[0]) {
      Element_GetNbOfNodes(EL) = 1 ;
      Element_GetDimension(EL) = 0 ;
      Element_GetNode(EL,0) = NO + 1 ;
    }
    if(N_EL > 1) {
      if(Node_GetCoordinate(NO + N_EL - 1)[0] == Node_GetCoordinate(NO + N_EL)[0]) {
        Element_GetNbOfNodes(EL + N_EL - 1) = 1 ;
        Element_GetDimension(EL + N_EL - 1) = 0 ;
        Element_GetNode(EL + N_EL - 1,0) = NO + N_EL - 1 ;
      }
    }
  }
}


void maillage(double* pt,int* ne,double dx_ini,int npt,Node_t* no)
/* Calcul des coordonnees (no) d'un maillage 1D */
{
  int    npt1 = npt - 1 ;
  double a,e,dx ;
  int    i ;

  Node_GetCoordinate(no)[0] = pt[0] ;

  dx = dx_ini ;
  for(i = 0 ; i < npt1 ; i++) {
    int j ;
    if(ne[i] > 1) e = log(fabs(pt[i + 1] - pt[i])/dx)/log((double) ne[i]) ; else e = 1. ;
    for(j = 1 ; j <= ne[i] ; j++) {
      a = ((double) j)/ne[i] ;
      Node_GetCoordinate(++no)[0] = pt[i] + (pt[i + 1] - pt[i])*pow(a,e) ;
    }
    a = fabs(Node_GetCoordinate(no)[0] - Node_GetCoordinate(no - 1)[0]) ;
    if(a > 0.) dx = a ;
  }
}


int mesh1d(double* point, double* l_c, int n_reg, Node_t** no)
/* Calcul des coordonnees (no) d'un maillage 1D
   et retourne le nombre d'elements */
{
#define MAX_NREG 100
  int    i,reg,n_el,n_no ;
  Node_t* p ;
  int    ne[MAX_NREG] ;

  if(n_reg > MAX_NREG) arret("mesh1d : trop de regions") ;

  n_el = 0 ;
  for(reg = 0 ; reg < n_reg ; reg++) {
    double l = fabs(point[reg + 1] - point[reg]) ;
    double b = l_c[reg + 1]/l_c[reg] ;
    double a = (l - l_c[reg])/(l - l_c[reg + 1]) ;
    int n ;

    if(l <= l_c[reg] || l <= l_c[reg + 1]) {
      n = 1 ;
    } else {
      if(a == 1.) n = floor(l/l_c[reg] + 0.5) ;
      else n = floor(log(b)/log(a) + 0.5) + 1 ;
    }

    ne[reg] = n ;

    n_el += n ;
  }

  n_no = n_el + 1 ;
  
  
  {
    Node_t* node = (Node_t*) malloc(n_no*sizeof(Node_t)) ;
    double* x = (double*) malloc(n_no*sizeof(double)) ;
    if(!node) arret("mesh1d(1) : impossible d\'allouer la memoire") ;
    if(!x) arret("mesh1d(2) : impossible d\'allouer la memoire") ;
  
    *no = node ;
    
    for(i = 0 ; i < n_no ; i++) {
      Node_GetCoordinate(node + i) = x + i ;
    }
  }


  p = *no ;
  for(reg = 0 ; reg < n_reg ; reg++) {
    double l = fabs(point[reg + 1] - point[reg]) ;
    double b = l_c[reg + 1]/l_c[reg] ;
    double a ;
    double l_i ;
    int n = ne[reg] ;

    /* a = (l - l_c[reg])/(l - l_c[reg+1]) ; */
    if(l_c[reg + 1] == l_c[reg]) a = 1. ;
    else if(n > 1) a = exp(log(b)/(n - 1)) ;
    else a = 0 ;
    
    /* l_i  = l_c[reg] ; */
    if(a == 1.) l_i = l/n ;
    else l_i = l*(1. - a)/(1. - a*b) ;

    Node_GetCoordinate(p)[0] = point[reg] ;
    Node_GetCoordinate(p + n)[0] = point[reg + 1] ;
    for(i = 1 ; i < n ; i++) {
      Node_GetCoordinate(p + i)[0] = Node_GetCoordinate(p + i - 1)[0] + l_i ;
      l_i *= a ;
    }
    p   += n ;
  }

  return n_el ;
}



void ecrit_mail_msh_1(Mesh_t* mesh,const char* nom)
/* fichier de maillage au format GMSH "nom.msh" */
{
  int    i ;
  char   nom_msh[Mesh_MaxLengthOfFileName] ;
  FILE*  fic_msh ;
  
  sprintf(nom_msh,"%s.msh",nom) ;
  fic_msh = fopen(nom_msh,"w") ;
  if(!fic_msh) {
    arret("erreur a l'ouverture du fichier %s",nom_msh) ;
  }

  /* Les noeuds */
  fprintf(fic_msh,"$NOD\n") ;
  /* nombre de noeuds */
  fprintf(fic_msh,"%d\n",N_NO) ;
  for(i=0;i<(int) N_NO;i++) {
    int    j ;
    /* le numero du noeud */
    fprintf(fic_msh,"%d",i+1) ;
    /* les coordonnees*/
    for(j=0;j<(int) DIM;j++) fprintf(fic_msh," %e",Node_GetCoordinate(NO + i)[j]) ;
    for(j=(int) DIM;j<3;j++) fprintf(fic_msh," 0") ;
    fprintf (fic_msh,"\n") ;
  }
  fprintf(fic_msh,"$ENDNOD\n") ;

  /* Les elements */
  fprintf(fic_msh,"$ELM\n") ;
  /* nombre d'elements lus */
  fprintf(fic_msh,"%d\n",N_EL) ;
  for(i=0;i<(int) N_EL;i++) {
    int    nn = Element_GetNbOfNodes(EL + i) ;
    int    imat = Element_GetMaterialIndex(EL + i) ;
    int    reg = Element_GetRegionIndex(EL + i) ;
    int    dim_el = Element_GetDimension(EL + i) ;
    int    j ;
    /* le numero d'element */
    fprintf(fic_msh,"%d",i+1) ;
    /* groupe et nb de noeuds par elements */
    fprintf(fic_msh," %d %d %d %d",gmsh_ElmType(dim_el,nn),imat+1,reg,nn) ;
    /* numerotation */
    for(j=0;j<nn;j++) {
      Node_t* node_j = Element_GetNode(EL + i,j) ;
      fprintf(fic_msh," %d",Node_GetNodeIndex(node_j) + 1) ;
    }
    fprintf(fic_msh,"\n") ;
  }
  fprintf(fic_msh,"$ENDELM\n") ;

  /* fermeture du ficher */
  fclose(fic_msh) ;
}


void ecrit_mail_msh_2(Mesh_t* mesh,const char* nom)
/* fichier de maillage au format GMSH "nom.msh" */
{
  int    i ;
  char   nom_msh[Mesh_MaxLengthOfFileName] ;
  FILE*  fic_msh ;
  
  /* Ouverture de fichier */
  sprintf(nom_msh,"%s.msh",nom) ;
  fic_msh = fopen(nom_msh,"w") ;
  if(!fic_msh) {
    arret("erreur a l'ouverture du fichier %s",nom_msh) ;
  }

  /* Le format */
  fprintf(fic_msh,"$MeshFormat\n") ;
  fprintf(fic_msh,"2 0 %lu\n",sizeof(double)) ;
  fprintf(fic_msh,"$EndMeshFormat\n") ;

  /* Les noeuds */
  fprintf(fic_msh,"$Nodes\n") ;
  /* nombre de noeuds */
  fprintf(fic_msh,"%d\n",N_NO) ;
  for(i=0;i<(int) N_NO;i++) {
    int    j ;
    double* x = Node_GetCoordinate(NO + i) ;
    /* le numero du noeud */
    fprintf(fic_msh,"%d",i+1) ;
    /* les coordonnees*/
    for(j=0;j<(int) DIM;j++) fprintf(fic_msh," %e",x[j]) ;
    for(j=(int) DIM;j<3;j++) fprintf(fic_msh," 0") ;
    fprintf (fic_msh,"\n") ;
  }
  fprintf(fic_msh,"$EndNoces\n") ;

  /* Les elements */
  fprintf(fic_msh,"$Elements\n") ;
  /* nombre d'elements lus */
  fprintf(fic_msh,"%d\n",N_EL) ;
  for(i=0;i<(int) N_EL;i++) {
    int    nn = Element_GetNbOfNodes(EL + i) ;
    int    imat = Element_GetMaterialIndex(EL + i) ;
    int    reg = Element_GetRegionIndex(EL + i) ;
    int    dim_el = Element_GetDimension(EL + i) ;
    int    j ;
    if(nn > 8) arret("ecrit_mail_msh_2 : trop de noeud") ;
    /* le numero d'element */
    fprintf(fic_msh,"%d",i+1) ;
    /* le type */
    fprintf(fic_msh," %d",gmsh_ElmType(dim_el,nn)) ;
    /* 3 tags */
    fprintf(fic_msh," 3") ;
    /* physical, elementary, partition */
    fprintf(fic_msh," %d %d 0",imat+1,reg) ;
    /* numerotation */
    for(j=0;j<nn;j++) {
      Node_t* node_j = Element_GetNode(EL + i,j) ;
      fprintf(fic_msh," %d",Node_GetNodeIndex(node_j) + 1) ;
    }
    fprintf(fic_msh,"\n") ;
  }
  fprintf(fic_msh,"$EndElements\n") ;

  /* fermeture du ficher */
  fclose(fic_msh) ;
}


int gmsh_ElmType(int dim,int nn)
/* N = i(Vertices) + j(Edges) + k(Faces) + l(Volume) */
{
  if(dim == 0) {
    switch (nn) {
    case 1 : return 15;             /* point */
    default: return 0;
    }
  } else if(dim == 1) {
    switch (nn) {
    case 2 : return 1;              /* line 1 */
    case 3 : return 8;              /* line 2 */
    default: return 0;
    }
  } else if(dim == 2) {
    switch (nn) {
    case 3 : return 2;              /* triangle 1 */
    case 6 : return 9;              /* triangle 2 */
    case 4 : return 3;              /* quadrangle 1 */
    case 9 : return 10;             /* quadrangle 2 */
    case 8 : return 16;             /* quadrangle 2 */
    default: return 0;
    }
  } else if(dim == 3) {
    switch (nn) {
    case 4 : return 4;              /* tetrahedron 1 */
    case 10: return 11;             /* tetrahedron 2 */
    case 8 : return 5;              /* hexahedron 1 */
    case 27: return 12;             /* hexahedron 2 */
    case 20: return 17;             /* hexahedron 2 */
    case 6 : return 6;              /* prism 1 */
    case 18: return 13;             /* prism 2 */
    case 15: return 18;             /* prism 2 */
    case 5 : return 7;              /* pyramid 1 */
    case 14: return 14;             /* pyramid 2 */
    case 13: return 19;             /* pyramid 2 */
    default: return 0;
    }
  } else {
    arret("gmsh_ElmType : dimension incorrect") ;
    return 0 ;
  }
}

int gmsh_NbNodes(int type)
/* N = i(Vertices) + j(Edges) + k(Faces) + l(Volume) */
{
  switch (type) {
  case 15: return 1;              /* point */
  case 1 : return 2;              /* line 1 */
  case 8 : return 2 + 1;          /* line 2 */
  case 2 : return 3;              /* triangle 1 */
  case 9 : return 3 + 3;          /* triangle 2 */
  case 3 : return 4;              /* quadrangle 1 */
  case 10: return 4 + 4 + 1;      /* quadrangle 2 */
  case 16: return 4 + 4;          /* quadrangle 2 */
  case 4 : return 4;              /* tetrahedron 1 */
  case 11: return 4 + 6;          /* tetrahedron 2 */
  case 5 : return 8;              /* hexahedron 1 */
  case 12: return 8 + 12 + 6 + 1; /* hexahedron 2 */
  case 17: return 8 + 12;         /* hexahedron 2 */
  case 6 : return 6;              /* prism 1 */
  case 13: return 6 + 9 + 3;      /* prism 2 */
  case 18: return 6 + 9;          /* prism 2 */
  case 7 : return 5;              /* pyramid 1 */
  case 14: return 5 + 8 + 1;      /* pyramid 2 */
  case 19: return 5 + 8;          /* pyramid 2 */
  default: return 0;
  }
}

int gmsh_DimElement(int type)
{
  switch (type) {
  case 15: return 0;              /* point */
  case 1 : return 1;              /* line 1 */
  case 8 : return 1;              /* line 2 */
  case 2 : return 2;              /* triangle 1 */
  case 9 : return 2;              /* triangle 2 */
  case 3 : return 2;              /* quadrangle 1 */
  case 10: return 2;              /* quadrangle 2 */
  case 16: return 2;              /* quadrangle 2 */
  case 4 : return 3;              /* tetrahedron 1 */
  case 11: return 3;              /* tetrahedron 2 */
  case 5 : return 3;              /* hexahedron 1 */
  case 12: return 3;              /* hexahedron 2 */
  case 17: return 3;              /* hexahedron 2 */
  case 6 : return 3;              /* prism 1 */
  case 13: return 3;              /* prism 2 */
  case 18: return 3;              /* prism 2 */
  case 7 : return 3;              /* pyramid 1 */
  case 14: return 3;              /* pyramid 2 */
  case 19: return 3;              /* pyramid 2 */
  default: return -1;
  }
}




void lit_mail_gmsh(Mesh_t* mesh,const char* nom_msh)
/* lit le maillage au format GMSH dans un fichier */
{
  char   line[Mesh_MaxLengthOfTextLine] ;
  FILE*  fic_msh ;

  /* ouverture du fichier */
  fic_msh = fopen(nom_msh,"r") ;
  if(!fic_msh) arret("lit_mail_gmsh (1) :impossible d\'ouvrir le fichier") ;

  do {
    fgets(line,sizeof(line),fic_msh) ;
  } while(line[0] != '$') ;

  /* fermeture du fichier */
  fclose(fic_msh) ;

  if(!strncmp(&line[1],"NOD",3)) { /* Version 1.0 */
    lit_mail_gmsh_1(mesh,nom_msh) ;
    return ;
  } else if(!strncmp(&line[1],"MeshFormat",10)) { /* Version 2.0 */
    lit_mail_gmsh_2(mesh,nom_msh) ;
    return ;
  }
  arret("lit_mail_gmsh (2) : non prevu") ;
}


void lit_mail_gmsh_1(Mesh_t* mesh,const char* nom_msh)
/* lit le maillage au format GMSH version 1.0 dans un fichier */
{
  int    i,n_c ;
  int    n_no_lu,n_el_lu ;
  char   mot[Mesh_MaxLengthOfKeyWord],line[Mesh_MaxLengthOfTextLine] ;
  Node_t** p_node ;
  FILE*  fic_msh ;

  /* ouverture du fichier */
  fic_msh = fopen(nom_msh,"r") ;
  if(fic_msh == NULL) arret("lit_mail_gmsh_1 (1) :impossible d\'ouvrir le fichier") ;

  /* analyse du fichier pour le calcul de n_no */
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$NOD") != 0) arret("lit_mail_gmsh_1 (1) : pas de $NOD") ;

  /* nombre de noeuds lus */
  fscanf(fic_msh,"%d",&n_no_lu) ;
  /* calcul de n_no */
  N_NO = 0 ;
  for(i = 0 ; i < n_no_lu ; i++) {
    int n ;
    /* le numero du noeud */
    fscanf(fic_msh,"%d",&n) ;
    if((int) N_NO < n) N_NO = n ;
    /* on lit le reste de la ligne */
    if(fgets(line,sizeof(line),fic_msh) == NULL) {
      fprintf(stdout,"erreur ou fin de fichier\n") ;
    }
  }
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$ENDNOD") != 0) arret("lit_mail_gmsh_1 (2) : pas de $ENDNOD") ;

  /* analyse du fichier pour le calcul de n_el */
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$ELM") != 0) arret("lit_mail_gmsh_1 (3) : pas de $ELM") ;

  /* nombre d'elements lus */
  fscanf(fic_msh,"%d",&n_el_lu) ;
  /* calcul de n_el et de n_c */
  N_EL = 0 ;
  n_c = 0 ;
  for(i = 0 ; i < n_el_lu ; i++) {
    int n ;
    unsigned short int   nn ;
    /* le numero d'element */
    fscanf(fic_msh,"%d",&n) ;
    if((int) N_EL < n) N_EL = n ;
    /* on lit le reste de la ligne */
    if(fgets(line,sizeof(line),fic_msh) == NULL) {
      fprintf(stdout,"erreur ou fin de fichier\n") ;
    }
    /* elm_type, identificateurs et nb de noeuds */
    sscanf(line,"%*d %*d %*d %hu",&nn) ;
    n_c += nn ;
  }
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$ENDELM") != 0) arret("lit_mail_gmsh_1 (4) : pas de $ENDELM") ;

  /* on recommence */
  rewind(fic_msh) ;

  /* NOEUDS */
  fscanf(fic_msh,"%s",mot) ;
  /* nombre de noeuds */
  fscanf(fic_msh,"%d",&n_no_lu) ;
  /* pointeurs */
  NO = (Node_t*) malloc(N_NO*sizeof(Node_t)) ;
  if(NO == NULL) arret("lit_mail_gmsh_1 (5) : impossible d\'allouer la memoire") ;
  
  {
    double* x = (double*) calloc(N_NO*DIM,sizeof(double)) ;
    if(!x) arret("lit_mail_gmsh_1 (6) : impossible d\'allouer la memoire") ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_GetCoordinate(NO + i) = x + i*DIM ;
      Node_GetNodeIndex(NO + i) = i ;
    }
  }

  /* les noeuds */
  for(i = 0 ; i < n_no_lu ; i++) {
    int j,n ;
    /* le numero du noeud */
    fscanf(fic_msh,"%d",&n) ; n -= 1 ;
    /* les coordonnees*/
    for(j = 0 ; j < DIM ; j++) {
      fscanf(fic_msh,"%le",Node_GetCoordinate(NO + n) + j) ;
    }
    /* on lit le reste de la ligne */
    if(fgets(line,sizeof(line),fic_msh) == NULL) {
      fprintf(stdout,"erreur ou fin de fichier\n") ;
    }
  }
  fscanf(fic_msh,"%s",mot) ;

  /* ELEMENTS */
  fscanf(fic_msh,"%s",mot) ;
  /* nombre d'elements */
  fscanf(fic_msh,"%d",&n_el_lu) ;
  /* pointeurs */
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(EL == NULL) {
    fprintf(stderr,"impossible d\'allouer la memoire (lit_mail_gmsh_1)") ;
    exit(EXIT_FAILURE) ;
  }
  
  p_node = (Node_t**) malloc(n_c*sizeof(Node_t*)) ;
  if(!p_node) arret("impossible d\'allouer la memoire (lit_mail_gmsh_1)") ;
  
  /* initialisation */
  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetNbOfNodes(EL + i) = 0 ;
    Element_GetMaterialIndex(EL + i) = -1 ;
    Element_GetRegionIndex(EL + i) = -1 ;
    Element_GetDimension(EL + i) = 0 ;
  }
  n_c = 0 ;
  /* les elements */
  for(i = 0 ; i < n_el_lu ; i++) {
    int j,n ;
    int imat, reg, nn ;
    int type ;
    /* le numero d'element */
    fscanf(fic_msh,"%d",&n) ; n -= 1 ;
    /* elm_type */
    fscanf(fic_msh,"%d",&type) ;
    /* identificateurs et nb de noeuds */
    fscanf(fic_msh,"%d %d %d",&imat,&reg,&nn) ;
    
    Element_GetElementIndex(EL + n) = n ;
    Element_GetMaterialIndex(EL + n) = imat - 1 ;
    Element_GetNbOfNodes(EL + n) = nn ;
    Element_GetRegionIndex(EL + n) = reg ;
    Element_GetDimension(EL + n) = gmsh_DimElement(type) ;
    if(Element_GetNbOfNodes(EL + n) > Element_MaxNbOfNodes) arret("trop de noeuds") ;
    /* numerotation */
    Element_GetPointerToNode(EL + n) = p_node + n_c ;
    n_c += nn ;
    for(j = 0 ; j < nn ; j++) {
      int nodeindex ;
      fscanf(fic_msh,"%d",&nodeindex) ;
      Element_GetNode(EL + n,j) = NO + nodeindex - 1 ;
    }
  }
  
  fscanf(fic_msh,"%s",mot) ;
  /* fermeture du fichier */
  fclose(fic_msh) ;
}




void lit_mail_gmsh_2(Mesh_t* mesh,const char* nom_msh)
/* lit le maillage au format GMSH version 2.0 dans un fichier */
{
  double version ;
  int    file_type,data_size ;
  int    i,n_c ;
  int    nb_nodes,nb_elements ;
  Node_t** p_node ;
  char   mot[Mesh_MaxLengthOfKeyWord],line[Mesh_MaxLengthOfTextLine] ;
  FILE*  fic_msh ;

  /* ouverture du fichier */
  fic_msh = fopen(nom_msh,"r") ;
  if(!fic_msh) arret("lit_mail_gmsh_2 (1) :impossible d\'ouvrir le fichier") ;

  /* analyse du fichier pour le format */
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$MeshFormat")) arret("lit_mail_gmsh_2 (2) : pas de $MeshFormat") ;
  fscanf(fic_msh, "%lf %d %d\n",&version,&file_type,&data_size) ;
  if(floor(version) != 2) arret("Error: Wrong msh file version") ;
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$EndMeshFormat")) arret("lit_mail_gmsh_2 (2) : pas de $EndFormat") ;


  /* analyse du fichier pour le calcul de n_no */
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$Nodes")) arret("lit_mail_gmsh_2 (2) : pas de $Nodes") ;

  /* nombre de noeuds lus */
  fscanf(fic_msh,"%d",&nb_nodes) ;
  /* calcul de n_no */
  N_NO = 0 ;
  for(i = 0 ; i < nb_nodes ; i++) {
    int n ;
    /* le numero du noeud */
    fscanf(fic_msh,"%d",&n) ;
    if((int) N_NO < n) N_NO = n ;
    /* on lit le reste de la ligne */
    if(!fgets(line,sizeof(line),fic_msh)) {
      arret("lit_mail_gmsh_2 (3) : erreur ou fin de fichier") ;
    }
  }
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$EndNodes")) arret("lit_mail_gmsh_2 (4) : pas de $EndNodes") ;

  /* analyse du fichier pour le calcul de n_el */
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$Elements")) arret("lit_mail_gmsh_2 (5) : pas de $Elements") ;

  /* nombre d'elements lus */
  fscanf(fic_msh,"%d",&nb_elements) ;
  /* calcul de n_el et de n_c */
  N_EL = 0 ;
  n_c = 0 ; /* nb de numeros de noeuds cumules */
  for(i = 0 ; i < nb_elements ; i++) {
    int n,elm_type ;
    /* le numero d'element */
    fscanf(fic_msh,"%d",&n) ;
    if((int) N_EL < n) N_EL = n ;
    /* on lit le reste de la ligne */
    if(!fgets(line,sizeof(line),fic_msh)) {
      arret("lit_mail_gmsh_2 (6) : erreur ou fin de fichier") ;
    }
    /* elm_type */
    sscanf(line,"%d",&elm_type) ;
    n_c += gmsh_NbNodes(elm_type) ;
  }
  fscanf(fic_msh,"%s",mot) ;
  if(strcmp(mot,"$EndElements")) arret("lit_mail_gmsh_2 (7) : pas de $EndElements") ;

  /* on recommence */
  rewind(fic_msh) ;

  /* FORMAT */
  fscanf(fic_msh,"%s",mot) ;
  fscanf(fic_msh, "%lf %d %d\n",&version,&file_type,&data_size) ;
  fscanf(fic_msh,"%s",mot) ;

  /* NOEUDS */
  fscanf(fic_msh,"%s",mot) ;
  /* nombre de noeuds */
  fscanf(fic_msh,"%d",&nb_nodes) ;
  /* pointeurs */
  NO = (Node_t*) malloc(N_NO*sizeof(Node_t)) ;
  if(NO == NULL) arret("lit_mail_gmsh_2 (8) : impossible d\'allouer la memoire") ;
  
  {
    double* x = (double*) calloc(N_NO*DIM,sizeof(double)) ;
    if(!x) arret("lit_mail_gmsh_2 (9) : impossible d\'allouer la memoire") ;
    
    for(i = 0 ; i < (int) N_NO ; i++) {
      Node_GetCoordinate(NO + i) = x + i*DIM ;
      Node_GetNodeIndex(NO + i) = i ;
    }
  }

  /* les noeuds */
  for(i = 0 ; i < nb_nodes ; i++) {
    int n,j ;
    /* le numero du noeud */
    fscanf(fic_msh,"%d",&n) ; n -= 1 ;
    /* les coordonnees*/
    for(j = 0 ; j < DIM ; j++) {
      fscanf(fic_msh,"%le",Node_GetCoordinate(NO + n) + j) ;
    }
    /* on lit le reste de la ligne */
    if(!fgets(line,sizeof(line),fic_msh)) {
      arret("lit_mail_gmsh_2 (10) : erreur ou fin de fichier") ;
    }
  }
  fscanf(fic_msh,"%s",mot) ;

  /* ELEMENTS */
  fscanf(fic_msh,"%s",mot) ;
  /* nombre d'elements */
  fscanf(fic_msh,"%d",&nb_elements) ;
  /* pointeurs */
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(EL == NULL) {
    arret("lit_mail_gmsh_2 (11) : impossible d\'allouer la memoire") ;
  }
  p_node = (Node_t**) malloc(n_c*sizeof(Node_t*)) ;
  if(!p_node) arret("lit_mail_gmsh_2 (12) : impossible d\'allouer la memoire") ;

  /* initialisation */
  for(i = 0 ; i < (int) N_EL ; i++) {
    Element_GetNbOfNodes(EL + i) = 0 ;
    Element_GetMaterialIndex(EL + i) = -1 ;
    Element_GetRegionIndex(EL + i) = -1 ;
    Element_GetDimension(EL + i) = 0 ;
  }

  n_c = 0 ; /* nb de numeros de noeuds cumules */
  /* les elements */
  for(i = 0 ; i < nb_elements ; i++) {
    int n,elm_type,nb_tags,j ;
    int physical,elementary,partition ;
    /* le numero d'element */
    fscanf(fic_msh,"%d",&n) ; n -= 1 ;
    /* elm_type */
    fscanf(fic_msh,"%d %d",&elm_type,&nb_tags) ;
    elementary = physical = partition = 1;
    for(j = 0; j < nb_tags; j++){
      int tag ;
      fscanf(fic_msh, "%d", &tag);	    
      if(j == 0)      physical   = tag ;
      else if(j == 1) elementary = tag ;
      else if(j == 2) partition  = tag ;
      /* ignore any other tags for now */
    }
    Element_GetElementIndex(EL + n) = n ;
    Element_GetMaterialIndex(EL + n) = physical - 1 ;
    Element_GetRegionIndex(EL + n) = elementary ;
    Element_GetNbOfNodes(EL + n) = gmsh_NbNodes(elm_type) ;
    Element_GetDimension(EL + n) = gmsh_DimElement(elm_type) ;
    if(!Element_GetNbOfNodes(EL + n)){
      arret("lit_mail_gmsh_2 (13) : Error: Unknown type for element"); 
    }
    if(Element_GetNbOfNodes(EL + n) > Element_MaxNbOfNodes) arret("lit_mail_gmsh_2 (14) : trop de noeuds") ;
    /* numerotation */
    Element_GetPointerToNode(EL + n) = p_node + n_c ;
    n_c += Element_GetNbOfNodes(EL + n) ;
    for(j = 0; j < Element_GetNbOfNodes(EL + n) ; j++) {
      int nodeindex ;
      fscanf(fic_msh,"%d",&nodeindex) ;
      Element_GetNode(EL + n,j) = NO + nodeindex - 1 ;
    }
  }
  
  fscanf(fic_msh,"%s",mot) ;
  /* fermeture du fichier */
  fclose(fic_msh) ;
}


void lit_mail_cesar(Mesh_t* mesh,const char* nom)
/* lit le maillage au format CESAR dans un fichier */
{
  int    dim_el[3][8]  = {{0,1,1,-1,-1,-1,-1,-1},{0,1,2,2,-1,2,-1,2},{0,1,2,3,-1,-1,-1,3}} ;
  double* x ;
  int    i,j ;
  /*
    char   type_el[3][8] = {{'P','L'},{'P','L','T','Q'},{'P','L','T','S','Y','I',' ','H'}} ; 
  */
  char   mot[Mesh_MaxLengthOfKeyWord] ;
  FILE*  fic_ces ;

  /* ouverture du fichier */
  fic_ces = fopen(nom,"r") ;
  if(fic_ces == NULL) {
    fprintf(stderr,"erreur a l'ouverture du fichier %s\n",nom) ;
    exit(EXIT_FAILURE) ;
  }

  /* les noeuds */
  fscanf(fic_ces,"%s",mot) ;
  if(strcmp(mot,"COOR") != 0) {
    fprintf(stdout,"pas de COOR !\n") ;
    exit(EXIT_SUCCESS) ;
  }
  
  /* nombre de noeuds */
  fscanf(fic_ces,"%*d %*d %u %*d",&N_NO) ;

  /* les noeuds */
  NO = (Node_t*) malloc(N_NO*sizeof(Node_t)) ;
  if(NO == NULL) arret("lit_mail_cesar(1) : impossible d\'allouer la memoire") ;

  x = (double*) malloc(N_NO*DIM*sizeof(double)) ;
  if(x == NULL) arret("lit_mail_cesar(2) : impossible d\'allouer la memoire") ;
  
  for(i = 0 ; i < (int) N_NO ; i++) {
    Node_GetCoordinate(NO + i) = x + i*DIM ;
    Node_GetNodeIndex(NO + i) = i ;
  }

  for(i = 0 ; i < (int) N_NO*DIM ; i++) fscanf(fic_ces,"%le",x + i) ;

  /* les elements */
  fscanf(fic_ces,"%s",mot) ;
  if(strcmp(mot,"ELEM") != 0) {
    fprintf(stdout,"pas de ELEM !\n") ;
    exit(EXIT_SUCCESS) ;
  }
  
  /* nombre d'elements */
  fscanf(fic_ces,"%*d %*d %u %*d",&N_EL) ;
  
  /* les elements */
  EL = (Element_t*) malloc(N_EL*sizeof(Element_t)) ;
  if(EL == NULL) {
    fprintf(stderr,"impossible d\'allouer la memoire") ;
    exit(EXIT_FAILURE) ;
  }
  
  /* nb de noeuds par element */
  fscanf(fic_ces,"%d",&j) ;
  for(i = 0 ; i < (int) N_EL ; i++) {
    short unsigned int nn ;
    fscanf(fic_ces,"%hu",&nn) ;
    nn -= j ;
    Element_GetElementIndex(EL + i) = i ;
    Element_GetNbOfNodes(EL + i) = nn ;
    Element_GetDimension(EL + i) = dim_el[DIM - 1][nn - 1] ;
    j += nn ;
  }
  
  for(i = 0 ; i < (int) N_EL ; i++) if(Element_GetNbOfNodes(EL + i) > Element_MaxNbOfNodes) {
    fprintf(stdout,"trop de noeuds\n") ;
    exit(EXIT_SUCCESS) ;
  }
  
  /* pointeurs */
  for(j = 0 , i = 0 ; i < (int) N_EL ; i++) j += Element_GetNbOfNodes(EL + i) ;
  {
    Node_t** no = (Node_t**) malloc(j*sizeof(Node_t*)) ;
    if(!no) arret("impossible d\'allouer la memoire") ;
    
    j = 0 ;
    for(i = 0 ; i < (int) N_EL ; i++) {
      Element_GetPointerToNode(EL + i) = no + j ;
      j += Element_GetNbOfNodes(EL + i) ;
    }
  }
  
  /* numerotation */
  for(i = 0 ; i < (int) N_EL ; i++) {
    int nn = Element_GetNbOfNodes(EL + i) ;
    for(j = 0 ; j < nn ; j++) {
      int nodeindex ;
      fscanf(fic_ces,"%d",&nodeindex) ;
      Element_GetNode(EL + i,j) = NO + nodeindex - 1 ;
    }
  }
  
  /* les noms -> les regions */
  for(i = 0 ; i < (int) N_EL ; i++) {
    int    n_type = 4 ;
    char   nom_el[][5] = {"OBS2","OBS3","OBT3","OBQ4"} ;
    fscanf(fic_ces,"%s",mot) ;
    Element_GetRegionIndex(EL + i) = 0 ;
    for(j = 0 ; j < n_type ; j++) {
      if(!strncmp(mot,nom_el[j],4)) {
        Element_GetRegionIndex(EL + i) = j + 1 ;
      }
    }
    if(Element_GetRegionIndex(EL + i)  == 0) arret("lit_mail_cesar : pas de region") ;
  }

  /* les groupes -> les materiaux */
  for(i = 0 ; i < (int) N_EL ; i++) {
    int imat ;
    fscanf(fic_ces,"%d",&imat) ;
    imat -= 1 ;
    Element_GetMaterialIndex(EL + i) = imat ;
    Element_GetRegionIndex(EL + i) = 10*Element_GetRegionIndex(EL + i) + imat ; /* formule pas satisfaisante ! */
  }

  fclose(fic_ces) ;
}




int*  Mesh_ReadInversePermutationOfNodes(DataFile_t* datafile,int n_no)
/** Read the inverse permutation vector of nodes in a file if it exists 
 *  or initialize it with identity function. Return a pointer to n_no int. 
 **/
{
  int    *perm = (int*) malloc(n_no*sizeof(int)) ;
  
  if(!perm) {
    arret("Mesh_ReadInversePermutationOfNodes(1)") ;
  }

  {
    char   nom_iperm[Mesh_MaxLengthOfFileName] ;
    
    {
      char*  filename = DataFile_GetFileName(datafile) ;
    
      if(strlen(filename) + 12 > Mesh_MaxLengthOfFileName) {
        arret("Mesh_ReadInversePermutationOfNodes(2)") ;
      }
    
      sprintf(nom_iperm,"%s.graph.iperm",filename) ;
    }
  
    {
      FILE* fic_iperm = fopen(nom_iperm,"r") ;
  
      if(!fic_iperm) {
        int  i ;
    
        for(i = 0 ; i < n_no ; i++) perm[i] = i ;
    
      } else {
        int  i ;
    
        for(i = 0 ; i < n_no ; i++) {
          int   j ;
      
          fscanf(fic_iperm,"%d",&j) ;
          perm[j] = i ;
        }
    
        fclose(fic_iperm) ;
      }
    }
  }
  
  return(perm) ;
}

