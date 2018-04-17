#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Tools/Math.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Message.h"
#include "Periodicities.h"
#include "Graph.h"


static Graph_t*  (Periodicities_ComputeGraph)(Mesh_t*) ;


Periodicities_t* (Periodicities_Create)(DataFile_t* datafile)
{
  Periodicities_t* periodicities = NULL ;
  
  {
    int n = DataFile_CountNbOfKeyWords(datafile,"PERIODICITIES,Periodicities",",") ;
    
    if(n == 0) {
      DataFile_CloseFile(datafile) ;
      return(NULL) ;
    } else if(n > 1) {
      arret("Periodicities_Create") ;
    }
    
    periodicities  = (Periodicities_t*) malloc(sizeof(Periodicities_t)) ;
  
    if(!periodicities) arret("Periodicities_Create") ;
  }
  
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"PERIODICITIES,Periodicities",",",1) ;
  
  Message_Direct("Enter in %s","Periodicities") ;
  Message_Direct("\n") ;
  
  
  {
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    int   n = atoi(line) ;
    
    Periodicities_GetNbOfPeriodicities(periodicities) = n ;
    if(n <= 0) return(periodicities) ;
  }


  {
    int n = Periodicities_GetNbOfPeriodicities(periodicities) ;
    Periodicity_t* periodicity = (Periodicity_t*) malloc(n*sizeof(Periodicity_t)) ;
    
    if(!periodicity) arret("Periodicities_Create(1)") ;
    Periodicities_GetPeriodicity(periodicities) = periodicity ;
  }
   
    
  {
    int n = Periodicities_GetNbOfPeriodicities(periodicities) ;
    double* vector = (double*) malloc(n*3*sizeof(double)) ;
    int i ;
    
    if(!vector) arret("Periodicities_Create(2)") ;
    
    for(i = 0 ; i < n ; i++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i ;
      
      Periodicity_GetPeriodVector(periodicity) = vector + 3*i ;
    }
  }


  /* Reading */
  {
    int n = Periodicities_GetNbOfPeriodicities(periodicities) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i ;
      char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
      char* pline ;

    
      /* Master region */
      if((pline = strstr(line,"MasterRegion"))) {
        pline = strchr(pline,'=') + 1 ;
        Periodicity_GetMasterRegion(periodicity) = atoi(pline) ;
      } else {
        arret("Periodicities_Create(3): no master region") ;
      }

    
      /* Slave region */
      if((pline = strstr(line,"SlaveRegion"))) {
        pline = strchr(pline,'=') + 1 ;
        Periodicity_GetSlaveRegion(periodicity) = atoi(pline) ;
      } else {
        arret("Periodicities_Create(4): no slave region") ;
      }
      
      
      /* Period vector */
      if((pline = strstr(line,"PeriodVector"))) {
        double* vector = Periodicity_GetPeriodVector(periodicity) ;
        int    j ;
        
        pline = strchr(pline,'=') + 1 ;
        
        for(j = 0 ; j < 3 ; j++) {
          vector[j] = strtod(pline,&pline) ;
        }
      } else {
        arret("Periodicities_Create(5): no period vector") ;
      }

    }
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(periodicities) ;
}





Graph_t*  (Periodicities_ComputeGraph)(Mesh_t* mesh)
/** Compute the graph matrix of nodes defined in periodicities.
 *  Return a pointer to Graph_t.
 **/
{
  Periodicities_t* periodicities = Mesh_GetPeriodicities(mesh) ;
  Graph_t* graph ;
  
  
  {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    /* Nb of connections per node (useful to size graph) */
    int*   nnz_no = (int*) calloc(n_no,sizeof(int)) ;
  
    if(!nnz_no) {
      arret("Periodicities_ComputeGraph(1): impossible d\'allouer la memoire") ;
    }
  
  
    /* An overestimation of nnz_no */
    {
      Element_t* elt = Mesh_GetElement(mesh) ;
      int n_elts  = Mesh_GetNbOfElements(mesh) ;
      int n_per   = Periodicities_GetNbOfPeriodicities(periodicities) ;
      int i_per ;
    
      for(i_per = 0 ; i_per < n_per ; i_per++) {
        Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i_per ;
        int masterreg = Periodicity_GetMasterRegion(periodicity) ;
        int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
        int    ie ;

        for(ie = 0 ; ie < n_elts ; ie++) {
          Element_t*  elt_i = elt + ie ;
          int eltreg_i = Element_GetRegionIndex(elt_i) ;
          
          if(eltreg_i == slavereg || eltreg_i == masterreg) {
            int nn  = Element_GetNbOfNodes(elt_i) ;
            int i ;
    
            for(i = 0 ; i < nn ; i++) {
              Node_t *node = Element_GetNode(elt_i,i) ;
              int k = Node_GetNodeIndex(node) ;
      
              nnz_no[k] += 1 ;
            }
          }
        }
      }
    }
    
    graph = Graph_Create(n_no,nnz_no) ;

    free(nnz_no) ;
  }
  
  
  /* Compute the graph for the nodes in the periodicities only */
  {
    int dim = Mesh_GetDimension(mesh) ;
    Elements_t* elts = Mesh_GetElements(mesh) ;
    Element_t* elt = Mesh_GetElement(mesh) ;
    int n_elts  = Mesh_GetNbOfElements(mesh) ;
    int n_per   = Periodicities_GetNbOfPeriodicities(periodicities) ;
    int i_per ;
    double hmin = Elements_GetMinimumSizeOfElements(elts) ;
    double tol = 0.01*fabs(hmin) ;
    
    
    for(i_per = 0 ; i_per < n_per ; i_per++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i_per ;
      int masterreg = Periodicity_GetMasterRegion(periodicity) ;
      int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
      double* periodvector = Periodicity_GetPeriodVector(periodicity) ;
      int    ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t*  elt_i = elt + ie ;
        Material_t* mat_i = Element_GetMaterial(elt_i) ;
        int eltreg_i = Element_GetRegionIndex(elt_i) ;
      
        if(!mat_i) continue ;
          
        if(eltreg_i == slavereg) {
          int nn_i  = Element_GetNbOfNodes(elt_i) ;
          int in_i ;
          
          /* Loop on slave nodes */
          for(in_i = 0 ; in_i < nn_i ; in_i++) {
            Node_t* node_i = Element_GetNode(elt_i,in_i) ;
            int n_unk = Node_GetNbOfUnknowns(node_i) ;
            int n_equ = Node_GetNbOfEquations(node_i) ;
            double* x_i = Node_GetCoordinate(node_i) ;
            Node_t* node_slave = node_i ;
            Node_t* node_master = NULL ;
            int je ;
            
            if(n_unk != n_equ) {
              arret("Periodicities_ComputeGraph(2)") ;
            }
          
            /* Find the associated master node */
            for(je = 0 ; je < n_elts ; je++) {
              Element_t*  elt_j = elt + je ;
              Material_t* mat_j = Element_GetMaterial(elt_j) ;
              int eltreg_j = Element_GetRegionIndex(elt_j) ;
      
              /* Skip material different from that of slave region */
              //if(mat_j != mat_i) continue ; /* Changed 29/02/2016 */
          
              if(eltreg_j == masterreg) {
                int nn_j  = Element_GetNbOfNodes(elt_j) ;
                int in_j ;
          
                for(in_j = 0 ; in_j < nn_j ; in_j++) {
                  Node_t* node_j = Element_GetNode(elt_j,in_j) ;
                  
                  /* The vector formed by the 2 nodes 
                   * should fit the period vector */
                  {
                    double* x_j = Node_GetCoordinate(node_j) ;
                    int notfit = 0 ;
                    int k ;
                  
                    for(k = 0 ; k < dim ; k++) {
                      double d = x_i[k] - x_j[k] - periodvector[k] ;
                      
                      if(fabs(d) > tol) {
                        notfit = 1 ; 
                        break ;
                      }
                    }
                    
                    if(notfit) continue ;
                    
                    node_master = node_j ;
                    goto Ifoundamaster ;
                  }
                }
              }
            }
            
            arret("Periodicities_ComputeGraph: no master node found") ;
            
            Ifoundamaster:
            
            {
              int i = Node_GetNodeIndex(node_slave) ;
              int j = Node_GetNodeIndex(node_master) ;
              int  degri = Graph_GetDegreeOfVertex(graph,i) ;
              int* listi = Graph_GetNeighborOfVertex(graph,i) ;
              int  degrj = Graph_GetDegreeOfVertex(graph,j) ;
              int* listj = Graph_GetNeighborOfVertex(graph,j) ;
        
	            if(i == j) continue ;
        
              /* Has j been already met? */
              {
                int met = 0 ;
                int k ;
            
                for(k = 0 ; k < degri ; k++) {
                  if(j == listi[k]) {
                    met = 1 ; 
                    break ;
                  }
                }
        
                if(met) continue ;
              }
        
              /* Not already met. So we increment with j and i */
              if(listi[degri] < 0) {
                Graph_GetDegreeOfVertex(graph,i) += 1 ;
	              listi[degri] = j ;
              } else {
                arret("Periodicities_ComputeGraph(3): not enough space") ;
              }
          
              if(listj[degrj] < 0) {
                Graph_GetDegreeOfVertex(graph,j) += 1 ;
	              listj[degrj] = i ;
              } else {
                arret("Periodicities_ComputeGraph(4): not enough space") ;
              }
            }
            
          }
        }
      }
    }
    
  }
  
  
  /* Nb of edges */
  Graph_UpdateTheNbOfEdges(graph) ;
  
  
  return(graph) ;
}



void  (Periodicities_UpdateGraph)(Mesh_t* mesh,Graph_t* graph)
/** Update the graph matrix of nodes by including 
 *  those defined in periodicities.
 **/
{
  int    n_no = Mesh_GetNbOfNodes(mesh) ;
  Graph_t*  pgraph = Periodicities_ComputeGraph(mesh) ;
  int in ;
    
  for(in = 0 ; in < n_no ; in++) {
    int  pdegrin = Graph_GetDegreeOfVertex(pgraph,in) ;
    int* plistin = Graph_GetNeighborOfVertex(pgraph,in) ;
    int* listin = Graph_GetNeighborOfVertex(graph,in) ;
    int j ;
        
    for(j = 0 ; j < pdegrin ; j++) {
      int jn = plistin[j] ;
      int  degrjn = Graph_GetDegreeOfVertex(graph,jn) ;
      int* listjn = Graph_GetNeighborOfVertex(graph,jn) ;
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
        arret("Periodicities_UpdateGraph(3): not enough space") ;
      }
          
      if(listjn[degrjn] < 0) {
        Graph_GetDegreeOfVertex(graph,jn) += 1 ;
        listjn[degrjn] = in ;
      } else {
        arret("Periodicities_UpdateGraph(4): not enough space") ;
      }

    }
  }
  
  Graph_Delete(&pgraph) ;
  
  
  /* Nb of edges */
  Graph_UpdateTheNbOfEdges(graph) ;
}



void  (Periodicities_UpdateMatrixPermutationNumbering)(Mesh_t* mesh)
/** Account for periodic mesh/BC. 
 *  Merge indexes of master and slave nodes.
 **/
{
  Periodicities_t* periodicities = Mesh_GetPeriodicities(mesh) ;
  
  if(!periodicities) return ;
  
  {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    Node_t* no  = Mesh_GetNode(mesh) ;
    Graph_t* graph = Periodicities_ComputeGraph(mesh) ;
    int in ;
    
    
    for(in = 0 ; in < n_no ; in++) {
      Node_t* node_i = no + in ;
      int n_unk = Node_GetNbOfUnknowns(node_i) ;
      int n_equ = Node_GetNbOfEquations(node_i) ;
      int  degrin = Graph_GetDegreeOfVertex(graph,in) ;
      int* listin = Graph_GetNeighborOfVertex(graph,in) ;
            
      if(n_unk != n_equ) {
        arret("Periodicities_UpdateMatrixPermutationNumbering(1)") ;
      }
      
      {
        int j ;
        
        for(j = 0 ; j < degrin ; j++) {
          int jn = listin[j] ;
          Node_t* node_j = no + jn ;

          {
            int k ;
              
            for(k = 0 ; k < n_unk ; k++) {
              int ki = Node_GetMatrixColumnIndex(node_i)[k] ;
              int kj = Node_GetMatrixColumnIndex(node_j)[k] ;

              /* The slave is i, the master is j */
              if(kj >= 0 && ki == -1) {
                Node_GetMatrixColumnIndex(node_i)[k] = kj ;
              }

              /* The slave is j, the master is i */
              if(ki >= 0 && kj == -1) {
                Node_GetMatrixColumnIndex(node_j)[k] = ki ;
              }
            }

            for(k = 0 ; k < n_equ ; k++) {
              int ki = Node_GetMatrixRowIndex(node_i)[k] ;
              int kj = Node_GetMatrixRowIndex(node_j)[k] ;
                
              if(kj >= 0 && ki == -1) {
                Node_GetMatrixRowIndex(node_i)[k] = kj ;
              }
                
              if(ki >= 0 && kj == -1) {
                Node_GetMatrixRowIndex(node_j)[k] = ki ;
              }
            }
          }
        }
      }
    }
  
    Graph_Delete(&graph) ;
  }
  
}


void  (Periodicities_ResetMatrixPermutationNumbering)(Mesh_t* mesh)
/** Set to a negative value (-1) the matrix row/column indexes 
 *  associated to the slave nodes of periodic mesh. 
 **/
{
  Periodicities_t* periodicities = Mesh_GetPeriodicities(mesh) ;
  
  if(!periodicities) return ;
  
  {
    Element_t* elt = Mesh_GetElement(mesh) ;
    int n_elts  = Mesh_GetNbOfElements(mesh) ;
    int n_per   = Periodicities_GetNbOfPeriodicities(periodicities) ;
    int i_per ;
    
    
    for(i_per = 0 ; i_per < n_per ; i_per++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i_per ;
      int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
      int    ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t*  elt_i = elt + ie ;
        Material_t* mat_i = Element_GetMaterial(elt_i) ;
        int eltreg_i = Element_GetRegionIndex(elt_i) ;
      
        if(!mat_i) continue ;
          
        if(eltreg_i == slavereg) {
          int    nn  = Element_GetNbOfNodes(elt_i) ;
          int    neq = Element_GetNbOfEquations(elt_i) ;
          int i ;

          for(i = 0 ; i < nn ; i++) {
            Node_t *node_i = Element_GetNode(elt_i,i) ;
            int   j ;

            for(j = 0 ; j < neq ; j++) {
              int ij = i*neq + j ;
              int jj ;
          
              /*  columns (unknowns) */
              if((jj = Element_GetUnknownPosition(elt_i)[ij]) >= 0) {
                Node_GetMatrixColumnIndex(node_i)[jj] = -1 ;
              }
          
              /*  rows (equations) */
              if((jj = Element_GetEquationPosition(elt_i)[ij]) >= 0) {
                Node_GetMatrixRowIndex(node_i)[jj] = -1 ;
              }
            }
          }
        }
      }
    }
  }
  
}




/* Functions not used */
#ifdef NOTDEFINED

static int**  (Periodicities_ComputeGraph1)(Mesh_t*) ;
extern void   (Periodicities_UpdateGraph1)(Mesh_t*,int**,int*) ;
static void   (Periodicities_UpdateMatrixPermutationNumbering1)(Mesh_t*) ;


int**  (Periodicities_ComputeGraph1)(Mesh_t* mesh)
/** Compute the graph matrix of nodes defined in periodicities.
 *  Return a pointer to an array of (n_no + 1) int* (graph).
 *  graph[i] is a pointer to int.
 *  graph[i][...] is the list of nodes connected to i.
 *  The nb of non zero terms of graph matrix
 *  is 2*nnz = (int) (graph[n_no] - graph[0]).
 **/
{
  Periodicities_t* periodicities = Mesh_GetPeriodicities(mesh) ;
  int    n_no = Mesh_GetNbOfNodes(mesh) ;
  int*   nnz_no ;
  int**  graph = (int**) malloc((n_no + 1)*sizeof(int*)) ;
  
  if(!graph) {
    arret("Periodicities_ComputeGraph(1)") ;
  }

  /* nnz_no is the nb of connections per node (useful to size graph) */
  nnz_no = (int*) calloc(n_no,sizeof(int)) ;
  
  if(!nnz_no) {
    arret("Periodicities_ComputeGraph(1): impossible d\'allouer la memoire") ;
  }
  
  
  /* An overestimation of nnz_no */
  {
    Element_t* elt = Mesh_GetElement(mesh) ;
    int n_elts  = Mesh_GetNbOfElements(mesh) ;
    int n_per   = Periodicities_GetNbOfPeriodicities(periodicities) ;
    int i_per ;
    
    for(i_per = 0 ; i_per < n_per ; i_per++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i_per ;
      int masterreg = Periodicity_GetMasterRegion(periodicity) ;
      int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
      int    ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t*  elt_i = elt + ie ;
        int eltreg_i = Element_GetRegionIndex(elt_i) ;
          
        if(eltreg_i == slavereg || eltreg_i == masterreg) {
          int nn  = Element_GetNbOfNodes(elt_i) ;
          int i ;
    
          for(i = 0 ; i < nn ; i++) {
            Node_t *node = Element_GetNode(elt_i,i) ;
            int k = Node_GetNodeIndex(node) ;
      
            nnz_no[k] += 1 ;
          }
        }
      }
    }
  }


  /* Allocate memory for graph (overestimation) */
  {
    int nnz ;
    int i ;
    
    nnz = 0 ;
    for(i = 0 ; i < n_no ; i++) {
      nnz += nnz_no[i] ;
    }
  
    graph[0] = (int*)  malloc(nnz*sizeof(int)) ;
  
    if(!graph[0]) {
      arret("Periodicities_ComputeGraph(2): not enough memory") ;
    }
  
    for(i = 0 ; i < n_no ; i++) {
      graph[i + 1] = graph[i] + nnz_no[i] ;
    }
  }
  
  
  /* Compute the graph for the nodes in the periodicities only */
  {
    int dim = Mesh_GetDimension(mesh) ;
    Elements_t* elts = Mesh_GetElements(mesh) ;
    Element_t* elt = Mesh_GetElement(mesh) ;
    int n_elts  = Mesh_GetNbOfElements(mesh) ;
    int n_per   = Periodicities_GetNbOfPeriodicities(periodicities) ;
    int i_per ;
    double hmin = Elements_GetMinimumSizeOfElements(elts) ;
    double tol = 0.01*hmin ;
    
    /* We re-initialize nnz_no to compute the actual space of graph */
    {
      int i ;
      
      for(i = 0 ; i < n_no ; i++) nnz_no[i] = 0 ;
    }
    
    for(i_per = 0 ; i_per < n_per ; i_per++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i_per ;
      int masterreg = Periodicity_GetMasterRegion(periodicity) ;
      int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
      double* periodvector = Periodicity_GetPeriodVector(periodicity) ;
      int    ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t*  elt_i = elt + ie ;
        Material_t* mat_i = Element_GetMaterial(elt_i) ;
        int eltreg_i = Element_GetRegionIndex(elt_i) ;
      
        if(!mat_i) continue ;
          
        if(eltreg_i == slavereg) {
          int nn_i  = Element_GetNbOfNodes(elt_i) ;
          int in_i ;
          
          /* Loop on slave nodes */
          for(in_i = 0 ; in_i < nn_i ; in_i++) {
            Node_t* node_i = Element_GetNode(elt_i,in_i) ;
            int n_unk = Node_GetNbOfUnknowns(node_i) ;
            int n_equ = Node_GetNbOfEquations(node_i) ;
            double* x_i = Node_GetCoordinate(node_i) ;
            Node_t* node_slave = node_i ;
            Node_t* node_master = NULL ;
            int je ;
            
            if(n_unk != n_equ) {
              arret("Periodicities_ComputeGraph(2)") ;
            }
          
            /* Find the associated master node */
            for(je = 0 ; je < n_elts ; je++) {
              Element_t*  elt_j = elt + je ;
              Material_t* mat_j = Element_GetMaterial(elt_j) ;
              int eltreg_j = Element_GetRegionIndex(elt_j) ;
      
              /* Skip material different from that of slave region */
              if(mat_j != mat_i) continue ;
          
              if(eltreg_j == masterreg) {
                int nn_j  = Element_GetNbOfNodes(elt_j) ;
                int in_j ;
          
                for(in_j = 0 ; in_j < nn_j ; in_j++) {
                  Node_t* node_j = Element_GetNode(elt_j,in_j) ;
                  
                  /* The vector formed by the 2 nodes 
                   * should fit the period vector */
                  {
                    double* x_j = Node_GetCoordinate(node_j) ;
                    int notfit = 0 ;
                    int k ;
                  
                    for(k = 0 ; k < dim ; k++) {
                      double d = x_i[k] - x_j[k] - periodvector[k] ;
                      
                      if(fabs(d) > tol) {notfit = 1 ; break ;}
                    }
                    
                    if(notfit) continue ;
                    
                    node_master = node_j ;
                    goto Ifoundamaster ;
                  }
                }
              }
            }
            
            arret("Periodicities_ComputeGraph: no master node found") ;
            
            Ifoundamaster:
            
            {
              int i = Node_GetNodeIndex(node_slave) ;
              int j = Node_GetNodeIndex(node_master) ;
        
	            if(i == j) continue ;
        
              /* Has j been already met? */
              {
                int met = 0 ;
                int k ;
            
                for(k = 0 ; k < nnz_no[i] ; k++) {
                  if(j == graph[i][k]) {met = 1 ; break ;}
                }
        
                if(met) continue ;
              }
        
              /* Not already met. So we increment with j and i */
              if(graph[i] + nnz_no[i] < graph[i + 1]) {
	              nnz_no[i] += 1 ;
	              graph[i][nnz_no[i] - 1] = j ;
              } else {
                arret("Periodicities_ComputeGraph(3): not enough space") ;
              }
          
              if(graph[j] + nnz_no[j] < graph[j + 1]) {
	              nnz_no[j] += 1 ;
	              graph[j][nnz_no[j] - 1] = i ;
              } else {
                arret("Periodicities_ComputeGraph(4): not enough space") ;
              }
            }
            
          }
        }
      }
    }
    
  }


  /* Compression of graph */
  {
    int i ;
    int j ;
    
    for(j = 0 , i = 0 ; i < n_no ; i++) {
      int k ;
      
      for(k = 0 ; k < nnz_no[i] ; k++) graph[0][j++] = graph[i][k] ;
    }
  
    for(i = 0 ; i < n_no ; i++) {
      graph[i + 1] = graph[i] + nnz_no[i] ;
    }
  }

  free(nnz_no) ;
    
  /* We reallocate the memory for graph[0] to save space */
  {
    int nnz2 = (int) (graph[n_no] - graph[0]) ;
    int* g0 = (int*) realloc(graph[0],nnz2*sizeof(int)) ;
      
    if(!g0) {
      arret("Periodicities_ComputeGraph(3): allocation impossible") ;
    }
  }
  
  return(graph) ;
}



void  (Periodicities_UpdateGraph1)(Mesh_t* mesh,int** graph,int* nnz_no)
/** Update the graph matrix of nodes by including 
 *  those defined in periodicities.
 **/
{
  int    n_no = Mesh_GetNbOfNodes(mesh) ;
  int**  pgraph = Periodicities_ComputeGraph1(mesh) ;
  int in ;
    
  for(in = 0 ; in < n_no ; in++) {
    int j ;
        
    for(j = 0 ; j < (int) (pgraph[in + 1] - pgraph[in]) ; j++) {
      int jn = pgraph[in][j] ;
        
      if(in == jn) continue ;
        
      /* Has jn been already met? */
      {
        int met = 0 ;
        int k ;
            
        for(k = 0 ; k < nnz_no[in] ; k++) {
          if(jn == graph[in][k]) {met = 1 ; break ;}
        }
        
        if(met) continue ;
      }
        
      /* Not already met. So we increment with jn and in */
      if(graph[in] + nnz_no[in] < graph[in + 1]) {
        nnz_no[in] += 1 ;
        graph[in][nnz_no[in] - 1] = jn ;
      } else {
        arret("Periodicities_UpdateGraph(3): not enough space") ;
      }
          
      if(graph[jn] + nnz_no[jn] < graph[jn + 1]) {
        nnz_no[jn] += 1 ;
        graph[jn][nnz_no[jn] - 1] = in ;
      } else {
        arret("Periodicities_UpdateGraph(4): not enough space") ;
      }

    }
  }
  
  free(pgraph[0]) ;
  free(pgraph) ;
}


void  (Periodicities_UpdateMatrixPermutationNumbering1)(Mesh_t* mesh)
/** Account for periodic mesh/BC. 
 *  Merge indexes of master and slave nodes.
 **/
{
  Periodicities_t* periodicities = Mesh_GetPeriodicities(mesh) ;
  
  if(!periodicities) return ;
  
  {
    int    n_no = Mesh_GetNbOfNodes(mesh) ;
    Node_t* no  = Mesh_GetNode(mesh) ;
    int** graph = Periodicities_ComputeGraph1(mesh) ;
    int in ;
    
    
    for(in = 0 ; in < n_no ; in++) {
      Node_t* node_i = no + in ;
      int n_unk = Node_GetNbOfUnknowns(node_i) ;
      int n_equ = Node_GetNbOfEquations(node_i) ;
            
      if(n_unk != n_equ) {
        arret("Periodicities_UpdateMatrixPermutationNumbering(1)") ;
      }
      
      {
        int j ;
        
        for(j = 0 ; j < (int) (graph[in + 1] - graph[in]) ; j++) {
          int jn = graph[in][j] ;
          Node_t* node_j = no + jn ;

          {
            int k ;
              
            for(k = 0 ; k < n_unk ; k++) {
              int ki = Node_GetMatrixColumnIndex(node_i)[k] ;
              int kj = Node_GetMatrixColumnIndex(node_j)[k] ;

              /* The slave is i, the master is j */
              if(kj >= 0 && ki == -1) {
                Node_GetMatrixColumnIndex(node_i)[k] = kj ;
              }

              /* The slave is j, the master is i */
              if(ki >= 0 && kj == -1) {
                Node_GetMatrixColumnIndex(node_j)[k] = ki ;
              }
            }

            for(k = 0 ; k < n_equ ; k++) {
              int ki = Node_GetMatrixRowIndex(node_i)[k] ;
              int kj = Node_GetMatrixRowIndex(node_j)[k] ;
                
              if(kj >= 0 && ki == -1) {
                Node_GetMatrixRowIndex(node_i)[k] = kj ;
              }
                
              if(ki >= 0 && kj == -1) {
                Node_GetMatrixRowIndex(node_j)[k] = ki ;
              }
            }
          }
        }
      }
    }
  
    free(graph[0]) ;
    free(graph) ;
  }
  
}
#endif
