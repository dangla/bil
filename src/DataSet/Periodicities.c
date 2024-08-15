#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Math_.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Message.h"
#include "Periodicities.h"
#include "Graph.h"
#include "Mry.h"


static Graph_t*  (Periodicities_ComputeGraph)(Mesh_t*) ;




Periodicities_t* (Periodicities_New)(const int n)
{
  Periodicities_t* periodicities = (Periodicities_t*) Mry_New(Periodicities_t) ;
  
  Periodicities_GetNbOfPeriodicities(periodicities) = 0 ;
  Periodicities_GetPeriodicity(periodicities) = NULL ;

  if(n > 0) {
    Periodicity_t* periodicity = (Periodicity_t*) Mry_New(Periodicity_t[n]) ;
    int i ;
    
    for(i = 0 ; i < n ; i++) {
      Periodicity_t* period = Periodicity_New() ;
      
      periodicity[i] = period[0] ;
      free(period) ;
    }
    
    Periodicities_GetNbOfPeriodicities(periodicities) = n ;
    Periodicities_GetPeriodicity(periodicities) = periodicity ;
  }
  
  return(periodicities) ;
}



Periodicities_t* (Periodicities_Create)(DataFile_t* datafile)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c = String_FindToken(filecontent,"PERIODICITIES,Periodicities",",") ;
  int   n_per = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Periodicities_t* periodicities = Periodicities_New(n_per) ;
  
  Message_Direct("Enter in %s","Periodicities") ;
  Message_Direct("\n") ;
  
  if(n_per <= 0) {
    return(periodicities) ;
  }


  /* Scan the datafile */
  {
    int i ;
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    for(i = 0 ; i < n_per ; i++) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) + i ;
      
      Message_Direct("Enter in %s %d","Periodicity",i+1) ;
      Message_Direct("\n") ;

      Periodicity_Scan(periodicity,datafile) ;
    }
  }
  
  return(periodicities) ;
}





void (Periodicities_Delete)(void* self)
{
  Periodicities_t* periodicities = (Periodicities_t*) self ;
  
  {
    int n = Periodicities_GetNbOfPeriodicities(periodicities) ;
    
    if(n > 0) {
      Periodicity_t* periodicity = Periodicities_GetPeriodicity(periodicities) ;
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        Periodicity_t* period = periodicity + i ;
        
        Periodicity_Delete(period) ;
      }
      
      free(periodicity) ;
    }
  }
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
        //int masterreg = Periodicity_GetMasterRegion(periodicity) ;
        //int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
        char* masterreg = Periodicity_GetMasterRegionName(periodicity) ;
        char* slavereg  = Periodicity_GetSlaveRegionName(periodicity) ;
        int    ie ;

        for(ie = 0 ; ie < n_elts ; ie++) {
          Element_t*  elt_i = elt + ie ;
          //int eltreg_i = Element_GetRegionTag(elt_i) ;
          char* eltreg_i = Element_GetRegionName(elt_i) ;
          
          if(String_Is(eltreg_i,slavereg) || String_Is(eltreg_i,masterreg)) {
            int nn  = Element_GetNbOfNodes(elt_i) ;
            int i ;
    
            for(i = 0 ; i < nn ; i++) {
              Node_t* node = Element_GetNode(elt_i,i) ;
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
      //int masterreg = Periodicity_GetMasterRegion(periodicity) ;
      //int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
      char* masterreg = Periodicity_GetMasterRegionName(periodicity) ;
      char* slavereg  = Periodicity_GetSlaveRegionName(periodicity) ;
      double* periodvector = Periodicity_GetPeriodVector(periodicity) ;
      int    ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t*  elt_i = elt + ie ;
        Material_t* mat_i = Element_GetMaterial(elt_i) ;
        //int eltreg_i = Element_GetRegionTag(elt_i) ;
        char* eltreg_i = Element_GetRegionName(elt_i) ;
      
        if(!mat_i) continue ;
          
        if(String_Is(eltreg_i,slavereg)) {
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
              //int eltreg_j = Element_GetRegionTag(elt_j) ;
              char* eltreg_j = Element_GetRegionName(elt_j) ;
      
              if(!mat_j) continue ;
              /* Skip material different from that of slave region */
              //if(mat_j != mat_i) continue ; /* Changed 29/02/2016 */
          
              if(String_Is(eltreg_j,masterreg)) {
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
  
  Graph_Delete(pgraph) ;
  free(pgraph) ;
  
  
  /* Nb of edges */
  Graph_UpdateTheNbOfEdges(graph) ;
}



void  (Periodicities_UpdateMatrixRowColumnIndexes)(Mesh_t* mesh)
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
        arret("Periodicities_UpdateMatrixRowColumnIndexes(1)") ;
      }
      
      {
        int j ;
        
        for(j = 0 ; j < degrin ; j++) {
          int jn = listin[j] ;
          Node_t* node_j = no + jn ;

          {
            int k ;
              
            for(k = 0 ; k < n_unk ; k++) {
              char* unk = Node_GetNameOfUnknown(node_i)[k] ;
              
              Node_UpdateMatrixColumnIndexForPeriodicity(node_i,node_j,unk) ;
            }

            for(k = 0 ; k < n_equ ; k++) {
              char* eqn = Node_GetNameOfEquation(node_i)[k] ;
              
              Node_UpdateMatrixRowIndexForPeriodicity(node_i,node_j,eqn) ;
            }
          }
        }
      }
    }
  
    Graph_Delete(graph) ;
    free(graph) ;
  }
  
}


void  (Periodicities_EliminateMatrixRowColumnIndexes)(Mesh_t* mesh)
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
      //int slavereg  = Periodicity_GetSlaveRegion(periodicity) ;
      char* slavereg  = Periodicity_GetSlaveRegionName(periodicity) ;
      int    ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t*  elt_i = elt + ie ;
        Material_t* mat_i = Element_GetMaterial(elt_i) ;
        //int eltreg_i = Element_GetRegionTag(elt_i) ;
        char* eltreg_i = Element_GetRegionName(elt_i) ;
      
        if(!mat_i) continue ;
          
        if(String_Is(eltreg_i,slavereg)) {
          int    nn  = Element_GetNbOfNodes(elt_i) ;
          int    neq = Element_GetNbOfEquations(elt_i) ;
          int i ;

          for(i = 0 ; i < nn ; i++) {
            Node_t* node_i = Element_GetNode(elt_i,i) ;
            int   j ;

            for(j = 0 ; j < neq ; j++) {
              //int ij = i*neq + j ;
          
              /*  columns (unknowns) */
              {
                char* unk = Element_GetNameOfUnknown(elt_i)[j] ;
                //int ii = Node_FindUnknownPositionIndex(node_i,unk) ;
                //int ii = Element_GetUnknownPosition(elt_i)[ij] ;
                
                Node_EliminateMatrixColumnIndexForPeriodicity(node_i,unk) ;
              }
          
              /*  rows (equations) */
              {
                char* eqn = Element_GetNameOfEquation(elt_i)[j] ;
                //int ii = Node_FindEquationPositionIndex(node_i,eqn) ;
                //int ii = Element_GetEquationPosition(elt_i)[ij] ;
                
                Node_EliminateMatrixRowIndexForPeriodicity(node_i,eqn) ;
              }
            }
          }
        }
      }
    }
  }
}

