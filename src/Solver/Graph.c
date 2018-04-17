#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "Message.h"
#include "Graph.h"
#include "AdjacencyList.h"



void*  Graph_Initialize(void* self,va_list) ;


Graph_t*  Graph_Create(int nvert,int* vert_nedges)
{
  Graph_t* graph = (Graph_t*) malloc(sizeof(Graph_t)) ;
  
  if(!graph) {
    arret("Graph_Create(1)") ;
  }
  
  
  {
    Graph_GetNbOfVertices(graph) = nvert ;
    Graph_GetNbOfEdges(graph) = 0 ;
  
  
    {
      AdjacencyList_t* adj = AdjacencyList_Create(nvert,vert_nedges) ;
    
      Graph_GetAdjacencyList(graph) = adj ;
    }
  }

  return(graph) ;
}



void*  Graph_Initialize(void* self,va_list ap)
{
  Graph_t* graph = (Graph_t*) self ;
  int nvert = va_arg(ap,int) ;
  int* vert_nedges = va_arg(ap,int*) ;
  

  Graph_GetNbOfVertices(graph) = nvert ;
  Graph_GetNbOfEdges(graph) = 0 ;
  
  
  {
    AdjacencyList_t* adj = AdjacencyList_Create(nvert,vert_nedges) ;
    
    Graph_GetAdjacencyList(graph) = adj ;
  }

  return(graph) ;
}



void Graph_Delete(Graph_t** graph)
{
  AdjacencyList_Delete(&Graph_GetAdjacencyList(*graph)) ;
  free(*graph) ;
  *graph = NULL ;
}

