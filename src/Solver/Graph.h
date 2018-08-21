#ifndef GRAPH_H
#define GRAPH_H


/* vacuous declarations and typedef names */
struct Graph_s          ; typedef struct Graph_s Graph_t ;


extern Graph_t*  Graph_Create(int,int*) ;
extern void      Graph_Delete(void*) ;


#define Graph_GetNbOfVertices(graph)              ((graph)->nvertices)
#define Graph_GetNbOfEdges(graph)                 ((graph)->nedges)
#define Graph_GetAdjacencyList(graph)             ((graph)->adj)


#define Graph_GetDegreeOfVertex(graph,i) \
        AdjacencyList_GetNbOfNeighbors(Graph_GetAdjacencyList(graph) + i)
        
#define Graph_GetNeighborOfVertex(graph,i) \
        AdjacencyList_GetNeighbor(Graph_GetAdjacencyList(graph) + i)


#define Graph_UpdateTheNbOfEdges(graph)               \
  {                                                   \
    int nvert = Graph_GetNbOfVertices(graph) ;        \
    int nedges = 0 ;                                  \
    int i ;                                           \
    for(i = 0 ; i < nvert ; i++) {                    \
      nedges += Graph_GetDegreeOfVertex(graph,i) ;    \
    }                                                 \
    Graph_GetNbOfEdges(graph) = nedges/2 ;              \
  }



#include "AdjacencyList.h"

struct Graph_s {              /* Graph */
  unsigned int  nvertices ;   /* Nb of vertices */
  unsigned int  nedges ;      /* Nb of edges */
  AdjacencyList_t* adj ;      /* Adjacency list */
} ;

#endif
