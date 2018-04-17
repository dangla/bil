#ifndef ADJACENCYLIST_H
#define ADJACENCYLIST_H


/* vacuous declarations and typedef names */
struct AdjacencyList_s  ; typedef struct AdjacencyList_s AdjacencyList_t ;


extern AdjacencyList_t* (AdjacencyList_Create)(int,int*) ;
extern void             (AdjacencyList_Delete)(AdjacencyList_t**) ;


#define AdjacencyList_GetNbOfNeighbors(adj)       ((adj)->ndest)
#define AdjacencyList_GetNeighbor(adj)            ((adj)->dest)



struct AdjacencyList_s {      /* Format */
  unsigned int  ndest ;       /* Nb of neighbors */
  int* dest ;                 /* Neighbors */
} ;

#endif
