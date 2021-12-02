#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include "Mry.h"
#include "Message.h"
#include "AdjacencyList.h"




AdjacencyList_t* (AdjacencyList_Create)(int nvert,int* vert_nedges)
{
  int nedges ;
  AdjacencyList_t* adj = (AdjacencyList_t*) Mry_New(AdjacencyList_t[nvert]) ;
  
  
  /* Nb of directed edges */
  {
    int i ;
    
    nedges = 0 ;
    for(i = 0 ; i < nvert ; i++) {
      nedges += vert_nedges[i] ;
    }
  }
  
  
  /* Allocate memory for the adjacency list */
  {
    int i ;
    int* list = (int*) Mry_New(int[nedges]) ;
    
    for(i = 0 ; i < nedges ; i++) {
      list[i] = -1 ;
    }
      
    {
      int j = 0 ;
      
      for(i = 0 ; i < nvert ; i++) {
        AdjacencyList_t* adji = adj + i ;
        
        AdjacencyList_GetNbOfNeighbors(adji) = 0 ;
        AdjacencyList_GetNeighbor(adji) = list + j ;
        j += vert_nedges[i] ;
      }
    }
  }

  return(adj) ;
}


void (AdjacencyList_Delete)(void* self)
{
  AdjacencyList_t* adj = (AdjacencyList_t*) self ;
  
  {
    int* list = AdjacencyList_GetNeighbor(adj) ;
    
    if(list) {
      free(list) ;
      AdjacencyList_GetNeighbor(adj) = NULL ;
    }
  }
}

