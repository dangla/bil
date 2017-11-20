#include <stdio.h>
#include "NodeSol.h"
#include "Message.h"


/* Extern functions */


NodeSol_t* (NodeSol_GetDeepNodeSol)(NodeSol_t* nodesol,unsigned int depth)
{
  while(depth--) nodesol = NodeSol_GetPreviousNodeSol(nodesol) ;
  return(nodesol) ;
}
