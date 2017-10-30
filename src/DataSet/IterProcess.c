#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Message.h"
#include "IterProcess.h"
#include "DataFile.h"
#include "DataSet.h"
#include "ObVals.h"
#include "Solver.h"
#include "Exception.h"



IterProcess_t*  IterProcess_Create(DataFile_t* datafile,ObVals_t* obvals)
{
  IterProcess_t* iterprocess = (IterProcess_t*) malloc(sizeof(IterProcess_t)) ;
  
  if(!iterprocess) arret("IterProcess_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"ALGO,ITER,Iterative Process",",",1) ;
  
  Message_Direct("Enter in %s","Iterative Process") ;
  Message_Direct("\n") ;
  
  
  DataFile_StoreFilePosition(datafile) ;


  /* Iterations */
  DataFile_MoveToStoredFilePosition(datafile) ;
  {
    char* line ;
    char* pline ;
    
    while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Iter"))) ;
    if((pline = strstr(line,"Iter"))) {
      pline = strchr(pline,'=') + 1 ;
      IterProcess_GetNbOfIterations(iterprocess) = atoi(pline) ;
    } else {
      arret("IterProcess_Create: no Iterations") ;
    }
    IterProcess_GetIterationIndex(iterprocess) = 0 ;
  }
    
  /* Tolerance */
  DataFile_MoveToStoredFilePosition(datafile) ;
  {
    char* line ;
    char* pline ;
    
    while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Tol"))) ;
    if((pline = strstr(line,"Tol"))) {
      pline = strchr(pline,'=') + 1 ;
      IterProcess_GetTolerance(iterprocess) = atof(pline) ;
    } else {
      arret("IterProcess_Create: no Tolerance") ;
    }
    IterProcess_GetCurrentError(iterprocess) = 0 ;
  }
    
  /* Repetition */
  DataFile_MoveToStoredFilePosition(datafile) ;
  {
    char* line ;
    char* pline ;
    
    while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Re"))) ;
    if((pline = strstr(line,"Re"))) {
      pline = strchr(pline,'=') + 1 ;
      IterProcess_GetNbOfRepetitions(iterprocess) = atoi(pline) ;
    } else {
      IterProcess_GetNbOfRepetitions(iterprocess) = 0 ;
    }
    IterProcess_GetRepetitionIndex(iterprocess) = 0 ;
  }
  
  /* Objective variations */
  IterProcess_GetObVals(iterprocess) = obvals ;
  
  DataFile_CloseFile(datafile) ;
  
  return(iterprocess) ;
}



void IterProcess_SetCurrentError(IterProcess_t* iterprocess,Nodes_t* nodes,Solver_t* solver)
{
  double*   x     = Solver_GetSolution(solver) ;
  int       nrows = Solver_GetNbOfRows(solver) ;
  ObVal_t*  obval = IterProcess_GetObVal(iterprocess) ;
  Node_t*   node  = Nodes_GetNode(nodes) ;
  unsigned int nb_nodes = Nodes_GetNbOfNodes(nodes) ;
  double err = 0. ;
  int    obvalindex = 0 ;
  int    nodeindex  = -1 ;
  unsigned int i ;
          
  for(i = 0 ; i < nb_nodes ; i++) {
    Node_t* nodi = node + i ;
    int nin = Node_GetNbOfUnknowns(nodi) ;
    int j ;
            
    for(j = 0 ; j < nin ; j++) {
      int   k = Node_GetMatrixColumnIndex(nodi)[j] ;
              
      if(k >= 0) {
        ObVal_t* obval_j = obval + Node_GetObValIndex(nodi)[j] ;
        double val = ObVal_GetValue(obval_j) ;
        double re = fabs(x[k])/val ;
                
        if(ObVal_GetType(obval_j) == 'r') {
          double* u_n = Node_GetPreviousUnknown(nodi) ;
                  
          if(fabs(u_n[j]) > 0.) re /= fabs(u_n[j]) ;
        }

        /* Sometimes re is strictly equal to zero */
        if(re >= err) {
          err = re ;
          obvalindex = Node_GetObValIndex(nodi)[j] ;
          nodeindex = i ;
        }
      }
    }
  }
  
  if(nrows > 0 && nodeindex < 0) {
    /* Raise an interrupt signal instead of exit */
    Message_Warning("IterProcess_SetCurrentError: can't compute error!") ;
    Exception_Interrupt ;
  }
          
  IterProcess_GetCurrentError(iterprocess) = err ;
  IterProcess_GetObValIndexOfCurrentError(iterprocess) = obvalindex ;
  IterProcess_GetNodeIndexOfCurrentError(iterprocess) = nodeindex ;
}


void IterProcess_PrintCurrentError(IterProcess_t* iterprocess)
{
  char*  name = IterProcess_GetNameOfTheCurrentError(iterprocess) ;
  double err  = IterProcess_GetCurrentError(iterprocess) ;
  int    iter = IterProcess_GetIterationIndex(iterprocess) ;
  int    inode = IterProcess_GetNodeIndexOfCurrentError(iterprocess) ;

  Message_Direct("  (%s[%d])Error = %4.2e (%d iters)\n",name,inode,err,iter) ;
  //Message_Direct("  (%s)Error = %4.2e (%d iters)\n",name,err,iter) ;
}
