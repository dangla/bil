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
#include "String_.h"
#include "Mry.h"




IterProcess_t*  (IterProcess_New)(void)
{
  IterProcess_t* iterprocess = (IterProcess_t*) Mry_New(IterProcess_t) ;
  
  /* Iterations */
  {
    IterProcess_GetNbOfIterations(iterprocess) = 5 ;
    IterProcess_GetIterationIndex(iterprocess) = 0 ;
  }
    
  /* Tolerance */
  {
    IterProcess_GetTolerance(iterprocess) = 1.e-4 ;
    IterProcess_GetCurrentError(iterprocess) = 0 ;
  }
    
  /* Repetition */
  {
    IterProcess_GetNbOfRepetitions(iterprocess) = 0 ;
    IterProcess_GetRepetitionIndex(iterprocess) = 0 ;
  }

  return(iterprocess) ;
}



IterProcess_t*  (IterProcess_Create)(DataFile_t* datafile,ObVals_t* obvals)
{
  IterProcess_t* iterprocess = IterProcess_New() ;
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"ALGO,ITER,Iterative Process",",") ;
  
  
  if(!c) {
    Message_FatalError("No Iterative Process") ;
  }
  
  
  Message_Direct("Enter in %s","Iterative Process") ;
  Message_Direct("\n") ;



  /* Objective variations */
  IterProcess_GetObVals(iterprocess) = obvals ;



  c = String_SkipLine(c) ;
      
  //DataFile_SetCurrentPositionInFileContent(datafile,c) ;
  

  /* Iterations */
  {
    int i ;
    int n = String_FindAndScanExp(c,"Iter",","," = %d",&i) ;
    
    if(n) {
      IterProcess_GetNbOfIterations(iterprocess) = i ;
    } else {
      arret("IterProcess_Create: no Iterations") ;
    }
  }
    
  /* Tolerance */
  {
    double tol ;
    int n = String_FindAndScanExp(c,"Tol",","," = %lf",&tol) ;
    
    if(n) {
      IterProcess_GetTolerance(iterprocess) = tol ;
    } else {
      arret("IterProcess_Create: no Tolerance") ;
    }
  }
    
  /* Repetition */
  {
    int i ;
    int n = String_FindAndScanExp(c,"Rep,Rec",","," = %d",&i) ;
    
    if(n) {
      IterProcess_GetNbOfRepetitions(iterprocess) = i ;
    }
  }


  return(iterprocess) ;
}



void  (IterProcess_Delete)(void* self)
{
  IterProcess_t* iterprocess = (IterProcess_t*) self ;
}



int (IterProcess_SetCurrentError)(IterProcess_t* iterprocess,Nodes_t* nodes,Solver_t* solver)
{
  unsigned int imatrix = Solver_GetMatrixIndex(solver) ;
  double*   x     = Solver_GetSolution(solver) ;
  int       nrows = Solver_GetNbOfRows(solver) ;
  ObVal_t*  obval = IterProcess_GetObVal(iterprocess) ;
  Node_t*   node  = Nodes_GetNode(nodes) ;
  unsigned int nb_nodes = Nodes_GetNbOfNodes(nodes) ;
  double err = 0. ;
  int    obvalindex = 0 ;
  int    nodeindex  = -1 ;
  
  
  if(nrows > 0) {
    unsigned int i ;
          
    for(i = 0 ; i < nb_nodes ; i++) {
      Node_t* nodi = node + i ;
      int nin = Node_GetNbOfUnknowns(nodi) ;
      int j ;
            
      for(j = 0 ; j < nin ; j++) {
        //int   k = Node_GetMatrixColumnIndex(nodi)[j] ;
        int   k = Node_GetSelectedMatrixColumnIndexOf(nodi,j,imatrix) ;
              
        if(k >= 0) {
          ObVal_t* obval_j = obval + Node_GetObValIndex(nodi)[j] ;
          double val = ObVal_GetValue(obval_j) ;
          double re = fabs(x[k])/val ;
                
          if(ObVal_IsRelativeValue(obval_j)) {
            double* u_n = Node_GetPreviousUnknown(nodi) ;
                  
            if(fabs(u_n[j]) > 0.) re /= fabs(u_n[j]) ;
          }

          /* Sometimes re is strictly equal to zero hence the >= */
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
      return(1) ;
      //Exception_Interrupt ;
    }
  }
          
  IterProcess_GetCurrentError(iterprocess) = err ;
  IterProcess_GetObValIndexOfCurrentError(iterprocess) = obvalindex ;
  IterProcess_GetNodeIndexOfCurrentError(iterprocess) = nodeindex ;
  
  return(0) ;
}



void (IterProcess_PrintCurrentError)(IterProcess_t* iterprocess)
{
  char*  name = IterProcess_GetNameOfTheCurrentError(iterprocess) ;
  double err  = IterProcess_GetCurrentError(iterprocess) ;
  int    iter = IterProcess_GetIterationIndex(iterprocess) ;
  int    inode = IterProcess_GetNodeIndexOfCurrentError(iterprocess) ;

  if(inode >= 0) {
    Message_Direct("  (%s[%d])Error = %4.2e (%d iters)\n",name,inode,err,iter) ;
    //Message_Direct("  (%s)Error = %4.2e (%d iters)\n",name,err,iter) ;
  } else {
    Message_Direct("  (none[-])Error = %4.2e (%d iters)\n",err,iter) ;
  }
}
