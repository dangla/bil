#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Message.h"
#include "Mry.h"
#include "TimeStep.h"
#include "DataFile.h"
#include "ObVals.h"


static TimeStep_t*  TimeStep_New(void) ;



TimeStep_t*  TimeStep_New(void)
{
  TimeStep_t* timestep = (TimeStep_t*) Mry_New(TimeStep_t) ;

  /* default values */
  TimeStep_GetInitialTimeStep(timestep)    = 0 ;
  TimeStep_GetMaximumTimeStep(timestep)    = 0 ;
  TimeStep_GetMinimumTimeStep(timestep)    = 0 ;
  TimeStep_GetReductionFactor(timestep)    = 0.5 ;
  TimeStep_GetMaximumCommonRatio(timestep) = 1.5 ;
  TimeStep_GetObVals(timestep)             = NULL ;
  TimeStep_GetLocation(timestep)           = 0 ;
  TimeStep_SetLocationAtBegin(timestep) ;
  TimeStep_GetSequentialIndex(timestep)    = 0 ;

  return(timestep) ;
}



TimeStep_t*  TimeStep_Create(DataFile_t* datafile,ObVals_t* obvals)
{
  TimeStep_t* timestep = TimeStep_New() ;
  char*   line ;
  char*   pline ;
  
  
  TimeStep_GetObVals(timestep) = obvals ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"ALGO,TIME,Time Steps",",",1) ;
  
  Message_Direct("Enter in %s","Time Steps") ;
  Message_Direct("\n") ;
  
  
  DataFile_StoreFilePosition(datafile) ;
    
  /* Dtini */
  DataFile_MoveToStoredFilePosition(datafile) ;
  while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Dtini"))) ;
  if(line && (pline = strstr(line,"Dtini"))) {
    pline = strchr(pline,'=') + 1 ;
    TimeStep_GetInitialTimeStep(timestep) = atof(pline) ;
  } else {
    arret("TimeStep_Create : no Dtini") ;
  }
    
  /* Dtmax */
  DataFile_MoveToStoredFilePosition(datafile) ;
  while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Dtmax"))) ;
  if(line && (pline = strstr(line,"Dtmax"))) {
    pline = strchr(pline,'=') + 1 ;
    TimeStep_GetMaximumTimeStep(timestep) = atof(pline) ;
  } else {
    arret("TimeStep_Create : no Dtmax") ;
  }
    
  /* Dtmin */
  DataFile_MoveToStoredFilePosition(datafile) ;
  while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Dtmin"))) ;
  if(line && (pline = strstr(line,"Dtmin"))) {
    pline = strchr(pline,'=') + 1 ;
    TimeStep_GetMinimumTimeStep(timestep) = atof(pline) ;
  } else {
    TimeStep_GetMinimumTimeStep(timestep) = 0. ;
  }
    
  /* Fraction */
  DataFile_MoveToStoredFilePosition(datafile) ;
  while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Reduction Factor"))) ;
  if(line && (pline = strstr(line,"Reduction Factor"))) {
    pline = strchr(pline,'=') + 1 ;
    TimeStep_GetReductionFactor(timestep) = atof(pline) ;
  } else {
    TimeStep_GetReductionFactor(timestep) = 0.5 ;
  }
    
  /* Raison */
  DataFile_MoveToStoredFilePosition(datafile) ;
  while((line = DataFile_ReadLineFromCurrentFilePosition(datafile)) && !(pline = strstr(line,"Common Ratio"))) ;
  if(line && (pline = strstr(line,"Common Ratio"))) {
    pline = strchr(pline,'=') + 1 ;
    TimeStep_GetMaximumCommonRatio(timestep) = atof(pline) ;
  } else {
    TimeStep_GetMaximumCommonRatio(timestep) = 1.5 ;
  }
  
  DataFile_CloseFile(datafile) ;
  
  return(timestep) ;
}




double TimeStep_ComputeTimeStep(TimeStep_t* timestep,Nodes_t* nodes,double tn,double dtn,double t1,double t2)
/** Return the current time increment dt so that the current time is t = tn + dt.
 *  The two times tn and t should lie within the range between t1 and t2.
 *  Also set the location of the time step dt within the range between t1 and t2 
 *  as being either at the beginning, in between or at the end of the  range [t1:t2].
 **/
{
  double dt ;
  double zero = 0. ;
  //static unsigned short int step ;
  //unsigned short int debut=0,entre=1,fin=2 ;
  
  
  /* Compute the time increment dt */
  if(dtn == zero) { /* for the first step */
    //step = debut ;
    TimeStep_SetLocationAtBegin(timestep) ;
    dt = TimeStep_GetInitialTimeStep(timestep) ;
    
//  } else if(tn == t1 && step == fin) { /* beginning of a new date */
  } else if(tn == t1 && TimeStep_IsLocatedAtEnd(timestep)) { /* beginning of a new date */
    TimeStep_SetLocationAtBegin(timestep) ;
    /* To get around the fact that for very tiny dt: (t1 - tn) = 0 ! */
    dt = TimeStep_GetInitialTimeStep(timestep) ;
    
  } else {
    unsigned int n_no = Nodes_GetNbOfNodes(nodes) ;
    Node_t* no = Nodes_GetNode(nodes) ;
    ObVal_t* obval = TimeStep_GetObVal(timestep) ;
    int sequentialindex = TimeStep_GetSequentialIndex(timestep) ;
    int    obvalindex = 0 ;
    int    nodeindex = 0 ;
    double varmax = zero ;
    double rs ;
    unsigned int i ;
    
    //step = entre ;
    TimeStep_SetLocationInBetween(timestep) ;
    
    for(i = 0 ; i < n_no ; i++) {
      Node_t* nodi = no + i ;
      int neq = Node_GetNbOfEquations(nodi) ;
      //double* u_1 = Node_GetCurrentUnknown(nodi) ;
      //double* u_n = Node_GetPreviousUnknown(nodi) ;
      double* u_n = Node_GetUnknownInDistantPast(nodi,2) ;
      double* u_1 = Node_GetPreviousUnknown(nodi) ;
      int*    node_seq_ind = Node_GetSequentialIndexOfUnknown(nodi) ;
      int j ;
      
      for(j = 0 ; j < neq ; j++) {
        if(node_seq_ind[j] == sequentialindex) {
          ObVal_t* obval_j = obval + Node_GetObValIndex(nodi)[j] ;
          double val = ObVal_GetValue(obval_j) ;
          double varrel = fabs(u_1[j] - u_n[j])/val ;
        
          if(ObVal_IsRelativeValue(obval_j)) {
            if(fabs(u_n[j]) > 0.) varrel /= fabs(u_n[j]) ;
          }
        
          if(varrel > varmax) {
            varmax = varrel ;
            obvalindex = Node_GetObValIndex(nodi)[j] ;
            nodeindex = i ;
          }
        }
      }
    }
    
    rs = TimeStep_GetMaximumCommonRatio(timestep) ;
    if(varmax > zero) rs = 1./varmax ;
    if(rs > TimeStep_GetMaximumCommonRatio(timestep)) {
      rs = TimeStep_GetMaximumCommonRatio(timestep) ;
    }
    
    dt = dtn * rs ;

    if(dt > TimeStep_GetMaximumTimeStep(timestep)) {
      dt = TimeStep_GetMaximumTimeStep(timestep) ;
    }
    
    if(dt < TimeStep_GetMinimumTimeStep(timestep) && dt > zero) {
      dt = TimeStep_GetMinimumTimeStep(timestep) ;
    }

    Message_Direct("(%s[%d])",ObVal_GetNameOfUnknown(TimeStep_GetObVal(timestep) + obvalindex),nodeindex) ;
    //Message_Direct("(%s)",ObVal_GetNameOfUnknown(TimeStep_GetObVal(timestep) + obvalindex)) ;
  }
  
  /* tn + dt should not be greater than the next date */
  if(tn + dt > t2) {
    dt = t2 - tn ;
  }

  //if(tn + dt == t2) step = fin ; 
  /* End of the current date step */
  if(tn + dt == t2) {
    TimeStep_SetLocationAtEnd(timestep) ;
  }

  if(dt < zero) arret("ComputeTimeStep: dt < 0") ;
  return(dt) ;
}
