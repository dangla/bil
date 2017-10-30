#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Message.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Elements.h"
#include "Fields.h"
#include "Functions.h"
#include "BConds.h"


static void     BConds_SetDefaultNameOfEquations(BConds_t*,Mesh_t*) ;

static BCond_t* BCond_Create(int) ;


BConds_t* BConds_Create(DataFile_t* datafile,Fields_t* fields,Functions_t* functions)
{
  BConds_t* bconds  = (BConds_t*) malloc(sizeof(BConds_t)) ;
  
  if(!bconds) arret("BConds_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"COND,Boundary Conditions",",",1) ;
  
  Message_Direct("Enter in %s","Boundary Conditions") ;
  Message_Direct("\n") ;
  
  
  /* Allocation of space for the boundary conditions */
  {
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    int n_bconds = atoi(line) ;
    
    BConds_GetNbOfBConds(bconds) = n_bconds ;
    
    if(n_bconds > 0) {
      BConds_GetBCond(bconds) = BCond_Create(n_bconds) ;
    }
  
    if(n_bconds <= 0) {
      DataFile_CloseFile(datafile) ;
      return(bconds) ;
    }
  }


  {
    int n_bconds = BConds_GetNbOfBConds(bconds) ;
    int i_cl ;
    
    for(i_cl = 0 ; i_cl < n_bconds ; i_cl++) {
      BCond_t* bcond = BConds_GetBCond(bconds) + i_cl ;
      char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
      char* pline ;

    
      /* Region */
      if((pline = strstr(line,"Reg"))) {
        pline = strchr(pline,'=') + 1 ;
        BCond_GetRegionIndex(bcond) = atoi(pline) ;
      } else {
        arret("BConds_Create: no Region") ;
      }
    
    
      /* Unknown */
      if((pline = strstr(line,"Unk")) || (pline = strstr(line,"Inc"))) {
        char   name_unk[BCond_MaxLengthOfKeyWord] ;
      
        pline = strchr(pline,'=') + 1 ;
      
        sscanf(pline,"%s",name_unk) ;
        strcpy(BCond_GetNameOfUnknown(bcond),name_unk) ;
      
        if(strlen(BCond_GetNameOfUnknown(bcond)) > BCond_MaxLengthOfKeyWord-1)  {
          arret("BConds_Create(6): mot trop long") ;
        } else if(isdigit(BCond_GetNameOfUnknown(bcond)[0])) {
          if(atoi(BCond_GetNameOfUnknown(bcond)) < 1) {
            arret("BConds_Create(7): numero non positif") ;
          }
        }
      } else {
        arret("BConds_Create: no Unknown") ;
      }
    
    
      /* Field */
      if((pline = strstr(line,"Field")) || (pline = strstr(line,"Champ"))) {
        int  n_fields = Fields_GetNbOfFields(fields) ;
        int  ich ;
      
        pline = strchr(pline,'=') + 1 ;
        sscanf(pline," %d",&ich) ;

        if(ich > n_fields) {
        
          arret("BConds_Create(10): champ non defini") ;
        
        } else if(ich > 0) {
          Field_t* field = Fields_GetField(fields) ;
        
          BCond_GetField(bcond) = field + ich - 1 ;
        
        } else {
        
          BCond_GetField(bcond) = NULL ;
        }
      } else {
        arret("BConds_Create: no Field") ;
      }
    
    
      /* Function */
      if((pline = strstr(line,"Func")) || (pline = strstr(line,"Fonc"))) {
        int  n_fcts = Functions_GetNbOfFunctions(functions) ;
        int  ifn ;
      
        pline = strchr(pline,'=') + 1 ;
        sscanf(pline," %d",&ifn) ;
      
        if(ifn > n_fcts) {
        
          arret("BConds_Create(11): fonction non definie") ;
        
        } else if(ifn > 0) {
          Function_t* fct = Functions_GetFunction(functions) ;
        
          BCond_GetFunction(bcond) = fct + ifn - 1 ;
        
        } else {
        
          BCond_GetFunction(bcond) = NULL ;
        
        }
      } else {
        arret("BConds_Create: no Function") ;
      }
    
    
      /* Equation (not mandatory) */
      strcpy(BCond_GetNameOfEquation(bcond)," ") ;
    
      if((pline = strstr(line,"Equ"))) {
        char   name_eqn[BCond_MaxLengthOfKeyWord] ;
      
        pline = strchr(pline,'=') + 1 ;
      
        sscanf(pline,"%s",name_eqn) ;
        strcpy(BCond_GetNameOfEquation(bcond),name_eqn) ;
      
        if(strlen(BCond_GetNameOfEquation(bcond)) > BCond_MaxLengthOfKeyWord)  {
          arret("BConds_Create(8): mot trop long") ;
        } else if(isdigit(BCond_GetNameOfEquation(bcond)[0])) {
          if(atoi(BCond_GetNameOfEquation(bcond)) < 1) {
            arret("BConds_Create(9): numero non positif") ;
          }
        }
      }
    }
  } 
  
  DataFile_CloseFile(datafile) ;
  
  return(bconds) ;
}





BCond_t* BCond_Create(int n_bconds)
{
  BCond_t* bcond = (BCond_t*) malloc(n_bconds*sizeof(BCond_t)) ;
    
  if(!bcond) arret("BCond_Create(1)") ;
    
    
  /* Allocation of space for the name of unknowns */
  {
    size_t sz = n_bconds*BCond_MaxLengthOfKeyWord*sizeof(char) ;
    char* name_unk = (char*) malloc(sz) ;
    int i ;
    
    if(!name_unk) arret("BCond_Create(2)") ;
  
    for(i = 0 ; i < n_bconds ; i++) {
      BCond_GetNameOfUnknown(bcond + i) = name_unk + i*BCond_MaxLengthOfKeyWord ;
    }
  }
    
    
  /* Allocation of space for the name of equations */
  {
    size_t sz = n_bconds*BCond_MaxLengthOfKeyWord*sizeof(char) ;
    char* name_eqn = (char*) malloc(sz) ;
    int i ;
    
    if(!name_eqn) arret("BCond_Create(3)") ;
    
    for(i = 0 ; i < n_bconds ; i++) {
      BCond_GetNameOfEquation(bcond + i) = name_eqn + i*BCond_MaxLengthOfKeyWord ;
    }
  }
  
  
  return(bcond) ;
}



void  BConds_SetDefaultNameOfEquations(BConds_t* bconds,Mesh_t* mesh)
/** Set the name of equations to be eliminated to their default names */
{
  int n_elts = Mesh_GetNbOfElements(mesh) ;
  Element_t* elt = Mesh_GetElement(mesh) ;
  BCond_t* bcond = BConds_GetBCond(bconds) ;
  int n_bconds = BConds_GetNbOfBConds(bconds) ;
  int    i_cl ;
  
  
  for(i_cl = 0 ; i_cl < n_bconds ; i_cl++) {
    BCond_t* bcond_i = bcond + i_cl ;
    int    reg_cl = BCond_GetRegionIndex(bcond_i) ;
    char*  inc_cl = BCond_GetNameOfUnknown(bcond_i) ;
    char*  eqn_cl = BCond_GetNameOfEquation(bcond_i) ;
    
    /* Set the name of equation to be eliminated */
    if(!strcmp(eqn_cl," ")) {
      int  ie ;

      for(ie = 0 ; ie < n_elts ; ie++) {
        Element_t* elt_i = elt + ie ;
        Material_t* mat = Element_GetMaterial(elt_i) ;
        int reg_el = Element_GetRegionIndex(elt_i) ;
      
        if(reg_el == reg_cl && mat != NULL) {
          /* Index of prescribed unknown */
          int jun = Element_FindUnknownPositionIndex(elt_i,inc_cl) ;
          
          if(jun < 0) {
            arret("BConds_SetDefaultNameOfEquations(1)") ;
          }
          
          {
            /* same j as that of inc_cl */
            int jeq = jun ;
        
            strcpy(eqn_cl,Element_GetNameOfEquation(elt_i)[jeq]) ;
          }
          
          break ;
        }
      }
    }
  }
}



void  BConds_ResetMatrixNumbering(BConds_t* bconds,Mesh_t* mesh)
/** Set to a negative value (-1) the matrix row/column 
 *  indexes which are prescribed by the BC */
{
  int n_elts = Mesh_GetNbOfElements(mesh) ;
  Element_t* elt = Mesh_GetElement(mesh) ;
  BCond_t* bcond = BConds_GetBCond(bconds) ;
  int n_bconds = BConds_GetNbOfBConds(bconds) ;
  int    i_cl ;
  
  
  BConds_SetDefaultNameOfEquations(bconds,mesh) ;

  
  /* Set to -1 the matrix row/column indexes.
   * So as to eliminate rows and columns due to boundary conditions */
  for(i_cl = 0 ; i_cl < n_bconds ; i_cl++) {
    BCond_t* bcond_i = bcond + i_cl ;
    int    reg_cl = BCond_GetRegionIndex(bcond_i) ;
    char*  inc_cl = BCond_GetNameOfUnknown(bcond_i) ;
    char*  eqn_cl = BCond_GetNameOfEquation(bcond_i) ;
    int    regionwasfound = 0 ;
    int    ie ;

    for(ie = 0 ; ie < n_elts ; ie++) {
      Element_t* elt_i = elt + ie ;
      Material_t* mat = Element_GetMaterial(elt_i) ;
      int reg_el = Element_GetRegionIndex(elt_i) ;
      
      if(reg_el == reg_cl && mat != NULL) {
        int    neq = Element_GetNbOfEquations(elt_i) ;
        int    nn  = Element_GetNbOfNodes(elt_i) ;
        int    i ;
        int    jun ;
        int    jeq ;

        /* Index of prescribed unknown */
        jun = Element_FindUnknownPositionIndex(elt_i,inc_cl) ;
        if(jun < 0) arret("BConds_ResetMatrixNumbering(1)") ;

        /* We eliminate the associated column */
        for(i = 0 ; i < nn ; i++) {
          Node_t* node_i = Element_GetNode(elt_i,i) ;
          int ij = i*neq + jun ;
          int jj = Element_GetUnknownPosition(elt_i)[ij] ;
          
          /* Set to a arbitrarily negative value (chosen to -1) */
          if(jj >= 0) Node_GetMatrixColumnIndex(node_i)[jj] = -1 ; 
        }

        /* Index of equation to be eliminated */
        if(!strcmp(eqn_cl," ")) { /* same j as inc_cl */
          jeq = jun ;
        } else {
          jeq = Element_FindEquationPositionIndex(elt_i,eqn_cl) ;
        }
        if(jeq < 0) arret("BConds_ResetMatrixNumbering(2)") ;

        /* We eliminate the associated row */
        for(i = 0 ; i < nn ; i++) {
          Node_t* node_i = Element_GetNode(elt_i,i) ;
          int ij = i*neq + jeq ;
          int jj = Element_GetEquationPosition(elt_i)[ij] ;
          
          /* Set to a arbitrarily negative value (chosen to -1) */
          if(jj >= 0) Node_GetMatrixRowIndex(node_i)[jj] = -1 ;
        }
        
        regionwasfound = 1 ;
      }
    }
    
    if(!regionwasfound) {
      arret("BConds_ResetMatrixNumbering(3): the region %d was not found",reg_cl) ;
    }
  }

}




void   BConds_AssignBoundaryConditions(BConds_t* bconds,Mesh_t* mesh,double t)
/** Assign the boundary conditions */
{
  unsigned short int dim = Mesh_GetDimension(mesh) ;
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int n_cl = BConds_GetNbOfBConds(bconds) ;
  BCond_t* cl = BConds_GetBCond(bconds) ;
  unsigned int    i_cl ;

  for(i_cl = 0 ; i_cl < n_cl ; i_cl++) {
    Function_t* fn = BCond_GetFunction(cl + i_cl) ;
    double ft = (fn) ? Function_ComputeValue(fn,t) : 1. ;
    int    reg_cl = BCond_GetRegionIndex(cl + i_cl) ;
    char*  inc_cl = BCond_GetNameOfUnknown(cl + i_cl) ;
    Field_t* ch_cl = BCond_GetField(cl + i_cl) ;
    unsigned int    ie ;

    for(ie = 0 ; ie < n_el ; ie++) {
      int  nn = Element_GetNbOfNodes(el + ie) ;
      Material_t* mat = Element_GetMaterial(el + ie) ;
      int reg_el = Element_GetRegionIndex(el + ie) ;
      
      if(reg_el == reg_cl && mat != NULL) {
        int    neq = Material_GetNbOfEquations(mat) ;
        int    i,j ;

        /* Index of prescribed unknown */
        j = Element_FindUnknownPositionIndex(el + ie,inc_cl) ;
        if(j < 0) arret("BConds_AssignBoundaryConditions(1)") ;
	
        /* We assign the prescribed value to the unknown */
        for(i = 0 ; i < nn ; i++) {
          int jj = Element_GetUnknownPosition(el + ie)[i*neq + j] ;
          
          if(jj >= 0) {
            double* u = Element_GetCurrentNodalUnknown(el + ie,i) ;
            
            if(ch_cl) {
              double* x = Element_GetNodeCoordinate(el + ie,i) ;
              
              u[jj] = ft*Field_ComputeValueAtPoint(ch_cl,x,dim) ;
              
            } else {
              
              u[jj] = 0. ;
              
            }
          }
        }
      }
    }
  }
}
