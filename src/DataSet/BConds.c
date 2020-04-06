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
#include "Mry.h"


static void     BConds_SetDefaultNameOfEquations(BConds_t*,Mesh_t*) ;


BConds_t* BConds_New(const int n_bconds)
{
  BConds_t* bconds  = (BConds_t*) Mry_New(BConds_t) ;
    
  BConds_GetNbOfBConds(bconds) = n_bconds ;
  
  
  /* Allocation of space for the boundary conditions */
  if(n_bconds > 0) {
    BCond_t* bcond  = (BCond_t*) Mry_New(BCond_t[n_bconds]) ;
    int i ;

    for(i = 0 ; i < n_bconds ; i++) {
      BCond_t* bc  = BCond_New() ;
      
      bcond[i] = bc[0] ;
    }

    BConds_GetBCond(bconds) = bcond ;
  }
  
  return(bconds) ;
}



BConds_t* BConds_Create(DataFile_t* datafile,Fields_t* fields,Functions_t* functions)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"COND,Boundary Conditions",",") ;
  int n_bconds = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  BConds_t* bconds = BConds_New(n_bconds) ;
  
  
  Message_Direct("Enter in %s","Boundary Conditions") ;
  Message_Direct("\n") ;
  
  
  if(n_bconds <= 0) {
    return(bconds) ;
  }
  



  /* Scan the datafile */
  {
    int i_cl ;
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    for(i_cl = 0 ; i_cl < n_bconds ; i_cl++) {
      BCond_t* bcond = BConds_GetBCond(bconds) + i_cl ;
    
      Message_Direct("Enter in %s %d","Boundary Condition",i_cl+1) ;
      Message_Direct("\n") ;
      
      BCond_GetFields(bcond) = fields ;
      BCond_GetFunctions(bcond) = functions ;
      
      BCond_Scan(bcond,datafile) ;
    }
    
  }
  
  return(bconds) ;
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



void  BConds_EliminateMatrixRowColumnIndexes(BConds_t* bconds,Mesh_t* mesh)
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
        if(jun < 0) arret("BConds_EliminateMatrixRowColumnIndexes(1)") ;

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
        if(jeq < 0) arret("BConds_EliminateMatrixRowColumnIndexes(2)") ;

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
      arret("BConds_EliminateMatrixRowColumnIndexes(3): the region %d was not found",reg_cl) ;
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
