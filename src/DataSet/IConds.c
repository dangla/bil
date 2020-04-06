#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "IConds.h"


static IConds_t* IConds_New(const int) ;



IConds_t* IConds_New(const int n_iconds)
{
  IConds_t* iconds  = (IConds_t*) Mry_New(IConds_t) ;
    
  IConds_GetNbOfIConds(iconds) = n_iconds ;
    
    
  /* Allocation of space for the name of file of nodal values */
  {
    char* filename = (char*) Mry_New(char[IConds_MaxLengthOfFileName]) ;
      
    IConds_GetFileNameOfNodalValues(iconds) = filename ;
    IConds_GetFileNameOfNodalValues(iconds)[0] = '\0' ;
  }
  
  
  /* Allocation of space for the boundary conditions */
  if(n_iconds > 0) {
    ICond_t* icond  = (ICond_t*) Mry_New(ICond_t[n_iconds]) ;
    int i ;

    for(i = 0 ; i < n_iconds ; i++) {
      ICond_t* ic  = ICond_New() ;
      
      icond[i] = ic[0] ;
    }

    IConds_GetICond(iconds) = icond ;
  }
  
  return(iconds) ;
}





IConds_t* IConds_Create(DataFile_t* datafile,Fields_t* fields,Functions_t* functions)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"INIT,Initialization,Initial Conditions",",") ;
  int n_iconds = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  IConds_t* iconds = IConds_New(n_iconds) ;
  
  
  Message_Direct("Enter in %s","Initial Conditions") ;
  Message_Direct("\n") ;
  
  
  
  if(n_iconds == 0) {
    return(iconds) ;
  }
  
  
  /* If n_conds < 0, the IC of the nodal unknowns of the entire mesh are read
   * from a file (see below in IConds_AssignInitialConditions) */
  if(n_iconds < 0) {
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    {
      char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
      char name[IConds_MaxLengthOfFileName] ;
      int n = String_FindAndScanExp(line,"File,Fichier",","," = %s",name) ;
        
      if(n) {
      
        if(strlen(name) > IConds_MaxLengthOfFileName-1)  {
          arret("IConds_Create: name too long") ;
        }
      
        strcpy(IConds_GetFileNameOfNodalValues(iconds),name) ;
      }
    }
    
    return(iconds) ;
  }
  
  
  
  {
    int  i_ic ;
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    for(i_ic = 0 ; i_ic < n_iconds ; i_ic++) {
      ICond_t* icond = IConds_GetICond(iconds) + i_ic ;
    
      Message_Direct("Enter in %s %d","Initial Condition",i_ic+1) ;
      Message_Direct("\n") ;
      
      ICond_GetFields(icond) = fields ;
      ICond_GetFunctions(icond) = functions ;
      
      ICond_Scan(icond,datafile) ;
      
    }
  }
  
  return(iconds) ;
}





void   IConds_AssignInitialConditions(IConds_t* iconds,Mesh_t* mesh,double t)
/** Assign the initial conditions */
{
  unsigned short int dim = Mesh_GetDimension(mesh) ;
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int n_ic = IConds_GetNbOfIConds(iconds) ;
  ICond_t* ic = IConds_GetICond(iconds) ;
  unsigned int    i_ic ;
  
  
  if(n_ic < 0) {
    char* nom = IConds_GetFileNameOfNodalValues(iconds) ;
    
    /* If a file name is given we read the nodal values */
    if(nom[0]) {
      int   n_nodes = Mesh_GetNbOfNodes(mesh) ;
      Node_t*  node = Mesh_GetNode(mesh) ;
      FILE*  fic_ini ;
      int i ;

      fic_ini = fopen(nom,"r") ;
      
      if(!fic_ini) {
        arret("IConds_AssignInitialConditions(10): can't open file") ;
      }
    
      /* We assign the prescribed values to unknowns */
      for(i = 0 ; i < n_nodes ; i++) {
        double* u = Node_GetCurrentUnknown(node + i) ;
        int n_unk = Node_GetNbOfUnknowns(node + i) ;
        int j ;
      
        for(j = 0 ; j < n_unk ; j++) {
          fscanf(fic_ini,"%le",u + j) ;
        }
      }
      
      fclose(fic_ini) ;
      
    } else {
      
        arret("IConds_AssignInitialConditions(5): no valid file name") ;
        
    }
    
    return ;
  }
  


  for(i_ic = 0 ; i_ic < n_ic ; i_ic++) {
    Function_t* fn = ICond_GetFunction(ic + i_ic) ;
    double ft = (fn) ? Function_ComputeValue(fn,t) : 1. ;
    int    reg_ic = ICond_GetRegionIndex(ic + i_ic) ;
    char*  inc_ic = ICond_GetNameOfUnknown(ic + i_ic) ;
    Field_t* ch_ic = ICond_GetField(ic + i_ic) ;
    char* nom = ICond_GetFileNameOfNodalValues(ic + i_ic) ;
    double* work = NULL ;
    unsigned int    ie ;
    
    
    /* If a file name is given we read the nodal values */
    if(nom[0] && !ch_ic) {
      int   n_nodes = Mesh_GetNbOfNodes(mesh) ;
      FILE*  fic_ini ;
      int i ;
      
      /* Work table */
      work = (double*) malloc(n_nodes*sizeof(double)) ;
  
      if(!work) {
        arret("IConds_AssignInitialConditions(1): not enough memory") ;
      }

      fic_ini = fopen(nom,"r") ;
      
      if(!fic_ini) {
        arret("IConds_AssignInitialConditions(2): can't open file %s",nom) ;
      }
      
      for(i = 0 ; i < n_nodes ; i++) {
        fscanf(fic_ini,"%le",work + i) ;
      }
      
      fclose(fic_ini) ;
    }
    
    
    for(ie = 0 ; ie < n_el ; ie++) {
      int  nn = Element_GetNbOfNodes(el + ie) ;
      Material_t* mat = Element_GetMaterial(el + ie) ;
      int reg_el = Element_GetRegionIndex(el + ie) ;
      
      if(reg_el == reg_ic && mat != NULL) {
        int    neq = Material_GetNbOfEquations(mat) ;
        int    i,j ;

        /* Index of prescribed unknown */
        j = Element_FindUnknownPositionIndex(el + ie,inc_ic) ;
        if(j < 0) arret("IConds_AssignInitialConditions(1)") ;
	
        /* We assign the prescribed value to the unknown */
        for(i = 0 ; i < nn ; i++) {
          int jj = Element_GetUnknownPosition(el + ie)[i*neq + j] ;
          
          if(jj >= 0) {
            double* u = Element_GetCurrentNodalUnknown(el + ie,i) ;
            
            if(ch_ic) {
              double* x = Element_GetNodeCoordinate(el + ie,i) ;
              
              u[jj] = ft*Field_ComputeValueAtPoint(ch_ic,x,dim) ;
              
            } else if(work) {
              Node_t* node = Element_GetNode(el + ie,i) ;
              int n = Node_GetNodeIndex(node) ;
              
              u[jj] = ft*work[n] ;
    
            } else {
              
              u[jj] = 0. ;
              
            }
          }
        }
      }
    }

    free(work) ;
  }
}
