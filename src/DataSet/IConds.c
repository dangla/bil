#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Message.h"
#include "DataFile.h"
#include "IConds.h"


IConds_t* IConds_Create(DataFile_t* datafile,Fields_t* fields,Functions_t* functions)
{
  IConds_t* iconds  = (IConds_t*) malloc(sizeof(IConds_t)) ;
  
  if(!iconds) arret("IConds_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"INIT,Initialization",",",1) ;
  
  Message_Direct("Enter in %s","Initialization") ;
  Message_Direct("\n") ;
    
    
  /* Allocation of space for the name of file of nodal values */
  {
    size_t sz = IConds_MaxLengthOfFileName*sizeof(char) ;
    char* filename = (char*) malloc(sz) ;
    
    if(!filename) arret("IConds_Create (3)") ;
      
    IConds_GetFileNameOfNodalValues(iconds) = filename ;
    IConds_GetFileNameOfNodalValues(iconds)[0] = '\0' ;
  }

  
  /* Allocation of space for the initial conditions */
  {
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    int n_iconds = atoi(line) ;
    
    IConds_GetNbOfIConds(iconds) = n_iconds ;
    
    if(n_iconds > 0) {
      IConds_GetICond(iconds) = ICond_Create(n_iconds) ;
    }

    if(n_iconds <= 0) {
      DataFile_CloseFile(datafile) ;
      return(iconds) ;
    }
  }
  
  
  {
    int  n_iconds = IConds_GetNbOfIConds(iconds) ;
    int  i_ic ;
    
    for(i_ic = 0 ; i_ic < n_iconds ; i_ic++) {
      ICond_t* icond = IConds_GetICond(iconds) + i_ic ;
      char*  line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
      char*  pline ;
      
      
      ICond_GetFields(icond) = fields ;
      ICond_GetFunctions(icond) = functions ;
    

      /* Region */
      if((pline = strstr(line,"Reg"))) {
        pline = strchr(pline,'=') + 1 ;
        ICond_GetRegionIndex(icond) = atoi(pline) ;
      } else {
        arret("IConds_Create: no Region") ;
      }
    
    
      /* Unknown */
      if((pline = strstr(line,"Unk")) || (pline = strstr(line,"Inc"))) {
        char   name_unk[IConds_MaxLengthOfKeyWord] ;
      
        pline = strchr(pline,'=') + 1 ;
        sscanf(pline,"%s",name_unk) ;
        strcpy(ICond_GetNameOfUnknown(icond),name_unk) ;
      
        if(strlen(ICond_GetNameOfUnknown(icond)) > ICond_MaxLengthOfKeyWord-1)  {
          arret("IConds_Create (6): name too long") ;
        } else if(isdigit(ICond_GetNameOfUnknown(icond)[0])) {
          if(atoi(ICond_GetNameOfUnknown(icond)) < 1) {
            arret("IConds_Create (7): index non positive") ;
          }
        }
      } else {
        arret("IConds_Create: no Unknown") ;
      }
    
    
      /* Field or File */
      ICond_GetField(icond) = NULL ;
      
      if((pline = strstr(line,"Field")) || (pline = strstr(line,"Champ"))) {
        int n_fields = Fields_GetNbOfFields(fields) ;
        int    ich ;
      
        pline = strchr(pline,'=') + 1 ;
        sscanf(pline," %d",&ich) ;
      
        if(ich > n_fields || ich < 0) {
          
          arret("IConds_Create: field not known") ;
          
        } else if(ich > 0) {
          Field_t* field = Fields_GetField(fields) ;
          
          ICond_GetField(icond) = field + ich - 1 ;
          
        } else {
          
          ICond_GetField(icond) = NULL ;
          
        }
      
      } else if((pline = strstr(line,"Fil")) || (pline = strstr(line,"Fic"))) {
        char   nom[IConds_MaxLengthOfFileName] ;
      
        pline = strchr(pline,'=') + 1 ;
        sscanf(pline,"%s",nom) ;
        strcpy(ICond_GetFileNameOfNodalValues(icond),nom) ;
      
        if(strlen(ICond_GetFileNameOfNodalValues(icond)) > ICond_MaxLengthOfFileName-1)  {
          arret("IConds_Create (6): name too long") ;
        }
    
      } else {
        arret("IConds_Create: no Field or File") ;
      }
      
    
      /* Function (not mandatory) */
      ICond_GetFunction(icond) = NULL ;
      
      if((pline = strstr(line,"Func")) || (pline = strstr(line,"Fonc"))) {
        int  n_fcts = Functions_GetNbOfFunctions(functions) ;
        int  ifn ;
      
        pline = strchr(pline,'=') + 1 ;
        sscanf(pline," %d",&ifn) ;
      
        if(ifn > n_fcts || ifn < 0) {
        
          arret("IConds_Create(11): undefined function") ;
        
        } else if(ifn > 0) {
          Function_t* fct = Functions_GetFunction(functions) ;
        
          ICond_GetFunction(icond) = fct + ifn - 1 ;
        
        } else {
        
          ICond_GetFunction(icond) = NULL ;
        
        }
      }
    }
  }
  
  DataFile_CloseFile(datafile) ;
  
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
