#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Models.h"
#include "String_.h"
#include "Mry.h"
#include "ObVals.h"






ObVals_t*  (ObVals_New)(const int n_obvals)
{
  ObVals_t* obvals = (ObVals_t*) Mry_New(ObVals_t) ;
  
  
  ObVals_GetNbOfObVals(obvals) = n_obvals ;
  
  
  /* Allocation of space for the objective values */
  if(n_obvals > 0) {
    ObVal_t* obval = (ObVal_t*) Mry_New(ObVal_t[n_obvals]) ;
    int i ;
    
    for(i = 0 ; i < n_obvals ; i++) {
      ObVal_t* ob = ObVal_New() ;
      
      obval[i] = ob[0] ;
      free(ob) ;
    }

    ObVals_GetObVal(obvals) = obval ;
  }
  
  return(obvals) ;
}



void  (ObVals_Delete)(void* self)
{
  ObVals_t* obvals = (ObVals_t*) self ;
  
  {
    int n_obvals = ObVals_GetNbOfObVals(obvals) ;
    ObVal_t* obval = ObVals_GetObVal(obvals) ;
    
    Mry_Delete(obval,n_obvals,ObVal_Delete) ;
    free(obval) ;
  }
}




ObVals_t*  (ObVals_Create)(DataFile_t* datafile,Mesh_t* mesh,Materials_t* mats)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"OBJE,Objective Variations",",") ;
  ObVals_t* obvals ;
  
  
  if(!c) {
    Message_FatalError("No Objective Variations") ;
  }
  
  
  Message_Direct("Enter in %s","Objective Variations") ;
  Message_Direct("\n") ;


  {
    Nodes_t* nodes = Mesh_GetNodes(mesh) ;
    int n_obvals = Nodes_ComputeNbOfUnknownFields(nodes) ;
    
    obvals = ObVals_New(n_obvals) ;
  }



  c = String_SkipLine(c) ;

  DataFile_SetCurrentPositionInFileContent(datafile,c) ;



  /* Scan the datafile for objective values */
  {
    ObVal_t* obval = ObVals_GetObVal(obvals) ;
    int n_obvals = ObVals_GetNbOfObVals(obvals) ;
    int i ;
  
    for(i = 0 ; i < n_obvals ; i++) {
      /* Check if a keyword is given twice */
      {
        char* line = DataFile_GetCurrentPositionInFileContent(datafile) ;
        char  name[ObVal_MaxLengthOfKeyWord] ;
        int    j ;
        
        String_Scan(line," %s",name) ;
    
        if(strlen(name) > ObVal_MaxLengthOfKeyWord) {
          arret("ObVals_Create: too long keyword") ;
        }
    
        for(j = 0 ; j < i ; j++) {
          if(!strcmp(name,ObVal_GetNameOfUnknown(obval + j))) {
            arret("ObVals_Create: keyword %s given twice",name) ;
          }
        }
      }
      
      ObVal_Scan(obval+i,datafile) ;
    }
  }
  
  


  /* Link up Nodes and ObVals
   * acces to objective values
   * initialize the objective value indexes at the nodes */
  {
    Nodes_t* nodes = Mesh_GetNodes(mesh) ;
    
    Nodes_GetObjectiveValues(nodes) = obvals ;
    Nodes_InitializeObValIndexes(nodes) ;
  }
  
  
  /* Copy objective values in those of models */
  {
    ObVal_t* obval = ObVals_GetObVal(obvals) ;
    int n_mats = Materials_GetNbOfMaterials(mats) ;
    Material_t* mat = Materials_GetMaterial(mats) ;
    int i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Model_t* model = Material_GetModel(mat + i) ;
      ObVal_t* model_obval = Model_GetObjectiveValue(model) ;
      char** name_unk = Model_GetNameOfUnknown(model) ;
      int nb_equ = Model_GetNbOfEquations(model) ;
      int j ;
      
      for(j = 0 ; j < nb_equ ; j++) {
        int k = ObVals_FindObValIndex(obvals,name_unk[j]) ;
        
        if(k >= 0) {
          
          model_obval[j] = obval[k] ;
          
        } else {
          
          arret("ObVals_Create: unknown %s not known",name_unk[j]) ;
          
        }
      }
    }
  }
  
  return(obvals) ;
}




      

int (ObVals_FindObValIndex)(ObVals_t* obvals,char* name)
{
  int n_obvals = ObVals_GetNbOfObVals(obvals) ;
  ObVal_t* obval = ObVals_GetObVal(obvals) ;
  
  {
    int i ;
      
    for(i = 0 ; i < n_obvals ; i++) {
      ObVal_t* obval_i = obval + i ;
      char* name_obval = ObVal_GetNameOfUnknown(obval_i) ;
          
      if(!strcmp(name,name_obval)) {
        return(i) ;
      }
    }
  }
    
  {
    arret("ObVals_FindObValIndex: unknown %s not found",name) ;
  }
      
  return(-1) ;
}
