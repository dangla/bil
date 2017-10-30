#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Models.h"
#include "ObVals.h"



ObVals_t*  ObVals_Create(DataFile_t* datafile,Mesh_t* mesh,Materials_t* mats)
{
  ObVals_t* obvals = (ObVals_t*) malloc(sizeof(ObVals_t)) ;
  
  if(!obvals) arret("ObVals_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"OBJE,Objective Variations",",",1) ;
  
  Message_Direct("Enter in %s","Objective Variations") ;
  Message_Direct("\n") ;
  
  
  /* Allocation of space for the objective values */
  {
    ObVal_t* obval = (ObVal_t*) malloc(Model_MaxNbOfEquations*sizeof(ObVal_t)) ;
    
    if(!obval) arret("ObVals_Create (7) : impossible d\'allouer la memoire") ;

    ObVals_GetObVal(obvals) = obval ;
  }
  
  
  /* Allocation of space in nodes for the index of objective values */
  {
    unsigned int n_dof = Nodes_GetNbOfDOF(Mesh_GetNodes(mesh)) ;
    unsigned short int* index = (unsigned short int*) malloc(n_dof*sizeof(unsigned short int)) ;
    
    if(!index) arret("ObVals_Create (1) : impossible d\'allouer la memoire") ;
    
    {
      Node_t* node = Mesh_GetNode(mesh) ;
      unsigned int n_nodes = Mesh_GetNbOfNodes(mesh) ;
      unsigned int i ;
      
      Node_GetObValIndex(node) = index ;
      
      for(i = 1 ; i < n_nodes ; i++) {
        Node_GetObValIndex(node + i) = Node_GetObValIndex(node + i - 1) + Node_GetNbOfEquations(node + i - 1) ;
      }
    }
    
    /* To access the objective values */
    {
      Nodes_t* nodes = Mesh_GetNodes(mesh) ;
    
      Nodes_GetObjectiveValues(nodes) = obvals ;
    }
  }


  /* Compute the nb of ObVals (ObVals_GetNbOfObVals), 
   * Initialize the names of unknown in ObVal (ObVal_GetNameOfUnknown)
   * Initialize the index of ObVal for each node unknown (Node_GetObValIndex)
   */
  {
    unsigned int n_nodes = Mesh_GetNbOfNodes(mesh) ;
    Node_t* node = Mesh_GetNode(mesh) ;
    ObVal_t* obval = ObVals_GetObVal(obvals) ;
    int n_obvals = 0 ;
    int i ;
    
    for(i = 0 ; i < (int) n_nodes ; i++) {
      int    ieq ;
      
      for(ieq = 0 ; ieq < Node_GetNbOfEquations(node + i) ; ieq++) {
        int    jeq ;
        
        for(jeq = 0 ; jeq < n_obvals ; jeq++) {
          if(!strcmp(ObVal_GetNameOfUnknown(obval + jeq),Node_GetNameOfUnknown(node + i)[ieq])) break ;
        }
        
        Node_GetObValIndex(node + i)[ieq] = jeq ;
        
        if(jeq == n_obvals) {
          n_obvals += 1 ;
          ObVal_GetNameOfUnknown(obval + jeq)  = Node_GetNameOfUnknown(node + i)[ieq] ;
          if(n_obvals > Model_MaxNbOfEquations) {
            arret("ObVals_Create (2) : trop d\'equations") ;
          }
        }
      }
    }
    ObVals_GetNbOfObVals(obvals) = n_obvals ;
  }



  /* Read the inputs of objective values */
  {
    ObVal_t* obval = ObVals_GetObVal(obvals) ;
    int n_obvals = ObVals_GetNbOfObVals(obvals) ;
    int i ;
    
    /* Default objective values */
    for(i = 0 ; i < n_obvals ; i++) {
      ObVal_GetValue(obval + i) = -1. ;
      ObVal_GetRelaxationFactor(obval + i) = 1. ;
    }
  
    for(i = 0 ; i < n_obvals ; i++) {
      char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
      int    ieq ;
      
      /* Read the name of unknown and identify its index ieq */
      {
        char  mot[ObVal_MaxLengthOfKeyWord] ;
      
        sscanf(line," %s",mot) ;
    
        //fscanf(ficd," %[^= ] =",mot) ;
        if(strlen(mot) > ObVal_MaxLengthOfKeyWord) {
          arret("ObVals_Create (3) : mot trop long") ;
        }
    
        for(ieq = 0 ; ieq < n_obvals ; ieq++) {
          if(!strcmp(mot,ObVal_GetNameOfUnknown(obval + ieq))) break ;
        }
    
        if(ieq == n_obvals) {
          arret("ObVals_Create (4) : mot_cle non connu") ;
        }
      }

      /* Read the inputs and initialize the obval */
      if(ObVal_GetValue(obval + ieq) < 0.) {
        char* pline = strchr(line,'=') + 1 ;
        
        /* la valeur */
        sscanf(pline,"%le",&ObVal_GetValue(obval + ieq)) ;
        //fscanf(ficd,"%le",&ObVal_GetValue(obval + ieq)) ;
        
        if(ObVal_GetValue(obval + ieq) < 0.) {
          arret("ObVals_Create (5) : valeur negative") ;
        }
      
        /* Default type */
        ObVal_GetType(obval + ieq) = 'a' ;
          
        /* Read the type "absolute" */
        if((pline = strstr(line,"Absolute"))) {
          ObVal_GetType(obval + ieq) = 'a' ;
        }
          
        /* Read the type "relative" */
        if((pline = strstr(line,"Relative"))) {
          ObVal_GetType(obval + ieq) = 'r' ;
        }
          
        /* Read the relaxation factor if any */
        if((pline = strstr(line,"Relaxation"))) {
          pline = strchr(pline,'=') + 1 ;
          sscanf(pline,"%le",&ObVal_GetRelaxationFactor(obval + ieq)) ;
        }
        
      } else {
        arret("ObVals_Create (6): keyword given twice") ;
      }
    }
  }
  
  DataFile_CloseFile(datafile) ;
  
  
  /* Copy objective values in those of models */
  {
    int n_obvals = ObVals_GetNbOfObVals(obvals) ;
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
        int ieq ;
      
        for(ieq = 0 ; ieq < n_obvals ; ieq++) {
          if(!strcmp(name_unk[j],ObVal_GetNameOfUnknown(obval + ieq))) break ;
        }
    
        if(ieq == n_obvals) {
          arret("ObVals_Create (7) : mot_cle non connu") ;
        }
        
        model_obval[j] = obval[ieq] ;
      }
    
    }
  }
  
  return(obvals) ;
}
