#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "Fields.h"

/*
static void   lit_grille(FieldGrid_t* ,int,char*) ;
*/
//static void           Field_ReadGrid(FieldGrid_t*,int,char*) ;



Fields_t* Fields_New(const int n_fields)
{
  Fields_t* fields = (Fields_t*) Mry_New(Fields_t) ;
  
  
  Fields_GetNbOfFields(fields) = n_fields ;
    
  {
    if(n_fields > 0) {
      Field_t* field = (Field_t*) Mry_New(Field_t[n_fields]) ;
      int i ;
      
      for(i = 0 ; i < n_fields ; i++) {
        Field_t* fld = Field_New() ;
        
        field[i] = fld[0] ;
      }
      
      Fields_GetField(fields) = field ;
    }
  }
  
  return(fields) ;
}



#if 0
//Fields_t* Fields_Create(DataFile_t* datafile,Materials_t* materials,Geometry_t* geometry)
Fields_t* Fields_Create(DataFile_t* datafile,Geometry_t* geometry)
{
  Fields_t* fields      = (Fields_t*) malloc(sizeof(Fields_t)) ;
  int dim = Geometry_GetDimension(geometry) ;
  int n_fields = 0 ;
  int    i ;
  
  if(!fields) arret("Fields_Create") ;

  #if 0
  {
    int n_mats = Materials_GetNbOfMaterials(materials) ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
      Material_GetFields(mat) = fields ;
    }
  }
  #endif
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"CHMP,FLDS,Fields",",",1) ;
  
  Message_Direct("Enter in %s","Fields") ;
  Message_Direct("\n") ;
  
  {
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
  
    n_fields = atoi(line) ;
    Fields_GetNbOfFields(fields) = n_fields ;
    
    if(n_fields > 0) {
      Fields_GetField(fields) = Field_Create(n_fields) ;
    }
  }
    
    
  if(n_fields <= 0) {
    DataFile_CloseFile(datafile) ;
    return(fields) ;
  }



  /* Read the fields */
  for(i = 0 ; i < n_fields ; i++) {
    Field_t* field = Fields_GetField(fields) + i ;
    char*   type = Field_GetType(field) ;
    char*   line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    char*   pline ;
    
    if((pline = strstr(line,"Type"))) {
      pline = strchr(pline,'=') + 1 ;
      sscanf(pline," %s",type) ;
    }

    /* Affine field */
    if(!strcmp(type,"affine")) {
      FieldAffine_t* affine = FieldAffine_Create() ;
      
      Field_GetFieldFormat(field) = affine ;
      
      /* Value */
      if((pline = strstr(line,"Val"))) {
        pline = strchr(pline,'=') + 1 ;
        FieldAffine_GetValue(affine) = atof(pline) ;
      } else {
        arret("Fields_Create(5): no Value") ;
      }
      
      /* Gradient */
      if((pline = strstr(line,"Grad"))) {
        double* grd = FieldAffine_GetGradient(affine) ;
        int    j ;
        
        pline = strchr(pline,'=') + 1 ;
        
        for(j = 0 ; j < dim ; j++) {
          grd[j] = strtod(pline,&pline) ;
        }
      } else {
        arret("Fields_Create(5): no Gradient") ;
      }
      
      /* Point */
      if((pline = strstr(line,"Poin"))) {
        double* x = FieldAffine_GetCoordinate(affine) ;
        int    j ;
        
        pline = strchr(pline,'=') + 1 ;
        
        for(j = 0 ; j < dim ; j++) {
          x[j] = strtod(pline,&pline) ;
        }
      } else {
        arret("Fields_Create(5): no Point") ;
      }
      
    /* Grid field */
    } else if(!strcmp(type,"grid")) {
 
#if 0
      /* Number of points in the grid */
      if((pline = strstr(line,"Numbers"))) {
        int    n_x = 1,n_y = 1,n_z = 1 ;
        
        pline = strchr(pline,'=') + 1 ;
        
        {
          if(dim == 1) {
            sscanf(pline,"%d",&n_x) ;
          } else if(dim == 2) {
            sscanf(pline,"%d %d",&n_x,&n_y) ;
          } else if(dim == 3) {
            sscanf(pline,"%d %d %d",&n_x,&n_y,&n_z) ;
          } else {
            arret("Fields_Create(2): dimension incompatible") ;
          }
        }
        
        {
          FieldGrid_t* grille = FieldGrid_Create(n_x,n_y,n_z) ;
      
          Field_GetFieldFormat(field) = grille ;
        }
        
      } else {
        arret("Fields_Create(5): no numbers") ;
      }
#endif
      
      /* File */
      if((pline = strstr(line,"File"))) {
        char name[Field_MaxLengthOfFileName] ;
        
        pline = strchr(pline,'=') + 1 ;
        
        if(strlen(pline) > Field_MaxLengthOfFileName) {
          arret("Fields_Create(5): too long file name") ;
        }
        
        sscanf(pline," %s",name) ;
      
        {
          FieldGrid_t* grid = FieldGrid_Create(name,dim) ;
        
          Field_GetFieldFormat(field) = grid ;
        }
  
      } else {
        arret("Fields_Create(5): no File") ;
      }
      
      
    /* Constant field */
    } else if(!strcmp(type,"constant")) {
      FieldConstant_t* cst = (FieldConstant_t*) malloc(sizeof(FieldConstant_t)) ;
      
      if(!cst) {
        arret("Fields_Create(6): allocation impossible") ;
      }
      
      Field_GetFieldFormat(field) = cst ;
      
      /* Value */
      if((pline = strstr(line,"Val"))) {
        pline = strchr(pline,'=') + 1 ;
        FieldConstant_GetValue(cst) = atof(pline) ;
      } else {
        arret("Fields_Create(6): no Value") ;
      }
      
    } else {
      arret("Fields_Create(5): unknown type") ;
    }
  }


  DataFile_CloseFile(datafile) ;
  
  return(fields) ;
}
#endif



#if 1
//Fields_t* Fields_Create(DataFile_t* datafile,Materials_t* materials,Geometry_t* geometry)
Fields_t* Fields_Create(DataFile_t* datafile,Geometry_t* geometry)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"CHMP,FLDS,Fields",",") ;
  int n_fields = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Fields_t* fields  = Fields_New(n_fields) ;
  
  Message_Direct("Enter in %s","Fields") ;
  Message_Direct("\n") ;
    
    
  if(n_fields <= 0) {
    return(fields) ;
  }

  #if 0
  {
    int n_mats = Materials_GetNbOfMaterials(materials) ;
    int i ;
    
    for(i = 0 ; i < n_mats ; i++) {
      Material_t* mat = Materials_GetMaterial(materials) + i ;
      
      Material_GetFields(mat) = fields ;
    }
  }
  #endif



  {
    int i ;
    
    c = String_SkipLine(c) ;
    
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    /* Read the fields */
    for(i = 0 ; i < n_fields ; i++) {
      Field_t* field = Fields_GetField(fields) + i ;
      
  
      Message_Direct("Enter in %s %d","Field",i+1) ;
      Message_Direct("\n") ;
      
      Field_Scan(field,datafile,geometry) ;
    }
  }
  
  return(fields) ;
}
#endif
