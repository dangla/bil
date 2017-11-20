#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "Fields.h"

/*
static void   lit_grille(FieldGrid_t* ,int,char*) ;
*/
static Field_t*       Field_Create(int) ;
static FieldAffine_t* FieldAffine_Create(void) ;
static FieldGrid_t*   FieldGrid_Create(char*,int) ;
//static void           Field_ReadGrid(FieldGrid_t*,int,char*) ;



Fields_t* Fields_Create(DataFile_t* datafile,Materials_t* materials,Geometry_t* geometry)
{
  Fields_t* fields      = (Fields_t*) malloc(sizeof(Fields_t)) ;
  int n_mats = Materials_GetNbOfMaterials(materials) ;
  int dim = Geometry_GetDimension(geometry) ;
  int n_fields = 0 ;
  int    i ;
  
  if(!fields) arret("Fields_Create") ;

  for(i = 0 ; i < n_mats ; i++) {
    Material_t* mat = Materials_GetMaterial(materials) + i ;
    Material_GetFields(mat) = fields ;
  }
  
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



Field_t* Field_Create(int n_fields)
{
  Field_t* field = (Field_t*) malloc(n_fields*sizeof(Field_t)) ;
    
  if(!field) {
    arret("Field_Create(1): allocation impossible") ;
  }

  /* Allocation of space for the type of field */
  {
    char* type = (char*) malloc(n_fields*Field_MaxLengthOfKeyWord*sizeof(char)) ;
    int i ;
    
    if(!type) {
      arret("Field_Create(2): allocation impossible") ;
    }
    
    for(i = 0 ; i < n_fields ; i++) {
      Field_GetType(field + i) = type + i*Field_MaxLengthOfKeyWord ;
      /* Default type */
      strcpy(Field_GetType(field + i),"affine") ;
    }
  }
  
  return(field) ;
}



FieldGrid_t* FieldGrid_Create(char* filename,int dim)
{
  FieldGrid_t* grid = (FieldGrid_t*) malloc(sizeof(FieldGrid_t)) ;
  int    n_x = 1,n_y = 1,n_z = 1 ;
  
  if(!grid) {
    arret("FieldGrid_Create(1): allocation impossible") ;
  }
        
  /* Read the numbers */
  {
    DataFile_t* dfile = DataFile_Create(filename) ;
    char*   line = DataFile_GetTextLine(dfile) ;
  
    DataFile_OpenFile(dfile,"r") ;
  
    DataFile_ReadLineFromCurrentFilePosition(dfile) ;

    if(dim == 1) {
      sscanf(line,"%d",&n_x) ;
    } else if(dim == 2) {
      sscanf(line,"%d %d",&n_x,&n_y) ;
    } else if(dim == 3) {
      sscanf(line,"%d %d %d",&n_x,&n_y,&n_z) ;
    } else {
      arret("FieldGrid_Create(2): dimension incompatible") ;
    }
  
    DataFile_CloseFile(dfile) ;
  
    DataFile_Delete(&dfile) ;
  }

  FieldGrid_GetNbOfPointsAlongX(grid) = n_x ;
  FieldGrid_GetNbOfPointsAlongY(grid) = n_y ;
  FieldGrid_GetNbOfPointsAlongZ(grid) = n_z ;
  
  /* Allocation of memory space for the file name */
  {
    size_t sz = Field_MaxLengthOfFileName*sizeof(char) ;
    char* name = (char*) malloc(sz) ;
    
    if(!name) {
      arret("FieldGrid_Create(2): allocation impossible") ;
    }
    
    FieldGrid_GetFileName(grid) = name ;
        
    if(strlen(filename) > Field_MaxLengthOfFileName) {
      arret("FieldGrid_Create(5): too long file name") ;
    }
    
    strcpy(name,filename) ;
  }
  
  /* Allocation of memory space for the coordinate */
  {
    size_t sz = (n_x + n_y + n_z)*sizeof(double) ;
    double* x = (double*) malloc(sz) ;
    double* y = x + n_x ;
    double* z = y + n_y ;
    
    if(!x) {
      arret("FieldGrid_Create(3): allocation impossible") ;
    }

    FieldGrid_GetCoordinateAlongX(grid) = x ;
    FieldGrid_GetCoordinateAlongY(grid) = y ;
    FieldGrid_GetCoordinateAlongZ(grid) = z ;

    /* Initialization of the coordinate */
    if(n_x > 0) FieldGrid_GetCoordinateAlongX(grid)[0] = 0. ;
    if(n_y > 0) FieldGrid_GetCoordinateAlongY(grid)[0] = 0. ;
    if(n_z > 0) FieldGrid_GetCoordinateAlongZ(grid)[0] = 0. ;
  }


  /* Allocation of memory space for the values */
  {
    size_t sz = n_x*n_y*n_z*sizeof(double) ;
    double* v = (double*) malloc(sz) ;
    
    if(!v) {
      arret("FieldGrid_Create(4): allocation impossible") ;
    }
    
    FieldGrid_GetValue(grid) = v ;
  }
  
  
  /* Read the grid */
  {
    DataFile_t* dfile = DataFile_Create(filename) ;
    char*   line = DataFile_GetTextLine(dfile) ;
  
    DataFile_OpenFile(dfile,"r") ;
  
    DataFile_ReadLineFromCurrentFilePosition(dfile) ;

    if(dim == 1) {
      sscanf(line,"%d",&n_x) ;
    } else if(dim == 2) {
      sscanf(line,"%d %d",&n_x,&n_y) ;
    } else if(dim == 3) {
      sscanf(line,"%d %d %d",&n_x,&n_y,&n_z) ;
    } else {
      arret("FieldGrid_Create(2): dimension incompatible") ;
    }
  
    if(dim >= 1) {
      double* x = FieldGrid_GetCoordinateAlongX(grid) ;
      
      DataFile_ReadDoublesFromCurrentFilePosition(dfile,x,n_x) ;
    }
  
    if(dim >= 2) {
      double* y = FieldGrid_GetCoordinateAlongY(grid) ;
      
      DataFile_ReadDoublesFromCurrentFilePosition(dfile,y,n_y) ;
    }
  
    if(dim >= 3) {
      double* z = FieldGrid_GetCoordinateAlongZ(grid) ;
      
      DataFile_ReadDoublesFromCurrentFilePosition(dfile,z,n_z) ;
    }

    /* Read the value of the field */
    {
      double* v = FieldGrid_GetValue(grid) ;
      
      DataFile_ReadDoublesFromCurrentFilePosition(dfile,v,n_x*n_y*n_z) ;
    }
  
    DataFile_CloseFile(dfile) ;
  
    DataFile_Delete(&dfile) ;
  }

  return(grid) ;
}



FieldAffine_t* FieldAffine_Create(void)
{
  FieldAffine_t* affine = (FieldAffine_t*) malloc(sizeof(FieldAffine_t)) ;
      
  if(!affine) {
    arret("FieldAffine_Create(1): allocation impossible") ;
  }

  /* Allocation of memory space for the gradient and the coordinate */
  {
    size_t sz = sizeof(double) ;
    double* grd = (double*) calloc(6,sz) ;
      
    if(!grd) {
      arret("FieldAffine_Create(2): allocation impossible") ;
    }
    
    FieldAffine_GetGradient(affine)   = grd ;
    FieldAffine_GetCoordinate(affine) = FieldAffine_GetGradient(affine) + 3 ;
  }
  
  return(affine) ;
}
