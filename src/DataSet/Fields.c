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
static double   champaffine(double*,int,FieldAffine_t) ;
static double   champgrille(double*,int,FieldGrid_t) ;



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



double Field_ComputeValueAtPoint(Field_t* ch,double* x,int dim)
{
  char*   type = Field_GetType(ch) ;
  double v ;

  if(!strcmp(type,"affine")) {
    FieldAffine_t* affine = (FieldAffine_t*) Field_GetFieldFormat(ch) ;
    v = champaffine(x,dim,*affine) ;
  } else if(!strcmp(type,"grid")) {
    FieldGrid_t* grille = (FieldGrid_t*) Field_GetFieldFormat(ch) ;
    v = champgrille(x,dim,*grille) ;
  } else if(!strcmp(type,"constant")) {
    FieldConstant_t* cst = (FieldConstant_t*) Field_GetFieldFormat(ch) ;
    v = FieldConstant_GetValue(cst) ;
  } else {
    arret("Field_ComputeValueAtPoint : type non prevu") ;
    return(0.) ;
  }

  return(v) ;
}


/* Intern functions */
#ifdef NOTDEFINED
void lit_grille(FieldGrid_t* grille,int dim,char* nom)
/* Lecture d'un champ sur une grille dans le fichier "nom" */
{
  char   line[Field_MaxLengthOfTextLine] ;
  FILE*   ficg ;
  int    n_x = 1,n_y = 1,n_z = 1 ;
  int    i ;

  ficg = fopen(nom,"r") ;
  if(!ficg) arret("lit_grille(1) : immpossible d ouvrir le fichier") ;

  do {
    fgets(line,sizeof(line),ficg) ;
  } while((char*) strtok(line," ") != NULL && line[0] == '#' && !feof(ficg)) ;

  if(dim == 1) {
    sscanf(line,"%d",&n_x) ;
  } else if(dim == 2) {
    sscanf(line,"%d %d",&n_x,&n_y) ;
  } else if(dim == 3) {
    sscanf(line,"%d %d %d",&n_x,&n_y,&n_z) ;
  } else {
    arret("lit_grille(2) : dimension incompatible") ;
  }

  FieldGrid_GetNbOfPointsAlongX(grille) = n_x ;
  FieldGrid_GetNbOfPointsAlongY(grille) = n_y ;
  FieldGrid_GetNbOfPointsAlongZ(grille) = n_z ;

  /* reservation de la memoire */
  FieldGrid_GetCoordinateAlongX(grille) = (double*) malloc((n_x+n_y+n_z)*sizeof(double)) ;
  if(!FieldGrid_GetCoordinateAlongX(grille)) arret("lit_grille(3) : impossible d\'allouer la memoire") ;

  FieldGrid_GetCoordinateAlongY(grille) = FieldGrid_GetCoordinateAlongX(grille) + n_x ;
  FieldGrid_GetCoordinateAlongZ(grille) = FieldGrid_GetCoordinateAlongY(grille) + n_y ;

  FieldGrid_GetValue(grille) = (double*)   malloc(n_x*n_y*n_z*sizeof(double)) ;
  if(!FieldGrid_GetValue(grille)) arret("lit_grille(4) : impossible d\'allouer la memoire") ;

  /* lecture de la grille */
  if(n_x > 0) FieldGrid_GetCoordinateAlongX(grille)[0] = 0. ;
  if(n_y > 0) FieldGrid_GetCoordinateAlongY(grille)[0] = 0. ;
  if(n_z > 0) FieldGrid_GetCoordinateAlongZ(grille)[0] = 0. ;
  
  if(dim >= 1) {
    for(i = 0 ; i < n_x ; i++) {
      fscanf(ficg,"%le",FieldGrid_GetCoordinateAlongX(grille) + i) ;
    }
  }
  
  if(dim >= 2) {
    for(i = 0 ; i < n_y ; i++) {
      fscanf(ficg,"%le",FieldGrid_GetCoordinateAlongY(grille) + i) ;
    }
  }
  
  if(dim >= 3) {
    for(i = 0 ; i < n_z ; i++) {
      fscanf(ficg,"%le",FieldGrid_GetCoordinateAlongZ(grille) + i) ;
    }
  }

  /* lecture du champ */
  for(i = 0 ; i < n_x*n_y*n_z ; i++) {
    fscanf(ficg,"%le",FieldGrid_GetValue(grille) + i) ;
  }

  fclose(ficg) ;
  return ;
}
#endif


#if 0
void Field_ReadGrid(FieldGrid_t* grid,int dim,char* filename)
/** Read a field as a grid in the file "filename" */
{
  DataFile_t* dfile = DataFile_Create(filename) ;
  char*   line = DataFile_GetTextLine(dfile) ;
  int    n_x = 1,n_y = 1,n_z = 1 ;
  
  DataFile_OpenFile(dfile,"r") ;
  
  DataFile_ReadLineFromCurrentFilePosition(dfile) ;

  if(dim >= 1) {
    int n = FieldGrid_GetNbOfPointsAlongX(grid) ;
    
    sscanf(line,"%d",&n_x) ;
    
    if(n_x <= n) {
      FieldGrid_GetNbOfPointsAlongX(grid) = n_x ;
    } else {
      arret("Field_ReadGrid: too great number") ;
    }
    
  }
  
  if(dim >= 2) {
    int n = FieldGrid_GetNbOfPointsAlongY(grid) ;
    
    sscanf(line,"%d %d",&n_x,&n_y) ;
    
    if(n_y <= n) {
      FieldGrid_GetNbOfPointsAlongY(grid) = n_y ;
    } else {
      arret("Field_ReadGrid: too great number") ;
    }
    
  }
  
  if(dim >= 3) {
    int n = FieldGrid_GetNbOfPointsAlongZ(grid) ;
    
    sscanf(line,"%d %d %d",&n_x,&n_y,&n_z) ;
    
    if(n_z <= n) {
      FieldGrid_GetNbOfPointsAlongZ(grid) = n_z ;
    } else {
      arret("Field_ReadGrid: too great number") ;
    }
    
  }
  
  if(dim >= 4) {
    arret("Field_ReadGrid(2) : dimension incompatible") ;
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
  
  return ;
}
#endif


double champaffine(double* x,int dim,chmpaffine_t ch)
{
  double v = ch.v,*x0 = ch.x,*grd = ch.g ;
  int    i ;
  
  for(i=0;i<dim;i++) v += grd[i]*(x[i] - x0[i]) ;
  
  return(v) ;
}


double champgrille(double* p,int dim,chmpgrille_t ch)
{
#define V(i,j,k) (ch.v[(i) + (j)*n_x + (k)*n_x*n_y])
  int    n_x = ch.n_x,n_y = ch.n_y,n_z = ch.n_z ;
  double* x = ch.x,*y = ch.y,*z = ch.z ;
  double x0 = p[0],y0 = p[1],z0 = p[2] ;
  double v ;
  int    ix1,iy1,iz1,ix2,iy2,iz2 ;

  if(dim > 0) {
    if(x0 <= x[0]) {
      ix1 = 0 ;
      ix2 = ix1 ;
    } else if(x0 >= x[n_x - 1]) {
      ix1 = n_x - 1 ;
      ix2 = ix1 ;
    } else {
      for(ix1=0;ix1<n_x-1;ix1++)  if(x0 < x[ix1+1]) break  ;
      ix2 = ix1 + 1 ;
    }
  }

  if(dim > 1) {
    if(y0 <= y[0]) {
      iy1 = 0 ;
      iy2 = iy1 ;
    } else if(y0 >= y[n_y - 1]) {
      iy1 = n_y - 1 ;
      iy2 = iy1 ;
    } else {
      for(iy1=0;iy1<n_y-1;iy1++)  if(y0 < y[iy1+1]) break  ;
      iy2 = iy1 + 1 ;
    }
  }

  if(dim > 2) {
    if(z0 <= z[0]) {
      iz1 = 0 ;
      iz2 = iz1 ;
    } else if(z0 >= z[n_z - 1]) {
      iz1 = n_z - 1 ;
      iz2 = iz1 ;
    } else {
      for(iz1=0;iz1<n_z-1;iz1++)  if(z0 < z[iz1+1]) break  ;
      iz2 = iz1 + 1 ;
    }
  }

  if(dim == 1) {
    double x1 = x[ix1],x2 = x[ix2] ;
    double v1 = V(ix1,0,0) ;
    double v2 = V(ix2,0,0) ;
    
    if(ix1 == ix2) v = v1 ;
    else v = v1 + (v2 - v1)*(x0 - x1)/(x2 - x1) ;

  } else if(dim == 2) {
    double x1 = x[ix1],x2 = x[ix2] ;
    double y1 = y[iy1],y2 = y[iy2] ;
    double d_x = x2 - x1,d_y = y2 - y1 ;
    double v11 = V(ix1,iy1,0) ;
    double v21 = V(ix2,iy1,0) ;
    double v12 = V(ix1,iy2,0) ;
    double v22 = V(ix2,iy2,0) ;
    
    if(ix1 == ix2 && iy1 == iy2) {
      v = v11 ;
      
    } else if(ix2 == ix1 + 1 && iy2 == iy1 + 1) {
      
      v = v11 + (v21 - v11)*(x0 - x1)/d_x + (v12 - v11)*(y0 - y1)/d_y \
      + (v22 + v11 - v21 - v12)*(x0 - x1)*(y0 - y1)/(d_x*d_y) ;
      
    } else if(ix1 == ix2) {
      
      v = v11 + (v12 - v11)*(y0 - y1)/d_y ;
      
    } else if(iy1 == iy2) {
      
      v = v11 + (v21 - v11)*(x0 - x1)/d_x ;
      
    } else {
      arret("champgrille : non prevu") ;
    }
    
    arret("champgrille : non prevu") ;
    
  } else if(dim == 3) {
    
    arret("champgrille : non prevu") ;
    
  } else {
    
    arret("champgrille : non prevu") ;
    
  }

  return(v) ;
#undef V
}
