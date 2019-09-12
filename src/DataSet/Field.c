#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "Mry.h"
#include "String.h"
#include "Field.h"

/*
static void   lit_grille(FieldGrid_t* ,int,char*) ;
*/
//static void           Field_ReadGrid(FieldGrid_t*,int,char*) ;
static double   champaffine(double*,int,FieldAffine_t) ;
static double   champgrille(double*,int,FieldGrid_t) ;





Field_t* Field_New(void)
{
  Field_t* field = (Field_t*) Mry_New(Field_t) ;


  /* Allocation of space for the type of field */
  {
    char* type = (char*) Mry_New(char[Field_MaxLengthOfKeyWord]) ;
    
    {
      Field_GetType(field) = type ;
      /* Default type */
      strcpy(Field_GetType(field),"affine") ;
    }
  }
  
  return(field) ;
}



void Field_Delete(void* self)
{
  Field_t** pfield = (Field_t**) self ;
  Field_t*   field = *pfield ;
  
  {
    char* type = Field_GetType(field) ;
    
    if(String_Is(type,"affine")) {
      void* fieldfmt = Field_GetFieldFormat(field) ;
      
      FieldAffine_Delete(&fieldfmt) ;
    }
  }
  
  free(Field_GetType(field)) ;
  free(field) ;
  
  *pfield = NULL ;
}



void Field_Scan(Field_t* field,DataFile_t* datafile,Geometry_t* geometry)
{
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
  int dim = Geometry_GetDimension(geometry) ;

  {
    char*   type = Field_GetType(field) ;
      
      
    /* Type (if given )*/
    {
      String_FindAndScanExp(line,"Type",","," = %s",type) ;
    }


    /* Affine field */
    if(String_Is(type,"affine")) {
      FieldAffine_t* affine = FieldAffine_Create() ;
      
      Field_GetFieldFormat(field) = affine ;
      
      /* Value */
      {
        double v ;
        int n = String_FindAndScanExp(line,"Val",","," = %lf",&v) ;
        
        if(n) {
          FieldAffine_GetValue(affine) = v ;
        } else {
          arret("Field_Scan: no value") ;
        }
      }
      
      /* Gradient */
      {
        int n = String_FindAndScanExp(line,"Grad",","," = ") ;
        
        if(n) {
          char* pline = String_GetAdvancedPosition ;
          double* grd = FieldAffine_GetGradient(affine) ;
          int   j ;
          
          /* If no conversion can be made, strtod return 0 */
          for(j = 0 ; j < 3 ; j++) {
            grd[j] = strtod(pline,&pline) ;
          }
        } else {
          arret("Field_Scan: no gradient") ;
        }
      }
      
      /* Point */
      {
        int n = String_FindAndScanExp(line,"Poin",","," = ") ;
        
        if(n) {
          char* pline = String_GetAdvancedPosition ;
          double* x = FieldAffine_GetCoordinate(affine) ;
          int   j ;
        
          for(j = 0 ; j < 3 ; j++) {
            x[j] = strtod(pline,&pline) ;
          }
        } else {
          arret("Field_Scan: no point") ;
        }
      }
      
      
    /* Grid field */
    } else if(String_Is(type,"grid")) {
      char name[Field_MaxLengthOfFileName] ;
      int n = String_FindAndScanExp(line,"File",","," = %s",name) ;
        
      if(strlen(name) > Field_MaxLengthOfFileName - 1) {
        arret("Field_Scan: too long file name") ;
      }
      
      /* File */
      if(n) {
        FieldGrid_t* grid = FieldGrid_Create(name,dim) ;
        
        Field_GetFieldFormat(field) = grid ;
      } else {
        arret("Field_Scan: no file") ;
      }
      
      
    /* Constant field */
    } else if(String_Is(type,"constant")) {
      FieldConstant_t* cst = (FieldConstant_t*) Mry_New(FieldConstant_t) ;
      
      Field_GetFieldFormat(field) = cst ;
      
      /* Value */
      {
        double v ;
        int n = String_FindAndScanExp(line,"Val",","," = %lf",&v) ;
        
        if(n) {
          FieldConstant_GetValue(cst) = v ;
        } else {
          arret("Field_Scan: no value") ;
        }
      }
      
    } else {
      arret("Field_Scan: unknown type") ;
    }
  }
}



Field_t* Field_Create(int n_fields)
{
  Field_t* field = (Field_t*) Mry_New(Field_t[n_fields]) ;

  /* Allocation of space for the type of field */
  {
    int i ;
      
    for(i = 0 ; i < n_fields ; i++) {
      Field_t* fld = Field_New() ;
        
      field[i] = fld[0] ;
    }
  }
  
  return(field) ;
}



#if 0
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
#endif



#if 1
FieldGrid_t* FieldGrid_Create(char* filename,int dim1)
{
  FieldGrid_t* grid = (FieldGrid_t*) Mry_New(FieldGrid_t) ;
  int    n_x = 1,n_y = 1,n_z = 1 ;


  /* Read the numbers */
  {
    DataFile_t* dfile = DataFile_Create(filename) ;
    char*   line = DataFile_ReadLineFromCurrentFilePositionInString(dfile) ;
    
    {
      char* c = line ;
      
      if((n_x = strtol(c,&c,10)) == 0) n_x = 1 ;
      if((n_y = strtol(c,&c,10)) == 0) n_y = 1 ;
      if((n_z = strtol(c,&c,10)) == 0) n_z = 1 ;
    }
  
    DataFile_Delete(&dfile) ;
  }

  FieldGrid_GetNbOfPointsAlongX(grid) = n_x ;
  FieldGrid_GetNbOfPointsAlongY(grid) = n_y ;
  FieldGrid_GetNbOfPointsAlongZ(grid) = n_z ;


  /* Allocation of memory space for the file name */
  {
    char* name = (char*) Mry_New(char[Field_MaxLengthOfFileName]) ;
    
    FieldGrid_GetFileName(grid) = name ;
        
    if(strlen(filename) > Field_MaxLengthOfFileName) {
      arret("FieldGrid_Create: too long file name") ;
    }
    
    strcpy(name,filename) ;
  }
  
  
  /* Allocation of memory space for the coordinates */
  {
    double* x = (double*) Mry_New(double[n_x + n_y + n_z]) ;
    double* y = x + n_x ;
    double* z = y + n_y ;

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
    double* v = (double*) Mry_New(double[n_x*n_y*n_z]) ;
    
    FieldGrid_GetValue(grid) = v ;
  }
  
  
  /* Read the grid */
  {
    DataFile_t* dfile = DataFile_Create(filename) ;
    char*   line = DataFile_ReadLineFromCurrentFilePositionInString(dfile) ;
    
    line = DataFile_GetCurrentPositionInFileContent(dfile) ;
    
    /* Read the coordinates of the grid */
    {
      int n = 0 ;
    
      if(n_x > 1) n += n_x ;
      if(n_y > 1) n += n_y ;
      if(n_z > 1) n += n_z ;
  
      if(n > 1) {
        double* x = FieldGrid_GetCoordinateAlongX(grid) ;
      
        String_ReadArray(line,n," %lf",x) ;
      
        line = String_GetAdvancedPosition ;
      }
    }

    /* Read the values at the grid points */
    {
      int n = n_x*n_y*n_z ;
      double* v = FieldGrid_GetValue(grid) ;
      
      String_ReadArray(line,n," %lf",v) ;
    }
  
    DataFile_Delete(&dfile) ;
  }

  return(grid) ;
}
#endif



FieldAffine_t* FieldAffine_Create(void)
{
  FieldAffine_t* affine = (FieldAffine_t*) Mry_New(FieldAffine_t) ;

  /* Allocation of memory space for the gradient and the coordinate */
  {
    double* grd = (double*) Mry_New(double[6]) ;
    
    FieldAffine_GetGradient(affine)   = grd ;
    FieldAffine_GetCoordinate(affine) = grd + 3 ;
  }
  
  return(affine) ;
}



void FieldAffine_Delete(void* self)
{
  FieldAffine_t** pfield = (FieldAffine_t**) self ;
  FieldAffine_t*   field = *pfield ;
  
  free(FieldAffine_GetGradient(field)) ;
  free(field) ;
  
  *pfield = NULL ;
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


double champaffine(double* x,int dim,FieldAffine_t ch)
{
  double v = ch.v,*x0 = ch.x,*grd = ch.g ;
  int    i ;
  
  for(i=0;i<dim;i++) v += grd[i]*(x[i] - x0[i]) ;
  
  return(v) ;
}


double champgrille(double* p,int dim,FieldGrid_t ch)
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



/* Not used */



/* Intern functions */
#if 0
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
