#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Mry.h"
#include "String_.h"
#include "Points.h"



Points_t*  (Points_New)(const int n_points)
{
  Points_t* points = (Points_t*) Mry_New(Points_t) ;


  /* Nb of points */
  {
    Points_GetNbOfPoints(points) = n_points ;
  }
  
  
  /* Pointer to point */
  {
    Point_t* point = (Point_t*) Mry_New(Point_t[n_points]) ;
    int i ;
    
    for(i = 0 ; i < n_points ; i++) {
      Point_t* pt = Point_New() ;
    
      point[i] = pt[0] ;
      free(pt) ;
    }
    
    Points_GetPoint(points) = point ;
  }
  
  return(points) ;
}



void  (Points_Delete)(void* self)
{
  Points_t* points = (Points_t*) self ;
  
  {
    int n_points = Points_GetNbOfPoints(points) ;
    Point_t* point = Points_GetPoint(points) ;
    
    Mry_Delete(point,n_points,Point_Delete) ;
    free(point) ;
  }
}




Points_t*  (Points_Create)(DataFile_t* datafile,Mesh_t* mesh)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"POIN,Points",",") ;
  int n_points = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Points_t* points = Points_New(n_points) ;
  
  
  Message_Direct("Enter in %s","Points") ;
  Message_Direct("\n") ;
  
  
  c = String_SkipLine(c) ;


  /* Read in the input data file */
  {
    Point_t* point = Points_GetPoint(points) ;
    char* c1 = String_CopyLine(c) ;
    
    
    /* If no token "Reg" is found in the line c1 */
    if(!String_FindToken(c1,"Reg")) {
      int dim = Mesh_GetDimension(mesh) ;
      double* x = (double*) Mry_New(double[3*n_points]) ;
    
      String_ScanArray(c,dim*n_points," %lf",x) ;
      
      /* The coordinates */
      {
        int i ;
        
        for(i = 0 ; i < n_points ; i++) {
          double* coor = Point_GetCoordinate(point + i) ;
          int j ;
      
          for(j = 0 ; j < dim ; j++) {
            coor[j] = x[dim*i + j] ;
          }
        }
      }
    
      free(x) ;

    /* If a token "Reg" is found in the line */
    } else {
      int i ;
    
      DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
      for(i = 0 ; i < n_points ; i++) {
        char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
  
  
        Message_Direct("Enter in %s %d","Point",i+1) ;
        Message_Direct("\n") ;
        
        Point_Scan(point + i,line) ;
    
      }
    }
  }
  
  
  /* The enclosing element */
  {
    Point_t* point = Points_GetPoint(points) ;
    int i ;
    
    for(i = 0 ; i < n_points ; i++) {
      Point_SetEnclosingElement(point + i,mesh) ;
    }
  }
  
  return(points) ;
}

