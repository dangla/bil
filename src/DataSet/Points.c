#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Points.h"


static void Point_SetEnclosingElement(Point_t*,Mesh_t*) ;


Points_t*  Points_Create(DataFile_t* datafile,Mesh_t* mesh)
{
  Points_t* points = (Points_t*) malloc(sizeof(Points_t)) ;
  int dim = Mesh_GetDimension(mesh) ;
  int n_points ;
  FILE *ficd ;
  int   i ;
  
  if(!points) arret("Points_Create") ;
  
  DataFile_OpenFile(datafile,"r") ;
  
  DataFile_SetFilePositionAfterKey(datafile,"POIN,Points",",",1) ;
  
  Message_Direct("Enter in %s","Points") ;
  Message_Direct("\n") ;

  
  ficd = DataFile_GetFileStream(datafile) ;


  /* Nb of points */
  {
    fscanf(ficd,"%d",&n_points) ;
    Points_GetNbOfPoints(points) = n_points ;
  }
  
  
  /* Memory space for points */
  {
    Point_t* point = (Point_t*) malloc(n_points*sizeof(Point_t)) ;
    
    if(!point) arret("Points_Create (2)") ;
    
    Points_GetPoint(points) = point ;
  }
  
  
  /* Memory space for coordinates */
  {
    double* coor = (double*) calloc(3*n_points,sizeof(double)) ;
    
    if(!coor) arret("Points_Create (1)") ;
    
    for(i = 0 ; i < n_points ; i++) {
      Point_t* point_i = Points_GetPoint(points) + i ;
      double  *x = coor + i*3 ;
      
      Point_GetCoordinate(point_i) = x ;
    }
  }
  
  
  /* Read in the input data file */
  {
    for(i = 0 ; i < n_points ; i++) {
      Point_t* point_i = Points_GetPoint(points) + i ;
      double  *x = Point_GetCoordinate(point_i) ;
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        fscanf(ficd,"%lf",x + j) ;
      }
    }
  }
  

  /* The enclosing element */
  for(i = 0 ; i < n_points ; i++) {
    Point_t* point_i = Points_GetPoint(points) + i ;
    
    Point_SetEnclosingElement(point_i,mesh) ;
  }
  
  
  DataFile_CloseFile(datafile) ;
  
  return(points) ;
}



void Point_SetEnclosingElement(Point_t* point,Mesh_t* mesh)
/** Set a pointer to the element which encloses the point */
{
  unsigned short int dim = Mesh_GetDimension(mesh) ;
  int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  double* pt = Point_GetCoordinate(point) ;
  double d0 = 0. ;
  int ie = -1 ;
  int    i ;
  

  for(i = 0 ; i < n_el ; i++) {
    int  nn = Element_GetNbOfNodes(el + i) ;
    Material_t* mat = Element_GetMaterial(el + i) ;
    double x_s[3] = {0.,0.,0.} ;
    double d = 0. ;
    int    j ;
    
    if(!mat) continue ;
    
    /* on prend celui dont le "centre" est le plus proche du point */
    for(j = 0 ; j < dim ; j++) {
      int in ;
      
      x_s[j] = 0. ;
      
      for(in = 0 ; in < nn ; in++) {
        x_s[j] += Element_GetNodeCoordinate(el + i,in)[j] ;
      }
      
      x_s[j] /= nn ;
      x_s[j] -= pt[j] ;
    }
    
    for(j = 0 ; j < dim ; j++) d += x_s[j]*x_s[j] ;
    
    if(ie < 0 || d < d0) {d0 = d ; ie = i ;}
  }
  
  Point_GetEnclosingElement(point) = (ie < 0) ? NULL : (el + ie) ;
    
  return ;
}
