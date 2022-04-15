#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Message.h"
#include "DataFile.h"
#include "Mesh.h"
#include "Mry.h"
#include "Point.h"






Point_t*  (Point_New)(void)
{
  Point_t* point = (Point_t*) Mry_New(Point_t) ;
  
  
  /* Memory space for coordinates */
  {
    double* coor = (double*) Mry_New(double[3]) ;

    Point_GetCoordinate(point) = coor ;
  }
  
  Point_GetRegionIndex(point) = 0 ;
  
  return(point) ;
}



void (Point_Delete)(void* self)
{
  Point_t* point = (Point_t*) self ;
  
  free(Point_GetCoordinate(point)) ;
}



void (Point_Scan)(Point_t* point,char* line)
{

  /* Coordinates */
  {
    int n = String_FindAndScanExp(line,"Coor",","," = ") ;
    
    if(n) {
      double* coor = Point_GetCoordinate(point) ;
      char* c = String_GetAdvancedPosition ;
    
      String_ScanArray(c,3," %lf",coor) ;
    } else {
      arret("Point_Scan: no coordinates") ;
    }
  }
  
  /* Region */
  {
    int i ;
    int n = String_FindAndScanExp(line,"Reg",","," = %d",&i) ;
    
    if(n) {
      Point_GetRegionIndex(point) = i ;
    } else {
      arret("Point_Scan: no region") ;
    }
  }
}



void (Point_SetEnclosingElement)(Point_t* point,Mesh_t* mesh)
/** Set a pointer to the element which encloses the point */
{
  unsigned short int dim = Mesh_GetDimension(mesh) ;
  int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  double* pt = Point_GetCoordinate(point) ;
  int reg = Point_GetRegionIndex(point) ;
  double d0 = 0. ;
  int ie = -1 ;
  int    i ;
  

  for(i = 0 ; i < n_el ; i++) {
    int  nn = Element_GetNbOfNodes(el + i) ;
    Material_t* mat = Element_GetMaterial(el + i) ;
    int reg_el = Element_GetRegionIndex(el + i) ;
    double x_s[3] = {0.,0.,0.} ;
    double d = 0. ;
    
    if((reg > 0) && (reg != reg_el)) continue ;
    
    if(!mat) continue ;
    
    /* Select the element whose center is the closest to the point */
    {
      int    j ;
      
      for(j = 0 ; j < dim ; j++) {
        int in ;
      
        x_s[j] = 0. ;
      
        for(in = 0 ; in < nn ; in++) {
          x_s[j] += Element_GetNodeCoordinate(el + i,in)[j] ;
        }
      
        x_s[j] /= nn ;
        x_s[j] -= pt[j] ;
      }
    
      for(j = 0 ; j < dim ; j++) {
        d += x_s[j]*x_s[j] ;
      }
    
      if(ie < 0 || d < d0) {
        d0 = d ;
        ie = i ;
      }
    }
  }
  
  Point_GetEnclosingElement(point) = (ie < 0) ? NULL : (el + ie) ;
    
  return ;
}
