#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <assert.h>
#include "Message.h"
#include "DataFile.h"
#include "Material.h"
#include "Curves.h"


/* Extern functions */


Material_t* (Material_Create)(int n_mat)
{
  int    i ;
  Material_t* material   = (Material_t*) malloc(n_mat*sizeof(Material_t)) ;
  
  if(!material) assert(material) ;
  
    
  /* Allocation of memory space for each material */
  for(i = 0 ; i < n_mat ; i++) {
    Material_t* mat = material + i ;
    
    
    /* Allocation of space for the code name of the model */
    {
      size_t sz = Material_MaxLengthOfKeyWord*sizeof(char) ;
      char* name = (char*) malloc(sz) ;
      
      if(!name) arret("Material_Create(1)") ;
      
      Material_GetCodeNameOfModel(mat) = name ;
    }
    
    /* The generic data */
    {
      Material_GetGenericData(mat) = NULL ;
    }
    
    /* The properties (also part of the generic data) */
    {
      double* pr = (double*) calloc(Material_MaxNbOfProperties,sizeof(double)) ;
    
      if(!pr) arret("Material_Create (2)") ;
    
      Material_GetNbOfProperties(mat) = 0 ;
      Material_GetProperty(mat) = pr ;
      
      Material_AppendData(mat,Material_MaxNbOfProperties,pr,double,"Parameters") ;
    }

    /* Curves */
    {
      Material_GetCurves(mat) = Curves_Create(Material_MaxNbOfCurves) ;
    }
    
    /* for compatibility with former coding of Material_t structure */
    mat->nc = Material_GetNbOfCurves(mat) ;
    mat->cb = Material_GetCurve(mat) ;
    
    /* The method */
    {
      char* meth = (char*) malloc(Material_MaxLengthOfKeyWord*sizeof(char)) ;
      
      if(!meth) arret("Material_Create (8)") ;
      
      Material_GetMethod(mat) = meth ;
    }
  }
  
  return(material) ;
}



int  (Material_ReadProperties)(Material_t* material,DataFile_t* datafile)
{
  Model_t* model = Material_GetModel(material) ;
  
  if(model) {
    Model_ReadMaterialProperties_t* readmatprop = Model_GetReadMaterialProperties(model) ;
    int n = readmatprop(material,datafile) ;
    
    if(n > Material_MaxNbOfProperties) {
      Message_RuntimeError("Material_ReadProperties: too many properties") ;
    }
    
    return(n) ;
  }
  
  Material_ScanProperties(material,datafile,NULL) ;
  
  return(0) ;
}



void (Material_ScanProperties)(Material_t* mat,DataFile_t* datafile,int (*pm)(const char*))
/** Read the material properties in the stream file ficd */
{
  //FILE *ficd = DataFile_GetFileStream(datafile) ;
  int    nd = Material_GetNbOfProperties(mat) ;
  short int    cont = 1 ;
  
  //if(!ficd) return ;

  while(cont) {
    char   mot[Material_MaxLengthOfKeyWord] = {'\n'} ;
    char*  line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    //char* equal = (line) ? strchr(line,'=') : NULL ;
    
    //if(!equal) continue ;
    
    if(!line) break ;
    
    sscanf(line," %[^= ]",mot) ;
    //sscanf(line,"%*[ ]%[^=]",mot) ;

    /* Reading some curves */
    if(!strncmp(mot,"Courbes",6) || !strncmp(mot,"Curves",5)) {
      Curves_t* curves = Material_GetCurves(mat) ;
      Curves_FreeBuffer(curves) ;
      Curves_ReadCurves(curves,line) ;
      
      if(Curves_GetNbOfCurves(curves) > Material_MaxNbOfCurves) {
        arret("Material_ScanProperties (2) : trop de courbes") ;
      }

    /* Reading the method */
    } else if(!strncmp(mot,"Method",6)) {
      char* p = strchr(line,'=') ;
      
      if(p) {
        char* cr = strchr(p,'\n') ;
        
        if(cr) *cr = '\0' ;
        
        p += strspn(p,"= ") ;
        //sscanf(p,"%s",Material_GetMethod(mat)) ;
        strcpy(Material_GetMethod(mat),p) ;
      }
      
    /* Reading the material properties and storing through pm */
    } else if(pm) {
      char   *p = strchr(line,'=') ;

      /* We assume that this is a property as long as "=" is found */
      if(p) {
        int i = (*pm)(mot) ;
        
        if(i >= 0) {
        
          sscanf(p+1,"%lf",Material_GetProperty(mat) + i) ;
          nd = (nd > i + 1) ? nd : i + 1 ;
        
        } else {
        
          Message_RuntimeError("%s is not known",mot) ;
          
        }
      
      /* ... otherwise we stop reading */
      } else {
        
        /* go out */
        cont = 0 ;
      }
      
    } else {
      break ;
    }
    
  }

  Material_GetNbOfProperties(mat) = nd ;
  return ;
}



void (Material_ScanProperties1)(Material_t* mat,FILE *ficd,int (*pm)(const char*),int nd)
/** Read the material properties in the stream file ficd */
{
  int    id = 0,ic = 0 ;
  
  if(!ficd) return ;
  
  Material_GetNbOfProperties(mat)  += nd ;
  /* mat->nc = 0 ; */

  while(id < nd) {
    char   mot[Material_MaxLengthOfKeyWord] ;
    char   line[Material_MaxLengthOfTextLine] ;
    int long pos = ftell(ficd) ;

    if(!fgets(line,sizeof(line),ficd)) arret("Material_ScanProperties1 (1) : erreur ou fin de fichier") ;
    sscanf(line," %[^= ] =",mot) ;

    if(!strncasecmp(mot,"courbes",6)) {
      do { /* pour lire plusieurs fois "Courbes" */
        Curves_t* curves = Material_GetCurves(mat) ;
        int i ;
        Curves_FreeBuffer(curves) ;
        i = Curves_ReadCurves(curves,line) ;
      
        if(Curves_GetNbOfCurves(curves) > Material_MaxNbOfCurves) {
          arret("Material_ScanProperties1 (2) : trop de courbes") ;
        }

        /* Material_GetNbOfCurves(mat) += i ; */
        ic      += i ;
        if(ic == i) id++ ; /* on incremente qu'une fois id */
        /* position dans le fichier */
        pos = ftell(ficd) ;
        if(!fgets(line,sizeof(line),ficd)) arret("Material_ScanProperties1 (2) : erreur ou fin de fichier") ;
        sscanf(line," %[^= ] =",mot) ;
      } while(!strncasecmp(mot,"courbes",6)) ;
      /* on retourne a la ligne precedente */
      fseek(ficd,pos,SEEK_SET) ;

    } else if(!strncasecmp(mot,"Method",6)) {
      char   *p = strchr(line,'=') + 1 ;
      sscanf(p,"%s",Material_GetMethod(mat)) ;
	
    } else {
      char   *p = strchr(line,'=') + 1 ;
      int    i = (*pm)(mot) ;

      if(i >= 0 && i < nd) {
        sscanf(p,"%lf",Material_GetProperty(mat) + i) ;
      } else if(i >= nd) {
        sprintf(line,"Material_ScanProperties1 (3) : \"%s\" ne peut etre stockee",mot) ; 
        arret(line) ;
      } else {
        /* on retourne a la ligne precedente */
        fseek(ficd,pos,SEEK_SET) ;
        return ;
      }
      id++ ;
    }
  }
}



void (Material_ScanProperties2)(Material_t* mat,FILE *ficd,int (*pm)(const char*),int nd,int nc)
/** Read the material properties in the stream file ficd */
{
  int    ic = 0,id = 0 ;
  
  if(!ficd) return ;
  
  Material_GetNbOfProperties(mat) += nd ;
  /* mat->nc = 0 ; */

  while(id < nd || ic < nc) {
    char   mot[Material_MaxLengthOfKeyWord] ;
    char   line[Material_MaxLengthOfTextLine] ;
    int long pos = ftell(ficd) ;

    if(!fgets(line,sizeof(line),ficd)) arret("Material_ScanProperties2 (1) : erreur ou fin de fichier") ;
    sscanf(line," %[^= ] =",mot) ;

    if(!strncasecmp(mot,"Courbes",6)) {
      Curves_t* curves = Material_GetCurves(mat) ;
      int i ;
      Curves_FreeBuffer(curves) ;
      i = Curves_ReadCurves(curves,line) ;
      
      if(Curves_GetNbOfCurves(curves) > Material_MaxNbOfCurves) {
        arret("Material_ScanProperties2 (2) : trop de courbes") ;
      }

      if(ic >= nc) arret("Material_ScanProperties2 (2) : trop de courbe donnees") ;

      /* Material_GetNbOfCurves(mat) += i ; */
      ic      += i ;

    } else if(!strncasecmp(mot,"Method",6)) {
      char   *p = strchr(line,'=') + 1 ;
      sscanf(p,"%s",Material_GetMethod(mat)) ;
      
    } else {
      char   *p = strchr(line,'=') + 1 ;
      int    i = (*pm)(mot) ;

      if(id >= nd) arret("Material_ScanProperties2 (3) : trop de proprietes donnees") ;

      if(i >= 0 && i < nd) sscanf(p,"%lf",Material_GetProperty(mat) + i) ;
      else if(i >= nd) {
        sprintf(line,"Material_ScanProperties2 (4) : \"%s\" ne peut etre stockee",mot) ; 
        arret(line) ;
      } else {
        if(ic < nc - 1) {
          sprintf(line,"Material_ScanProperties2 (5) : il y a des courbes non lues") ; 
          arret(line) ;
        }
        /* on retourne a la ligne precedente */
        fseek(ficd,pos,SEEK_SET) ;
        return ;
      }
      id++ ;
    }
  }
}
