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
#include "Mry.h"
#include "String_.h"


/* Extern functions */


Material_t* (Material_New)(void)
{
  Material_t* material   = (Material_t*) Mry_New(Material_t) ;
  
    
  /* Allocation of memory space for the material */
  {
    Material_t* mat = material ;
    
    
    /* Allocation of space for the code name of the model */
    {
      char* name = (char*) Mry_New(char[Material_MaxLengthOfKeyWord]) ;
      
      Material_GetCodeNameOfModel(mat) = name ;
    }
    
    /* The generic data */
    {
      Material_GetGenericData(mat) = NULL ;
    }
    
    /* The properties (also part of the generic data) */
    {
      double* pr = (double*) Mry_New(double[Material_MaxNbOfProperties]) ;
    
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
      char* meth = (char*) Mry_New(char[Material_MaxLengthOfKeyWord]) ;
      
      Material_GetMethod(mat) = meth ;
    }
  }
  
  return(material) ;
}



void (Material_Delete)(void* self)
{
  Material_t* material = (Material_t*) self ;
  
  {
    char* name = Material_GetCodeNameOfModel(material) ;
    
    if(name) {
      free(name) ;
      Material_GetCodeNameOfModel(material) = NULL ;
    }
  }
  
  {
    Curves_t* curves = Material_GetCurves(material) ;
    
    if(curves) {
      Curves_Delete(curves) ;
      free(curves) ;
      Material_GetCurves(material) = NULL ;
    }
  }
  
  {
    GenericData_t* gdat = Material_GetGenericData(material) ;
    
    if(gdat) {
      GenericData_Delete(gdat) ;
      free(gdat) ;
      Material_GetGenericData(material) = NULL ;
    }
  }
  
  {
    char* meth = Material_GetMethod(material) ;

    if(meth) {
      free(meth) ;
      Material_GetMethod(material) = NULL ;
    }
  }
}



void (Material_Scan)(Material_t* mat,DataFile_t* datafile,Geometry_t* geom)
{
  
    {
      char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
        
      /* Read and store the code name of the model */
      {
        char   codename[Material_MaxLengthOfKeyWord] ;
        char*  code = codename + 1 ;
        //int n = String_FindAndScanExp(line,"Model =,",","," %s",codename+1) ;
      
        if(String_Is(line,"Model",5)) {
          char* c = String_FindAndSkipToken(line,"=") ;
          
          String_ScanStringUntil(c,codename + 1,"(" String_SpaceChars) ;
          //String_Scan(line,"%*s = %s",codename + 1) ;
        } else {
          String_Scan(line,"%s",codename + 1) ;
        }
      
        if(isdigit(codename[1])) {
          codename[0] = 'm' ;
          code = codename ;
        }
      
        /* Code name of the model */
        strcpy(Material_GetCodeNameOfModel(mat),code) ;
      }
      
      
      /* Find or append a model and point to it */
      {
        char*  code = Material_GetCodeNameOfModel(mat) ;
        Models_t* usedmodels = Material_GetUsedModels(mat) ;
        Model_t* matmodel = Models_FindOrAppendModel(usedmodels,code,geom,datafile) ;
        int modind  = Models_FindModelIndex(usedmodels,code) ;
        
        Material_GetModel(mat) = matmodel ;
        Material_GetModelIndex(mat) = modind ;
      }
      
      
      /* Read and store the user-defined name of unknowns */
      {
        char* c = String_FindChar(line,'(') ;
        
        if(c) {
          Model_t* matmodel = Material_GetModel(mat) ;
          int neq = Model_GetNbOfEquations(matmodel) ;
          int i ;
          
          c[0] = ' ' ;
          for(i = 0 ; i < neq ; i++) {
            char name[Material_MaxLengthOfKeyWord] ;
            char* c1 = String_FindAnyChar(c,",)\n") ;
            
            if(c1) {
              c1[0] = ' ' ;
              String_Scan(c,"%s",name) ;
              Model_CopyNameOfUnknown(matmodel,i,name) ;
              c = c1 ;
            } else {
              Message_FatalError("Material_Scan") ;
            }
          }
        }
      }
      /* Read and store the user-defined name of equations */
      {
        char* c = String_FindChar(line,'(') ;
        c = String_FindChar(c,'(') ;
        
        if(c) {
          Model_t* matmodel = Material_GetModel(mat) ;
          int neq = Model_GetNbOfEquations(matmodel) ;
          int i ;
          
          c[0] = ' ' ;
          for(i = 0 ; i < neq ; i++) {
            char name[Material_MaxLengthOfKeyWord] ;
            char* c1 = String_FindAnyChar(c,",)\n") ;
            
            if(c1) {
              c1[0] = ' ' ;
              String_Scan(c,"%s",name) ;
              Model_CopyNameOfEquation(matmodel,i,name) ;
              c = c1 ;
            } else {
              Message_FatalError("Material_Scan") ;
            }
          }
        }
      }


      /* for compatibility with old version */
      {
        if(Material_GetModel(mat)) {
          mat->eqn = Material_GetNameOfEquation(mat) ;
          mat->inc = Material_GetNameOfUnknown(mat) ;
        }
      }


      /* Input material data */
      /* A model pointing to a null pointer serves to build curves only */
      {
        Material_GetNbOfProperties(mat) = Material_ReadProperties(mat,datafile) ;
      }
    
    
      /* for compatibility with old version */
      {
        if(Material_GetModel(mat)) {
          if(Material_GetNbOfEquations(mat) == 0) {
            Material_GetNbOfEquations(mat) = mat->neq ;
          }
        }
      }
      mat->nc = Material_GetNbOfCurves(mat) ;
    
    
      if(!Material_GetModel(mat)) {
        //Message_Warning("Material_Scan: Model not known") ;
        Message_FatalError("Material_Scan: Model not known") ;
      }

    }
}



int  (Material_ReadProperties)(Material_t* material,DataFile_t* datafile)
{
  Model_t* model = Material_GetModel(material) ;
  
  if(model) {
    Model_ReadMaterialProperties_t* readmatprop = Model_GetReadMaterialProperties(model) ;
    
    if(readmatprop) {
      int n = readmatprop(material,datafile) ;
    
      if(n > Material_MaxNbOfProperties) {
        Message_RuntimeError("Material_ReadProperties: too many properties") ;
      }
    
      return(n) ;
    }
  }
  
  Material_ScanProperties(material,datafile,NULL) ;
  
  return(0) ;
}



#if 0
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
#endif



#if 1
void (Material_ScanProperties)(Material_t* mat,DataFile_t* datafile,int (*pm)(const char*))
/** Read the material properties in the string of the file content */
{
  int    nd = Material_GetNbOfProperties(mat) ;
  short int    cont = 1 ;
  

  while(cont) {
    char   mot[Material_MaxLengthOfKeyWord] = {'\n'} ;
    char*  line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
    
    if(!line) break ;
    
    {
      //String_ScanStringUntil(line,mot,"=") ;
      //String_Scan(line,"%*[ ]%[^=]",mot) ;
      //sscanf(line," %*[ ]%[^=]",mot) ;
      int n = String_Scan(line," %[^= ]",mot) ;
      //String_Scan(line,"%*[ ]%[^=]",mot) ;
      
      if(n > Material_MaxLengthOfKeyWord) {
        arret("Material_ScanProperties: too many characters") ;
      }
    }

    /* Reading some curves */
    if(String_Is(mot,"Courbes",6) || String_Is(mot,"Curves",5)) {
      Curves_t* curves = Material_GetCurves(mat) ;
      
      Curves_ReadCurves(curves,line) ;
      
      if(Curves_GetNbOfCurves(curves) > Material_MaxNbOfCurves) {
        arret("Material_ScanProperties: too many curves") ;
      }

    /* Reading the method */
    } else if(String_Is(mot,"Method",6)) {
      char* p = String_FindChar(line,'=') ;
      
      if(p) {
        char* cr = String_FindChar(p,'\n') ;
        
        if(cr) *cr = '\0' ;
        
        p = String_FindAndSkipToken(p,"=") ;
        p = String_SkipBlankChars(p) ;
        //sscanf(p,"%s",Material_GetMethod(mat)) ;
        strcpy(Material_GetMethod(mat),p) ;
      }
      
    /* Reading the material properties and storing through pm */
    } else if(pm) {
      char* p = String_FindChar(line,'=') ;

      /* We assume that this is a property as long as "=" is found */
      if(p) {
        int i = (*pm)(mot) ;
        
        if(i >= 0) {
        
          String_Scan(p+1,"%lf",Material_GetProperty(mat) + i) ;
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
#endif



#if 0
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
        int i = Curves_ReadCurves(curves,line) ;
      
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
      int i = Curves_ReadCurves(curves,line) ;
      
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
#endif
