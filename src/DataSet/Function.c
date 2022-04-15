#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Mry.h"
#include "String_.h"
#include "DataFile.h"
#include "Curves.h"
#include "Function.h"



static int Function_ReadInFile(Function_t*,char*) ;


Function_t*  (Function_New)(const int n)
{
  Function_t* function = (Function_t*) Mry_New(Function_t) ;
    
  Function_GetNbOfPoints(function) = n ;

  {
    double* x = (double*) Mry_New(double[n]) ;
    
    Function_GetXValue(function) = x ;
  }

  {
    double* f = (double*) Mry_New(double[n]) ;
    
    Function_GetFValue(function) = f ;
  }
  
  return(function) ;
}



void (Function_Delete)(void* self)
{
  Function_t* function = (Function_t*) self ;
  
  {
    double* x = Function_GetXValue(function) ;
    
    if(x) {
      free(x) ;
      Function_GetXValue(function) = NULL ;
    }
  }
  
  {
    double* f = Function_GetFValue(function) ;
    
    if(f) {
      free(f) ;
      Function_GetFValue(function) = NULL ;
    }
  }
}



int (Function_Scan)(Function_t* function,DataFile_t* datafile)
{
  char* cur = DataFile_GetCurrentPositionInFileContent(datafile) ;
  char* line = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;

  line = String_SkipBlankChars(line) ;
  
  
    /* Read direct */
    if(String_Is(line,"Ntimes",1)) {
      int n_tm = 0 ;
      int n = String_FindAndScanExp(line,"N",","," = %d",&n_tm) ;
      
      if(n) {
        Function_t* fct = Function_New(n_tm) ;
        
        function[0] = fct[0] ;
        free(fct) ;
      }
      
      
      /* Read the F(T) */
      if(n_tm > 0) {
        double* t = Function_GetXValue(function) ;
        double* f = Function_GetFValue(function) ;
        //char* c = line ;
        char* c = cur ;
        
        c = String_FindToken(c,"F") ;
        
        if(c) {
          
          String_ScanArrays(c,n_tm," F( %lf ) = %lf",t,f) ;
          
          c = String_GetAdvancedPosition ;
          
          DataFile_SetCurrentPositionInFileContent(datafile,c) ;
          
        } else {
          
          arret("Function_Scan: no key F() found") ;
          
        }
        
        #if 0
        {
          int j ;
        
          for(j = 0 ; j < n_tm ; j++) {
          
            do {
              n = String_FindAndScanExp(c,"F",",","(%lf)",t + j) ;
              String_FindAndScanExp(c," =",","," %lf",f + j) ;
            
              if(!n) {
                c = DataFile_ReadLineFromCurrentFilePositionInString(datafile) ;
              }
            } while(!n && c) ;
          
            if(!c) {
              arret("Function_Scan: not enough data") ;
            }
          
          }
        }
        #endif
      }
      
      return(1) ;
    }



    /* Read in a file */
    if(String_Is(line,"File",2)) {
      char name[Function_MaxLengthOfFileName] ;
      int n_fn ;
      
      String_FindAndScanExp(line,"Fi",","," = %s",name) ;
      
      if(strlen(name) > Function_MaxLengthOfFileName) {
        arret("Function_Scan: too long file name") ;
      }
      
      {
        char* c = String_FindAndSkipToken(line,"=") ;
        
        n_fn = Function_ReadInFile(function,c) ;
      }
      
      return(n_fn) ;
    }



    {
      arret("Function_Scan: keyword not known") ;
    }
    
  return(0) ;
}



double (Function_ComputeValue)(Function_t* fn,double t)
{
  if(fn) {
    int   nb_points = Function_GetNbOfPoints(fn) ;
    double* tm = Function_GetXValue(fn) ;
    double* ft = Function_GetFValue(fn) ;
    
    /*
      Cas t < t[0] 
    */
    if(t <= tm[0]) {
      return(ft[0]) ;
      
    /*
      Cas t > t[n-1] 
    */
    } else if(t >= tm[nb_points - 1]) {
      return(ft[nb_points - 1]) ;
      
    /*
      Cas t[i1] <= t <= t[i2]
      Calcul des deux points i1=Min(i) et i2=Max(i) du tableau fn
      correspondant au plus petit intervalle de temps [t1;t2]
      contenant t. Deux cas de figure se presentent:
      1) t2 > t1 et i2 = i1+1;
      2) t2 = t1 et i2 >= i1 avec i2-i1 maximum.
    */
    } else {
      int    i1 = 0, i2 = nb_points - 1 ;
      double t1, t2, f1, f2 ;
      int    i ;
      
      for(i = 0 ; i < nb_points ; i++) {
        if(t >= tm[i]) i1 = i ;
        if(t <= tm[i]) break  ;
      }
      
      for(i = i1 ; i < nb_points ; i++) if(t < tm[i]) break ;
      
      for( ; i >= 0 ; i--) {
        if(t <= tm[i]) i2 = i ;
        if(t >= tm[i]) break  ;
      }
      
      t1 = tm[i1] ;
      f1 = ft[i1] ;
      t2 = tm[i2] ;
      f2 = ft[i2] ;
      
      /* 1) Intervalle non nul */
      if(t2 > t1) {
        return (f1 + (f2 - f1)*(t - t1)/(t2 - t1)) ;
        
      /* 2) Intervalle nul */
      } else {
        /* un seul point */
        if(i1 == i2) {
          return (f1) ;
        /* plusieurs points */
        } else {
          arret("plusieurs pas de temps nuls (fonction)") ;
          return(0.) ;
        }
      }
    }
  }
  
  return(0.) ;
}




int (Function_ReadInFile)(Function_t* fn,char* line1)
/* Lecture des fonctions du temps dans le fichier "nom"
   retourne le nb de fonctions lues */
{
  int    n_points,n_fonctions ;
  char   line[Function_MaxLengthOfTextLine],*c ;
  FILE   *fict ;
  int long pos ;
  char nom[Function_MaxLengthOfFileName] ;
  
  String_Scan(line1," %s",nom) ;
      
  if(strlen(nom) > Function_MaxLengthOfFileName) {
    arret("Function_ReadInFile: too long file name") ;
  }

  fict = fopen(nom,"r") ;
  
  if(fict) {
    /* nb de fonctions */
    n_fonctions = 0 ;
    
    do {
      fgets(line,sizeof(line),fict) ;
      c = (char*) strtok(line," \n") ;
      /* } while(c != NULL && *c == '#' && !feof(fict)) ; */
    } while((c == NULL || *c == '#') && !feof(fict)) ;
    
    while((char*) strtok(NULL," \n") != NULL) n_fonctions++ ;

    /* nb de points */
    n_points = 1 ;
    
    while(!feof(fict) && fgets(line,sizeof(line),fict)) {
      c = (char*) strtok(line," \n") ;
      if(c != NULL && *c != '#') n_points++ ;
    }

    /* reservation de la memoire */
    {
      int i ;
    
      for(i = 0 ; i < n_fonctions ;i++) {
        Function_t* fct = Function_New(n_points) ;
        
        fn[i] = fct[0] ;
        free(fct) ;
      }
    }
  
    {
      /* positionnement dans le fichier */
      rewind(fict) ;
    
      do {
        pos = ftell(fict) ;
        fgets(line,sizeof(line),fict) ;
        c = (char*) strtok(line," \n") ;
        /* } while(c != NULL && *c == '#' && !feof(fict)) ; */
      } while((c == NULL || *c == '#') && !feof(fict)) ;
    
      fseek(fict,pos,SEEK_SET) ;
  
  
      /* Read time and function values */
      {
        int i ;
        
        for(i = 0 ; i < n_points ; i++) {
          double x ;
          int    j ;
    
          fscanf(fict,"%le",&x) ;
    
          for(j = 0 ; j < n_fonctions ; j++) {
            double* t = Function_GetXValue(fn + j) ;
            double* f = Function_GetFValue(fn + j) ;
            
            t[i] = x ;
            fscanf(fict,"%le",f + i) ;
          }
        }
      }
    }
    
  } else {
    int nbofcurves = 10 ;
    Curves_t* curves = Curves_Create(nbofcurves) ;
    Curve_t*  curve  = Curves_GetCurve(curves) ;
    
    Curves_ReadCurves(curves,line1) ;
    
    n_fonctions = Curves_GetNbOfCurves(curves) ;
    
    n_points = Curve_GetNbOfPoints(curve) ;
    
    //arret("Function_ReadInFile: unable to open the file") ;

    /* reservation de la memoire */
    {
      int i ;
    
      for(i = 0 ; i < n_fonctions ;i++) {
        Function_t* fct = Function_New(n_points) ;
        
        fn[i] = fct[0] ;
        free(fct) ;
      }
    }
  
    {
      double* t = Curve_CreateSamplingOfX(curve) ;
      int i ;
      
      /* Read time and function values */
      for(i = 0 ; i < n_points ; i++) {
        int    j ;
    
        for(j = 0 ; j < n_fonctions ; j++) {
          Curve_t* curve_j = curve + j ;
          double* y = Curve_GetYValue(curve_j) ;
          double* x = Function_GetXValue(fn + j) ;
          double* f = Function_GetFValue(fn + j) ;
          
          x[i] = t[i] ;
          f[i] = y[i] ;
        }
      }
      
      free(t) ;
    }
    
    Curves_Delete(curves) ;
    free(curves) ;
  }

  fclose(fict) ;

  return(n_fonctions) ;
}
