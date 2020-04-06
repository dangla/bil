#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "DataFile.h"
#include "String.h"
#include "Mry.h"
#include "Functions.h"


static int    lit_fonction(Function_t*,char*) ;



Functions_t* (Functions_New)(const int n_fncts)
{
  Functions_t* functions   = (Functions_t*) Mry_New(Functions_t) ;
  
  
  Functions_GetNbOfFunctions(functions) = n_fncts ;
  Functions_GetFunction(functions) = NULL ;

  if(n_fncts > 0) {
    Function_t* function = (Function_t*) Mry_New(Function_t[n_fncts]) ;
    
    Functions_GetFunction(functions) = function ;
  }
  
  return(functions) ;
}



Functions_t* (Functions_Create)(DataFile_t* datafile)
{
  char* filecontent = DataFile_GetFileContent(datafile) ;
  char* c  = String_FindToken(filecontent,"FONC,FUNC,Functions",",") ;
  int n_fncts = (c = String_SkipLine(c)) ? atoi(c) : 0 ;
  Functions_t* functions = Functions_New(n_fncts) ;
  
  
  Message_Direct("Enter in %s","Functions") ;
  Message_Direct("\n") ;
  
  
  if(n_fncts <= 0) {
    return(functions) ;
  }


  {
    int i_fn ;
    
    c = String_SkipLine(c) ;
      
    DataFile_SetCurrentPositionInFileContent(datafile,c) ;
    
    for(i_fn = 0 ; i_fn < n_fncts ; i_fn++) {
      Function_t* function = Functions_GetFunction(functions) + i_fn ;
  
  
      Message_Direct("Enter in %s %d","Function",i_fn + 1) ;
      Message_Direct("\n") ;
      
      i_fn += Function_Scan(function,datafile) - 1 ;
      
      if(i_fn >= n_fncts) {
        arret("Functions_Create: too many functions") ;
      }
    }
  }
  
  return(functions) ;
}



/* Intern Functions */

int lit_fonction(Function_t* fn,char* nom)
/* Lecture des fonctions du temps dans le fichier "nom"
   retourne le nb de fonctions lues */
{
  int    n_points,n_fonctions ;
  char   line[Function_MaxLengthOfTextLine],*c ;
  FILE   *fict ;
  int    i ;
  int long pos ;

  fict = fopen(nom,"r") ;
  if(!fict) arret("lit_fonction(1) : immpossible d ouvrir le fichier") ;

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
    double* t = (double*) malloc(n_points*sizeof(double)) ;
    double* f = (double*) malloc(n_fonctions*n_points*sizeof(double)) ;
    if(!t) arret("lit_fonction(2) : impossible d\'allouer la memoire") ;
    if(!f) arret("lit_fonction(3) : impossible d\'allouer la memoire") ;
    
    for(i = 0 ; i < n_fonctions ;i++) {
      Function_GetNbOfPoints(fn + i) = n_points ;
      Function_GetXValue(fn + i) = t ; /* meme temps pour les fonctions */
      Function_GetFValue(fn + i) = f + i*n_points ;
    }
  }

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
  for(i = 0 ; i < n_points ; i++) {
    double* x = Function_GetXValue(fn) ;
    int    j ;
    
    fscanf(fict,"%le",x + i) ;
    
    for(j = 0 ; j < n_fonctions ; j++) {
      double* f = Function_GetFValue(fn + j) ;
      
      fscanf(fict,"%le",f + i) ;
    }
  }

  fclose(fict) ;

  return(n_fonctions) ;
}
