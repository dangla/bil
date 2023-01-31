#ifndef MATERIAL_H
#define MATERIAL_H


/* vacuous declarations and typedef names */

/* class-like structures */
struct Material_s     ; typedef struct Material_s     Material_t ;


#include <stdlib.h>
#include "Model.h"
#include "DataFile.h"
#include "Geometry.h"


extern Material_t* (Material_New)             (void) ;
extern void        (Material_Delete)          (void*) ;
extern void        (Material_Scan)            (Material_t*,DataFile_t*,Geometry_t*) ;
extern int         (Material_ReadProperties)  (Material_t*,DataFile_t*) ;
extern void        (Material_ScanProperties)  (Material_t*,DataFile_t*,Model_ComputePropertyIndex_t*) ;
extern void        (Material_ScanProperties1) (Material_t*,FILE*,Model_ComputePropertyIndex_t*,int) ;
extern void        (Material_ScanProperties2) (Material_t*,FILE*,Model_ComputePropertyIndex_t*,int,int) ;



#define Material_MaxLengthOfKeyWord            (30)
#define Material_MaxLengthOfTextLine           (500)

#define Material_MaxNbOfCurves                 (20)     /* Max nb of curves per mat */
#define Material_MaxNbOfProperties             (200)    /* Max nb of scalar inputs */



#define Material_GetNbOfProperties(MAT)   ((MAT)->n)
#define Material_GetProperty(MAT)         ((MAT)->pr)
#define Material_GetCurves(MAT)           ((MAT)->curves)
#define Material_GetFields(MAT)           ((MAT)->fields)
#define Material_GetFunctions(MAT)        ((MAT)->functions)
#define Material_GetModel(MAT)            ((MAT)->model)
#define Material_GetMethod(MAT)           ((MAT)->method)
#define Material_GetCodeNameOfModel(MAT)  ((MAT)->codenameofmodel)
#define Material_GetGenericData(MAT)      ((MAT)->genericdata)
#define Material_GetUsedModels(MAT)       ((MAT)->models)
#define Material_GetModelIndex(MAT)       ((MAT)->modelindex)



/* Material properties */
#define Material_GetNbOfCurves(MAT) \
        Curves_GetNbOfCurves(Material_GetCurves(MAT))

#define Material_GetCurve(MAT) \
        Curves_GetCurve(Material_GetCurves(MAT))

#define Material_GetNbOfFields(MAT) \
        Fields_GetNbOfFields(Material_GetFields(MAT))

#define Material_GetField(MAT) \
        Fields_GetField(Material_GetFields(MAT))

#define Material_GetNbOfFunctions(MAT) \
        Functions_GetNbOfFunctions(Material_GetFunctions(MAT))

#define Material_GetFunction(MAT) \
        Functions_GetFunction(Material_GetFunctions(MAT))

#define Material_GetDimension(MAT) \
        Geometry_GetDimension(Model_GetGeometry(Material_GetModel(MAT)))

#define Material_FindCurve(MAT,S) \
        Curves_FindCurve(Material_GetCurves(MAT),S)
        
#define Material_GetPropertyValue(MAT,S) \
        (Material_GetProperty(MAT) + Model_GetComputePropertyIndex(Material_GetModel(MAT))(S))[0]

#define Material_SetPropertiesToZero(MAT,N) \
        do { \
          int Material_i ; \
          for(Material_i = 0 ; Material_i < N ; Material_i++) { \
            Material_GetProperty(MAT)[Material_i] = 0 ; \
          } \
        } while(0)

/*
** #define Material_ReadProperties(MAT,datafile) \
*          Model_GetReadMaterialProperties(Material_GetModel(MAT))(MAT,datafile)
*/

/* Equations/unknowns */
#define Material_GetNbOfEquations(MAT) \
        Model_GetNbOfEquations(Material_GetModel(MAT))

#define Material_GetNameOfEquation(MAT) \
        Model_GetNameOfEquation(Material_GetModel(MAT))

#define Material_GetNameOfUnknown(MAT) \
        Model_GetNameOfUnknown(Material_GetModel(MAT))
        
#define Material_CopyNameOfEquation(MAT,index,name) \
        (strcpy(Material_GetNameOfEquation(MAT)[index],name))

#define Material_CopyNameOfUnknown(MAT,index,name) \
        (strcpy(Material_GetNameOfUnknown(MAT)[index],name))

#define Material_GetObjectiveValue(MAT) \
        Model_GetObjectiveValue(Material_GetModel(MAT))

#define Material_GetSequentialIndexOfUnknown(MAT) \
        Model_GetSequentialIndexOfUnknown(Material_GetModel(MAT))
        
        
#include "TypeId.h"

/* GenericData */
#define Material_AppendGenericData(MAT,GD) \
        do { \
          if(Material_GetGenericData(MAT)) { \
            GenericData_Append(Material_GetGenericData(MAT),GD) ; \
          } else { \
            Material_GetGenericData(MAT) = GD ; \
          } \
        } while(0)
        
#define Material_FindGenericData(MAT,...) \
        GenericData_Find(Material_GetGenericData(MAT),__VA_ARGS__)
        
#define Material_FindData(MAT,...) \
        GenericData_FindData(Material_GetGenericData(MAT),__VA_ARGS__)
        
#define Material_FindNbOfData(MAT,...) \
        GenericData_FindNbOfData(Material_GetGenericData(MAT),__VA_ARGS__)

#define Material_AppendData(MAT,...) \
        Material_AppendGenericData(MAT,GenericData_Create(__VA_ARGS__))





#include "Fields.h"
#include "Functions.h"
#include "Curves.h"
#include "GenericData.h"
#include "Models.h"

struct Material_s {           /* material */
  char*   codenameofmodel ;   /**< Code name of the model */
  char*   method ;            /**< Characterize a method */
  int     n ;                 /**< Nb of properties */
  double* pr ;                /**< The properties */
  GenericData_t* genericdata ;
  Curves_t* curves ;          /**< Curves */
  Fields_t* fields ;          /**< Fields */
  Functions_t* functions ;    /**< Time functions */
  Models_t* models ;          /**< Used models */
  Model_t* model ;            /**< Model */
  int modelindex ;            /**< Model index */
  
  /* for compatibility with former version (should be eliminated) */
  unsigned short int neq ;    /**< nombre d'equations du modele */
  char**   eqn ;              /**< nom des equations */
  char**   inc ;              /**< nom des inconnues */
  int      nc ;               /**< nb of curves */
  Curve_t* cb ;               /**< curves */
  
#ifdef NOTDEFINED             /* NON UTILISE POUR LE MOMENT */
  int      nfd ;              /**< nombre de donnees formelles */
  int      fdl ;              /**< longueur en caractere des donnees formelles */
  char**   fd ;               /**< les donnees formelles en mode char */
#endif
} ;


/* Old notations which should be eliminated */
#define mate_t                 Material_t
#define dmat                   Material_ScanProperties1
#define lit_mate               Material_ScanProperties2

#endif
