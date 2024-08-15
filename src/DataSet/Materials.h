#ifndef MATERIALS_H
#define MATERIALS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structures */
struct Materials_s    ; typedef struct Materials_s    Materials_t ;


#include "DataFile.h"
#include "Geometry.h"
#include "Fields.h"
#include "Functions.h"
#include "Models.h"


extern Materials_t* (Materials_Create)(DataFile_t*,Geometry_t*,Fields_t*,Functions_t*,Models_t*) ;
extern Materials_t* (Materials_New)   (const int,Models_t*) ;
extern void         (Materials_Delete)(void*) ;


#define Materials_GetNbOfMaterials(MATS)  ((MATS)->n_mat)
#define Materials_GetMaterial(MATS)       ((MATS)->mat)
#define Materials_GetUsedModels(MATS)     ((MATS)->models)


#define Materials_GetNbOfUsedModels(MATS) \
        Models_GetNbOfModels(Materials_GetUsedModels(MATS))
        
#define Materials_GetUsedModel(MATS) \
        Models_GetModel(Materials_GetUsedModels(MATS))



#include "Material.h"


struct Materials_s {          /* materials */
  unsigned int n_mat ;        /**< Nb of materials */
  Material_t* mat ;           /**< Material */
  Models_t* models ;          /**< Used models */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
