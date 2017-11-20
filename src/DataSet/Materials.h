#ifndef MATERIALS_H
#define MATERIALS_H


/* vacuous declarations and typedef names */

/* class-like structures */
struct Materials_s    ; typedef struct Materials_s    Materials_t ;


#include "DataFile.h"
#include "Geometry.h"


extern Materials_t* (Materials_Create)(DataFile_t*,Geometry_t*) ;


#define Materials_GetNbOfMaterials(MATS)  (MATS->n_mat)
#define Materials_GetMaterial(MATS)       (MATS->mat)
#define Materials_GetUsedModels(MATS)     (MATS->models)


#define Materials_GetNbOfUsedModels(MATS) \
        Models_GetNbOfModels(Materials_GetUsedModels(MATS))



#include "Models.h"
#include "Material.h"


struct Materials_s {          /* materials */
  unsigned int n_mat ;        /**< Nb of materials */
  Material_t* mat ;           /**< Material */
  Models_t* models ;          /**< Used models */
} ;

#endif
