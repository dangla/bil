#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "Options.h"
#include "Mesh.h"
#include "Message.h"
#include "BilLib.h"
#include "SuperLUformat.h"



/* Extern functions */


#ifdef SLU_DIR

SuperLUformat_t* SuperLUformat_Create(Mesh_t* mesh)
/* Create a matrix in SuperLUformat format */
{
  int n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
  NCformat_t*    asluNC = NCformat_Create(mesh) ;
  SuperLUformat_t* aslu = (SuperLUformat_t*) malloc(sizeof(SuperLUformat_t)) ;

  assert(aslu) ;
  
  SuperLUformat_GetStorageType(aslu) = SLU_NC ;
  SuperLUformat_GetDataType(aslu)    = SLU_D ;
  SuperLUformat_GetMatrixType(aslu)  = SLU_GE ;
  SuperLUformat_GetNbOfColumns(aslu) = n_col ;
  SuperLUformat_GetNbOfRows(aslu)    = n_col ;
  SuperLUformat_GetStorage(aslu)     = (void*) asluNC ;
  
  return(aslu) ;
}



void SuperLUformat_Delete(void* self)
{
  SuperLUformat_t** aslu = (SuperLUformat_t**) self ;
  NCformat_t* asluNC = (NCformat_t*) SuperLUformat_GetStorage(*aslu) ;
  
  NCformat_Delete(&asluNC) ;
  
  free(*aslu) ;
  *aslu = NULL ;
}

#endif
