#ifndef VECTORSTORAGEFORMAT_H
#define VECTORSTORAGEFORMAT_H



enum VectorStorageFormat_e {        /* format of the vector to be stored */
  VectorStorageFormat_Array,        /* Array of doubles */
  VectorStorageFormat_PETSCVEC,     /* PETSc Vec format */
  VectorStorageFormat_NULL
} ;


/* class-like structure "VectorStorageFormat_t" and attributes */

/* vacuous declarations and typedef names */
typedef enum VectorStorageFormat_e VectorStorageFormat_e ;
struct VectorStorageFormat_s   ; typedef struct VectorStorageFormat_s  VectorStorageFormat_t ;



#include "Options.h"

extern VectorStorageFormat_t* (VectorStorageFormat_Create)(Options_t*) ;
extern void                   (VectorStorageFormat_Delete)(void*) ;


/** The getters */
#define VectorStorageFormat_GetType(MSF)              ((MSF)->type)
#define VectorStorageFormat_GetOptions(MSF)           ((MSF)->options)


#include "Utils.h"

#define VectorStorageFormat_Is(MSF,KEY) \
        (VectorStorageFormat_GetType(MSF) == Utils_CAT(VectorStorageFormat_,KEY))

#define VectorStorageFormat_Type(KEY) \
        ((VectorStorageFormat_e) Utils_CAT(VectorStorageFormat_,KEY))



/* complete the structure types by using the typedef */

struct VectorStorageFormat_s {
  VectorStorageFormat_e type ;
  Options_t* options ;
} ;

#endif
