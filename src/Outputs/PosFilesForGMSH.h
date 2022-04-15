#ifndef POSFILESFORGMSH_H
#define POSFILESFORGMSH_H


/* class-like structures "OutputFiles_t" */

/* vacuous declarations and typedef names */
struct PosFilesForGMSH_s ; typedef struct PosFilesForGMSH_s PosFilesForGMSH_t ;



/* Declaration of Macros, Methods and Structures */

#include "DataSet.h"

extern PosFilesForGMSH_t*  (PosFilesForGMSH_Create)(DataSet_t*) ;
extern void                (PosFilesForGMSH_Delete)(void*) ;
extern void                (PosFilesForGMSH_ParsedFileFormat)(PosFilesForGMSH_t*) ;
extern void                (PosFilesForGMSH_ASCIIFileFormat)(PosFilesForGMSH_t*) ;



#define PosFilesForGMSH_GetDataSet(PFG)                ((PFG)->dataset)
#define PosFilesForGMSH_GetOutputFiles(PFG)            ((PFG)->outputfiles)
#define PosFilesForGMSH_GetBilVersion(PFG)             ((PFG)->version)



#include "OutputFiles.h"

/* complete the structure types by using the typedef */
struct PosFilesForGMSH_s {
  double version ;
  DataSet_t* dataset ;
  OutputFiles_t* outputfiles ;  /* The output files */
} ;


#endif
