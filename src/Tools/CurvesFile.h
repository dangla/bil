#ifndef CURVESFILE_H
#define CURVESFILE_H



/* vacuous declarations and typedef names */

/* class-like structures */
struct CurvesFile_s       ; typedef struct CurvesFile_s       CurvesFile_t ;


extern CurvesFile_t*   (CurvesFile_Create)(void) ;
extern void            (CurvesFile_Delete)(CurvesFile_t**) ;
/* extern void            (CurvesFile_MoveToFilePositionStartingInputData)(CurvesFile_t*) ; */
extern int             (CurvesFile_Initialize)(CurvesFile_t*,const char *) ;
extern int             (CurvesFile_WriteCurves)(CurvesFile_t*) ;



#include "TextFile.h"

#define CurvesFile_MaxLengthOfTextLine      (500)
#define CurvesFile_MaxNbOfCurves            (20)
#define CurvesFile_SizeOfBuffer             (CurvesFile_MaxLengthOfTextLine*sizeof(char))
#define CurvesFile_MaxLengthOfKeyWord       (30)
#define CurvesFile_MaxLengthOfFileName      (TextFile_MaxLengthOfFileName)
#define CurvesFile_MaxLengthOfCurveName     (30)

#define CurvesFile_FieldDelimiters           " ,\t"



#define CurvesFile_GetTextFile(curvesfile)          ((curvesfile)->textfile)
#define CurvesFile_GetFilePositionStartingInputData(curvesfile)      ((curvesfile)->inputpos)
#define CurvesFile_GetTextLine(curvesfile)          ((curvesfile)->line)
#define CurvesFile_GetNbOfCurves(curvesfile)        ((curvesfile)->n_curves)
#define CurvesFile_GetNbOfPoints(curvesfile)        ((curvesfile)->n_points)
#define CurvesFile_GetScaleType(curvesfile)         ((curvesfile)->scale)
#define CurvesFile_GetBuffer(curvesfile)            ((curvesfile)->buffer)
#define CurvesFile_GetCommandLine(curvesfile)       ((curvesfile)->cmdline)
#define CurvesFile_GetCurrentPositionInTheCommandLine(curvesfile)       ((curvesfile)->pcmdline)
#define CurvesFile_GetCurves(curvesfile)            ((curvesfile)->readcurves)



#define CurvesFile_GetFileName(curvesfile)          (TextFile_GetFileName(CurvesFile_GetTextFile(curvesfile)))
#define CurvesFile_GetFileStream(curvesfile)        (TextFile_GetFileStream(CurvesFile_GetTextFile(curvesfile)))
#define CurvesFile_GetFilePosition(curvesfile)      (TextFile_GetFilePosition(CurvesFile_GetTextFile(curvesfile)))



/*
 *  Function-like macros
 */
#define CurvesFile_OpenFile(curvesfile,mode)   (TextFile_OpenFile(CurvesFile_GetTextFile(curvesfile),mode))

#define CurvesFile_CloseFile(curvesfile)       (TextFile_CloseFile(CurvesFile_GetTextFile(curvesfile)))

#define CurvesFile_Exists(curvesfile)          (TextFile_Exists(CurvesFile_GetTextFile(curvesfile)))

#define CurvesFile_DoesNotExist(curvesfile)    (!CurvesFile_Exists(curvesfile))

#define CurvesFile_StoreFilePosition(curvesfile) (TextFile_StoreFilePosition(CurvesFile_GetTextFile(curvesfile)))

#define CurvesFile_MoveToStoredFilePosition(curvesfile)  (TextFile_MoveToStoredFilePosition(CurvesFile_GetTextFile(curvesfile)))

#define CurvesFile_ReadLineFromCurrentFilePosition(curvesfile)  (TextFile_ReadLineFromCurrentFilePosition(CurvesFile_GetTextFile(curvesfile),CurvesFile_GetTextLine(curvesfile),CurvesFile_MaxLengthOfTextLine))

#define CurvesFile_AllocateInBuffer(curvesfile,sz)  (Buffer_Allocate(CurvesFile_GetBuffer(curvesfile),(sz)))

#define CurvesFile_FreeBuffer(curvesfile)           (Buffer_Free(CurvesFile_GetBuffer(curvesfile)))

#define CurvesFile_FileCopy(curvesfile)        (TextFile_FileCopy(CurvesFile_GetTextFile(curvesfile)))

#define CurvesFile_FileStreamCopy(curvesfile)  (TextFile_FileStreamCopy(CurvesFile_GetTextFile(curvesfile)))



#include "Buffer.h"
#include "Curves.h"

struct CurvesFile_s {         /* File of discretized curves */
  const char* cmdline ;       /* Command line used to build the curves */
  const char* pcmdline ;      /* Current position in the command line */
  TextFile_t* textfile ;      /* Text file */
  /* fpos_t *inputpos ; */          /* File stream position which starts the input */
  unsigned int n_curves ;     /* Nb of curves */
  unsigned int n_points ;     /* Nb of points */
  char   scale ;              /* Scale = n(ormal-scale) or l(og-scale) */
  Buffer_t *buffer ;          /* Buffer */
  Curves_t *readcurves ;      /* Already read curves */
  char* line ;                /* Pointer to text lines */
} ;


#endif
