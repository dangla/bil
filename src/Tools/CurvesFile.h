#ifndef CURVESFILE_H
#define CURVESFILE_H



/* vacuous declarations and typedef names */

/* class-like structures */
struct CurvesFile_s       ; typedef struct CurvesFile_s       CurvesFile_t ;


extern CurvesFile_t*   (CurvesFile_Create)(void) ;
extern void            (CurvesFile_Delete)(void*) ;
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



#define CurvesFile_GetTextFile(CF)                         ((CF)->textfile)
#define CurvesFile_GetFilePositionStartingInputData(CF)    ((CF)->inputpos)
#define CurvesFile_GetTextLine(CF)                         ((CF)->line)
#define CurvesFile_GetNbOfCurves(CF)                       ((CF)->n_curves)
#define CurvesFile_GetNbOfPoints(CF)                       ((CF)->n_points)
#define CurvesFile_GetScaleType(CF)                        ((CF)->scale)
#define CurvesFile_GetCommandLine(CF)                      ((CF)->cmdline)
#define CurvesFile_GetCurrentPositionInTheCommandLine(CF)  ((CF)->pcmdline)
#define CurvesFile_GetCurves(CF)                           ((CF)->readcurves)




/*
 *  Function-like macros
 */
#define CurvesFile_GetFileName(CF) \
        TextFile_GetFileName(CurvesFile_GetTextFile(CF))
        
#define CurvesFile_GetFileStream(CF) \
        TextFile_GetFileStream(CurvesFile_GetTextFile(CF))

#define CurvesFile_GetFilePosition(CF) \
        TextFile_GetFilePosition(CurvesFile_GetTextFile(CF))

#define CurvesFile_OpenFile(CF,mode) \
        TextFile_OpenFile(CurvesFile_GetTextFile(CF),mode)

#define CurvesFile_CloseFile(CF) \
        TextFile_CloseFile(CurvesFile_GetTextFile(CF))

#define CurvesFile_Exists(CF) \
        TextFile_Exists(CurvesFile_GetTextFile(CF))

#define CurvesFile_DoesNotExist(CF) \
        (!CurvesFile_Exists(CF))

#define CurvesFile_StoreFilePosition(CF) \
        TextFile_StoreFilePosition(CurvesFile_GetTextFile(CF))

#define CurvesFile_MoveToStoredFilePosition(CF) \
        TextFile_MoveToStoredFilePosition(CurvesFile_GetTextFile(CF))

#define CurvesFile_ReadLineFromCurrentFilePosition(CF) \
        TextFile_ReadLineFromCurrentFilePosition(CurvesFile_GetTextFile(CF),CurvesFile_GetTextLine(CF),CurvesFile_MaxLengthOfTextLine)

#define CurvesFile_FileCopy(CF) \
        TextFile_FileCopy(CurvesFile_GetTextFile(CF))

#define CurvesFile_FileStreamCopy(CF) \
        TextFile_FileStreamCopy(CurvesFile_GetTextFile(CF))



#include "Buffer.h"
#include "Curves.h"

struct CurvesFile_s {         /* File of discretized curves */
  const char* cmdline ;       /* Command line used to build the curves */
  const char* pcmdline ;      /* Current position in the command line */
  TextFile_t* textfile ;      /* Text file */
  /* fpos_t* inputpos ; */          /* File stream position which starts the input */
  unsigned int n_curves ;     /* Nb of curves */
  unsigned int n_points ;     /* Nb of points */
  char   scale ;              /* Scale = n(ormal-scale) or l(og-scale) */
  Curves_t* readcurves ;      /* Already read curves */
  char* line ;                /* Pointer to text lines */
} ;


#endif
