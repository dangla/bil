#ifndef OUTPUTFILES_H
#define OUTPUTFILES_H


/* class-like structures "OutputFiles_t" */

/* vacuous declarations and typedef names */
struct OutputFiles_s ; typedef struct OutputFiles_s OutputFiles_t ;



/* Declaration of Macros, Methods and Structures */

#include "DataSet.h"

extern OutputFiles_t*   (OutputFiles_Create)(char*,int,int) ;
extern void    (OutputFiles_Delete)(OutputFiles_t**) ;
extern void    (OutputFiles_PostProcessForGmshParsedFileFormat)(OutputFiles_t*,DataSet_t*) ;
extern void    (OutputFiles_PostProcessForGmshASCIIFileFormat)(OutputFiles_t*,DataSet_t*) ;
extern void    (OutputFiles_BackupSolutionAtTime)(OutputFiles_t*,DataSet_t*,double,int) ;
extern void    (OutputFiles_BackupSolutionAtPoint)(OutputFiles_t*,DataSet_t*,double,double) ;


/* Function-like macros */
#define OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) \
        TextFile_ReadLineFromCurrentFilePosition(textfile,OutputFiles_GetTextLine(outputfiles),OutputFiles_MaxLengthOfTextLine)


#include "Views.h"
#include "TextFile.h"

#define OutputFiles_MaxNbOfViews           (Views_MaxNbOfViews)
#define OutputFiles_MaxLengthOfViewName    (View_MaxLengthOfViewName)
#define OutputFiles_MaxLengthOfFileName    (TextFile_MaxLengthOfFileName)
#define OutputFiles_MaxLengthOfTextLine    ((OutputFiles_RecordNumberLength)*(OutputFiles_MaxNbOfViews*9 + 3))

#define OutputFiles_RecordFieldWidth       (14)
#define OutputFiles_RecordPrecision        (6)
#define OutputFiles_RecordNumberFormat     "% -14.6e "
#define OutputFiles_RecordNumberLength     (14)


#define OutputFiles_GetDataFileName(outputfiles)            ((outputfiles)->filename)
#define OutputFiles_GetNbOfDateFiles(outputfiles)           ((outputfiles)->n_dates)
#define OutputFiles_GetNbOfPointFiles(outputfiles)          ((outputfiles)->n_points)
//#define OutputFiles_GetDateFile1(outputfiles)                ((outputfiles)->datefile)
//#define OutputFiles_GetPointFile1(outputfiles)               ((outputfiles)->pointfile)
#define OutputFiles_GetTextLine(outputfiles)                ((outputfiles)->line)
#define OutputFiles_GetDateOutputFile(outputfiles)          ((outputfiles)->dateoutputfile)
#define OutputFiles_GetPointOutputFile(outputfiles)         ((outputfiles)->pointoutputfile)

//#define OutputFiles_GetDateFile(outputfiles)                ((outputfiles)->datefile)
//#define OutputFiles_GetDateFile(outputfiles)                (OutputFile_GetTextFile(OutputFiles_GetDateOutputFile(outputfiles)))
//#define OutputFiles_GetPointFile(outputfiles)               ((outputfiles)->pointfile))
//#define OutputFiles_GetPointFile(outputfiles)                (OutputFile_GetTextFile(OutputFiles_GetPointOutputFile(outputfiles)))



#include "OutputFile.h"

/* complete the structure types by using the typedef */
struct OutputFiles_s {            /* Output files */
  char*  filename ;               /* name of the data file */
  int    n_dates ;                /* Nb of dates */
  int    n_points ;               /* Nb of points */
  /* char**  datefilename ;         *//* names of the date output files */
  /* char**  pointfilename ;        *//* names of the point output files */
  /* FILE**  datefilestream ;       *//* results at specified times */
  /* FILE**  pointfilestream ;      *//* results at specified points */
  //TextFile_t* datefile ;          /* The date files */
  //TextFile_t* pointfile ;         /* The point files */
  OutputFile_t* dateoutputfile ;  /* The date output files */
  OutputFile_t* pointoutputfile ; /* The point output files */
  char* line ;                    /* Pointer to text lines */
} ;


#endif
