#ifndef OUTPUTFILES_H
#define OUTPUTFILES_H


/* class-like structures "OutputFiles_t" */

/* vacuous declarations and typedef names */
struct OutputFiles_s ; typedef struct OutputFiles_s OutputFiles_t ;



/* Declaration of Macros, Methods and Structures */

#include "DataSet.h"

extern OutputFiles_t*   (OutputFiles_Create)(char*,int,int) ;
extern void    (OutputFiles_Delete)(void*) ;
extern void    (OutputFiles_PostProcessForGmshParsedFileFormat)(OutputFiles_t*,DataSet_t*) ;
extern void    (OutputFiles_PostProcessForGmshASCIIFileFormat)(OutputFiles_t*,DataSet_t*) ;
extern void    (OutputFiles_BackupSolutionAtTime_)(OutputFiles_t*,DataSet_t*,double,int) ;
extern void    (OutputFiles_BackupSolutionAtPoint_)(OutputFiles_t*,DataSet_t*,double,double) ;


/* Function-like macros */
#define OutputFiles_ReadLineFromCurrentFilePosition(OFS,textfile) \
        TextFile_ReadLineFromCurrentFilePosition(textfile,OutputFiles_GetTextLine(OFS),OutputFiles_MaxLengthOfTextLine)
        

#define OutputFiles_BackupSolutionAtTime(OFS,...) \
        if(OFS) OutputFiles_BackupSolutionAtTime_(OFS,__VA_ARGS__)


#define OutputFiles_BackupSolutionAtPoint(OFS,...) \
        if(OFS) OutputFiles_BackupSolutionAtPoint_(OFS,__VA_ARGS__)


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


#define OutputFiles_GetDataFileName(OFS)            ((OFS)->filename)
#define OutputFiles_GetNbOfDateFiles(OFS)           ((OFS)->n_dates)
#define OutputFiles_GetNbOfPointFiles(OFS)          ((OFS)->n_points)
//#define OutputFiles_GetDateFile1(OFS)                ((OFS)->datefile)
//#define OutputFiles_GetPointFile1(OFS)               ((OFS)->pointfile)
#define OutputFiles_GetTextLine(OFS)                ((OFS)->line)
#define OutputFiles_GetDateOutputFile(OFS)          ((OFS)->dateoutputfile)
#define OutputFiles_GetPointOutputFile(OFS)         ((OFS)->pointoutputfile)
#define OutputFiles_GetResults(OFS)                 ((OFS)->results)

//#define OutputFiles_GetDateFile(OFS)                ((OFS)->datefile)
//#define OutputFiles_GetDateFile(OFS)                (OutputFile_GetTextFile(OutputFiles_GetDateOutputFile(OFS)))
//#define OutputFiles_GetPointFile(OFS)               ((OFS)->pointfile))
//#define OutputFiles_GetPointFile(OFS)                (OutputFile_GetTextFile(OutputFiles_GetPointOutputFile(OFS)))





#include "OutputFile.h"
#include "Results.h"

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
  Results_t*    results ;         /* Allocated space for the results */
  char* line ;                    /* Pointer to text lines */
} ;


#endif
