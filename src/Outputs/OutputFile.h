#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H


/* class-like structures "OutputFiles_t" */

/* vacuous declarations and typedef names */
struct OutputFile_s  ; typedef struct OutputFile_s  OutputFile_t ;



/* Declaration of Macros, Methods and Structures */

extern char    (OutputFile_TypeOfCurrentFile) ;


#define OutputFile_GetTextFile(outputfile)            ((outputfile)->textfile)


#define OutputFile_IsPointType      (OutputFile_TypeOfCurrentFile == 'p')
#define OutputFile_IsTimeType       (OutputFile_TypeOfCurrentFile == 't')




struct OutputFile_s {             /* Output file */
  TextFile_t* textfile ;          /* The file */
} ;

#endif
