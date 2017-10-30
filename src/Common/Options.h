#ifndef OPTIONS_H
#define OPTIONS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Options_s      ; typedef struct Options_s      Options_t ;


extern Options_t*  Options_Create(void) ;


#define Options_MaxLengthOfKeyWord               (30)

#define Options_GetPrintData(options)              ((options)->debug)
#define Options_GetResolutionMethod(options)       ((options)->method)
#define Options_GetPrintLevel(options)             ((options)->level)
#define Options_GetModule(options)                 ((options)->module)
#define Options_GetGraphMethod(options)            ((options)->graph)
#define Options_GetInputFileName(options)          ((options)->inputfile)
#define Options_GetElementOrderingMethod(options)  ((options)->eordering)
#define Options_GetNodalOrderingMethod(options)    ((options)->nordering)
#define Options_GetPostProcessingMethod(options)   ((options)->postprocess)

#define Options_GetDebug(options)                (Options_GetPrintData(options))


struct Options_s {            /* options */
  char   *debug ;             /* data to be printed */
  char   *method ;            /* resolution method */
  char   *level ;             /* print level */
  char   *module ;            /* module */
  char   *graph ;             /* Graph method */
  char   *inputfile ;         /* Input data file */
  char   *eordering ;         /* Element ordering method */
  char   *nordering ;         /* Nodal ordering method */
  char   *postprocess ;       /* Post-processing method */
} ;


#endif
