#ifndef COMMANDLINE_H
#define COMMANDLINE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct CommandLine_s  ; typedef struct CommandLine_s  CommandLine_t ;


extern CommandLine_t*    (CommandLine_Create)(int,char**) ;
extern void              (CommandLine_Delete)(void*) ;


#define CommandLine_GetNbOfArg(cmd)            ((cmd)->argc)
#define CommandLine_GetArg(cmd)                ((cmd)->argv)



struct CommandLine_s {        /* Command line */
  int    argc ;               /* Nb of command line arguments */
  char** argv ;               /* Command line arguments */
} ;


#endif
