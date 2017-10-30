#ifndef COMMANDLINE_H
#define COMMANDLINE_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct CommandLine_s  ; typedef struct CommandLine_s  CommandLine_t ;


extern void              (CommandLine_Initialize)(int,char**) ;
extern CommandLine_t*    (CommandLine_GetInstance)(void) ;


#define CommandLine_GetNbOfArg(cmd)            ((cmd)->argc)
#define CommandLine_GetArg(cmd)                ((cmd)->argv)



#define CommandLine_NbOfArg \
        CommandLine_GetNbOfArg(CommandLine_GetInstance())
        
#define CommandLine_Arg \
        CommandLine_GetArg(CommandLine_GetInstance())


struct CommandLine_s {        /* Command line */
  int    argc ;               /* Nb of command line arguments */
  char** argv ;               /* Command line arguments */
} ;


#endif
