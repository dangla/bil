#ifndef EXCEPTION_H
#define EXCEPTION_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Exception_s     ; typedef struct Exception_s     Exception_t ;


extern Exception_t*    (Exception_GetInstance)(void) ;


#define Exception_GetEnvironment(exception)        ((exception)->env)
#define Exception_GetExceptionType(exception)      ((exception)->type)


#include <signal.h>
#include <setjmp.h>

#define Exception_Environment \
        Exception_GetEnvironment(Exception_GetInstance())
      
#define Exception_ExceptionType \
        Exception_GetExceptionType(Exception_GetInstance())

#define Exception_SaveEnvironment \
        (Exception_ExceptionType = setjmp(Exception_Environment))

#define Exception_RestoreEnvironment(val) \
        (longjmp(Exception_Environment,val))

#define Exception_Interrupt \
        (raise(SIGINT))

#define Exception_IsCaught \
        (Exception_ExceptionType)

#define Exception_IsNotCaught \
        (!Exception_IsCaught)


struct Exception_s {        /* Exception handler */
  jmp_buf env ;             /* Environment */
  int type ;                /* Signal sent by longjmp */
} ;

#endif
