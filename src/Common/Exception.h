#ifndef EXCEPTION_H
#define EXCEPTION_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct Exception_s     ; typedef struct Exception_s     Exception_t ;


extern Exception_t*    (Exception_GetInstance)(void) ;
extern void            (Exception_Delete)(void*) ;


#define Exception_GetEnvironment(exception)        ((exception)->env)
#define Exception_GetExceptionType(exception)      ((exception)->type)
#define Exception_GetDelete(exception)             ((exception)->Delete)


#include <signal.h>
#include <setjmp.h>
#include "Message.h"

#define Exception_Environment \
        Exception_GetEnvironment(Exception_GetInstance())
      
#define Exception_ExceptionType \
        Exception_GetExceptionType(Exception_GetInstance())

#define Exception_SaveEnvironment \
        (Exception_ExceptionType = setjmp(Exception_Environment))

#define Exception_RestoreEnvironment(val) \
        (longjmp(Exception_Environment,val))

#define Exception_Interrupt \
        do { \
          if(raise(SIGINT)) { \
            Message_FatalError("Error raising the signal") ; \
          } \
        } while (0)

#define Exception_IsCaught \
        (Exception_ExceptionType)

#define Exception_IsNotCaught \
        (!Exception_IsCaught)
        

/* Orders given through the exception mechanism */

/* Order to backup */
#define Exception_BackupAndTerminate \
        Exception_RestoreEnvironment(1)
        
#define Exception_OrderToBackupAndTerminate \
        (Exception_ExceptionType == 1)
        

/* Order to reduce the time step */
#define Exception_ReiterateWithSmallerTimeStep \
        Exception_RestoreEnvironment(2)

#define Exception_OrderToReiterateWithSmallerTimeStep \
        (Exception_ExceptionType == 2)
        

/* Order to initialize the time step */
#define Exception_ReiterateWithInitialTimeStep \
        Exception_RestoreEnvironment(3)
        
#define Exception_OrderToReiterateWithInitialTimeStep \
        (Exception_ExceptionType == 3)


#include "GenericObject.h"

struct Exception_s {        /* Exception handler */
  jmp_buf env ;             /* Environment */
  int type ;                /* Signal read by longjmp (1,2,3) and return by setjmp after a long jump */
  GenericObject_Delete_t* Delete ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
