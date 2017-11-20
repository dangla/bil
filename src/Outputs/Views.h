#ifndef VIEWS_H
#define VIEWS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Views_s      ; typedef struct Views_s      Views_t ;


extern Views_t*   Views_Create(int) ;
extern void       Views_Delete(Views_t**) ;


#define Views_GetNbOfViews(views)      ((views)->nbofviews)
#define Views_GetView(views)           ((views)->view)


#define Views_MaxNbOfViews           (100)


#include "View.h"


struct Views_s {            /* Views */
  unsigned int nbofviews ;  /* nb of views */
  View_t *view ;            /* view */
} ;

#endif
