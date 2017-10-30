#ifndef VIEWS_H
#define VIEWS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Views_s      ; typedef struct Views_s      Views_t ;
/*     1. Views_t attributes */
struct View_s       ; typedef struct View_s       View_t ;


/* 1. Views_t */
extern Views_t*   Views_Create(int) ;
extern void       Views_Delete(Views_t**) ;


#define Views_GetNbOfViews(views)      ((views)->nbofviews)
#define Views_GetView(views)           ((views)->view)


#define Views_MaxNbOfViews           (100)


struct Views_s {            /* Views */
  unsigned int nbofviews ;  /* nb of views */
  View_t *view ;            /* view */
} ;



/* 2. View_t */
extern View_t*    View_Create(int) ;
extern void       View_Delete(View_t**) ;


#define View_MaxLengthOfViewName    (50)      /* Max length of view name */


#define View_GetNbOfComponents(view)  ((view)->n)
#define View_GetNameOfView(view)      ((view)->name)
#define View_GetGlobalIndex(view)     ((view)->index)



struct View_s {               /* View (scalar, vector, tensor) */
  short int n ;               /* Nb of components (1,3,9) */
  char*   name ;              /* Name of the view */
  int     index ;             /* Index of the view in the global set of views */
} ;


#endif
