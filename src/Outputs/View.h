#ifndef VIEW_H
#define VIEW_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct View_s       ; typedef struct View_s       View_t ;


extern View_t*    View_Create(int) ;
extern void       View_Delete(void*) ;


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
