#ifndef NODE_H
#define NODE_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Node_s         ; typedef struct Node_s         Node_t ;
typedef Node_t*   Node_tt ;


#include "Utils.h"


//extern Node_t*  (Node_Create)(const int) ;
//extern void     (Node_Delete)(void*) ;
extern int      (Node_FindUnknownPositionIndex)(const Node_t*,const char*) ;
extern int      (Node_FindEquationPositionIndex)(const Node_t*,const char*) ;
extern void     (Node_EliminateMatrixColumnIndexForOverlappingNodes)(const Node_t*,const char*) ;
extern void     (Node_EliminateMatrixRowIndexForOverlappingNodes)(const Node_t*,const char*) ;
extern void     (Node_UpdateMatrixColumnIndexForOverlappingNodes)(const Node_t*,const char*) ;
extern void     (Node_UpdateMatrixRowIndexForOverlappingNodes)(const Node_t*,const char*) ;
extern void     (Node_MakeUnknownContinuousAtOverlappingNodes)(const Node_t*,const char*) ;
extern void     (Node_MakeEquationContinuousAtOverlappingNodes)(const Node_t*,const char*) ;
extern Node_t*  (Node_OverlappingNodes3)(const Node_t*,int*,Node_t*) ;


#define Node_OverlappingNodes(...) \
        Utils_CAT_NARG(Node_OverlappingNodes,__VA_ARGS__)(__VA_ARGS__)

#define Node_OverlappingNodes2(...) \
        Node_OverlappingNodes3(__VA_ARGS__,NULL)




#include "Model.h"

/* Some constants */

#define Node_MaxNbOfEquations \
        Model_MaxNbOfEquations

#define Node_SizeOfBuffer \
        (100*sizeof(Node_t))


/* Accessors */
#define Node_GetNodeIndex(NOD)            ((NOD)->index)
#define Node_GetCoordinate(NOD)           ((NOD)->x)
#define Node_GetNbOfEquations(NOD)        ((NOD)->neq)
#define Node_GetNbOfUnknowns(NOD)         ((NOD)->nin)
#define Node_GetNameOfEquation(NOD)       ((NOD)->eqn)
#define Node_GetNameOfUnknown(NOD)        ((NOD)->inc)
#define Node_GetSequentialIndexOfUnknown(NOD)        ((NOD)->sequentialindex)
#define Node_GetObValIndex(NOD)           ((NOD)->i_obj)
#define Node_GetMatrixColumnIndex(NOD)    ((NOD)->colindex)
#define Node_GetMatrixRowIndex(NOD)       ((NOD)->rowindex)
#define Node_GetNodeSol(NOD)              ((NOD)->sol)
#define Node_GetNbOfElements(NOD)         ((NOD)->nel)
#define Node_GetPointerToElement(NOD)     ((NOD)->element)
#define Node_GetBuffers(NOD)              ((NOD)->buffers)



/* Buffer */
#define Node_GetBuffer(NOD) \
        Buffers_GetBufferOfCurrentThread(Node_GetBuffers(NOD))


/* Access to elements
 * ------------------ */
#define Node_GetElement(NOD,i) \
        (Node_GetPointerToElement(NOD)[i])



/* Access to NodeSol
 * ----------------- */
#define Node_GetCurrentNodeSol(NOD) \
        Node_GetNodeSol(NOD)
        
#define Node_GetPreviousNodeSol(NOD) \
        NodeSol_GetPreviousNodeSol(Node_GetNodeSol(NOD))
        
#define Node_GetNextNodeSol(NOD) \
        NodeSol_GetNextNodeSol(Node_GetNodeSol(NOD))

#define Node_GetNodeSolInDistantPast(NOD,i) \
        NodeSol_GetNodeSolInDistantPast(Node_GetNodeSol(NOD),i)

#define Node_GetNodeSolInDistantFuture(NOD,i) \
        NodeSol_GetNodeSolInDistantFuture(Node_GetNodeSol(NOD),i)



/* Access to unknowns
 * ------------------ */
#define Node_GetCurrentUnknown(NOD) \
        NodeSol_GetUnknown(Node_GetNodeSol(NOD))

#define Node_GetPreviousUnknown(NOD) \
        NodeSol_GetUnknown(Node_GetPreviousNodeSol(NOD))

#define Node_GetNextUnknown(NOD) \
        NodeSol_GetUnknown(Node_GetNextNodeSol(NOD))

#define Node_GetUnknownInDistantPast(NOD,i) \
        NodeSol_GetUnknown(Node_GetNodeSolInDistantPast(NOD,i))

#define Node_GetUnknownInDistantFuture(NOD,i) \
        NodeSol_GetUnknown(Node_GetNodeSolInDistantFuture(NOD,i))



/* Operations on buffer
 * -------------------- */
#define Node_AllocateInBuffer(NOD,sz) \
        Buffer_Allocate(Node_GetBuffer(NOD),(sz))
        
#define Node_FreeBuffer(NOD)  \
        Buffer_Free(Node_GetBuffer(NOD))
        
#define Node_FreeBufferFrom(NOD,p) \
        Buffer_FreeFrom(Node_GetBuffer(NOD),(char*) (p))
        
        
        
/* Operations on matrix row/column indexes
 * for unknown at specific nodes (index < 0)
 * ----------------------------------------- */
#define Node_IndexForBCond             (-1)
#define Node_IndexForPeriodicity       (-2)
#define Node_IndexForOverlappingNode   (-3)

/* 1. node of boundary conditions */
#define Node_EliminateMatrixColumnIndexForBCond(NOD,name) \
        do { \
          int idof = Node_FindUnknownPositionIndex(NOD,name) ; \
          if(Node_GetMatrixColumnIndex(NOD) && idof >= 0) { \
            Node_GetMatrixColumnIndex(NOD)[idof] = Node_IndexForBCond ; \
          } \
        } while(0)

#define Node_EliminateMatrixRowIndexForBCond(NOD,name) \
        do { \
          int idof = Node_FindEquationPositionIndex(NOD,name) ; \
          if(Node_GetMatrixRowIndex(NOD) && idof >= 0) { \
            Node_GetMatrixRowIndex(NOD)[idof] = Node_IndexForBCond ; \
          } \
        } while(0)
        
        
/* 2. slave node of periodic mesh */
#define Node_EliminateMatrixColumnIndexForPeriodicity(NOD,name) \
        do { \
          int idof = Node_FindUnknownPositionIndex(NOD,name) ; \
          if(Node_GetMatrixColumnIndex(NOD) && idof >= 0) { \
            Node_GetMatrixColumnIndex(NOD)[idof] = Node_IndexForPeriodicity ; \
          } \
        } while(0)

#define Node_EliminateMatrixRowIndexForPeriodicity(NOD,name) \
        do { \
          int idof = Node_FindEquationPositionIndex(NOD,name) ; \
          if(Node_GetMatrixRowIndex(NOD) && idof >= 0) { \
            Node_GetMatrixRowIndex(NOD)[idof] = Node_IndexForPeriodicity ; \
          } \
        } while(0)

#define Node_UpdateMatrixColumnIndexForPeriodicity(NODI,NODJ,name) \
        do { \
          if(Node_GetMatrixColumnIndex(NODI) && Node_GetMatrixColumnIndex(NODJ)) { \
            int idof = Node_FindUnknownPositionIndex(NODI,name) ; \
            int jdof = Node_FindUnknownPositionIndex(NODJ,name) ; \
            if(idof >= 0 && jdof >= 0) { \
              int ki = Node_GetMatrixColumnIndex(NODI)[idof] ; \
              int kj = Node_GetMatrixColumnIndex(NODJ)[jdof] ; \
              /* The slave is i, the master is j */ \
              if(kj >= 0 && ki == Node_IndexForPeriodicity) { \
                Node_GetMatrixColumnIndex(NODI)[idof] = kj ; \
              } \
              /* The slave is j, the master is i */ \
              if(ki >= 0 && kj == Node_IndexForPeriodicity) { \
                Node_GetMatrixColumnIndex(NODJ)[jdof] = ki ; \
              } \
            } \
          } \
        } while(0)
            
#define Node_UpdateMatrixRowIndexForPeriodicity(NODI,NODJ,name) \
        do { \
          if(Node_GetMatrixRowIndex(NODI) && Node_GetMatrixRowIndex(NODJ)) { \
            int idof = Node_FindEquationPositionIndex(NODI,name) ; \
            int jdof = Node_FindEquationPositionIndex(NODJ,name) ; \
            if(idof >= 0 && jdof >= 0) { \
              int ki = Node_GetMatrixRowIndex(NODI)[idof] ; \
              int kj = Node_GetMatrixRowIndex(NODJ)[jdof] ; \
              /* The slave is i, the master is j */ \
              if(kj >= 0 && ki == Node_IndexForPeriodicity) { \
                Node_GetMatrixRowIndex(NODI)[idof] = kj ; \
              } \
              /* The slave is j, the master is i */ \
              if(ki >= 0 && kj == Node_IndexForPeriodicity) { \
                Node_GetMatrixRowIndex(NODJ)[jdof] = ki ; \
              } \
            } \
          } \
        } while(0)


/* 3. overlapping nodes of zero-thickness elements */


/* Matrix row/column indexes of selected matrix */

#define Node_IndexForUnselectedMatrix   (-10)

#define Node_GetSelectedMatrixColumnIndexOf(NOD,i,ind) \
        ((i >= 0) && (Node_GetSequentialIndexOfUnknown(NOD)[i] == ind)) ? Node_GetMatrixColumnIndex(NOD)[i] : Node_IndexForUnselectedMatrix

#define Node_GetSelectedMatrixRowIndexOf(NOD,i,ind) \
        ((i >= 0) && (Node_GetSequentialIndexOfUnknown(NOD)[i] == ind)) ? Node_GetMatrixRowIndex(NOD)[i] : Node_IndexForUnselectedMatrix
        



#include "Element.h"
#include "NodeSol.h"
#include "Buffers.h"


struct Node_s {               /* noeud */
  unsigned int index ;        /* node index */
  unsigned short int nel ;    /* nb of elements */
  Element_tt* element ;       /* pointers to elements */
  double* x ;                 /* coordonnees */
  unsigned short int neq ;    /* nombre d'equations au noeud */
  unsigned short int nin ;    /* nombre d'inconnues au noeud */
  char**   inc ;              /* nom des inconnues */
  char**   eqn ;              /* nom des equations */
  int*     sequentialindex ;  /* Sequential indexes of unknowns/equations */
  unsigned short int* i_obj ; /* indices des inconnues dans obj */
  int*    colindex ;          /* column index (unknowns) */
  int*    rowindex ;          /* row index (equations) */
  NodeSol_t* sol ;            /* Nodal solution */
  Buffers_t* buffers ;        /* Buffers */
} ;

#endif
