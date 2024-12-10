#ifndef NODE_H
#define NODE_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct Node_s         ; typedef struct Node_s         Node_t ;
typedef Node_t*   Node_tt ;


#include "Utils.h"
#include "Solutions.h"
#include "Buffers.h"


extern Node_t*  (Node_New)(const int) ;
extern void     (Node_CreateMore)(Node_t*,Buffers_t*) ;
extern void     (Node_Delete)(void*) ;
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
#define Node_GetNodeIndex(NOD)            ((NOD)->NodeIndex)
#define Node_GetCoordinate(NOD)           ((NOD)->Coordinate)
#define Node_GetNbOfEquations(NOD)        ((NOD)->NbOfEquations)
#define Node_GetNbOfUnknowns(NOD)         ((NOD)->NbOfUnknowns)
#define Node_GetNameOfEquation(NOD)       ((NOD)->NameOfEquation)
#define Node_GetNameOfUnknown(NOD)        ((NOD)->NameOfUnknown)
#define Node_GetSequentialIndexOfUnknown(NOD)        ((NOD)->SequentialIndexOfUnknown)
#define Node_GetObValIndex(NOD)           ((NOD)->ObValIndex)
#define Node_GetMatrixColumnIndex(NOD)    ((NOD)->MatrixColumnIndex)
#define Node_GetMatrixRowIndex(NOD)       ((NOD)->MatrixRowIndex)
#define Node_GetNbOfElements(NOD)         ((NOD)->NbOfElements)
#define Node_GetPointerToElement(NOD)     ((NOD)->PointerToElement)
#define Node_GetBuffers(NOD)              ((NOD)->Buffers)
#define Node_GetSolutions(NOD)            ((NOD)->Solutions)



/* Buffer */
#define Node_GetBuffer(NOD) \
        Buffers_GetBufferOfCurrentThread(Node_GetBuffers(NOD))


/* Access to elements
 * ------------------ */
#define Node_GetElement(NOD,i) \
        (Node_GetPointerToElement(NOD)[i])



/* Access to solution */
#define Node_GetSolution(NOD) \
        Solutions_GetSolution(Node_GetSolutions(NOD))
        
#define Node_GetPreviousSolution(NOD) \
        Solutions_GetPreviousSolution(Node_GetSolutions(NOD))
        
#define Node_GetNextSolution(NOD) \
        Solutions_GetNextSolution(Node_GetSolutions(NOD))
        
#define Node_GetSolutionInDistantPast(NOD,i) \
        Solutions_GetSolutionInDistantPast(Node_GetSolutions(NOD),i)


/* Access to NodeSol
 * ----------------- */
#define Node_GetNodeSol(NOD) \
        (Solution_GetNodeSol(Node_GetSolution(NOD)) + Node_GetNodeIndex(NOD))

#define Node_GetCurrentNodeSol(NOD) \
        Node_GetNodeSol(NOD)
        
#define Node_GetPreviousNodeSol(NOD) \
        (Solution_GetNodeSol(Node_GetPreviousSolution(NOD)) + Node_GetNodeIndex(NOD))
        
#define Node_GetNextNodeSol(NOD) \
        (Solution_GetNodeSol(Node_GetNextSolution(NOD)) + Node_GetNodeIndex(NOD))

#define Node_GetNodeSolInDistantPast(NOD,i) \
        (Solution_GetNodeSol(Node_GetSolutionInDistantPast(NOD,i)) + Node_GetNodeIndex(NOD))



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
//#include "NodeSol.h"
#include "Buffers.h"
#include "Solutions.h"


struct Node_s {
  Element_tt* PointerToElement ;
  Buffers_t* Buffers ;
  Solutions_t* Solutions ;             /* Pointer to the global solutions */
  double* Coordinate ;
  char**   NameOfUnknown ;
  char**   NameOfEquation ;
  int*     SequentialIndexOfUnknown ;  /* Sequential indexes of unknowns/equations */
  unsigned int* ObValIndex ;     /* indices des inconnues dans obj */
  int*    MatrixColumnIndex ;          /* column index (unknowns) */
  int*    MatrixRowIndex ;             /* row index (equations) */
  unsigned int NodeIndex ;             /* node index */
  unsigned int NbOfElements ;
  unsigned int NbOfEquations ;   /* nombre d'equations au noeud */
  unsigned int NbOfUnknowns ;    /* nombre d'inconnues au noeud */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
