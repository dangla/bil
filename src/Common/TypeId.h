#ifndef TYPEID_H
#define TYPEID_H


/* Type identifiers */
enum TypeId_e {
  /* From C */
  TypeId_undefined,
  TypeId_unsigned_char,
  TypeId_char,
  TypeId_unsigned_int,
  TypeId_short_int,
  TypeId_int,
  TypeId_unsigned_long,
  TypeId_long_int,
  TypeId_float,
  TypeId_double,
  TypeId_long_double,
  /* From Bil */
  TypeId_BCond_t,
  TypeId_BConds_t,
  TypeId_Buffer_t,
  TypeId_CommandLine_t,
  TypeId_Context_t,
  TypeId_Curve_t,
  TypeId_Curves_t,
  TypeId_CurvesFile_t,
  TypeId_DataFile_t,
  TypeId_DataSet_t,
  TypeId_Date_t,
  TypeId_Dates_t,
  TypeId_Element_t,
  TypeId_Elements_t,
  TypeId_ElementSol_t,
  TypeId_ElementsSol_t,
  TypeId_Exception_t,
  TypeId_FEM_t,
  TypeId_Field_t,
  TypeId_Fields_t,
  TypeId_Function_t,
  TypeId_Functions_t,
  TypeId_FVM_t,
  TypeId_GenericData_t,
  TypeId_Geometry_t,
  TypeId_Graph_t,
  TypeId_ICond_t,
  TypeId_IConds_t,
  TypeId_IntFct_t,
  TypeId_IntFcts_t,
  TypeId_IterProcess_t,
  TypeId_Load_t,
  TypeId_Loads_t,
  TypeId_Material_t,
  TypeId_Materials_t,
  TypeId_Math_t,
  TypeId_Matrix_t,
  TypeId_Mesh_t,
  TypeId_Message_t,
  TypeId_Model_t,
  TypeId_Models_t,
  TypeId_Module_t,
  TypeId_Modules_t,
  TypeId_Node_t,
  TypeId_Nodes_t,
  TypeId_NodeSol_t,
  TypeId_NodesSol_t,
  TypeId_ObVals_t,
  TypeId_ObVal_t,
  TypeId_Options_t,
  TypeId_OutputFile_t,
  TypeId_OutputFiles_t,
  TypeId_Periodicity_t,
  TypeId_Periodicities_t,
  TypeId_Point_t,
  TypeId_Points_t,
  TypeId_Result_t,
  TypeId_Results_t,
  TypeId_ShapeFct_t,
  TypeId_ShapeFcts_t,
  TypeId_Solution_t,
  TypeId_Solutions_t,
  TypeId_Solver_t,
  TypeId_TextFile_t,
  TypeId_TimeStep_t,
  TypeId_Unit_t,
  TypeId_Units_t,
  TypeId_View_t,
  TypeId_Views_t,
  TypeId_last
} ;

typedef enum TypeId_e     TypeId_t ;
        

#include "Utils.h"
#include "Tuple.h"
#include "Algos.h"


/* T  stands for a real type */
/* ID stands for a TypeId */


/* Create a TypeId */
#define TypeId_Create(T) \
        Tuple_SEQ(TypeId_Create_(T))


/* Test the type */
#define TypeId_Is(ID,T) \
        (ID == TypeId_Create(T))
        
#define TypeId_SetTo(ID,T) \
        do {ID = TypeId_Create(T) ;} while(0)


/* Implementation */

#define TypeId_Create_(...) \
        (Utils_CAT(TypeId_,__VA_ARGS__))
        
#define TypeId_long \
        TypeId_long_(
        
#define TypeId_long_(...) \
        Utils_CAT(TypeId_long_,__VA_ARGS__))
        
#define TypeId_short \
        TypeId_short_(
        
#define TypeId_short_(...) \
        Utils_CAT(TypeId_short_,__VA_ARGS__))
        
#define TypeId_unsigned \
        TypeId_unsigned_(
        
#define TypeId_unsigned_(...) \
        Utils_CAT(TypeId_unsigned_,__VA_ARGS__))



#endif
