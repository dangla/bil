#ifndef TYPEID_H
#define TYPEID_H

#ifdef __CPLUSPLUS
extern "C" {
#endif




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
  TypeId_Damage_t,
  TypeId_DataFile_t,
  TypeId_DataSet_t,
  TypeId_Date_t,
  TypeId_Dates_t,
  TypeId_Elasticity_t,
  TypeId_Element_t,
  TypeId_Elements_t,
  TypeId_ElementSol_t,
  TypeId_ElementsSol_t,
  TypeId_Exception_t,
  TypeId_FEM_t,
  TypeId_FEM2_t,
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
  TypeId_InternationalSystemOfUnits_t,
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
  TypeId_Plasticity_t,
  TypeId_Point_t,
  TypeId_Points_t,
  TypeId_Result_t,
  TypeId_Results_t,
  TypeId_Session_t,
  TypeId_ShapeFct_t,
  TypeId_ShapeFcts_t,
  TypeId_Solution_t,
  TypeId_Solutions_t,
  TypeId_Solver_t,
  TypeId_Solvers_t,
  TypeId_TextFile_t,
  TypeId_TimeStep_t,
  TypeId_Unit_t,
  TypeId_Units_t,
  TypeId_View_t,
  TypeId_Views_t,
  /* From SuperLU_DIST */
  TypeId_dScalePermstruct_t,
  TypeId_dLUstruct_t,
  TypeId_gridinfo_t,
  /* From Petsc */
  TypeId_KSP,
  TypeId_PC,
  /* The end */
  TypeId_last
} ;

typedef enum TypeId_e     TypeId_e ;
struct TypeId_s ; typedef struct TypeId_s    TypeId_t ;


#define TypeId_Create(T) \
        TypeId_Create_(TypeId_IdNumber(T),sizeof(T))
        
        
extern TypeId_t* (TypeId_Create_)(const TypeId_e,const size_t) ;
extern void      (TypeId_Delete)(void*) ;



#define TypeId_GetIdNumber(TID)      ((TID)->id)
#define TypeId_GetSize(TID)          ((TID)->size)

        

#include "Utils.h"
#include "Tuple.h"
#include "Algos.h"
#include "Arith.h"


/* T  stands for a real type */
/* TID stands for a TypeId */


/* id number */
#define TypeId_IdNumber(T) \
        Tuple_SEQ(TypeId_IdNumber_(T))


/* Test the type */
#if 0
#define TypeId_Is(TID,T) \
        (TypeId_GetIdNumber(TID) == TypeId_IdNumber(T))
        
#define TypeId_SetTo(TID,T) \
        do {TypeId_GetIdNumber(TID) = TypeId_IdNumber(T) ;} while(0)
#endif
        
        
/* Code */
#if 0
#define TypeId_Code(T) \
        Utils_CAT(TypeId_IdNumber(T),_c)
        

/* Return a type from its CODE */
#define TypeId_Type(CODE) \
        Tuple_ELEM(Arith_INCR(CODE),Tuple_TUPLE(TypeId_List))
#endif
        
        


/* Implementation */

#define TypeId_IdNumber_(...) \
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









#define TypeId_List \
        undefined, \
        unsigned char, \
        char, \
        unsigned int, \
        short int, \
        int, \
        unsigned long, \
        long int, \
        float, \
        double, \
        long double, \
        BCond_t, \
        BConds_t, \
        Buffer_t, \
        CommandLine_t, \
        Context_t, \
        Curve_t, \
        Curves_t, \
        CurvesFile_t, \
        Damage_t, \
        DataFile_t, \
        DataSet_t, \
        Date_t, \
        Dates_t, \
        Elasticity_t, \
        Element_t, \
        Elements_t, \
        ElementSol_t, \
        ElementsSol_t, \
        Exception_t, \
        FEM_t, \
        Field_t, \
        Fields_t, \
        Function_t, \
        Functions_t, \
        FVM_t, \
        GenericData_t, \
        Geometry_t, \
        Graph_t, \
        ICond_t, \
        IConds_t, \
        InternationalSystemOfUnits_t, \
        IntFct_t, \
        IntFcts_t, \
        IterProcess_t, \
        Load_t, \
        Loads_t, \
        Material_t, \
        Materials_t, \
        Math_t, \
        Matrix_t, \
        Mesh_t, \
        Message_t, \
        Model_t, \
        Models_t, \
        Module_t, \
        Modules_t, \
        Node_t, \
        Nodes_t, \
        NodeSol_t, \
        NodesSol_t, \
        ObVals_t, \
        ObVal_t, \
        Options_t, \
        OutputFile_t, \
        OutputFiles_t, \
        Periodicity_t, \
        Periodicities_t, \
        Plasticity_t, \
        Point_t, \
        Points_t, \
        Result_t, \
        Results_t, \
        Session_t, \
        ShapeFct_t, \
        ShapeFcts_t, \
        Solution_t, \
        Solutions_t, \
        Solver_t, \
        Solvers_t, \
        TextFile_t, \
        TimeStep_t, \
        Unit_t, \
        Units_t, \
        View_t, \
        Views_t, \
        last



/* Codes */
  /* From C */
#define  TypeId_undefined_c                      0
#define  TypeId_unsigned_char_c                  1
#define  TypeId_char_c                           2
#define  TypeId_unsigned_int_c                   3
#define  TypeId_short_int_c                      4
#define  TypeId_int_c                            5
#define  TypeId_unsigned_long_c                  6
#define  TypeId_long_int_c                       7
#define  TypeId_float_c                          8
#define  TypeId_double_c                         9
#define  TypeId_long_double_c                    10
  /* From Bil */
#define  TypeId_BCond_t_c                        11
#define  TypeId_BConds_t_c                       12
#define  TypeId_Buffer_t_c                       13
#define  TypeId_CommandLine_t_c                  14
#define  TypeId_Context_t_c                      15
#define  TypeId_Curve_t_c                        16
#define  TypeId_Curves_t_c                       17
#define  TypeId_CurvesFile_t_c                   18
#define  TypeId_Damage_t_c                       19
#define  TypeId_DataFile_t_c                     20
#define  TypeId_DataSet_t_c                      21
#define  TypeId_Date_t_c                         22
#define  TypeId_Dates_t_c                        23
#define  TypeId_Elasticity_t_c                   24
#define  TypeId_Element_t_c                      25
#define  TypeId_Elements_t_c                     26
#define  TypeId_ElementSol_t_c                   27
#define  TypeId_ElementsSol_t_c                  28
#define  TypeId_Exception_t_c                    29
#define  TypeId_FEM_t_c                          30
#define  TypeId_FEM2_t_c
#define  TypeId_Field_t_c                        31
#define  TypeId_Fields_t_c                       32
#define  TypeId_Function_t_c                     33
#define  TypeId_Functions_t_c                    34
#define  TypeId_FVM_t_c                          35
#define  TypeId_GenericData_t_c                  36
#define  TypeId_Geometry_t_c                     37
#define  TypeId_Graph_t_c                        38
#define  TypeId_ICond_t_c                        39
#define  TypeId_IConds_t_c                       40
#define  TypeId_InternationalSystemOfUnits_t_c   41
#define  TypeId_IntFct_t_c                       42
#define  TypeId_IntFcts_t_c                      43
#define  TypeId_IterProcess_t_c                  44  
#define  TypeId_Load_t_c                         45
#define  TypeId_Loads_t_c                        46
#define  TypeId_Material_t_c                     47
#define  TypeId_Materials_t_c                    48
#define  TypeId_Math_t_c                         49
#define  TypeId_Matrix_t_c                       50
#define  TypeId_Mesh_t_c                         51
#define  TypeId_Message_t_c                      52
#define  TypeId_Model_t_c                        53
#define  TypeId_Models_t_c                       54
#define  TypeId_Module_t_c                       55
#define  TypeId_Modules_t_c                      56
#define  TypeId_Node_t_c                         57
#define  TypeId_Nodes_t_c                        58
#define  TypeId_NodeSol_t_c                      59
#define  TypeId_NodesSol_t_c                     60
#define  TypeId_ObVals_t_c                       61
#define  TypeId_ObVal_t_c                        62
#define  TypeId_Options_t_c                      63
#define  TypeId_OutputFile_t_c                   64
#define  TypeId_OutputFiles_t_c                  65
#define  TypeId_Periodicity_t_c                  66
#define  TypeId_Periodicities_t_c                67
#define  TypeId_Plasticity_t_c                   68
#define  TypeId_Point_t_c                        69
#define  TypeId_Points_t_c                       70
#define  TypeId_Result_t_c                       71
#define  TypeId_Results_t_c                      72
#define  TypeId_Session_t_c                      73
#define  TypeId_ShapeFct_t_c                     74
#define  TypeId_ShapeFcts_t_c                    75
#define  TypeId_Solution_t_c                     76
#define  TypeId_Solutions_t_c                    77
#define  TypeId_Solver_t_c                       78
#define  TypeId_Solvers_t_c                      79
#define  TypeId_TextFile_t_c                     80
#define  TypeId_TimeStep_t_c                     81
#define  TypeId_Unit_t_c                         82
#define  TypeId_Units_t_c                        83
#define  TypeId_View_t_c                         84
#define  TypeId_Views_t_c                        85
        


struct TypeId_s {
  TypeId_e id ;           /* Id number of the data */
  size_t size ;           /* Size of elementary data */
} ;




#ifdef __CPLUSPLUS
}
#endif
#endif
