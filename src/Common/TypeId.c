#include "TypeId.h"
#include "BCond.h"
#include "GenericData.h"
#include "Message.h"
#include "Methods/FEM.h"
#include "Methods/FVM.h"
#include "Exception.h"
#include "InternationalSystemOfUnits.h"
#include <stdio.h>
#include <stdlib.h>





void  (TypeId_Delete)(TypeId_t typ,void* self)
{
  switch(typ) {
    case TypeId_undefined       : break ;
    case TypeId_unsigned_char   : return(free(self));
    case TypeId_char            : return(free(self));
    case TypeId_unsigned_int    : return(free(self));
    case TypeId_short_int       : return(free(self));
    case TypeId_int             : return(free(self));
    case TypeId_unsigned_long   : return(free(self));
    case TypeId_long_int        : return(free(self));
    case TypeId_float           : return(free(self));
    case TypeId_double          : return(free(self));
    case TypeId_long_double     : return(free(self));
    case TypeId_BCond_t         : break ;
    case TypeId_BConds_t        : break ;
    case TypeId_Buffer_t        : break ;
    case TypeId_CommandLine_t   : break ;
    case TypeId_Context_t       : break ;
    case TypeId_Curve_t         : break ;
    case TypeId_Curves_t        : break ;
    case TypeId_CurvesFile_t    : break ;
    case TypeId_DataFile_t      : break ;
    case TypeId_DataSet_t       : break ;
    case TypeId_Date_t          : break ;
    case TypeId_Dates_t         : break ;
    case TypeId_Element_t       : break ;
    case TypeId_Elements_t      : break ;
    case TypeId_ElementSol_t    : break ;
    case TypeId_ElementsSol_t   : break ;
    case TypeId_Exception_t     : return(Exception_Delete(self)) ;
    case TypeId_FEM_t           : return(FEM_Delete(self));
    case TypeId_Field_t         : break ;
    case TypeId_Fields_t        : break ;
    case TypeId_Function_t      : break ;
    case TypeId_Functions_t     : break ;
    case TypeId_FVM_t           : return(FVM_Delete(self));
    case TypeId_GenericData_t   : return(GenericData_Delete(self));
    case TypeId_Geometry_t      : break ;
    case TypeId_Graph_t         : break ;
    case TypeId_ICond_t         : break ;
    case TypeId_IConds_t        : break ;
    case TypeId_InternationalSystemOfUnits_t : return(InternationalSystemOfUnits_Delete(self)) ;
    case TypeId_IntFct_t        : break ;
    case TypeId_IntFcts_t       : break ;
    case TypeId_IterProcess_t   : break ;
    case TypeId_Load_t          : break ;
    case TypeId_Loads_t         : break ;
    case TypeId_Material_t      : break ;
    case TypeId_Materials_t     : break ;
    case TypeId_Math_t          : break ;
    case TypeId_Matrix_t        : break ;
    case TypeId_Mesh_t          : break ;
    case TypeId_Message_t       : return(Message_Delete(self));
    case TypeId_Model_t         : break ;
    case TypeId_Models_t        : break ;
    case TypeId_Module_t        : break ;
    case TypeId_Modules_t       : break ;
    case TypeId_Node_t          : break ;
    case TypeId_Nodes_t         : break ;
    case TypeId_NodeSol_t       : break ;
    case TypeId_NodesSol_t      : break ;
    case TypeId_ObVals_t        : break ;
    case TypeId_ObVal_t         : break ;
    case TypeId_Options_t       : break ;
    case TypeId_OutputFile_t    : break ;
    case TypeId_OutputFiles_t   : break ;
    case TypeId_Periodicity_t   : break ;
    case TypeId_Periodicities_t : break ;
    case TypeId_Point_t         : break ;
    case TypeId_Points_t        : break ;
    case TypeId_Result_t        : break ;
    case TypeId_Results_t       : break ;
    case TypeId_Session_t       : break ;
    case TypeId_ShapeFct_t      : break ;
    case TypeId_ShapeFcts_t     : break ;
    case TypeId_Solution_t      : break ;
    case TypeId_Solutions_t     : break ;
    case TypeId_Solver_t        : break ;
    case TypeId_TextFile_t      : break ;
    case TypeId_TimeStep_t      : break ;
    case TypeId_Unit_t          : break ;
    case TypeId_Units_t         : break ;
    case TypeId_View_t          : break ;
    case TypeId_Views_t         : break ;
    default                     : break ;
  }
  
  Message_FatalError("TypeId_Delete: unknown type") ;
  
  return ;
}
