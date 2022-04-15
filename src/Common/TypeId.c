#include <stdio.h>
#include <stdlib.h>

#include "TypeId.h"
#include "BCond.h"
#include "Damage.h"
#include "Elasticity.h"
#include "Plasticity.h"
#include "ElementsSol.h"
#include "GenericData.h"
#include "Message.h"
#include "FEM.h"
#include "FVM.h"
#include "Exception.h"
#include "InternationalSystemOfUnits.h"
#include "Math_.h"





void  (TypeId_Delete)(TypeId_t typ,void* self)
{
  switch(typ) {
    case TypeId_undefined       : break ;
    case TypeId_unsigned_char   : return ;
    case TypeId_char            : return ;
    case TypeId_double          : return ;
    case TypeId_long_double     : return ;
    case TypeId_float           : return ;
    case TypeId_unsigned_int    : return ;
    case TypeId_short_int       : return ;
    case TypeId_int             : return ;
    case TypeId_unsigned_long   : return ;
    case TypeId_long_int        : return ;
    case TypeId_BCond_t         : break ;
    case TypeId_BConds_t        : break ;
    case TypeId_Buffer_t        : break ;
    case TypeId_CommandLine_t   : break ;
    case TypeId_Context_t       : break ;
    case TypeId_Curve_t         : break ;
    case TypeId_Curves_t        : break ;
    case TypeId_CurvesFile_t    : break ;
    case TypeId_Damage_t        : Damage_Delete(self); return ;
    case TypeId_DataFile_t      : break ;
    case TypeId_DataSet_t       : DataSet_Delete(self); return ;
    case TypeId_Date_t          : break ;
    case TypeId_Dates_t         : break ;
    case TypeId_Elasticity_t    : Elasticity_Delete(self); return ;
    case TypeId_Element_t       : break ;
    case TypeId_Elements_t      : break ;
    case TypeId_ElementSol_t    : break ;
    case TypeId_ElementsSol_t   : ElementsSol_Delete(self); return ;
    case TypeId_Exception_t     : Exception_Delete(self); return ;
    case TypeId_FEM_t           : FEM_Delete(self); return ;
    case TypeId_Field_t         : break ;
    case TypeId_Fields_t        : break ;
    case TypeId_Function_t      : break ;
    case TypeId_Functions_t     : break ;
    case TypeId_FVM_t           : FVM_Delete(self); return ;
    case TypeId_GenericData_t   : GenericData_Delete(self); return ;
    case TypeId_Geometry_t      : break ;
    case TypeId_Graph_t         : break ;
    case TypeId_ICond_t         : break ;
    case TypeId_IConds_t        : break ;
    case TypeId_InternationalSystemOfUnits_t : InternationalSystemOfUnits_Delete(self); return ;
    case TypeId_IntFct_t        : break ;
    case TypeId_IntFcts_t       : break ;
    case TypeId_IterProcess_t   : break ;
    case TypeId_Load_t          : break ;
    case TypeId_Loads_t         : break ;
    case TypeId_Material_t      : break ;
    case TypeId_Materials_t     : break ;
    case TypeId_Math_t          : Math_Delete(self); return ;
    case TypeId_Matrix_t        : break ;
    case TypeId_Mesh_t          : break ;
    case TypeId_Message_t       : Message_Delete(self); return ;
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
    case TypeId_Options_t       : Options_Delete(self); return ;
    case TypeId_OutputFile_t    : break ;
    case TypeId_OutputFiles_t   : break ;
    case TypeId_Periodicity_t   : break ;
    case TypeId_Periodicities_t : break ;
    case TypeId_Plasticity_t    : Plasticity_Delete(self); return ;
    case TypeId_Point_t         : break ;
    case TypeId_Points_t        : break ;
    case TypeId_Result_t        : break ;
    case TypeId_Results_t       : break ;
    case TypeId_Session_t       : break ;
    case TypeId_ShapeFct_t      : break ;
    case TypeId_ShapeFcts_t     : break ;
    case TypeId_Solution_t      : break ;
    case TypeId_Solutions_t     : Solutions_Delete(self); return ;
    case TypeId_Solver_t        : Solver_Delete(self); return ;
    case TypeId_Solvers_t       : Solvers_Delete(self); return ;
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
