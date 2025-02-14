#include "FEM2.h"
#include "FEM.h"
#include "Message.h"
#include "Matrix.h"
#include "Residu.h"
#include "Session.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <assert.h>



static FEM2_t*  (FEM2_Create)(void) ;

/* 
   Extern Functions 
*/
FEM2_t* (FEM2_Create)(void)
{
  FEM2_t*  fem2    = (FEM2_t*) Mry_New(FEM2_t) ;
  
  
  /* Space allocation for buffer */
  {
    Buffers_t* buf = Buffers_Create(FEM2_SizeOfBuffer) ;
    
    FEM2_GetBuffers(fem2) = buf ;
  }
  
  return(fem2) ;
}



void (FEM2_Delete)(void* self)
{
  FEM2_t* fem2 = (FEM2_t*) self ;
  
  {
    Buffers_t* buf = FEM2_GetBuffers(fem2) ;
    
    if(buf) {
      Buffers_Delete(buf)  ;
      free(buf) ;
      FEM2_GetBuffers(fem2) = NULL ;
    }
  }
}



FEM2_t*  (FEM2_GetInstance)(DataSet_t* dataset,Solvers_t* solvers,Solutions_t* sols_n,Solutions_t* sols)
{
  GenericData_t* gdat = Session_FindGenericData(FEM2_t,"FEM2") ;
  
  if(!gdat) {
    FEM2_t* fem2 = FEM2_Create() ;
    
    gdat = GenericData_Create(1,fem2,FEM2_t,"FEM2") ;
    
    Session_AddGenericData(gdat) ;
    
    assert(gdat == Session_FindGenericData(FEM2_t,"FEM2")) ;
  }
  
  {
    FEM2_t* fem2 = (FEM2_t*) GenericData_GetData(gdat) ;
  
    FEM2_GetDataSet(fem2) = dataset ;
    FEM2_GetSolvers(fem2) = solvers ;
    FEM2_GetCurrentSolutions(fem2) = sols ;
    FEM2_GetPreviousSolutions(fem2) = sols_n ;
    FEM2_FreeBuffer(fem2) ;
  
    return(fem2) ;
  }
}





int (FEM2_HomogenizeTangentStiffnessTensor)(FEM2_t* fem2,double t,double dt,double* c)
/** Compute the homogenized tangent stiffness tensor of a microstructure.
 *  Note: the method is presented in "NumericalHomogenization.pdf" 
 *  which can be found under the same directory.
 *  Inputs:
 *  - fem2: should contain
 *         - mesh: the mesh of the microstructure
 *         - solver: the solver created from mesh
 *  - t, dt: the time and the time step
 *  Output:
 *  - c: the stiffness tensor as an array of 81 doubles.
 *  Return 0 if it succeeds or else if it fails.
 */
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  DataSet_t* dataset = FEM2_GetDataSet(fem2) ;
  Solver_t* solver = FEM2_GetSolver(fem2) ;
  Solutions_t* sols = FEM2_GetCurrentSolutions(fem2) ;
  Mesh_t* mesh = DataSet_GetMesh(dataset) ;
  int dim = Mesh_GetDimension(mesh) ;
  Matrix_t* a = Solver_GetMatrix(solver) ;
  Residu_t* residu = Solver_GetResidu(solver) ;
  double*   b = (double*) Residu_GetRHS(residu) ;
  double*   u = (double*) Residu_GetSolution(residu) ;
  int nrhs = Residu_GetNbOfRHS(residu) ;
  int ncol = Residu_GetLengthOfRHS(residu) ;
  double*  pb[9] = {b,b+ncol,b+2*ncol,b+ncol,b+3*ncol,b+4*ncol,b+2*ncol,b+4*ncol,b+5*ncol} ;
  double*  pu[9] = {u,u+ncol,u+2*ncol,u+ncol,u+3*ncol,u+4*ncol,u+2*ncol,u+4*ncol,u+5*ncol} ;
  double    E[9] = {1,0.5,0.5,0.5,1,0.5,0.5,0.5,1} ;
  
  if(nrhs < 6) {
    Message_FatalError("FEM2_HomogenizeStiffnessTensor: the nb of rhs is %d whereas it must be at least 6",nrhs) ;
  }
  
    
  Matrix_SetValuesToZero(a) ;
  Residu_SetValuesToZero(residu) ;
  
  Solutions_InitializeMeshPointers(sols,mesh) ;


  /* Initializations */
  {
    int j ;
    
    for(j = 0 ; j < 81 ; j++) c[j] = 0 ;
  }
  

  /* The matrix and the r.h.s. */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
#define NE (Element_MaxNbOfNodes*Model_MaxNbOfEquations)
    double ke[NE*NE] ;
#undef NE
    int ie ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Material_t* mat = Element_GetMaterial(el + ie) ;
    
      if(mat) {
      
        Element_FreeBuffer(el + ie) ;
        {
          int i = Element_ComputeMatrix(el + ie,t,dt,ke) ;
        
          if(i != 0) return(i) ;
        }
      
        Matrix_AssembleElementMatrix(a,el+ie,ke) ;
      
        {
          int  nn = Element_GetNbOfNodes(el + ie) ;
          int neq = Element_GetNbOfEquations(el + ie) ;
          int ndof = nn*neq ;
          int n ;
                
          for(n = 0 ; n < nn ; n++) {
            Node_t* node_n = Element_GetNode(el + ie,n) ;
            double* x_n = Node_GetCoordinate(node_n) ;
            int m ;
                
            for(m = 0 ; m < nn ; m++) {
              Node_t* node_m = Element_GetNode(el + ie,m) ;
              double* x_m = Node_GetCoordinate(node_m) ;
              int i ;
        
              for(i = 0 ; i < dim ; i++) {
                int ni = n*neq + i ;
                int jj_row = Element_GetEquationPosition(el + ie)[ni] ;
                
                /*  The r.h.s. stored in pb */
                if(jj_row >= 0) {
                  int row_i = Node_GetMatrixRowIndex(node_n)[jj_row] ;
                  int j ;
          
                  if(row_i >= 0) {
                    for(j = 0 ; j < dim ; j++) {
                      int mj = m*neq + j ;
                      int ij = ni * ndof + mj ;
                      int k ;
                      
                      for(k = j ; k < dim ; k++) {
                        int      jk = 3 * j + k ;
                        double* bjk = pb[jk] ;

                        bjk[row_i] -= ke[ij] * E[jk] * x_m[k] ;
                      }
                    }
                  }
                }
                
                /*  First part of C */
                {
                  int j ;
          
                  for(j = 0 ; j < dim ; j++) {
                    int mj = m*neq + j ;
                    int ij = ni * ndof + mj ;
                    int k ;
                      
                    for(k = 0 ; k < dim ; k++) {
                      int l ;
          
                      for(l = 0 ; l < dim ; l++) {
                        C(k,i,l,j) += x_n[k] * ke[ij] * x_m[l] ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  /* The solutions */
  {
    int k ;

    for(k = 0 ; k < dim ; k++) {
      int l ;
          
      for(l = k ; l < dim ; l++) {
        Solver_GetRHS(solver) = pb[3 * k + l] ;
        Solver_GetSolution(solver) = pu[3 * k + l] ;
        Solver_Solve(solver) ;
      }
    }
    
    Solver_GetRHS(solver) = b ;
    Solver_GetSolution(solver) = u ;
  }


  /* The macroscopic tangent stiffness matrix */
  {
    int i ;
        
    for(i = 0 ; i < dim ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        int k ;
          
        for(k = 0 ; k < dim ; k++) {
          int l ;
          
          for(l = 0 ; l < dim ; l++) {
            double  ub = 0. ;
            
            {
              double* uik = pu[3 * i + k] ;
              double* blj = pb[3 * l + j] ;
              int p ;
              
              for(p = 0 ; p < ncol ; p++) {
                ub += uik[p] * blj[p] ;
              }
            }
            
            {
              double  Eik =  E[3 * i + k] ;
              double  Elj =  E[3 * l + j] ;

              C(k,i,j,l) -= ub / (Eik * Elj) ;
            }
          }
        }
      }
    }
  }

  {
    double vol = FEM_ComputeVolume(mesh) ;
    int j ;
    
    for(j = 0 ; j < 81 ; j++) c[j] /= vol ;
  }
  
  return(0) ;
#undef C
}




int (FEM2_ComputeHomogenizedStressTensor)(FEM2_t* fem2,double t,double dt,double* dstrain,double* stress)
/** Compute the homogenized stress tensor of a microstructure
 *  subjected to an increment of macro-strain "dstrain".
 *  The stress at time "t" is computed by starting from the
 *  solution stored in "sols_n" obtained at time "t-dt".
 *  The number of time steps from t-dt to t is determined by the
 *  objective variations given in the dataset.
 *  Inputs:
 *  - fem2: should contain
 *          - dataset:
 *          - solver:
 *          - sols_n, sols:
 *  - t, dt:
 *  - dstrain: strain increment as an array of 9 doubles
 *  Output:
 *  - stress: stress tensor as an array of 9 doubles.
 */
{
  DataSet_t* dataset = FEM2_GetDataSet(fem2) ;
  Solver_t* solver = FEM2_GetSolver(fem2) ;
  Solutions_t* sols_n = FEM2_GetPreviousSolutions(fem2) ;
  Solutions_t* sols = FEM2_GetCurrentSolutions(fem2) ;
  
  /* Set input data of the microstructure */
  {
    /* Update the time function */
    {
      Functions_t* fcts = DataSet_GetFunctions(dataset) ;
      int nfcts = 9 ;
      int j ;
      
      for(j = 0 ; j < nfcts ; j++) {
        Function_t* func = Functions_GetFunction(fcts) + j ;
        double* x = Function_GetXValue(func) ;
        double* f = Function_GetFValue(func) ;

        x[0] = t - dt ;
        x[1] = t ;
        f[0] = 0 ;
        f[1] = dstrain[j] ;
      }
    }
  }
    
  /* Compute the microstructure */
  {
    Module_t* module = DataSet_GetModule(dataset) ;
    Dates_t*   dates  = DataSet_GetDates(dataset) ;
    Date_t*    date   = Dates_GetDate(dates) ;
    TimeStep_t*  timestep  = DataSet_GetTimeStep(dataset) ;
    double dtini = TimeStep_GetInitialTimeStep(timestep) ;
    double t_n = Solutions_GetTime(sols_n) ;
    
    {
      /* t should be equal to t_n + dt. */
      if(fabs(t - dt - t_n) > 1.e-4*dtini) {
        Message_FatalError("FEM2_ComputeStressTensor: t_n = %e ; t - dt = %e",t_n,t-dt) ;
      }
      
      /* There are 2 dates */
      Date_GetTime(date) = t_n ;
      Date_GetTime(date + 1) = t ;
    
      {
        Solution_t* sol   = Solutions_GetSolution(sols) ;
        Solution_t* sol_n = Solutions_GetSolution(sols_n) ;
        
        Solution_Copy(sol,sol_n) ;
      }
      
      {
        //Mesh_t* mesh = DataSet_GetMesh(dataset) ;
        int i ;
        
        /* Both method work */
        #if 0
        i = Module_SolveProblem(module,dataset,sols,solver,NULL) ;
        #else
        {
          Module_InitializeProblem(module,dataset,sols) ;
          i = Module_Increment(module,dataset,sols,solver,NULL,t_n,t) ;
        }
        #endif
        
        if(i < 0) {
          return(i) ;
        }
      }
    }
  }

  /* Backup stresses as averaged stresses */
  {
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;

    Solutions_InitializeMeshPointers(sols,mesh) ;
    FEM_AverageStresses(mesh,stress) ;
  }
  
  return(0) ;
}



void (FEM2_InitializeMicrostructureDataSet)(FEM2_t* fem2)
{
  DataSet_t* dataset = FEM2_GetDataSet(fem2) ;
  
  /* Set input data of the microstructure */
  {
    /* Update the macro-gradient and the macro-fctindex */
    {
      Materials_t* mats = DataSet_GetMaterials(dataset) ;
      int nmats = Materials_GetNbOfMaterials(mats) ;
      int j ;
    
      for(j = 0 ; j < nmats ; j++) {
        Material_t* mat = Materials_GetMaterial(mats) + j ;
        Model_t* model = Material_GetModel(mat) ;
        Model_ComputePropertyIndex_t* pidx = Model_GetComputePropertyIndex(model) ;
        
        if(!pidx) {
          arret("FEM2_InitializeMicrostructureDataSet: Model_GetComputePropertyIndex(model) undefined") ;
        }
        
        {
          double* grd = Material_GetProperty(mat) + pidx("macro-gradient") ;
          double* fid = Material_GetProperty(mat) + pidx("macro-fctindex") ;
          int i ;
    
          for(i = 0 ; i < 9 ; i++) {
            grd[i] = 1 ;
            fid[i] = i + 1 ;
          }
        }
      }
    }
    
    
    /* Check and update the function of time */
    {
      Functions_t* fcts = DataSet_GetFunctions(dataset) ;
      int nfcts = Functions_GetNbOfFunctions(fcts) ;
      int j ;
      
      if(nfcts < 9) {
        Message_FatalError("FEM2_InitializeMicrostructureDataSet: the nb of functions is %d but it must be 9 at least",nfcts) ;
      }
      
      nfcts = 9 ;
      for(j = 0 ; j < nfcts ; j++) {
        Function_t* func = Functions_GetFunction(fcts) + j ;
        int npts = Function_GetNbOfPoints(func) ;
          
        if(npts < 2) {
          Message_FatalError("FEM2_InitializeMicrostructureDataSet: the nb of points is %d but it must be 2 at least",npts) ;
        }
          
        {
          double* t = Function_GetXValue(func) ;
          double* f = Function_GetFValue(func) ;
            
          Function_GetNbOfPoints(func) = 2 ;
          t[0] = 0 ;
          t[1] = 1 ;
          f[0] = 0 ;
          f[1] = 1 ;
        }
      }
    }
    
    /* The dates */
    {
      {
        Dates_t* dates = DataSet_GetDates(dataset) ;
        int     nbofdates  = Dates_GetNbOfDates(dates) ;
          
        if(nbofdates < 2) {
          Message_FatalError("FEM2_InitializeMicrostructureDataSet: the nb of dates is %d but it must be 2 at least",nbofdates) ;
        }
      
        Dates_GetNbOfDates(dates) = 2 ;
      }
    }
  }
}




#if 0
int (FEM2_HomogenizeTangentStiffnessTensor1)(Mesh_t* mesh,Solver_t* solver,double t,double dt,double* c)
/** Compute the homogenized tangent stiffness tensor of a microstructure
 *  by the finite element method.
 *  Inputs:
 *  - mesh: contained the mesh of the microstructure
 *  - solver: the solver created from mesh
 *  - t, dt: the time and the time step
 *  Output:
 *  - c: the stiffness tensor as an array of 81 doubles.
 *  Return 0 if succeeds or else if fails.
 */
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  int dim = Mesh_GetDimension(mesh) ;
  Matrix_t* a = Solver_GetMatrix(solver) ;
  Residu_t* residu = Solver_GetResidu(solver) ;
  double*   b = (double*) Residu_GetRHS(residu) ;
  double*   u = (double*) Residu_GetSolution(residu) ;
  int nrhs = Residu_GetNbOfRHS(residu) ;
  int ncol = Residu_GetLengthOfRHS(residu) ;
  double*  pb[9] = {b,b+ncol,b+2*ncol,b+ncol,b+3*ncol,b+4*ncol,b+2*ncol,b+4*ncol,b+5*ncol} ;
  double*  pu[9] = {u,u+ncol,u+2*ncol,u+ncol,u+3*ncol,u+4*ncol,u+2*ncol,u+4*ncol,u+5*ncol} ;
  double    E[9] = {1,0.5,0.5,0.5,1,0.5,0.5,0.5,1} ;
  
  if(nrhs < 6) {
    arret("FEM2_HomogenizeStiffnessTensor") ;
  }


  /* Initializations */
  {
    int j ;
    
    for(j = 0 ; j < 81 ; j++) c[j] = 0 ;
  }
  
    
  Matrix_SetValuesToZero(a) ;
  Residu_SetValuesToZero(residu) ;
  

  /* The matrix and the r.h.s. */
  {
    int n_el = Mesh_GetNbOfElements(mesh) ;
    Element_t* el = Mesh_GetElement(mesh) ;
#define NE (Element_MaxNbOfNodes*Model_MaxNbOfEquations)
    double ke[NE*NE] ;
#undef NE
    int ie ;
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Material_t* mat = Element_GetMaterial(el + ie) ;
    
      if(mat) {
      
        Element_FreeBuffer(el + ie) ;
        {
          int i = Element_ComputeMatrix(el + ie,t,dt,ke) ;
        
          if(i != 0) return(i) ;
        }
      
        Matrix_AssembleElementMatrix(a,el+ie,ke) ;
      
        {
          int  nn = Element_GetNbOfNodes(el + ie) ;
          int neq = Element_GetNbOfEquations(el + ie) ;
          int ndof = nn*neq ;
          int n ;
                
          for(n = 0 ; n < nn ; n++) {
            Node_t* node_n = Element_GetNode(el + ie,n) ;
            double* x_n = Node_GetCoordinate(node_n) ;
            int m ;
                
            for(m = 0 ; m < nn ; m++) {
              Node_t* node_m = Element_GetNode(el + ie,m) ;
              double* x_m = Node_GetCoordinate(node_m) ;
              int i ;
        
              for(i = 0 ; i < dim ; i++) {
                int ni = n*neq + i ;
                int jj_row = Element_GetEquationPosition(el + ie)[ni] ;
                
                /*  The r.h.s. stored in pb */
                if(jj_row >= 0) {
                  int row_i = Node_GetMatrixRowIndex(node_n)[jj_row] ;
                  int j ;
          
                  if(row_i >= 0) {
                    for(j = 0 ; j < dim ; j++) {
                      int mj = m*neq + j ;
                      int ij = ni * ndof + mj ;
                      int k ;
                      
                      for(k = j ; k < dim ; k++) {
                        int      jk = 3 * j + k ;
                        double* bjk = pb[jk] ;

                        bjk[row_i] -= ke[ij] * E[jk] * x_m[k] ;
                      }
                    }
                  }
                }
                
                /*  First part of C */
                {
                  int j ;
          
                  for(j = 0 ; j < dim ; j++) {
                    int mj = m*neq + j ;
                    int ij = ni * ndof + mj ;
                    int k ;
                      
                    for(k = 0 ; k < dim ; k++) {
                      int l ;
          
                      for(l = 0 ; l < dim ; l++) {
                        C(k,i,l,j) += x_n[k] * ke[ij] * x_m[l] ;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  /* The solutions */
  {
    int k ;

    for(k = 0 ; k < dim ; k++) {
      int l ;
          
      for(l = k ; l < dim ; l++) {
        Solver_GetRHS(solver) = pb[3 * k + l] ;
        Solver_GetSolution(solver) = pu[3 * k + l] ;
        Solver_Solve(solver) ;
      }
    }
    
    Solver_GetRHS(solver) = b ;
    Solver_GetSolution(solver) = u ;
  }


  /* The macroscopic tangent stiffness matrix */
  {
    int i ;
        
    for(i = 0 ; i < dim ; i++) {
      int j ;
      
      for(j = 0 ; j < dim ; j++) {
        int k ;
          
        for(k = 0 ; k < dim ; k++) {
          int l ;
          
          for(l = 0 ; l < dim ; l++) {
            double  ub = 0. ;
            
            {
              double* uik = pu[3 * i + k] ;
              double* blj = pb[3 * l + j] ;
              int p ;
              
              for(p = 0 ; p < ncol ; p++) {
                ub += uik[p] * blj[p] ;
              }
            }
            
            {
              double  Eik =  E[3 * i + k] ;
              double  Elj =  E[3 * l + j] ;

              C(k,i,j,l) -= ub / (Eik * Elj) ;
            }
          }
        }
      }
    }
  }

  {
    double vol = FEM_ComputeVolume(mesh) ;
    int j ;
    
    for(j = 0 ; j < 81 ; j++) c[j] /= vol ;
  }
  
  return(0) ;
#undef C
}






int (FEM2_ComputeHomogenizedStressTensor1)(DataSet_t* dataset,Solver_t* solver,double t,double dt,Solutions_t* sols_n,Solutions_t* sols,double* dstrain,double* stress)
/** Compute the homogenized stress tensor of a microstructure.
 *  Inputs:
 *  - dataset:
 *  - solver:
 *  - t, dt:
 *  - sols_n, sols:
 *  - dstrain: strain increment as an array of 9 doubles
 *  Output:
 *  - stress: stress tensor as an array of 9 doubles.
 */
{
  /* Set input data of the microstructure */
  {
    /* Update the time function */
    {
      Functions_t* fcts = DataSet_GetFunctions(dataset) ;
      int nfcts = 9 ;
      int j ;
      
      for(j = 0 ; j < nfcts ; j++) {
        Function_t* func = Functions_GetFunction(fcts) + j ;
        double* x = Function_GetXValue(func) ;
        double* f = Function_GetFValue(func) ;

        x[0] = t - dt ;
        x[1] = t ;
        f[0] = 0 ;
        f[1] = dstrain[j] ;
      }
    }
  }
    
  /* Compute the microstructure */
  {
    Module_t* module = DataSet_GetModule(dataset) ;
    Dates_t*   dates  = DataSet_GetDates(dataset) ;
    Date_t*    date   = Dates_GetDate(dates) ;
    TimeStep_t*  timestep  = DataSet_GetTimeStep(dataset) ;
    double dtini = TimeStep_GetInitialTimeStep(timestep) ;
    double t_n = Solutions_GetTime(sols_n) ;
    
    {
      /* t should be equal to t_n + dt. */
      if(fabs(t - dt - t_n) > 1.e-4*dtini) {
        Message_FatalError("FEM2_ComputeStressTensor: t_n = %e ; t - dt = %e",t_n,t-dt) ;
      }
      
      /* There are 2 dates */
      Date_GetTime(date) = t_n ;
      Date_GetTime(date + 1) = t ;
    
      {
        Solution_t* sol   = Solutions_GetSolution(sols) ;
        Solution_t* sol_n = Solutions_GetSolution(sols_n) ;
        
        Solution_Copy(sol,sol_n) ;
      }
      
      {
        int i ;
        
        #if 1
        i = Module_SolveProblem(module,dataset,sols,solver,NULL) ;
        #else
        {
          Module_InitializeProblem(module,dataset,sols) ;
          i = Module_Increment(module,dataset,sols,solver,NULL,t_n,t) ;
        }
        #endif
        
        if(i < 0) {
          return(i) ;
          //Message_FatalError("FEM2_ComputeStressTensor: something went wrong") ;
          //Exception_Interrupt ;
        }
      }
    }
  }

  /* Backup stresses as averaged stresses */
  {
    Mesh_t* mesh = DataSet_GetMesh(dataset) ;
    
    FEM_AverageStresses(mesh,stress) ;
  }
  
  return(0) ;
}
#endif
