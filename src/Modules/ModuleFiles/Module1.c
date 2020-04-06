#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "Context.h"
#include "CommonModule.h"


#define AUTHORS  "Dangla"
#define TITLE    "Fully Coupled Approach"

#include "PredefinedMethods.h"

static Module_ComputeProblem_t   calcul ;
static Module_SolveProblem_t     Algorithm ;
static int    ComputeInitialState(Mesh_t*,double) ;
static int    ComputeExplicitTerms(Mesh_t*,double) ;
static int    ComputeMatrix(Mesh_t*,double,double,Matrix_t*) ;
static void   ComputeResidu(Mesh_t*,double,double,double*,Loads_t*) ;
static int    ComputeImplicitTerms(Mesh_t*,double,double) ;



/*
  Extern functions
*/

int SetModuleProp(Module_t* module)
{
  Module_CopyShortTitle(module,TITLE) ;
  Module_CopyNameOfAuthors(module,AUTHORS) ;
  Module_GetComputeProblem(module) = calcul ;
  Module_GetSolveProblem(module) = Algorithm ;
  return(0) ;
}

/*
  Intern functions
*/

int calcul(DataSet_t* jdd)
{
  Mesh_t* mesh = DataSet_GetMesh(jdd) ;
  const int n_sol = 2 ; /* Must be 2 at minimum but works with more */
  Solutions_t* sols = Solutions_Create(mesh,n_sol) ;

  /* Execute this line to set only one allocation of space for explicit terms. */
  /* This is not mandatory except in some models where constant terms are saved as 
   * explicit terms and updated only once during initialization. 
   * It is then necessary to merge explicit terms. Otherwise it is not mandatory.
   * Should be eliminated in the future. */
  Solutions_MergeExplicitTerms(sols) ;
  /* This is done 11/05/2015 */
  //Message_Warning("Explicit terms are not merged anymore in this version.") ;
  
  
  {
    DataFile_t* datafile = DataSet_GetDataFile(jdd) ;
    int i = 0 ;
  
  /* Set up the system of equations */
    {
      BConds_t* bconds = DataSet_GetBConds(jdd) ;
      
      //Mesh_SetMatrixRowColumnIndexes(mesh,bconds) ;
    }
    
  /* 1. Initial time */
    {
      Solution_t* sol = Solutions_GetSolution(sols) ;
      Dates_t*  dates = DataSet_GetDates(jdd) ;
      Date_t* date    = Dates_GetDate(dates) ;
      double t0       = Date_GetTime(date) ;
      
      Solution_GetTime(sol) = t0 ;
    }
    
  /* 2. Calculation */
    {
      char*   filename = DataFile_GetFileName(datafile) ;
      Dates_t*  dates    = DataSet_GetDates(jdd) ;
      int     nbofdates  = Dates_GetNbOfDates(dates) ;
      Points_t* points   = DataSet_GetPoints(jdd) ;
      int     n_points   = Points_GetNbOfPoints(points) ;
      Options_t* options = DataSet_GetOptions(jdd) ;
      Solver_t* solver = Solver_Create(mesh,options,1) ;
      OutputFiles_t* outputfiles = OutputFiles_Create(filename,nbofdates,n_points) ;
      
      i = Algorithm(jdd,sols,solver,outputfiles) ;
      
      Solver_Delete(&solver) ;
      OutputFiles_Delete(&outputfiles) ;
    }
      
  /* 3. Store for future resumption */
    {
      Solution_t* sol = Solutions_GetSolution(sols) ;
      double t =  Solution_GetTime(sol) ;
      
      Mesh_InitializeSolutionPointers(mesh,sols) ;
      Mesh_StoreCurrentSolution(mesh,datafile,t) ;
    }
    
    return(i) ;
  }
}


/*
  Intern functions
*/


static int   Algorithm(DataSet_t* jdd,Solutions_t* sols,Solver_t* solver,OutputFiles_t* outputfiles)
/** On input sols should point to the initial solution
 *  except if the context tells that initialization should be performed.
 *  On output sols points to the last converged solution.
 *  Return 0 if convergence has met, -1 otherwise.
 */
{
#define SOL_1     Solutions_GetSolution(sols)
#define SOL_n     Solution_GetPreviousSolution(SOL_1)

#define T_n       Solution_GetTime(SOL_n)
#define DT_n      Solution_GetTimeStep(SOL_n)
#define STEP_n    Solution_GetStepIndex(SOL_n)

#define T_1       Solution_GetTime(SOL_1)
#define DT_1      Solution_GetTimeStep(SOL_1)
#define STEP_1    Solution_GetStepIndex(SOL_1)

  DataFile_t*    datafile    = DataSet_GetDataFile(jdd) ;
  Options_t*     options     = DataSet_GetOptions(jdd) ;
  Mesh_t*        mesh        = DataSet_GetMesh(jdd) ;
  BConds_t*      bconds      = DataSet_GetBConds(jdd) ;
  Loads_t*       loads       = DataSet_GetLoads(jdd) ;
  Dates_t*       dates       = DataSet_GetDates(jdd) ;
  TimeStep_t*    timestep    = DataSet_GetTimeStep(jdd) ;
  IterProcess_t* iterprocess = DataSet_GetIterProcess(jdd) ;
  
  Nodes_t*       nodes       = Mesh_GetNodes(mesh) ;
  unsigned int   nbofdates   = Dates_GetNbOfDates(dates) ;
  Date_t*        date        = Dates_GetDate(dates) ;

  unsigned int   idate ;
  double t_0 ;
  
  
  /*
   * 1. Initialization
   */
  Mesh_InitializeSolutionPointers(mesh,sols) ;

  
  {
    int i = Mesh_LoadCurrentSolution(mesh,datafile,&T_1) ;
    
    idate = 0 ;
    
    if(i) {
      while(idate + 1 < nbofdates && T_1 >= Date_GetTime(date + idate + 1)) idate++ ;
      
      Message_Direct("Continuation ") ;
      
      if(DataFile_ContextIsFullInitialization(datafile)) {
        Message_Direct("(full initialization) ") ;
      } else if(DataFile_ContextIsPartialInitialization(datafile)) {
        Message_Direct("(partial initialization) ") ;
      } else if(DataFile_ContextIsNoInitialization(datafile)) {
        Message_Direct("(no initialization) ") ;
      }
      
      Message_Direct("at t = %e (between steps %d and %d)\n",T_1,idate,idate+1) ;
    }
    
    if(DataFile_ContextIsInitialization(datafile)) {
      IConds_t* iconds = DataSet_GetIConds(jdd) ;
    
      IConds_AssignInitialConditions(iconds,mesh,T_1) ;

      ComputeInitialState(mesh,T_1) ;
    }
  }
  
  
  /*
   * 2. Backup
   */
  t_0 = T_1 ;
  OutputFiles_BackupSolutionAtPoint(outputfiles,jdd,T_1,t_0) ;
  OutputFiles_BackupSolutionAtTime(outputfiles,jdd,T_1,idate) ;
  
  
  /*
   * 3. Loop on dates
   */
  for(; idate < nbofdates - 1 ; idate++) {
    Date_t* date_i = date + idate ;
    
    /*
     * 3.1 Loop on time steps
     */
    do {
      /*
       * 3.1.1 Looking for a new solution at t + dt
       * We step forward (point to the next solution) 
       */
      Solutions_StepForward(sols) ;
      Mesh_InitializeSolutionPointers(mesh,sols) ;
      
      /*
       * 3.1.1b Save the environment. 
       * That means that this is where the environment
       * is restored after a nonlocal jump.
       */
      Exception_SaveEnvironment ;
      
      /*
       * 3.1.1c Backup the previous solution:
       * if the saved environment was restored after a nonlocal jump
       * and 
       * if the exception mechanism orders to do it.
       */
      {
        if(Exception_OrderToBackupAndTerminate) {
          backupandreturn :
          Solutions_StepBackward(sols) ;
          Mesh_InitializeSolutionPointers(mesh,sols) ;
          OutputFiles_BackupSolutionAtTime(outputfiles,jdd,T_1,idate+1) ;
          return(-1) ;
        }
      }
      
      /*
       * 3.1.2 Compute the explicit terms with the previous solution
       */
      {
        int i = ComputeExplicitTerms(mesh,T_n) ;
        
        if(i != 0) {
          Message_Direct("\n") ;
          Message_Direct("Algorithm(1): undefined explicit terms\n") ;
          /* Backup the previous solution */
          if(T_n > t_0) {
            goto backupandreturn ;
          }
          return(-1) ;
        }
      }
        
      /*
       * 3.1.3 Compute and set the time step
       */
      {
        double t1 = Date_GetTime(date_i) ;
        double t2 = Date_GetTime(date_i + 1) ;
        double dt = TimeStep_ComputeTimeStep(timestep,nodes,T_n,DT_n,t1,t2) ;
        
        DT_1 = dt ;
        STEP_1 = STEP_n + 1 ;
      }
      
      /*
       * 3.1.3b Initialize the repetition index
       */
      IterProcess_GetRepetitionIndex(iterprocess) = 0 ;
      
      
      /*
       * 3.1.3c Reduce the time step 
       * if the exception mechanism orders to do it.
       */
      {
        if(Exception_OrderToReiterateWithSmallerTimeStep) {
          repeatwithreducedtimestep :
          
          IterProcess_IncrementRepetitionIndex(iterprocess) ;
          DT_1 *= TimeStep_GetReductionFactor(timestep) ;
          
        } else if(Exception_OrderToReiterateWithInitialTimeStep) {
          repeatwithinitialtimestep :
          
          IterProcess_IncrementRepetitionIndex(iterprocess) ;
          DT_1 *= TimeStep_GetReductionFactor(timestep) ;
          {
            double t_ini = TimeStep_GetInitialTimeStep(timestep) ;
              
            if(DT_1 > t_ini) DT_1 = t_ini ;
          }
        }
      }

      
      /*
       * 3.1.3d The time at which we compute
       */
      {
        int irecom = IterProcess_GetRepetitionIndex(iterprocess) ;
        
        if(irecom > 0) Message_Direct("Repetition no %d\n",irecom) ;
      }
      T_1 = T_n + DT_1 ;
      Message_Direct("Step %d  t = %e (dt = %4.2e)",STEP_1,T_1,DT_1) ;
      
      /*
       * 3.1.4 Initialize the unknowns
       */
      Mesh_SetCurrentUnknownsWithBoundaryConditions(mesh,bconds,T_1) ;
      
      /*
       * 3.1.5 Loop on iterations
       */
      IterProcess_GetIterationIndex(iterprocess) = 0 ;
      while(IterProcess_LastIterationIsNotReached(iterprocess)) {
        IterProcess_IncrementIterationIndex(iterprocess) ;
        
        /*
         * 3.1.5.1 The implicit terms (constitutive equations)
         */
        {
          int i = ComputeImplicitTerms(mesh,T_1,DT_1) ;
          
          if(i != 0) {
            if(IterProcess_LastRepetitionIsNotReached(iterprocess)) {
              goto repeatwithinitialtimestep ;
            } else {
              int iter = IterProcess_GetIterationIndex(iterprocess) ;
              
              Message_Direct("\n") ;
              Message_Direct("Algorithm(2): undefined implicit terms at iteration %d\n",iter) ;
              goto backupandreturn ;
            }
          }
        }
        
        /*
         * 3.1.5.2 The residu
         */
        {
          double*  rhs = Solver_GetRHS(solver) ;
          
          ComputeResidu(mesh,T_1,DT_1,rhs,loads) ;
          
          {
            char*  debug = Options_GetPrintedInfos(options) ;
            
            if(!strcmp(debug,"residu")) {
              Solver_Print(solver,debug) ;
            }
          }
        }
        
        /*
         * 3.1.5.3 The matrix
         */
        {
          Matrix_t*  a = Solver_GetMatrix(solver) ;
          int i = ComputeMatrix(mesh,T_1,DT_1,a) ;
          
          if(i != 0) {
            if(IterProcess_LastRepetitionIsNotReached(iterprocess)) {
              goto repeatwithinitialtimestep ;
            } else {
              int iter = IterProcess_GetIterationIndex(iterprocess) ;
              
              Message_Direct("\n") ;
              Message_Direct("Algorithm(3): undefined matrix at iteration %d\n",iter) ;
              goto backupandreturn ;
            }
          }
          
          {
            char*  debug = Options_GetPrintedInfos(options) ;
            
            if(!strncmp(debug,"matrix",4)) {
              Solver_Print(solver,debug) ;
            }
          }
        }
        
        /*
         * 3.1.5.4 Resolution
         */
        {
          int i = Solver_Solve(solver) ;
          
          if(i != 0) {
            if(IterProcess_LastRepetitionIsNotReached(iterprocess)) {
              goto repeatwithinitialtimestep ;
            } else {
              int iter = IterProcess_GetIterationIndex(iterprocess) ;
              
              Message_Direct("\n") ;
              Message_Direct("Algorithm(4): unable to solve at iteration %d\n",iter) ;
              goto backupandreturn ;
            }
          }
        }
        
        /*
         * 3.1.5.5 Update the unknowns
         */
        Mesh_UpdateCurrentUnknowns(mesh,solver) ;
        
        /*
         * 3.1.5.6 The error
         */
        {
          int i = IterProcess_SetCurrentError(iterprocess,nodes,solver) ;
          
          if(i != 0) {
            if(IterProcess_LastRepetitionIsNotReached(iterprocess)) {
              goto repeatwithinitialtimestep ;
            } else {
              int iter = IterProcess_GetIterationIndex(iterprocess) ;
              
              Message_Direct("\n") ;
              Message_Direct("Algorithm(5): unable to compute error at iteration %d\n",iter) ;
              goto backupandreturn ;
            }
          }
        }
        
        /*
         * 3.1.5.7 We get out if convergence is met
         */
        if(IterProcess_ConvergenceIsMet(iterprocess)) break ;
        
        {
          if(Options_IsToPrintOutAtEachIteration(options)) {
            if(IterProcess_LastIterationIsNotReached(iterprocess)) {
              IterProcess_PrintCurrentError(iterprocess) ;
            }
          }
        }
      }
      
      {
        IterProcess_PrintCurrentError(iterprocess) ;
      }
      
      /*
       * 3.1.6 Back to 3.1.3 with a smaller time step
       */
      if(IterProcess_ConvergenceIsNotMet(iterprocess)) {
        if(IterProcess_LastRepetitionIsNotReached(iterprocess)) {
          goto repeatwithreducedtimestep ;
        }
      }
      
      /*
       * 3.1.7 Backup for specific points
       */
      OutputFiles_BackupSolutionAtPoint(outputfiles,jdd,T_1,t_0) ;
      /*
       * 3.1.8 Go to 3.2 if convergence was not met
       */
      if(IterProcess_ConvergenceIsNotMet(iterprocess)) break ;
    } while(T_1 < Date_GetTime(date_i + 1)) ;
    
    /*
     * 3.2 Backup for this time
     */
    OutputFiles_BackupSolutionAtTime(outputfiles,jdd,T_1,idate+1) ;
    
    /*
     * 3.3 Go to 4. if convergence was not met
     */
    if(IterProcess_ConvergenceIsNotMet(iterprocess)) break ;
  }
  
  /*
   * 4. Step backward if convergence was not met
   */
  if(IterProcess_ConvergenceIsNotMet(iterprocess)) {
    Solutions_StepBackward(sols) ;
    return(-1) ;
  }
  
  return(0) ;

#undef T_n
#undef DT_n
#undef STEP_n
#undef T_1
#undef DT_1
#undef STEP_1
#undef SOL_n
#undef SOL_1
}


int ComputeInitialState(Mesh_t* mesh,double t)
{
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int    ie ;

  for(ie = 0 ; ie < n_el ; ie++) {
    Material_t* mat = Element_GetMaterial(el + ie) ;
    
    if(mat) {
      int i ;
      
      Element_FreeBuffer(el + ie) ;
      i = Element_ComputeInitialState(el + ie,t) ;
      if(i != 0) return(i) ;
    }
  }
  
  return(0) ;
}


int ComputeExplicitTerms(Mesh_t* mesh,double t)
{
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int    ie ;

  for(ie = 0 ; ie < n_el ; ie++) {
    Material_t* mat = Element_GetMaterial(el + ie) ;
    
    if(mat) {
      int    i ;
      
      Element_FreeBuffer(el + ie) ;
      i = Element_ComputeExplicitTerms(el + ie,t) ;
      if(i != 0) return(i) ;
    }
  }
  
  return(0) ;
}


int ComputeMatrix(Mesh_t* mesh,double t,double dt,Matrix_t* a)
{
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int    ie ;

  Matrix_SetValuesToZero(a) ;
  
  for(ie = 0 ; ie < n_el ; ie++) {
    Material_t* mat = Element_GetMaterial(el + ie) ;
    
    if(mat) {
#define NE (Element_MaxNbOfNodes*Model_MaxNbOfEquations)
      double ke[NE*NE] ;
#undef NE
      int    i ;
      
      Element_FreeBuffer(el + ie) ;
      i = Element_ComputeMatrix(el + ie,t,dt,ke) ;
      if(i != 0) return(i) ;
      
      Matrix_AssembleElementMatrix(a,el+ie,ke) ;
    }
  }
  
  return(0) ;
}




void ComputeResidu(Mesh_t* mesh,double t,double dt,double* r,Loads_t* loads)
{
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int n_cg = Loads_GetNbOfLoads(loads) ;
  Load_t* cg = Loads_GetLoad(loads) ;
  unsigned int    ie,i_cg ;
#define NE (Element_MaxNbOfNodes*Model_MaxNbOfEquations)
  double re[NE] ;
#undef NE
  double zero = 0. ;
  
  {
    unsigned int    n_col = Mesh_GetNbOfMatrixColumns(mesh) ;
    unsigned int    j ;
    
    for(j = 0 ; j < n_col ; j++) r[j] = zero ;
  }
  
  /* Residu */
  for(ie = 0 ; ie < n_el ; ie++) {
    int  nn = Element_GetNbOfNodes(el + ie) ;
    Material_t* mat = Element_GetMaterial(el + ie) ;
    
    if(mat) {
      int    neq = Material_GetNbOfEquations(mat) ;
      int i ;
      
      Element_FreeBuffer(el + ie) ;
      Element_ComputeResidu(el + ie,t,dt,re) ;
      
      /* assembling */
      for(i = 0 ; i < nn ; i++) {
        Node_t* node_i = Element_GetNode(el + ie,i) ;
        int    j ;
        for(j = 0 ; j < neq ; j++) {
          int ij = i*neq + j ;
          int jj = Element_GetUnknownPosition(el + ie)[ij] ;
          if(jj >= 0) {
            int k = Node_GetMatrixColumnIndex(node_i)[jj] ;
            if(k >= 0) r[k] += re[ij] ;
          }
        }
      }
    }
  }
  
  /* Loads */
  for(i_cg = 0 ; i_cg < n_cg ; i_cg++) {
    int reg_cg = Load_GetRegionIndex(cg + i_cg) ;
    for(ie = 0 ; ie < n_el ; ie++) if(Element_GetRegionIndex(el + ie) == reg_cg) {
      int  nn = Element_GetNbOfNodes(el + ie) ;
      Material_t* mat = Element_GetMaterial(el + ie) ;
    
      if(mat) {
        int    neq = Material_GetNbOfEquations(mat) ;
        int i ;
        
        Element_FreeBuffer(el + ie) ;
        Element_ComputeLoads(el + ie,t,dt,cg + i_cg,re) ;
        
        /* assembling */
        for(i = 0 ; i < nn ; i++) {
          Node_t* node_i = Element_GetNode(el + ie,i) ;
          int    j ;
          for(j = 0 ; j < neq ; j++) {
            int ij = i*neq + j ;
            int jj = Element_GetUnknownPosition(el + ie)[ij] ;
            if(jj >= 0) {
              int k = Node_GetMatrixColumnIndex(node_i)[jj] ;
              if(k >= 0) r[k] += re[ij] ;
            }
          }
        }
      }
    }
  }
}


int ComputeImplicitTerms(Mesh_t* mesh,double t,double dt)
{
  unsigned int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  unsigned int    ie ;

  for(ie = 0 ; ie < n_el ; ie++) {
    Material_t* mat = Element_GetMaterial(el + ie) ;
    
    if(mat) {
      int    i ;
      
      Element_FreeBuffer(el + ie) ;
      i = Element_ComputeImplicitTerms(el + ie,t,dt) ;
      if(i != 0) return(i) ;
    }
  }
  
  return(0) ;
}
