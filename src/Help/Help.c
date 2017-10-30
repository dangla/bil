#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Models.h"
#include "Help.h"


static void   Help_Geometry(void) ;
static void   Help_Mesh(void) ;
static void   Help_Material(void) ;
static void   Help_Fields(void) ;
static void   Help_Initialization(void) ;
static void   Help_Functions(void) ;
static void   Help_BoundaryConditions(void) ;
static void   Help_Loads(void) ;
static void   Help_Points(void) ;
static void   Help_Dates(void) ;
static void   Help_ObjectiveValues(void) ;
static void   Help_IterativeProcess(void) ;
static void   Help_TimeStep(void) ;
static void   Help_BuiltinCurves(void) ;


/* Extern functions */

void Help_HelpOnline(void)
{
  char mot[5] ;
  
  debut:
  
  Message_Direct("\n") ;
  Message_Direct("Help topics available:") ;
  Message_Direct("\n") ;
  
  
/*01234567890123456789012345678901234567890123456789012345678901234567890*/
  Message_Direct("Mandatory fields\n") ;
  Message_Direct("\
  Geometry         Mesh             Material              Fields\n") ;
  Message_Direct("\
  Initialization   Functions        Boundary Conditions   Loads\n") ;
  Message_Direct("\
  Points           Dates            Objective Variations  Iterative Process\n") ;
  Message_Direct("\
  Time Steps\n") ;
  Message_Direct("\n") ;
  Message_Direct("\n") ;
  
  
  Message_Direct("Optional fields\n") ;
  Message_Direct("\
  Periodicities         Units\n") ;
  Message_Direct("\n") ;
  Message_Direct("\n") ;
  
  
  Message_Direct("Other topics\n") ;
  Message_Direct("\
  Builtin Curves\n") ;
  Message_Direct("\n") ;
  Message_Direct("\n") ;
  
  
  Message_Direct("Help topic (q=quit): ") ;
  fscanf(stdin,"%s",mot) ;
  
  if(!strcmp(mot,"q")) return ;
  
  else if(!strncasecmp(mot,"Geometry"            ,2)) Help_Geometry() ;
  else if(!strncasecmp(mot,"Mesh"                ,2)) Help_Mesh() ;
  else if(!strncasecmp(mot,"Material"            ,2)) Help_Material() ;
  else if(!strncasecmp(mot,"Fields"              ,2)) Help_Fields() ;
  else if(!strncasecmp(mot,"Initialization"      ,2)) Help_Initialization() ;
  else if(!strncasecmp(mot,"Functions"           ,2)) Help_Functions() ;
  else if(!strncasecmp(mot,"Boundary Conditions" ,2)) Help_BoundaryConditions() ;
  else if(!strncasecmp(mot,"Loads"               ,2)) Help_Loads() ;
  else if(!strncasecmp(mot,"Points"              ,2)) Help_Points() ;
  else if(!strncasecmp(mot,"Dates"               ,2)) Help_Dates() ;
  else if(!strncasecmp(mot,"Objective Variations",2)) Help_ObjectiveValues() ;
  else if(!strncasecmp(mot,"Iterative Process"   ,2)) Help_IterativeProcess() ;
  else if(!strncasecmp(mot,"Time Steps"          ,2)) Help_TimeStep() ;
  else if(!strncasecmp(mot,"Builtin Curves"      ,2)) Help_BuiltinCurves() ;
  else {
    Message_Direct("Topic %s not found\n",mot) ;
  }
  
  goto debut ;
}


void Help_WriteData(char *nom)
/* Creation du fichier de donnees */
{
  int    n_mat = 1,n_el = 0,dim ;
  char   mot[100] ;
  int    i,n ;
  FILE   *ficd ;
  
  Message_Direct("Creation of %s and writing the input data\n",nom) ;
  
  ficd = fopen(nom,"w") ;
  if(!ficd) {
    arret("error in opening file %s",nom) ;
  }

  Message_Direct("1. Geometry:\n") ;
  fprintf(ficd,"Geometry\n") ;
  {
    Message_Direct("dimension: ") ;
    fscanf(stdin,"%d",&dim) ; fprintf(ficd,"%d",dim) ;
    if(dim < 3) {
      Message_Direct("geometry (plan/axis/sphe): ") ;
      fscanf(stdin,"%s",mot) ; fprintf(ficd," %s",mot) ;
    }
    fprintf(ficd,"\n") ;
  }

  Message_Direct("2. Mesh:\n") ;
  fprintf(ficd,"Mesh\n") ;
  if(dim == 1) {
    int    j ;
    double dx,pt ;
    Message_Direct("number of points (nb of regions + 1): ") ;
    fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d",n) ;
    Message_Direct("point coordinates: ") ;
    for(i=0;i<n;i++) {
      fscanf(stdin,"%lf",&pt) ; fprintf(ficd," %g",pt) ;
    }
    fprintf(ficd,"\n") ;
    Message_Direct("length of the first element starting from the first node: ") ;
    fscanf(stdin,"%lf",&dx) ; fprintf(ficd,"%g\n",dx) ;
    Message_Direct("number of elements per region: ") ;
    n_el = 0 ;
    for(i=0;i<n-1;i++) {
      fscanf(stdin,"%d",&j) ;
      fprintf(ficd,"%d ",j) ;
      n_el += j ;
    }
    fprintf(ficd,"\n") ;
    n_mat = 0 ;
    for(i=0;i<n-1;i++) {
      Message_Direct("material index of region %d: ",i+1) ;
      fscanf(stdin,"%d",&j) ;
      fprintf(ficd,"%d ",j) ;
      if(j > n_mat) n_mat = j ;
    }
    fprintf(ficd,"\n") ;
  } else {
    char   nom1[100] ;
    Message_Direct("file name: ") ;
    fscanf(stdin,"%s",nom1) ; fprintf(ficd,"%s\n",nom1) ;
  }

  for(i=0;i<n_mat;i++) {
    Message_Direct("3.%d Material %d:\n",i+1,i+1) ;
    fprintf(ficd,"%s\n","Material") ;
    Message_Direct("models ? ") ;
    fscanf(stdin,"%s",mot) ;
    fprintf(ficd,"Model = %s\n",mot) ;
    Models_Print(mot,ficd) ;
  }

  Message_Direct("4. Fields:\n") ;
  fprintf(ficd,"Fields\n") ;
  Message_Direct("A field v(x) is expressed as:\n\
 \tv(x) = v_o + g_o*(x - x_o)\n\
 where\n\
 x_o is the coordinate vector of a point,\n\
 v_o is the value v(x_o),\n\
 g_o is the gradient of v.\n") ;
  Message_Direct("number of fields to be defined (0 = no one): ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    int    j ;
    double v,g[3],z[3] ;
    Message_Direct("value: ") ;
    fscanf(stdin,"Value = %lf",&v) ; fprintf(ficd,"Region = %g",v) ;
    Message_Direct("gradient g(i = 1,%d): ",dim) ;
    for(j=0;j<dim;j++) fscanf(stdin,"%lf",g+j) ;
    fprintf(ficd,"Gradient =") ;
    for(j=0;j<dim;j++) fprintf(ficd," %g",g[j]) ;
    Message_Direct("point x(i = 1,%d): ",dim) ;
    for(j=0;j<dim;j++) fscanf(stdin,"%lf",z+j) ;
    fprintf(ficd,"Point =") ;
    for(j=0;j<dim;j++) fprintf(ficd," %g",z[j]) ;
    fprintf(ficd,"\n") ;
  }

  Message_Direct("5. Initial Conditions:\n") ;
  fprintf(ficd,"Initial Conditions\n") ;
  Message_Direct("number of initialisation (0  = no one): ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    int    j ;
    Message_Direct("number of regions: ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd,"Reg = %d",j) ;
    Message_Direct("name of the parameter to be initialzed: ") ;
    fscanf(stdin,"%s",mot) ; fprintf(ficd," Unknown = %s",mot) ;
    Message_Direct("field index: ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd," Field = %d",j) ;
    fprintf(ficd,"\n") ;
  }

  Message_Direct("6. Functions (of time):\n") ;
  fprintf(ficd,"Functions\n") ;
  Message_Direct("number of functions to be defined (0 = no one): ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    int    j,k ;
    double t,f ;
    Message_Direct("number of points for the %dth function: ",i+1) ;
    fscanf(stdin,"%d",&k) ; fprintf(ficd,"%d\n",k) ;
    for(j=0;j<k;j++) {
      Message_Direct("%d give (t,f): ",j+1) ;
      fscanf(stdin,"%lf %lf",&t,&f) ; fprintf(ficd,"F(%g) = %g\n",t,f) ;
    }
  }

  Message_Direct("7. Boundary Conditions:\n") ;
  fprintf(ficd,"Boundary Conditions\n") ;
  Message_Direct("number of boundary conditions (0 = no one): ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    int    j ;
    Message_Direct("region index: ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd,"Region = %d",j) ;
    Message_Direct("name of the parameter to be prescribed: ") ;
    fscanf(stdin,"%s",mot) ; fprintf(ficd," Unknown = %s",mot) ;
    Message_Direct("field index (0 = default null field): ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd," Field = %d",j) ;
    Message_Direct("function index \"f(t)\" (0 = unit function \"f(t)=1\"): ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd," Function = %d",j) ;
    fprintf(ficd,"\n") ;
  }

  Message_Direct("8. Loads:\n") ;
  fprintf(ficd,"Loads\n") ;
  Message_Direct("number of loads (0 = no one): ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    int    j ;
    Message_Direct("region index: ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd,"Region = %d",j) ;
    Message_Direct("name of de equation which the load is applied to: ") ;
    fscanf(stdin,"%s",mot) ; fprintf(ficd," Equation = %s",mot) ;
    Message_Direct("type of load (force, flux, pressure etc...): ") ;
    fscanf(stdin,"%s",mot) ; fprintf(ficd," Type = %s",mot) ;
    Message_Direct("field index (0 = default null field): ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd," Field = %d",j) ;
    Message_Direct("function index \"f(t)\" (0 = unit function \"f(t)=1\"): ") ;
    fscanf(stdin,"%d",&j) ; fprintf(ficd," Function = %d",j) ;
    fprintf(ficd,"\n") ;
  }

  Message_Direct("9. Points (outputs at some points):\n") ;
  fprintf(ficd,"Points\n") ;
  Message_Direct("number of points (0 = no one): ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    int    j ;
    double x ;
    Message_Direct("coordinates of the %dth point: ",i+1) ;
    for(j=0;j<dim;j++) {fscanf(stdin,"%lf",&x) ; fprintf(ficd,"%g ",x) ;}
    fprintf(ficd,"\n") ;
  }

  Message_Direct("10. Dates:\n") ;
  fprintf(ficd,"Dates\n") ;
  Message_Direct("number of dates + 1: ") ;
  fscanf(stdin,"%d",&n) ; fprintf(ficd,"%d\n",n) ;
  for(i=0;i<n;i++) {
    double t ;
    Message_Direct("time %d: ",i) ;
    fscanf(stdin,"%lf",&t) ; fprintf(ficd,"%g\n",t) ;
  }

  Message_Direct("11. Objective Variations:\n") ;
  fprintf(ficd,"Objective Variations\n") ;
  Message_Direct("number of objective values: ") ;
  fscanf(stdin,"%d",&n) ;
  for(i=0;i<n;i++) {
    double v ;
    Message_Direct("name of the parameter %d: ",i+1) ;
    fscanf(stdin,"%s",mot) ; fprintf(ficd,"%s = ",mot) ;
    Message_Direct("objective variation of %s: ",mot) ;
    fscanf(stdin,"%lf",&v) ; fprintf(ficd,"%g\n",v) ;
  }

  Message_Direct("12. Iterative Process:\n") ;
  fprintf(ficd,"Iterative Process\n") ;
  {
    double v ;
    Message_Direct("max number of iterations: ") ;
    fscanf(stdin,"%d",&n) ; fprintf(ficd,"Iteration = %d ",n) ;
    Message_Direct("tolerance: ") ;
    fscanf(stdin,"%lf",&v) ; fprintf(ficd,"Tolerance = %g ",v) ;
    Message_Direct("max number of repetitions: ") ;
    fscanf(stdin,"%d",&n) ; fprintf(ficd,"Repetition = %d ",n) ;
    fprintf(ficd,"\n") ;
  }

  Message_Direct("13. Time steps:\n") ;
  fprintf(ficd,"Time Steps\n") ;
  {
    double v ;
    Message_Direct("initial time increment: ") ;
    fscanf(stdin,"%lf",&v) ; fprintf(ficd,"Dtini = %g ",v) ;
    Message_Direct("max time increment: ") ;
    fscanf(stdin,"%lf",&v) ; fprintf(ficd,"Dtmax = %g\n",v) ;
    fprintf(ficd,"\n") ;
  }

  /* END */
  Message_Direct("End\n") ;
  fclose(ficd) ;
}

/* Intern functions */

void Help_Geometry(void)
{
  Message_Direct("\n") ;
  Message_Direct("Geometry: ") ;
  Message_Direct("dimension and symmetry\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t dim\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
  \t dim = dimension (1,2,3)\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
  \t For dim = 1,2:\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
->\t sym\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  \t sym = symmetry (\"plan\",\"axis\",\"sphe\")\n") ;

  Message_Direct("\n") ;
  Message_Direct("Example: ") ;
  Message_Direct("\n") ;
  Message_Direct("Geometry\n") ;
  Message_Direct("2 plan\n") ;
}


void Help_Mesh(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of the mesh\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t filename\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  \t filename      = name of the mesh file\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  \t Accepted filenames are:\n") ;

  Message_Direct("\
  \t *.msh \t(GMSH format);\n\
  \t *.ces \t(CESAR format);\n\
  \t *.m1d \t(only for 1D problem);\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  \t File format \"*.m1d\" is given as follows:\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t np\n\
->\t (pt,lc)[]\n\
->\t mat[]\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  \t np                   = nb of points. There are np-1 regions.\n\
  \t pt[]  (i = 1...np)   = x coordinate of points.\n\
  \t lc[]  (i = 1...np)   = mesh element size at the points.\n\
  \t mat[] (i = 1...np-1) = index number of the materials in the regions.\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\n\
  Instead of \"filename\" and for 1D problem only, we can give the following\n\
  input data directly after the key-word MESH:\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t np\n\
->\t pt[]\n\
->\t dxini\n\
->\t ne[]\n\
->\t mat[]\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  \t np                   = nb of points\n\
  \t pt[]  (i = 1...np)   = x coordinate of points.\n\
  \t dxini                = length of the first element attached to the first point\n\
  \t ne[]  (i = 1...np-1) = nb of elements in the regions\n\
  \t mat[] (i = 1...np-1) = index number of the materials in the regions.\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  NB: in 1D problem, to define a 1 node element give the same x-coordinate \
  to 2 points. These 2 points will generate only one node and one element.\n") ;

  Message_Direct("\n") ;
  Message_Direct("Example: ") ;
  Message_Direct("\n") ;
  Message_Direct("Mesh\n") ;
  Message_Direct("mymesh.msh\n") ;
}


void Help_Material(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of the material properties\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->(opt)\t \"Method = \",method\n\
->     \t \"Model  = \",model\n") ;

  Message_Direct("\
  \t method        = string defining an existing method\n\
  \t model         = code name of the model (basename of its file)\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  Enter the material properties.\n") ;
  
  Message_Direct("\
  Give as many times as there are material properties.\n") ;

  Message_Direct("\
  for(i = 0 , i < Nb of Properties , i++) {\n") ;

  Message_Direct("\
->\t key,\"=\",v\n") ;

  Message_Direct("\
  }\n") ;

  Message_Direct("\
  \t key           = string for the name of the property\n\
  \t v             = value or string assigned to the property\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  To display the available models, type\n") ;
  
  Message_Direct("\
  \t \"bil -m\"\n") ;

  Message_Direct("\
  To display an example for the input data of the model I, type\n") ;

  Message_Direct("\
  \t \"bil -m\" I\n") ;

  Message_Direct("Example:") ;
  Message_Direct("\n") ;
  Message_Direct("Material\n") ;
  Message_Direct("Model = XXX\n") ;
  Message_Direct("porosity = 0.2\n") ;
  Message_Direct("Curves  = my_file\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  In this example \"Curves\" is a reserved key-word.\n\
  \"my_file\" is the name of a file containing some curves\n\
  given as a series of columns.\n\
  The first column is the abcissa sampled with a constant step\n\
  (only the first and the last values are recorded).\n\
  The other columns are the values of curves at the specified abcissa.\n\
  The number of lines (points) is determined automatically.\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
  If \"my_file\" is not an existing file, Bil keeps on reading in the line\n\
  so as to find the informations to build the file as follows:\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  ...my_file p_c = Range{x1 = 0 , x2 = 1.e8 , n = 1001}\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  \"p_c\" is the name given to the abcissa. \"Range\" is a reserved key-word.\n\
  The curves are sampled in the range 0:1.e8 with a step equal to 1.e5\n\
  (i.e. (x2-x1)/(n-1)). The reading continues as follows:\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\n\
  ...n = 1001} S_l = Van-Genuchten(1){a = 1.e6 , m = 0.6}\n") ;

  Message_Direct("\
  \"S_l\" is the name given to the first curve which is calculated with the\n\
  model called Van-Genuchten and the given input data. The integer 1 between\n\
  parenthesis indicates that the curve is calculated by taking the first column\n\
  as input.\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  The reading can continue in the same line to build another curve as\n\
  follows:\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  ...m = 0.6} tau_l = Pow(2){lam = 0.4}\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  This new curve called \"tau_l\" is a power law of \"S_l\" (column 2 as given\n\
  in parenthesis). We can continue as the following general rule:\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  ...lam = 0.4} Y = Model(i){...}\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  \"Y\" is the name of the curve built with column \"i\" as input and with\n\
  the model \"Model\" and the parameters between braces.\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  See the \"BuiltinCurves\" topic to display all the built-in curves.\n") ;
}


void Help_Fields(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of fields\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t nch\n\
  \t nch = number of fields\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  Enter as many times as nch:\n") ;

  Message_Direct("\
  for(i = 0 ; i < nfn ; i++) {\n\
->(opt)\t \"Type =\",type\n\
->\t Enter the input data for the type (see below) \n\
  }\n\
  \t type          = string displaying the type of field\n\
  \t                 (affine,grid,etc..)\n\
  ") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  By default \"Type = affine\" if not given\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
  If \"Type = affine\": uniform gradient field\n") ;

  Message_Direct("\
->\t \"Value =\",v,\"Gradient =\",g[],\"Point =\",z[]\n") ;

  Message_Direct("\
  \t v                  = value of prameter at point A\n\
  \t g[i] (i = 1...dim) = parameter gradient vector at point A\n\
  \t z[i] (i = 1...dim) = coordinates of point A\n") ;

  
  Message_Direct("\n") ;

  Message_Direct("\
  If \"Type = grid\": field defined in a grid\n") ;

  Message_Direct("\
->\"File =\",name\n") ;

  Message_Direct("\
  \t name          = name of the file where the following\n\
  \t                 data are recorded\n") ;

  Message_Direct("\
  if dim = 1 (n_y = n_z = 1)\n\
->\t n_x\n\
->\t x[]\n\
  if dim = 2 (n_z = 1)\n\
->\t n_x,n_y\n\
->\t x[],y[]\n\
  if dim = 3\n\
->\t n_x,n_y,n_z\n\
->\t x[],y[],z[]\n\
  \t n_x,n_y,n_z   = number of points along the x-axis, y-axis and z-axis\n\
  \t                 respectively\n") ;

  Message_Direct("\
  \t x[i] (i = 1...n_x) = point coordinates along the x-axis\n\
  \t y[i] (i = 1...n_y) = point coordinates along the y-axis\n\
  \t z[i] (i = 1...n_z) = point coordinates along the z-axis\n") ;

  Message_Direct("\
->\t v[]\n\
  \t v[i] (i = 1...N) = values at the nodes of the grid\n\
  \t                    (N = n_x*n_y*n_z), by giving\n\
  \t                    first those of the x-axis, then\n\
  \t                    y-axis and z-axis\n") ;
}


void Help_Initialization(void)
{
  Message_Direct("\n") ;
  Message_Direct("Initialization of the primary unknowns\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t ninit\n\
  \t ninit = number of initialization\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  If ninit = 0, by default the values are equal to zero\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  If ninit > 0, give as many times as ninit:\n\
  for(i = 0 ; i < ninit ; i++) {\n\
->\t \"Region =\",reg,\"Unknown =\",cle,\"Field =\",ich\n\
  \t or\n\
->\t \"Region =\",reg,\"Unknown =\",cle,\"File =\",filename\n\
  }\n") ;

  Message_Direct("\
  \t reg           = index of the region\n\
  \t cle           = string displaying the parameter\n\
  \t ich           = index of the field\n\
  \t filename      = file name where the values of the primary unknowns\n\
  \t                 at the node of the mesh are stored\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
  If ninit < 0, reading in file\n\
->\t name\n\
  \t name          = file name where the whole primary unknowns\n\
  \t                 at the nodes of the mesh are stored\n") ;
}


void Help_Functions(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of time functions\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t nfn\n\
  \t nfn  = number of functions\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  Give as many times as nfn \n\
  for(i = 0 ; i < nfn ; i++) {\n\
->\t \"N =\",np\n\
  \t np       = number of points for the function\n\
->\t Enter the input data for function (see below)\n\
  \t or\n\
->\t \"File =\",name\n\
  \t name     = name of the file where the function is defined.\n\
  }\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  If \"N\" was given:\n") ;
  
  Message_Direct("\
    give as many values as np\n\
    for(j = 0 ; j < np ; j++) {\n\
->\t \"F(\",t,\") =\",f\n\
    } \n\
  } \n\
  \t t        = time\n\
  \t f        = value f(t)\n\
  ") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  If \"File\" was given:\n") ;
  
  Message_Direct("\
    The file \"name\" should contain 2 (or more) space-separated columns:\n\
    t  F_1(t)  F_2(t)  ...\n\
    if more than 2 columns are given the number of functions is incremented accordingly\n\
    ") ;
    
  Message_Direct("\n") ;
}

void Help_BoundaryConditions(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of boundary conditions\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t ncl\n\
  \t ncl = number of boundary conditions\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  Give as many times as ncl\n\
  for(i = 0 ; i < ncl ; i++) {\n\
->\t \"Region =\",reg,\"Unknown  =\",cle,\n\
  \t \"Field  =\",ich,\"Function =\",ifn \n\
  }\n\
  ") ;
  
  Message_Direct("\
  \t reg           = index of the region\n") ;

  Message_Direct("\
  \t cle           = string of characters representing the unknown\n") ;

  Message_Direct("\
  \t ich           = index of the field connected to this B.C.\n") ;

  Message_Direct("\
  \t ifn           = index of the function f(t) connected to this B.C.\n\
  The boundary condition at nodes is then given by f(t)*field(x)\n\
  (by default f(t) = 1)\n") ;

  Message_Direct("\n") ;
  Message_Direct("Example: ") ;
  Message_Direct("\n") ;
  Message_Direct("Boundary Conditions\n") ;
  Message_Direct("1\n") ;
  Message_Direct("Region = 100  ") ;
  Message_Direct("Unknown = p_co2  ") ;
  Message_Direct("Field = 2  ") ;
  Message_Direct("Function = 1  ") ;
  Message_Direct("\n") ;
}


void Help_Loads(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of loads\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t nch\n\
  \t nch  = Nb of loads\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  Give as many times as nch\n\
  for(i = 0 ; i < nch ; i++) {\n\
->\t \"Region =\",reg,\"Equation =\",cle,\n\
->\t \"Type   =\",typ,\"Field    =\",ich,\"Function =\",ifn\n\
  }\n\
  ") ;
  
  Message_Direct("\
  \t reg           = Index of the region\n") ;

  Message_Direct("\
  \t cle           = String of characters representing the equation\n") ;

  Message_Direct("\
  \t typ           = String of characters representing the loading type\n\
  \t                 (pression,flux,etc...)\n") ;

  Message_Direct("\
  \t ich           = index of the field connected to this load\n") ;

  Message_Direct("\
  \t ifn           = index of the function f(t) connected to this load\n") ;
  
  Message_Direct("\n") ;

  Message_Direct("\
  The load at the boundary is then given by:\n\
  \t - for a type \"flux\": df/dt(t)*field(x) (by default f(t) = t)\n\
  \t - for a type \"pressure\": f(t)*field(x) (by default f(t) = 1)\n") ;
}


void Help_Points(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of points\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t np\n\
->\t z[][]\n\
  \t np      = Nb of points\n\
  \t z[i][j] = Coordinate j (j = 1 , dim) of point i (i = 1 , np)\n") ;
}


void Help_Dates(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of dates\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t n,t[]\n\
  \t n   = Nb of dates\n\
  \t t[] = The dates\n") ;
}


void Help_ObjectiveValues(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of objective variations of primary unknowns.\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  Give as many times as the number of types of primary unknowns:\n\
  for(i = 0 ; i < Nb of unknowns ; i++) {\n\
->\t key,\"=\",V_obj\n\
  }\n\
  \t key                 = string representing the unknown \n\
  \t V_obj               = objective variation of the unknown \n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  The resolution of the set of non linear algebraic equations is performed\n\
  by an iteration procedure (e.g. Newton's method) until a convergence \n\
  criterion is met under the form\n\
  \n\
  (U(n+1) - U(n))/V_obj < small number \n\
  ") ;
}


void Help_IterativeProcess(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of convergence criteria\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t \"Iterations =\",niter,\"Tolerance =\",tol,\"Repetition =\",nrep\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
\t niter               = maximum number of iterations\n") ;
  Message_Direct("\
\t tol                 = tolerance\n") ;
  Message_Direct("\
\t nrep                = maximum number of repetitions\n") ;

  Message_Direct("\n") ;
  
  Message_Direct("\
  If the process doesn't converge the time step is divided by 2 as \n\
  many times as the number of repetitions defined by \"nrep\" \n\
  except if \"nrep\" is zero.\n") ;
}


void Help_TimeStep(void)
{
  Message_Direct("\n") ;
  Message_Direct("Definition of parameters required to compute the time steps\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  These parameters allow to compute the time steps as follows:\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  From the beginning of each date the first time step is initialized \n\
  with \"dtini\". \n") ;
  
  Message_Direct("\
  The following time steps are then determined in such a way that \n\
  the solution increment is of the order of magnitude of the value \n\
  defined by the key-word OBJE.\n") ;
  Message_Direct("\
  The reducing factor is limited to 1.5 and the time step \n\
  itself to \"dtmax\".\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
->\t \"Dtini =\",dtini,\"Dtmax =\",dtmax\n") ;

  Message_Direct("\n") ;

  Message_Direct("\
\t dtini               = initial time step\n") ;
  Message_Direct("\
\t dtmax               = maximum time step\n") ;
}

void Help_BuiltinCurves(void)
{
  Message_Direct("\n") ;
  Message_Direct("Built-in curves\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Range:\n") ;
  
  Message_Direct("\
  x = Range{x1 = %%lf , x2 = %%lf , n = %%d}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Langmuir: y0 * x / (x0 + x)\n") ;
  
  Message_Direct("\
  y = Langmuir(i){y0 = %%lf , x0 = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("LangmuirN: y0 * (x/x0)**n / (1 + (x/x0)**n)\n") ;
  
  Message_Direct("\
  y = Langmuir(i){y0 = %%lf , x0 = %%lf , n = %%d}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Freundlich: a * (x**b)\n") ;
  
  Message_Direct("\
  y = Freundlich(i){a = %%lf , b = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Mualem_wet: sqrt(x)*(1 - (1 - x**(1/m))**m)**2\n") ;
  
  Message_Direct("\
  y = Mualem_liq(i){m = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Mualem_liq: sqrt(x)*(1 - (1 - x**(1/m))**m)**2\n") ;
  
  Message_Direct("\
  y = Mualem_liq(i){m = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Van-Genuchten: (1 + (x/x0)**(1/(1-m)))**(-m)\n") ;
  
  Message_Direct("\
  y = Van-Genuchten(i){x0 = %%lf , m = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Nav-Genuchten: y0 * (x**(-1/m) - 1)**(1-m)\n") ;
  
  Message_Direct("\
  y = Nav-Genuchten(i){y0 = %%lf , m = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Mualem_dry: Mualem's model for drying paths\n") ;
  
  Message_Direct("\
  y = Mualem_dry(i){a_w = %%lf, m_w = %%lf , a_d = %%lf, m_d = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Mualem_gas: sqrt(1-x)*(1 - x**(1/m))**(2*m)\n") ;
  
  Message_Direct("\
  y = Mualem_gas(i){m = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Van-Genuchten_gas: ((1-x)**p)*(1 - x**(1/m))**(q*m)\n") ;
  
  Message_Direct("\
  y = Van-Genuchten_gas(i){m = %%lf , p = %%lf , q = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Monlouis-Bonnaire: ((1-x)**(5.5))*(1 - x**(1/m))**(2*m)\n") ;
  
  Message_Direct("\
  y = Monlouis-Bonnaire(i){m = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Millington: (1 - x)**b\n") ;
  
  Message_Direct("\
  y = Millington(i){b = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("CSH3Poles: solid solution with 3 end-members\n") ;
  
  Message_Direct("Write 4 curves: C/S , H/S , V_csh , S_sh \n") ;
  
  Message_Direct("\
  y = CSH3Poles(i){y_Tob = %%lf , y_Jen = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("CSHLangmuirN:\n") ;
  
  Message_Direct("For x:") ;
  Message_Direct("y1 * (x/x1)**n1 / (1 + (x/x1)**n1)") ;
  Message_Direct(" + ") ;
  Message_Direct("y2 * (x/x2)**n2 / (1 + (x/x2)**n2)") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  x = CSHLangmuirN(i){y1 = %%lf , x1 = %%lf , n1 = %%lf , y2 = %%lf , x2 = %%lf , n2 = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("For SHt:") ;
  Message_Direct("(1 + (x/x1)**n1))**(-y1/n1) * (1 + (x/x2)**n2))**(-y2/n2)") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("\
  SHt = CSHLangmuirN(i){y1 = %%lf , x1 = %%lf , n1 = %%lf , y2 = %%lf , x2 = %%lf , n2 = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("MolarDensityOfPerfectGas: x/RT\n") ;
  
  Message_Direct("\
  y = MolarDensityOfPerfectGas(i){T = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("MolarDensityOfCO2: Redlich-Kwong's model\n") ;
  
  Message_Direct("\
  y = MolarDensityOfCO2(i){T = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("RedlichKwongCO2: Redlich-Kwong's model (in mole/L)\n") ;
  
  Message_Direct("\
  y = RedlichKwongCO2(i){T = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("ViscosityOfCO2: Fenghour's model\n") ;
  
  Message_Direct("\
  y = ViscosityOfCO2(i){T = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Affine: a*x + b\n") ;
  
  Message_Direct("\
  y = Affine(i){a = %%lf , b = %%lf}\n") ;
  
  Message_Direct("\n") ;
  
  Message_Direct("Expressions: F(a,b,c,...,x)\n") ;
  
  Message_Direct("\
  y = Expressions(i){a = %%lf ; b = %%lf ; c = %%f ; ... ; y = F(a,b,c,...,x)}\n") ;
  
  Message_Direct("\n") ;
}
