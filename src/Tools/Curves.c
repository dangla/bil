#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>

#include "Message.h"
#include "Tools/Math.h"
#include "Buffer.h"
#include "CurvesFile.h"
#include "Curves.h"


static double courbe_nor(double,Curve_t*) ;
static double dcourbe_nor(double,Curve_t*) ;
static double icourbe_nor(double,Curve_t*) ;
static double courbe_log(double,Curve_t*) ;
static double dcourbe_log(double,Curve_t*) ;
static double icourbe_log(double,Curve_t*) ;


/* Extern functions */

Curves_t* Curves_Create(unsigned int n_curves)
{
  Curves_t *curves   = (Curves_t*) malloc(sizeof(Curves_t)) ;
  
  if(!curves) arret("Curves_Create (0)") ;
  
  Curves_GetNbOfAllocatedCurves(curves) = n_curves ;
  Curves_GetNbOfCurves(curves) = 0 ; /* Important initialization */
  
  {
    Curve_t* cv = (Curve_t*) malloc(n_curves*sizeof(Curve_t)) ;
    
    if(!cv) arret("Curves_Create (1)") ;
    
    Curves_GetCurve(curves) = cv ;
  }
    
  /* Space allocation for buffer */
  Curves_GetBuffer(curves) = Buffer_Create(Curves_SizeOfBuffer) ;
  
  return(curves) ;
}


void Curves_Delete(Curves_t** curves)
{
  free(Curves_GetCurve(*curves)) ;
  free(Curves_GetBuffer(*curves)) ;
  free(*curves) ;
}


int Curves_Append(Curves_t* curves, Curve_t* cv)
/** Append curves with cv 
 *  Return the index of the appended curve */
{
    int NbOfCurves = Curves_GetNbOfCurves(curves) ;
    Curve_t* curve = Curves_GetCurve(curves) ;
    
    if(Curves_CannotAppendCurves(curves,1)) {
      arret("Curves_Append") ;
    }
    
    curve[NbOfCurves] = *cv ;
    
    Curves_GetNbOfCurves(curves) += 1 ;
    
    return(NbOfCurves) ;
}


int    Curves_FindCurveIndex(Curves_t* curves,const char* yname)
/** Find the curve position index whose y-axis name is pointed to by yname.
 *  Return -1 if the name is not found. */
{
  int n = Curves_GetNbOfCurves(curves) ;
  int    i ;

  if(isdigit(yname[0])) { /* numeric characters */
    i  = atoi(yname) - 1 ;
    
    if(i >= n) i = -1 ;
    
  } else {                /* alphanumeric characters */
    Curve_t* curve = Curves_GetCurve(curves) ;
    size_t ylen = strlen(yname) ;
    
    for(i = 0 ; i < n ; i++) {
      char* yyname = Curve_GetNameOfYAxis(curve + i) ;
      if(!strncmp(yname,yyname,ylen)) break ;
    }
    
    if(i == n) i = -1 ;
  }

  //if(i < 0) arret("Curves_FindCurveIndex(2): position not known") ;
  return(i) ;
}



Curve_t* Curves_FindCurve(Curves_t* curves,const char* label)
{
  int i = Curves_FindCurveIndex(curves,label) ;
  Curve_t* curve ;
    
  if(i < 0) {
    curve = NULL ;
  } else {
    curve = Curves_GetCurve(curves) + i ;
  }
  
  return(curve) ;
}



Curve_t* Curve_Create(unsigned int n_points)
{
  Curve_t* curve   = (Curve_t*) malloc(sizeof(Curve_t)) ;
  
  if(!curve) arret("Curve_Create (0)") ;

  Curve_GetNbOfPoints(curve) = n_points ;
  
  
  /* Allocate memory space for the values */
  {
    double* mry = (double *) malloc((n_points + 2)*sizeof(double)) ;
    
    if(!mry) arret("Curve_Create (1) : not enough memory") ;
    
    Curve_GetXRange(curve) = mry ;
    Curve_GetYValue(curve) = mry + 2 ;
  }
  
  
  /* Allocate memory space for the names of axis */
  {
    char* name = (char*) malloc(2*Curve_MaxLengthOfCurveName*sizeof(char)) ;
    
    if(!name) arret("Curve_Create (2) : not enough memory") ;
    
    Curve_GetNameOfXAxis(curve) = name ;
    Curve_GetNameOfYAxis(curve) = name + Curve_MaxLengthOfCurveName ;
    
    name[0] = '\0' ;
    name[Curve_MaxLengthOfCurveName] = '\0' ;
  }
  
  return(curve) ;
}


void Curve_Delete(Curve_t** curve)
{
  free(Curve_GetXRange(*curve)) ;
  free(Curve_GetNameOfXAxis(*curve)) ;
  free(*curve) ;
}


Curve_t* Curve_CreateDerivative(Curve_t* curve)
{
  int n_points    = Curve_GetNbOfPoints(curve) ;
  Curve_t* dcurve = Curve_Create(n_points) ;
  double* dy      = Curve_GetYValue(dcurve) ;
  double* x       = Curve_CreateSamplingOfX(curve) ;
  int i ;
  
  Curve_GetNbOfPoints(dcurve) = Curve_GetNbOfPoints(curve) ;
  Curve_GetXRange(dcurve)[0]  = Curve_GetXRange(curve)[0] ;
  Curve_GetXRange(dcurve)[1]  = Curve_GetXRange(curve)[1] ;
  Curve_GetScaleType(dcurve)  = Curve_GetScaleType(curve) ;
  
  for(i = 0 ; i < n_points ; i++) {
    dy[i] = Curve_ComputeDerivative(curve,x[i]) ;
  }
  
  free(x) ;

  return(dcurve) ;
}


Curve_t* Curve_CreateIntegral(Curve_t* curve)
{
  int n_points    = Curve_GetNbOfPoints(curve) ;
  Curve_t* icurve = Curve_Create(n_points) ;
  double* intydx  = Curve_GetYValue(icurve) ;
  double* x       = Curve_CreateSamplingOfX(curve) ;
  int i ;
  
  Curve_GetNbOfPoints(icurve) = Curve_GetNbOfPoints(curve) ;
  Curve_GetXRange(icurve)[0]  = Curve_GetXRange(curve)[0] ;
  Curve_GetXRange(icurve)[1]  = Curve_GetXRange(curve)[1] ;
  Curve_GetScaleType(icurve)  = Curve_GetScaleType(curve) ;
  
  for(i = 0 ; i < n_points ; i++) {
    intydx[i] = Curve_ComputeIntegral(curve,x[i]) ;
  }
  
  free(x) ;

  return(icurve) ;
}


Curve_t* Curve_CreateInverse(Curve_t* curve,const char scale)
/** Create the inverse of a function assumed as monotoneous.
 *  The sampling is defined by scale.
 *  Return a pointer to Curve_t.
 */
{
  int     n = Curve_GetNbOfPoints(curve) ;
  Curve_t* evruc = Curve_Create(n) ;
  
  Curve_GetNbOfPoints(evruc) = Curve_GetNbOfPoints(curve) ;
  Curve_GetScaleType(evruc)  = scale ;
  Curve_GetXRange(evruc)[0]  = Curve_GetYValue(curve)[0] ;
  Curve_GetXRange(evruc)[1]  = Curve_GetYValue(curve)[n - 1] ;
  
  {
    double *xo  = Curve_CreateSamplingOfX(curve) ;
    double *yo  = Curve_GetYValue(curve) ;
    double *yi  = Curve_CreateSamplingOfX(evruc) ;
    double *xi  = Curve_GetYValue(evruc) ;
    double ymax = MAX(fabs(yo[0]),fabs(yo[n - 1])) ;
    double err ;
    double tol = 1.e-10*ymax ;
    int i ;
    
    /* Compute the xi's */
    xi[0    ] = Curve_GetXRange(curve)[0] ;
    xi[n - 1] = Curve_GetXRange(curve)[1] ;
    
    for(i = 1 ; i < n - 1 ; i++) {
      double y = yi[i] ;
      double x = xi[i - 1] ;
      
      /* Initialize x more efficiently */
      {
        int j = 1 ;
        
        while( (y - yo[j - 1])*(y - yo[j]) > 0) j++ ;
        
        x = xo[j - 1] ;
      }
      
      /* Given y, find x = Inv[f](y) i.e. y = f(x) */
      {
        int    j = 0 ;
    
        do {
          double f  = Curve_ComputeValue(curve,x) ;
          double df = Curve_ComputeDerivative(curve,x) ;
          double dx ;
      
          if(df == 0) {
            arret("Curve_CreateInverse: inversion impossible") ;
          }
          
          dx = - (f - y)/df ;
      
          x += dx ;
      
          err = fabs(f - y) ;
      
          if(j++ > 50) {
            printf("f    = %e\n",f) ;
            printf("df   = %e\n",df) ;
            printf("dx   = %e\n",dx) ;
            arret("Curve_CreateInverse: no convergence") ;
          }
        } while(err > tol) ;
        
        xi[i] = x ;
      }
    }
    
    free(yi) ;
    free(xo) ;
  }
  
  return(evruc) ;
}


double* Curve_CreateSamplingOfX(Curve_t *curve)
{
  int n_points = Curve_GetNbOfPoints(curve) ;
  double *x    = Curve_GetXRange(curve) ;
  char scale   = Curve_GetScaleType(curve) ;
  int ni    = n_points - 1 ;
  double x1 = x[0] ;
  double x2 = x[1] ;
  double *xs = (double*) malloc(n_points*sizeof(double)) ;
  
  if(!xs) arret("Curve_CreateSamplingOfX: not enough memory") ;
  
  if(scale == 'n') {
    double dx = (x2 - x1) ;
    int i ;
    
    for(i = 0 ; i < n_points ; i++) {
      double n = ((double) i)/ni ;
      
      xs[i] = x1 + n*dx ;
    }
    
  } else if(scale == 'l') {
    double ratio = x2/x1 ;
    int i ;
    
    if(x1 <= 0. || x2 <= 0.) {
      arret("Curve_CreateSamplingOfX(1)") ;
    }
    
    for(i = 0 ; i < n_points ; i++) {
      double n = ((double) i)/ni ;
      
      xs[i] = x1*pow(ratio,n) ;
    }
    
  } else {
    arret("Curve_CreateSamplingOfX(2)") ;
  }

  return(xs) ;
}


double Curve_ComputeValue(Curve_t *cb,double a)
/** Return the value at a */
{
  if(Curve_GetScaleType(cb) == 'n') return(courbe_nor(a,cb)) ;
  else if(Curve_GetScaleType(cb) == 'l') return(courbe_log(a,cb)) ;
  else arret("Curve_ComputeValue: option non prevue") ;
  return(0.) ;
}

double Curve_ComputeDerivative(Curve_t *cb,double a)
/** Return the derivative at a */
{
  if(Curve_GetScaleType(cb) == 'n') return(dcourbe_nor(a,cb)) ;
  else if(Curve_GetScaleType(cb) == 'l') return(dcourbe_log(a,cb)) ;
  else arret("Curve_ComputeDerivative: option non prevue") ;
  return(0.) ;
}

double Curve_ComputeIntegral(Curve_t *cb,double a)
/** Return the integral from begin to a */
{
  if(Curve_GetScaleType(cb) == 'n') return(icourbe_nor(a,cb)) ;
  else if(Curve_GetScaleType(cb) == 'l') return(icourbe_log(a,cb)) ;
  else arret("Curve_ComputeIntegral: option non prevue") ;
  return(0.) ;
}

char* Curve_PrintInFile(Curve_t* curve)
{
  int     n = Curve_GetNbOfPoints(curve) ;
  double* x = Curve_CreateSamplingOfX(curve) ;
  double* y = Curve_GetYValue(curve) ;
  /* Create a temporary filename */
  /* char*   targetfile = (char*) tempnam(".",NULL) ; *//* no portability */
  char*   targetfile = (char*) tmpnam(NULL) ;
  FILE*   target     = fopen(targetfile,"w") ;
  int i ;
  
  for(i = 0 ; i < n ; i++) {
    fprintf(target,"%e %e\n",x[i],y[i]) ;
  }
  
  fclose(target) ;
  
  free(x) ;
  
  return(targetfile) ;
}




int   Curves_ReadCurves(Curves_t* curves,const char* dline)
/** Read new curves as defined in the filename found in dline 
 *  and append "curves" accordingly.
 *  Return the nb of curves that has been read */
{
  int    n_curves ;
  int    n_points ;
  CurvesFile_t* curvesfile = CurvesFile_Create() ;
  int    i ;

  /* Is there a file name? */
  if(!CurvesFile_Initialize(curvesfile,dline)) {
    /* Create it if needed. */
    CurvesFile_GetCurves(curvesfile) = curves ;
    CurvesFile_WriteCurves(curvesfile) ;
    /* We reinitialize again for the file position starting input data*/
    /* CurvesFile_Initialize(curvesfile,dline) ; */
  }


  /* Nb of curves */
  n_curves = CurvesFile_GetNbOfCurves(curvesfile) ;


  /* Nb of points */
  n_points = CurvesFile_GetNbOfPoints(curvesfile) ;


  {    
    if(Curves_CannotAppendCurves(curves,n_curves)) {
      arret("Curves_ReadCurves (4) : trop de courbes") ;
    }
  }
  
  
  /* Allocate memory for the curves */
  for(i = 0 ; i < n_curves ; i++) {
    Curve_t *cb_i = Curves_GetCurve(curves) + Curves_GetNbOfCurves(curves) + i ;
    Curve_t *cb   = Curve_Create(n_points) ;
    
    *cb_i = *cb ;
  }
  
  
  /* Read and store the names of x-axis and y-axis */
  {
    char *line ;
    
    CurvesFile_OpenFile(curvesfile,"r") ;
    
    do {
      
      line = CurvesFile_ReadLineFromCurrentFilePosition(curvesfile) ;
      
    } while((line) && (line[0] == '#') && (strncmp(line,"# Labels:",9) != 0)) ;
    
    if(strncmp(line,"# Labels:",9) == 0) {
      
      line += 9 ;
      
      /* Read the label of x-axis */
      line = strtok(line," ") ;
      
      if(line) {
        char xlabel[Curve_MaxLengthOfCurveName] ;
        char* c ;
        
        /* Save the label of x-axis */
        strcpy(xlabel,line) ;
            
        if((c = strchr(xlabel,'('))) c[0] = '\0' ;
      
        for(i = 0 ; i < n_curves ; i++) {
          Curve_t *cb_i = Curves_GetCurve(curves) + Curves_GetNbOfCurves(curves) + i ;
          char* xname = Curve_GetNameOfXAxis(cb_i) ;
          char* yname = Curve_GetNameOfYAxis(cb_i) ;
        
          /* Store the label of x-axis */
          strcpy(xname,xlabel) ;
        
          /* Read the label of y-axis */
          line = strtok(NULL," ") ;
        
          /* Store the label of y-axis */
          if(line) {
            strcpy(yname,line) ;
            
            if((c = strchr(yname,'('))) c[0] = '\0' ;
          }
        }
      }
    }
    
    CurvesFile_CloseFile(curvesfile) ;
  }


  /* Read the x-values and y-values */
  {
    double a_1,a_2 ;
    char *line ;
    FILE   *fict = CurvesFile_OpenFile(curvesfile,"r") ;
    
    do {
      CurvesFile_StoreFilePosition(curvesfile) ;
      
      line = CurvesFile_ReadLineFromCurrentFilePosition(curvesfile) ;
      
    } while((line) && (line[0] == '#')) ;
      
    CurvesFile_MoveToStoredFilePosition(curvesfile) ;
  
    /* Set position in the file */
    /* CurvesFile_MoveToFilePositionStartingInputData(curvesfile) ; */

  
    /* Read x and y */
    for(i = 0 ; i < n_points ; i++) {
      double x ;
      int j ;
      
      fscanf(fict,"%le",&x) ;
    
      if(i == 0) {
        a_1 = x ;
      } else {
        a_2 = x ;
      }
    
      for(j = 0 ; j < n_curves ; j++) {
        Curve_t *cb_j = Curves_GetCurve(curves) + Curves_GetNbOfCurves(curves) + j ;
        double *y = Curve_GetYValue(cb_j) + i ;
        char fmt[] = "%*["CurvesFile_FieldDelimiters"] %le" ;
      
        fscanf(fict,fmt,y) ;
      }
    }
  
    /* The same first and last points for all the curves */
    for(i = 0 ; i < n_curves ; i++) {
      Curve_t *cb_i = Curves_GetCurve(curves) + Curves_GetNbOfCurves(curves) + i ;
    
      Curve_GetXRange(cb_i)[0] = a_1 ;
      Curve_GetXRange(cb_i)[1] = a_2 ;
    }
  
    CurvesFile_CloseFile(curvesfile) ;
  }


  /* Scale */
  {
    char   scale = CurvesFile_GetScaleType(curvesfile) ;

    for(i = 0 ; i < n_curves ; i++) {
      Curve_t *cb_i = Curves_GetCurve(curves) + Curves_GetNbOfCurves(curves) + i ;
      
      Curve_GetScaleType(cb_i) = scale ;
    }
  }
  
  Curves_GetNbOfCurves(curves) += n_curves ;
  
  CurvesFile_Delete(&curvesfile) ;

  return(n_curves) ;
}



/* Intern functions */

double courbe_nor(double a,Curve_t *cb)
{
  /* Retourne la valeur de la courbe en a */
  int    ni = Curve_GetNbOfPoints(cb) - 1 ;
  double a1 = Curve_GetXRange(cb)[0] ;
  double a2 = Curve_GetXRange(cb)[1] ;

  if(a < a1) return(Curve_GetYValue(cb)[0]) ;
  else if(a > a2) return(Curve_GetYValue(cb)[ni]) ;
  else {
    double da = (a2 - a1)/ni ;
    double r  = (a - a1)/da ;
    int i  = floor(r) ;
    double a0 = a1 + i*da ;
    double dv = (Curve_GetYValue(cb)[i+1] - Curve_GetYValue(cb)[i])/da ;
    double v  = Curve_GetYValue(cb)[i] + dv*(a - a0) ;
    return(v) ;
  }
}

double dcourbe_nor(double a,Curve_t *cb)
{
  /* Retourne la derivee de la courbe en a */
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double da  = (Curve_GetXRange(cb)[1] - Curve_GetXRange(cb)[0])/n_i ;

  return((courbe_nor(a + da,cb) - courbe_nor(a - da,cb))*0.5/da) ;
}


double icourbe_nor(double a,Curve_t *cb)
/* Return the integral computed from cb */ 
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double *x = Curve_GetXRange(cb) ;
  double *y = Curve_GetYValue(cb) ;
  double x1 = x[0] ;
  double x2 = x[1] ;
  double dx = (x2 - x1)/n_i ;
  double fa = Curve_ComputeValue(cb,a) ;
  double intydx = 0 ;
  int    i ;
  
  for(i = 0 ; x1 + (i + 1)*dx < a ; i++) {
    intydx += y[i] + y[i + 1] ;
  }
  
  intydx *= dx*0.5 ;
  
  intydx += (y[i] + fa)*(a - x1 - i*dx)*0.5 ;
  
  return(intydx) ;
}

double courbe_log(double a,Curve_t *cb)
/* Retourne la valeur en a de la courbe echantillonnee en base log10 */
{
  int    ni = Curve_GetNbOfPoints(cb) - 1 ;
  double a1 = Curve_GetXRange(cb)[0] ;
  double a2 = Curve_GetXRange(cb)[1] ;

  if(a < a1) return(Curve_GetYValue(cb)[0]) ;
  else if(a > a2) return(Curve_GetYValue(cb)[ni]) ;
  else {
    double loga1 = log10(a1) ;
    double loga2 = log10(a2) ;
    double dloga = (loga2 - loga1)/ni ;
    double loga  = log10(a) ;
    double r  = (loga - loga1)/dloga ;
    int    i  = floor(r) ;
    double loga0 = loga1 + i*dloga ;
    double dv = (Curve_GetYValue(cb)[i+1] - Curve_GetYValue(cb)[i])/dloga ;
    double v  = Curve_GetYValue(cb)[i] + dv*(loga - loga0) ;
    return(v) ;
  }
}

double dcourbe_log(double a,Curve_t *cb)
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double a1 = Curve_GetXRange(cb)[0],a2 = Curve_GetXRange(cb)[1] ;
  double loga1 = log10(a1),loga2 = log10(a2) ;
  double dloga = (loga2 - loga1)/n_i ;
  double ada = a*pow(10.,dloga),da = ada - a ;

  if(a < a1) return(0.) ; /* pour le cas a = 0 ! */
  else return((courbe_log(a + da,cb) - courbe_log(a - da,cb))*0.5/da) ;
}

double icourbe_log(double a,Curve_t *cb)
/* Return the integral curve computed from cb */ 
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double *x = Curve_GetXRange(cb) ;
  double *y = Curve_GetYValue(cb) ;
  double x1 = x[0] ;
  double x2 = x[1] ;
  double logx1 = log10(x1) ;
  double logx2 = log10(x2) ;
  double dlogx = (logx2 - logx1)/n_i ;
  double loga  = log10(a) ;
  double ln10  = log(10) ;
  double ya = Curve_ComputeValue(cb,a) ;
  double intydx = 0 ;
  int    i ;
  
  for(i = 0 ; logx1 + (i + 1)*dlogx < loga ; i++) {
    double logxi = logx1 + i*dlogx ;
    double xi    = pow(10,logxi) ;
    double logxj = logxi + dlogx ;
    double xj    = pow(10,logxj) ;
    
    intydx += y[i]*xi + y[i + 1]*xj ;
  }
  
  intydx *= dlogx*0.5*ln10 ;
  
  {
    double logxi = logx1 + i*dlogx ;
    double xi    = pow(10,logxi) ;
    
    intydx += (y[i]*xi + ya*a)*(loga - logxi)*0.5*ln10 ;
  }
  
  return(intydx) ;
}



/* Not used from here */

#ifdef NOTDEFINED
int   (Curves_WriteCurves1)(char *dline)
/** Write curves in a file as discrete data */
{
#define MAX_COLX  5
#define FCOLX(i) (frmt + 4*(MAX_COLX + 1 - i))
  char   frmt[] = "%*s %*s %*s %*s %*s %lf" ;
  int    n_courbes = 0 ;
  FILE   *fic ;
  char   filename[Curve_MaxLengthOfFileName] ;
  char   model[Curve_MaxLengthOfKeyWord],label[Curve_MaxLengthOfKeyWord] ;
  char   *line,scale = 'n' ;
  int    colx ;

  /* Scale ? */
  sscanf(dline," %s",model) ;
  if(strstr(model,"_log")) scale = 'l' ;


  /* File name */
  line = strchr(dline,'=') + 1 ;
  sscanf(line," %s",filename) ;
  line  = strstr(line,filename) + strlen(filename) ;


  /* Write the first column x-axis */
  {
    int    n_points ;
    double x_1,x_2 ;
    int i ;
    
    /* Read the range and the nb of points */
    if(sscanf(line," %s = %[^{]",label,model) == 2) {
      if(!strcmp(model,"Range")) {
        line = strchr(line,'{') ;
        sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %d }",&x_1,&x_2,&n_points) ;
        line = strchr(line,'}') + 1 ;
      } else {
        arret("le mot-cle n\'est pas \"Range\"") ;
        return(0) ;
      }
    } else {
      arret("impossible de construire la courbe \"%s\"",filename) ;
      return(0) ;
    }


    fic = fopen(filename,"w") ;
    if(!fic) arret("Curves_WriteCurves(1) : impossible d ouvrir le fichier") ;

    /* The first lines must contain some informations:
     * - The name of models used to compute the columns 
     * - The name used for the columns */
    fprintf(fic,"# Models: X-axis(%d)\n",n_courbes+1) ;
    fprintf(fic,"# Labels: %s(%d)\n",label,n_courbes+1) ;

    /* Write the column x-axis */
    for(i = 0 ; i < n_points ; i++){
      int    n1 = n_points - 1 ;
      double x,dx = (x_2 - x_1)/n1 ;
      
      if(scale == 'l') x = x_1*pow(x_2/x_1,((double) i)/n1) ;
      else x = x_1 + i*dx ;
      fprintf(fic,"%e\n",x) ;   
    }
    fclose(fic) ;
  }


  /* Write other columns y-axis */
  /* Read the key-word of the curve */
  while(sscanf(line," %s = %[^(](%d)",label,model,&colx) == 3) {
    char   ltmp[Curve_MaxLengthOfTextLine] ;
    FILE   *ftmp = tmpfile() ; /* temporary file */
    
    n_courbes += 1 ;
    
    if(colx > MAX_COLX || colx > n_courbes) {
      arret("Curves_WriteCurves(1) : colonne impossible") ;
    }
    
    line = strchr(line,'{') ;

    /* Copy the file "filename" into the temporary file ftmp and
     * add the name of the model and the label of the y-axis */
    fic = fopen(filename,"r") ;
    if(!fic) arret("Curves_WriteCurves(1) : impossible d ouvrir le fichier") ;
    
    while(fgets(ltmp,sizeof(ltmp),fic)) {
      /* We add the name of the model referred to by "model" */
      if(!strncmp(ltmp,"# Models",8)) {
        sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n_courbes+1) ;
        
      /* We add the name of the curve referred to by "label" */
      } else if(!strncmp(ltmp,"# Labels",7)) {
        sprintf(strchr(ltmp,'\n')," %s(%d)\n",label,n_courbes+1) ;
      }
      
      fprintf(ftmp,"%s",ltmp) ;
    }
    
    fclose(fic) ;

    rewind(ftmp) ;

    fic = fopen(filename,"w") ;
    if(!fic) arret("Curves_WriteCurves(1) : impossible d ouvrir le fichier") ;

    /* Write a new column depending on the name of the model */
    if(!strcmp(model,"Freundlich")) {
      double alpha,beta ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&alpha,&beta) ;

      AppendCurves(fic,ftmp,colx,Freundlich,alpha,beta) ;
      
      
    } else if(!strcmp(model,"Langmuir")) {
      double c0,ca_max ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&ca_max,&c0) ;

      AppendCurves(fic,ftmp,colx,Langmuir,ca_max,c0) ;

      
    } else if(!strcmp(model,"Mualem_wet") || !strcmp(model,"Mualem_liq")){
      double m ;
      
      sscanf(line,"{ %*s = %lf }",&m) ;
      
      AppendCurves(fic,ftmp,colx,MualemVanGenuchten,m) ;
      
      
    } else if(!strcmp(model,"Mualem_dry")){
      double m_w,m_d,a_d,a_w ;
      char   ltmp1[Curve_MaxLengthOfTextLine] ;
      FILE   *ftmp1 = tmpfile() ;
      double kh_max ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&a_w,&m_w,&a_d,&m_d) ;

      /* calcul de kh dans un fichier temporaire */
      {
        double kh = 0. ;
        double p ;
        do {
          fgets(ltmp,sizeof(ltmp),ftmp) ;
        } while(ltmp[0] == '#') ;
	
        sscanf(ltmp,FCOLX(colx),&p) ;
        fprintf(ftmp1,"%e %e\n",p,kh) ;
	
        while(fgets(ltmp,sizeof(ltmp),ftmp)) {
          if(ltmp[0] != '#') {
            double pdp = p ;
            sscanf(ltmp,FCOLX(colx),&p) ;
            {
              double dp  = p - pdp ;
              double s_d = vangenuchten(pdp/a_d,m_d) ;
              double s_w = vangenuchten(pdp/a_w,m_w) ;
              double hdh = (s_d - s_w)/(1. - s_w) ;
              double h,dh ;
	      
              s_d = vangenuchten(p/a_d,m_d) ;
              s_w = vangenuchten(p/a_w,m_w) ;
              h   = (s_d - s_w)/(1. - s_w) ;
	      
              dh  = h - hdh ;
              kh += dh/(p - 0.5*dp) ;
	      
              fprintf(ftmp1,"%e %e\n",p,kh) ;
            }
          }
        }
        kh_max = kh ;
      }
      
      rewind(ftmp1) ;
      rewind(ftmp) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double kh ;
          double p ;
	  
          fgets(ltmp1,sizeof(ltmp1),ftmp1) ;
          sscanf(ltmp1,FCOLX(2),&kh) ;
          kh = 1. - kh/kh_max ;
	  
          sscanf(ltmp,FCOLX(colx),&p) ;
	  
          {
            double s_w  = vangenuchten(p/a_w,m_w) ;
            double kl   = 1 - pow(1 - pow(s_w,1./m_w),m_w) ;
	    
            double s_d  = vangenuchten(p/a_d,m_d) ;
	    
            double k_rd = sqrt(s_d)*(kl + (1. - kl)*kh) ;
	    
            sprintf(strchr(ltmp,'\n')," %e\n",k_rd) ;
          }
        }
        fprintf(fic,"%s",ltmp) ;
      }
      fclose(ftmp1) ;
      
    } else if(!strcmp(model,"Mualem_gas")){
      double m ;
      
      sscanf(line,"{ %*s = %lf }",&m) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,k_rg ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          k_rg = mualem_gas(s_w,m) ;
          sprintf(strchr(ltmp,'\n')," %e\n",k_rg) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Van-Genuchten")){
      double a,m ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&a,&m) ;
      
      AppendCurves(fic,ftmp,colx,VanGenuchten,a,m) ;
      
    } else if(!strcmp(model,"Nav-Genuchten")){
      double a,m ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&a,&m) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,p_c ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          p_c = a*navgenuchten(s_w,m) ;
          sprintf(strchr(ltmp,'\n')," %e\n",p_c) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Baroghel")){
      double a,b,v,mu ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&a,&b,&v,&mu) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double hr,s_w,p_c ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          hr = a*pow(s_w,v) + (1 - a)*(1 + pow(b,mu))*pow(s_w,mu)/(pow(s_w,mu) + pow(b,mu)) ;
          p_c = -1.372e8*log(hr) ;
          sprintf(strchr(ltmp,'\n')," %e\n",p_c) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Francy_liq")){
      double s_r ;
      
      sscanf(line,"{ %*s = %lf }",&s_r) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,k_rw ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          if(s_w > s_r) k_rw = pow((s_w - s_r)/(1 - s_r),3.) ;
          else k_rw = 0. ;
          sprintf(strchr(ltmp,'\n')," %e\n",k_rw) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Francy_ion")){
      double s_r ;
      
      sscanf(line,"{ %*s = %lf }",&s_r) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,k_re ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          if(s_w > s_r) k_re = 1./(1 + pow(-3*log((s_w - s_r)/(1 - s_r)),2.)) ;
          else k_re = 0. ;
          sprintf(strchr(ltmp,'\n')," %e\n",k_re) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Pow")){
      double a ;
      
      sscanf(line,"{ %*s = %lf }",&a) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,k_re ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          k_re = pow(s_w,a) ; 
          sprintf(strchr(ltmp,'\n')," %e\n",k_re) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Barry")){
      double a ;
      
      sscanf(line,"{ %*s = %lf }",&a) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,k_sl ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          k_sl = 1./(1 + a*pow(1.-s_w,4.)) ; 
          sprintf(strchr(ltmp,'\n')," %e\n",k_sl) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Monlouis")){
      double m ;
      
      sscanf(line,"{ %*s = %lf }",&m) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,k_rl ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          k_rl = monlouis_bonnaire(s_w,m) ;
          sprintf(strchr(ltmp,'\n')," %e\n",k_rl) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Millington")){
      double b ;
      
      sscanf(line,"{ %*s = %lf }",&b) ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,tau_g ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          tau_g = millington(s_w,b) ;
          sprintf(strchr(ltmp,'\n')," %e\n",tau_g) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Bazant")){

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double s_w,tau_l ;
          sscanf(ltmp,FCOLX(colx),&s_w) ;
          tau_l = bazant(s_w) ; 
          sprintf(strchr(ltmp,'\n')," %e\n",tau_l) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Integral")){
      int    coly ;
      double x,y,int_ydx = 0. ;
      double x0,y0 ;
      int    ligne = 0 ;
      
      sscanf(line,"{ %*s = %d }",&coly) ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          sscanf(ltmp,FCOLX(colx),&x) ;
          sscanf(ltmp,FCOLX(coly),&y) ;
          if(ligne++ > 0) {
            double ydx = 0.5*(y + y0)*(x - x0) ;
            int_ydx += ydx ;
          }
          sprintf(strchr(ltmp,'\n')," %e\n",int_ydx) ;
          x0 = x ;
          y0 = y ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Surface_stress")){
      int    coly ;
      double RT ;
      double x,y,int_ydx = 0. ;
      double x0,y0 ;
      int    ligne = 0 ;
      
      sscanf(line,"{ %*s = %d , %*s = %lf }",&coly,&RT) ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          sscanf(ltmp,FCOLX(colx),&x) ;
          sscanf(ltmp,FCOLX(coly),&y) ;
          y *= RT/x ;
          if(ligne++ > 0) {
            double ydx = 0.5*(y + y0)*(x - x0) ;
            int_ydx += ydx ;
          }
          sprintf(strchr(ltmp,'\n')," %e\n",int_ydx) ;
          x0 = x ;
          y0 = y ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"CSH2Poles")){
      double k_CH  = 6.456e-6 ;
      /* Amorpheous Silica */
      double k_SH = 1.93642e-3 ;
      double x_SH = 0 ;
      double y_SH = 1 ;
      double z_SH = 2 ;
      double v_SH = 43.e-3 ;
      /* Jennite */
      double k_Jen0 = 4.39e-18,k_Jen = 2.39e-16 ;
      double x_Jen0 = 1.66666 ,x_Jen = 1.5     ;
      double                   y_Jen = 0.9     ;
      double z_Jen0 = 2.66666 ,z_Jen = 2.4     ;
      double v_Jen0 = 81.8e-3 ,v_Jen = 73.6e-3 ;
      
      double q_SH = 1. ;
      
      /* 
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf }",&k_CH,&k_SH,&k_Jen) ;
      */
      sscanf(line,"{ %*s = %lf }",&y_Jen) ;
      
      x_Jen   = x_Jen0*y_Jen ;
      z_Jen   = z_Jen0*y_Jen ;
      k_Jen   = pow(k_Jen0,y_Jen) ;
      v_Jen   = v_Jen0*y_Jen ;
      
      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(!strncmp(ltmp,"# Models",8)) {
          int j ;
          for(j=0;j<3;j++) sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n_courbes+2+j) ;
          fprintf(fic,"%s",ltmp) ;
          fgets(ltmp,sizeof(ltmp),ftmp) ;
          sprintf(strchr(ltmp,'\n')," y(%d)\n",n_courbes+2) ;
          sprintf(strchr(ltmp,'\n')," V(%d)\n",n_courbes+3) ;
          sprintf(strchr(ltmp,'\n')," q_SH(%d)\n",n_courbes+4) ;
        } else if(ltmp[0] != '#') {
          double q_CH ;
          double q_CSH[2] ;
          double a_CSH[2] ;
          double x_CSH[2] ;
          double y_CSH[2] ;
          double z_CSH[2] ;
          double k_CSH[2] ;
          double v_CSH[2] ;
          int    j ;
          k_CSH[0] = k_SH   ; 
          x_CSH[0] = x_SH   ; y_CSH[0] = y_SH   ; z_CSH[0] = z_SH ;
          k_CSH[1] = k_Jen  ; 
          x_CSH[1] = x_Jen  ; y_CSH[1] = y_Jen  ; z_CSH[1] = z_Jen ;
          v_CSH[0] = v_SH ;
          v_CSH[1] = v_Jen ;
          sscanf(ltmp,FCOLX(colx),&q_CH) ;
          for(j=0;j<2;j++) {
            a_CSH[j] = pow(k_CH*q_CH,x_CSH[j])*pow(k_SH,y_CSH[j])/k_CSH[j] ;
          }
          q_SH = fraction_Silica(a_CSH,y_CSH,2,q_SH) ;
          for(j=0;j<2;j++) {
            q_CSH[j] = a_CSH[j]*pow(q_SH,y_CSH[j]) ;
          }
          {
            double x_m = 0,y_m = 0,z_m = 0,v_m = 0 ;
            for(j=0;j<2;j++) {
              x_m += q_CSH[j]*x_CSH[j] ;
              y_m += q_CSH[j]*y_CSH[j] ;
              z_m += q_CSH[j]*z_CSH[j] ;
              v_m += q_CSH[j]*v_CSH[j] ;
            }
            x_m /= y_m ; z_m /= y_m ; v_m /= y_m ;
            sprintf(strchr(ltmp,'\n')," %e %e %e %e\n",x_m,z_m,v_m,q_SH) ;
          }
        }
        fprintf(fic,"%s",ltmp) ;
      }
      /* n_courbes += 4 ; */
      n_courbes += 3 ;
      
    } else if(!strcmp(model,"CSH3Poles")){
      double k_CH  = 6.456e-6 ;
      /* Amorpheous Silica */
      double k_SH = 1.93642e-3 ;
      double x_SH = 0 ;
      double y_SH = 1 ;
      double z_SH = 2 ;
      double v_SH = 43.e-3 ;
      /* Jennite */
      double k_Jen0 = 4.39e-18,k_Jen = 2.39e-16 ;
      double x_Jen0 = 1.66666 ,x_Jen = 1.5     ;
      double                   y_Jen = 0.9     ;
      double z_Jen0 = 2.66666 ,z_Jen = 2.4     ;
      double v_Jen0 = 81.8e-3 ,v_Jen = 73.6e-3 ;
      /* Tobermorite */
      double k_Tob0 = 1.684e-12,k_Tob = 6.42e-22 ;
      double x_Tob0 = 0.83333  ,x_Tob = 1.5 ;
      double                    y_Tob = 1.8 ;
      double z_Tob0 = 1.83333  ,z_Tob = 3.3 ;
      double v_Tob0 = 54.8e-3  ,v_Tob = 98.6e-3 ;
      
      double q_SH = 1. ;
      
      /*
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&k_CH,&k_SH,&k_Tob,&k_Jen) ;
      */
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&y_Tob,&y_Jen) ;
      
      x_Jen   = x_Jen0*y_Jen ;
      z_Jen   = z_Jen0*y_Jen ;
      k_Jen   = pow(k_Jen0,y_Jen) ;
      v_Jen   = v_Jen0*y_Jen ;
      x_Tob   = x_Tob0*y_Tob ;
      z_Tob   = z_Tob0*y_Tob ;
      k_Tob   = pow(k_Tob0,y_Tob) ;
      v_Tob   = v_Tob0*y_Tob ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(!strncmp(ltmp,"# Models",8)) {
          int j ;
          for(j=0;j<3;j++) sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n_courbes+2+j) ;
          fprintf(fic,"%s",ltmp) ;
          fgets(ltmp,sizeof(ltmp),ftmp) ;
          sprintf(strchr(ltmp,'\n')," z(%d)\n",n_courbes+2) ;
          sprintf(strchr(ltmp,'\n')," V(%d)\n",n_courbes+3) ;
          sprintf(strchr(ltmp,'\n')," q_SH(%d)\n",n_courbes+4) ;
          /*
          sprintf(strchr(ltmp,'\n')," q_Tob(%d)\n",n_courbes+5) ;
          sprintf(strchr(ltmp,'\n')," q_Jen(%d)\n",n_courbes+6) ;
          */
        } else if(ltmp[0] != '#') {
          double q_CH ;
          double q_CSH[3] ;
          double a_CSH[3] ;
          double x_CSH[3] ;
          double y_CSH[3] ;
          double z_CSH[3] ;
          double k_CSH[3] ;
          double v_CSH[3] ;
          int    j ;
          k_CSH[0] = k_SH   ; 
          x_CSH[0] = x_SH   ; y_CSH[0] = y_SH   ; z_CSH[0] = z_SH ;
          k_CSH[1] = k_Tob  ; 
          x_CSH[1] = x_Tob  ; y_CSH[1] = y_Tob  ; z_CSH[1] = z_Tob ;
          k_CSH[2] = k_Jen  ; 
          x_CSH[2] = x_Jen  ; y_CSH[2] = y_Jen  ; z_CSH[2] = z_Jen ;
          v_CSH[0] = v_SH ;
          v_CSH[1] = v_Tob ;
          v_CSH[2] = v_Jen ;
          sscanf(ltmp,FCOLX(colx),&q_CH) ;
          for(j=0;j<3;j++) {
            a_CSH[j] = pow(k_CH*q_CH,x_CSH[j])*pow(k_SH,y_CSH[j])/k_CSH[j] ;
          }
          q_SH = fraction_Silica(a_CSH,y_CSH,3,q_SH) ;
          for(j=0;j<3;j++) {
            q_CSH[j] = a_CSH[j]*pow(q_SH,y_CSH[j]) ;
          }
          {
            double x_m = 0,y_m = 0,z_m = 0 ;
            double v_m = 0 ;
            for(j=0;j<3;j++) {
              x_m += q_CSH[j]*x_CSH[j] ;
              y_m += q_CSH[j]*y_CSH[j] ;
              z_m += q_CSH[j]*z_CSH[j] ;
              v_m += q_CSH[j]*v_CSH[j] ;
            }
            x_m /= y_m ; z_m /= y_m ; v_m /= y_m ;
            sprintf(strchr(ltmp,'\n')," %e %e %e %e\n",x_m,z_m,v_m,q_SH) ;
          }
          /*
          for(j=0;j<3;j++) {
            sprintf(strchr(ltmp,'\n')," %e\n",q_CSH[j]) ;
          }
          */
        }
        fprintf(fic,"%s",ltmp) ;
      }
      /* n_courbes += 5 ; */
      n_courbes += 3 ;
      
    } else if(!strcmp(model,"CSH4Poles")){
      /* Portlandite */
      double k_CH   = 6.456e-6 ;
      /* Amorpheous Silica */
      double k_SH = 1.93642e-3 ;
      double x_SH = 0 ;
      double y_SH = 1 ;
      double z_SH = 0 ;
      double v_SH = 29.e-3 ;
      /* Jennite */
      double k_Jen0 = 4.4e-18 ,k_Jen = 2.39e-16 ;
      double x_Jen0 = 1.66666 ,x_Jen = 1.5     ;
      double                   y_Jen = 0.9     ;
      double z_Jen0 = 2.1     ,z_Jen = 2.4     ;
      double v_Jen0 = 78.e-3  ,v_Jen = 73.6e-3 ;
      /* Tobermorite */
      double k_Tob0 = 1.684e-12,k_Tob1 = 5.53e-29 ,k_Tob2 = 6.42e-22 ;
      double x_Tob0 = 0.83333  ,x_Tob1 = 2        ,x_Tob2 = 1.5 ;
      double                    y_Tob1 = 2.4      ,y_Tob2 = 1.8 ;
      double z_Tob0 = 1.83333  ,z_Tob1 = 4.4      ,z_Tob2 = 3.3 ;
      double v_Tob0 = 53.5e-3  ,v_Tob1 = 131.47e-3,v_Tob2 = 98.6e-3 ;
      
      double q_SH = 1. ;
      
      /*
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&k_CH,&k_SH,&k_Tob1,&k_Tob2,&k_Jen) ;
      */
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&y_Tob1,&y_Tob2,&y_Jen,&z_SH) ;
      
      /* Amorpheous silica */
      v_SH 	  = 29.e-3 + 7.e-3*z_SH;
      /* Jennite */
      x_Jen   = x_Jen0*y_Jen ;
      z_Jen   = z_Jen0*y_Jen ;
      k_Jen   = pow(k_Jen0,y_Jen) ;
      v_Jen   = v_Jen0*y_Jen ;
      /* Tobermorite 0 */
      v_Tob0 = x_Tob0/x_Jen0*v_Jen0 + (1 - x_Tob0/x_Jen0)*v_SH ;
      z_Tob0 = x_Tob0/x_Jen0*z_Jen0 + (1 - x_Tob0/x_Jen0)*z_SH ;
      /* Tobermorite I */
      x_Tob1  = x_Tob0*y_Tob1 ;
      z_Tob1  = z_Tob0*y_Tob1 ;
      k_Tob1  = pow(k_Tob0,y_Tob1) ;
      v_Tob1  = v_Tob0*y_Tob1 ;
      /* Tobermorite II */
      x_Tob2  = x_Tob0*y_Tob2 ;
      z_Tob2  = z_Tob0*y_Tob2 ;
      k_Tob2  = pow(k_Tob0,y_Tob2) ;
      v_Tob2  = v_Tob0*y_Tob2 ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(!strncmp(ltmp,"# Models",8)) {
          int j ;
          for(j=0;j<3;j++) sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n_courbes+2+j) ;
          fprintf(fic,"%s",ltmp) ;
          fgets(ltmp,sizeof(ltmp),ftmp) ;
          sprintf(strchr(ltmp,'\n')," z(%d)\n",n_courbes+2) ;
          sprintf(strchr(ltmp,'\n')," V(%d)\n",n_courbes+3) ;
          sprintf(strchr(ltmp,'\n')," q_SH(%d)\n",n_courbes+4) ;
          /*
          sprintf(strchr(ltmp,'\n')," q_Tob1(%d)\n",n_courbes+5) ;
          sprintf(strchr(ltmp,'\n')," q_Tob2(%d)\n",n_courbes+6) ;
          sprintf(strchr(ltmp,'\n')," q_Jen(%d)\n",n_courbes+7) ;
          */
        } else if(ltmp[0] != '#') {
          double q_CH ;
          double q_CSH[4] ;
          double a_CSH[4] ;
          double x_CSH[4] ;
          double y_CSH[4] ;
          double z_CSH[4] ;
          double k_CSH[4] ;
          double v_CSH[4] ;
          int    j ;
          k_CSH[0] = k_SH   ; 
          x_CSH[0] = x_SH   ; y_CSH[0] = y_SH   ; z_CSH[0] = z_SH ;
          k_CSH[1] = k_Tob1 ; 
          x_CSH[1] = x_Tob1 ; y_CSH[1] = y_Tob1 ; z_CSH[1] = z_Tob1 ;
          k_CSH[2] = k_Tob2 ; 
          x_CSH[2] = x_Tob2 ; y_CSH[2] = y_Tob2 ; z_CSH[2] = z_Tob2 ;
          k_CSH[3] = k_Jen  ; 
          x_CSH[3] = x_Jen  ; y_CSH[3] = y_Jen  ; z_CSH[3] = z_Jen ;
          v_CSH[0] = v_SH ;
          v_CSH[1] = v_Tob1 ;
          v_CSH[2] = v_Tob2 ;
          v_CSH[3] = v_Jen ;
          sscanf(ltmp,FCOLX(colx),&q_CH) ;
          for(j = 0 ; j < 4 ; j++) {
            a_CSH[j] = pow(k_CH*q_CH,x_CSH[j])*pow(k_SH,y_CSH[j])/k_CSH[j] ;
          }
          q_SH = fraction_Silica(a_CSH,y_CSH,4,q_SH) ;
          for(j = 0 ; j < 4 ; j++) {
            q_CSH[j] = a_CSH[j]*pow(q_SH,y_CSH[j]) ;
          }
          {
            double x_m = 0,y_m = 0,z_m = 0 ;
            double v_m = 0 ;
            for(j = 0 ; j < 4 ; j++) {
              x_m += q_CSH[j]*x_CSH[j] ;
              y_m += q_CSH[j]*y_CSH[j] ;
              z_m += q_CSH[j]*z_CSH[j] ;
              v_m += q_CSH[j]*v_CSH[j] ;
            }
            x_m /= y_m ; z_m /= y_m ; v_m /= y_m ;
            sprintf(strchr(ltmp,'\n')," %e %e %e %e\n",x_m,z_m,v_m,q_SH) ;
          }
          /*
          for(j=1;j<4;j++) {
            sprintf(strchr(ltmp,'\n')," %e\n",q_CSH[j]) ;
          }
          */
        }
        fprintf(fic,"%s",ltmp) ;
      }
      /* n_courbes += 6 ; */
      n_courbes += 3 ;
      
    } else if(!strcmp(model,"CSHnPoles")){
      /* Portlandite */
      double k_CH   = 6.456e-6 ;
      /* Amorpheous Silica */
      double k_SH = 1.93642e-3 ;
      double x_SH = 0 ;
      double z_SH = 2 ;
      double v_SH = 43.e-3     ;
      /* Jennite */
      double k_Jen = 4.39e-18  ;
      double x_Jen = 1.66666   ;
      double z_Jen = x_Jen     ;
      double v_Jen = 81.8e-3   ;
      /* Tobermorite */
      double k_Tob = 1.684e-12 ;
      double x_Tob = 0.83333   ;
      double z_Tob = 1.        ;
      double v_Tob = 54.8e-3   ;
      
      double q_SH = 1. ;
      int    n_poles = 0 ;
      
      /*
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&k_CH,&k_SH,&k_Tob1,&k_Tob2,&k_Jen) ;
      */
      sscanf(line,"{ %*s = %d }",&n_poles) ;
      
      if(n_poles < 2 || n_poles > 3) arret("Curves_WriteCurves : CSHnpoles") ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(!strncmp(ltmp,"# Models",8)) {
          int j ;
          for(j=0;j<3;j++) sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n_courbes+2+j) ;
          fprintf(fic,"%s",ltmp) ;
          fgets(ltmp,sizeof(ltmp),ftmp) ;
          sprintf(strchr(ltmp,'\n')," z(%d)\n",n_courbes+2) ;
          sprintf(strchr(ltmp,'\n')," V(%d)\n",n_courbes+3) ;
          sprintf(strchr(ltmp,'\n')," q_SH(%d)\n",n_courbes+4) ;
        } else if(ltmp[0] != '#') {
          double q_CH ;
          double q_CSH[3] ;
          double a_CSH[3] ;
          double x_CSH[3] ;
          double z_CSH[3] ;
          double k_CSH[3] ;
          double v_CSH[3] ;
          int    j ;
          double sum_a = 0 ;
          k_CSH[0] = k_SH   ; 
          x_CSH[0] = x_SH   ; z_CSH[0] = z_SH ;
          k_CSH[1] = k_Jen  ; 
          x_CSH[1] = x_Jen  ; z_CSH[1] = z_Jen ;
          k_CSH[2] = k_Tob  ; 
          x_CSH[2] = x_Tob  ; z_CSH[2] = z_Tob ;
          v_CSH[0] = v_SH ;
          v_CSH[1] = v_Jen ;
          v_CSH[2] = v_Tob ;
          sscanf(ltmp,FCOLX(colx),&q_CH) ;
          for(j=0;j<n_poles;j++) {
            a_CSH[j] = pow(k_CH*q_CH,x_CSH[j])*k_SH/k_CSH[j] ;
            sum_a += a_CSH[j] ;
          }
          /* q_SH = fraction_Silica(a_CSH,y_CSH,n_poles,q_SH) ; */
          q_SH = 1/sum_a ;
          for(j=0;j<n_poles;j++) {
            q_CSH[j] = a_CSH[j]*q_SH ;
          }
          {
            double x_m = 0,z_m = 0 ;
            double v_m = 0 ;
            for(j=0;j<n_poles;j++) {
              x_m += q_CSH[j]*x_CSH[j] ;
              z_m += q_CSH[j]*z_CSH[j] ;
              v_m += q_CSH[j]*v_CSH[j] ;
            }
            sprintf(strchr(ltmp,'\n')," %e %e %e %e\n",x_m,z_m,v_m,q_SH) ;
          }
        }
        fprintf(fic,"%s",ltmp) ;
      }
      n_courbes += 3 ;
      
    } else if(!strcmp(model,"fCSH4Poles")){
      /* Portlandite */
      double k_CH   = 6.456e-6 ;
      /* Amorpheous Silica */
      double k_SH = 1.93642e-3 ;
      double x_SH = 0 ;
      double y_SH = 1 ;
      double z_SH = 2 ;
      double v_SH = 43.e-3 ;
      /* Jennite */
      double k_Jen0 = 4.39e-18,k_Jen = 2.39e-16 ;
      double x_Jen0 = 1.66666 ,x_Jen = 1.5     ;
      double                   y_Jen = 0.9     ;
      double z_Jen0 = 2.66666 ,z_Jen = 2.4     ;
      double v_Jen0 = 81.8e-3 ,v_Jen = 73.6e-3 ;
      /* Tobermorite */
      double k_Tob0 = 1.684e-12,k_Tob1 = 5.53e-29 ,k_Tob2 = 6.42e-22 ;
      double x_Tob0 = 0.83333  ,x_Tob1 = 2        ,x_Tob2 = 1.5 ;
      double                    y_Tob1 = 2.4      ,y_Tob2 = 1.8 ;
      double z_Tob0 = 1.83333  ,z_Tob1 = 4.4      ,z_Tob2 = 3.3 ;
      double v_Tob0 = 54.8e-3  ,v_Tob1 = 131.47e-3,v_Tob2 = 98.6e-3 ;
      
      double q_SH = 1 ;
      double q_CH = 0 ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf }",&y_Tob1,&y_Tob2,&y_Jen) ;
      
      x_Jen   = x_Jen0*y_Jen ;
      z_Jen   = z_Jen0*y_Jen ;
      k_Jen   = pow(k_Jen0,y_Jen) ;
      v_Jen   = v_Jen0*y_Jen ;
      x_Tob1  = x_Tob0*y_Tob1 ;
      z_Tob1  = z_Tob0*y_Tob1 ;
      k_Tob1  = pow(k_Tob0,y_Tob1) ;
      v_Tob1  = v_Tob0*y_Tob1 ;
      x_Tob2  = x_Tob0*y_Tob2 ;
      z_Tob2  = z_Tob0*y_Tob2 ;
      k_Tob2  = pow(k_Tob0,y_Tob2) ;
      v_Tob2  = v_Tob0*y_Tob2 ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(!strncmp(ltmp,"# Models",8)) {
          int j ;
          for(j=0;j<3;j++) sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n_courbes+2+j) ;
          fprintf(fic,"%s",ltmp) ;
          fgets(ltmp,sizeof(ltmp),ftmp) ;
          sprintf(strchr(ltmp,'\n')," z(%d)\n",n_courbes+2) ;
          sprintf(strchr(ltmp,'\n')," V(%d)\n",n_courbes+3) ;
          sprintf(strchr(ltmp,'\n')," g(%d)\n",n_courbes+4) ;
        } else if(ltmp[0] != '#') {
          double q_CSH[4] ;
          double a_CSH[4] ;
          double x_CSH[4] ;
          double y_CSH[4] ;
          double z_CSH[4] ;
          double k_CSH[4] ;
          double v_CSH[4] ;
          double x_m,x_m1 = 0 ;
          int    j ;
          k_CSH[0] = k_SH   ; 
          x_CSH[0] = x_SH   ; y_CSH[0] = y_SH   ; z_CSH[0] = z_SH ;
          k_CSH[1] = k_Tob1 ; 
          x_CSH[1] = x_Tob1 ; y_CSH[1] = y_Tob1 ; z_CSH[1] = z_Tob1 ;
          k_CSH[2] = k_Tob2 ; 
          x_CSH[2] = x_Tob2 ; y_CSH[2] = y_Tob2 ; z_CSH[2] = z_Tob2 ;
          k_CSH[3] = k_Jen  ; 
          x_CSH[3] = x_Jen  ; y_CSH[3] = y_Jen  ; z_CSH[3] = z_Jen ;
          v_CSH[0] = v_SH ;
          v_CSH[1] = v_Tob1 ;
          v_CSH[2] = v_Tob2 ;
          v_CSH[3] = v_Jen ;
          sscanf(ltmp,FCOLX(colx),&x_m) ;
          
          do {
            double q_CHn = q_CH ;
            double q_SHn = q_SH ;
            double dq_CH = (q_CH < 1.e-8) ? 1.e-8 : 1.e-5*q_CH ;
            double dq_SH ;
            double x_max ;
            
            q_CH = q_CHn + dq_CH ;
            
            if(q_CH < 1) {
              for(j=0;j<4;j++) {
                a_CSH[j] = pow(k_CH*q_CH,x_CSH[j])*pow(k_SH,y_CSH[j])/k_CSH[j] ;
              }
              q_SH  = fraction_Silica(a_CSH,y_CSH,4,q_SHn) ;
              dq_SH = q_SH - q_SHn ;
              x_m1  = - dq_SH/q_SH*q_CH/dq_CH ;
              x_max = x_m1 ;

            } else {
              /* Arbitraire */
              for(j=0;j<4;j++) {
                a_CSH[j] = pow(k_CH/q_CH,x_CSH[j])*pow(k_SH,y_CSH[j])/k_CSH[j] ;
              }
              q_SH  = fraction_Silica(a_CSH,y_CSH,4,q_SHn) ;
              dq_SH = q_SH - q_SHn ;
              x_m1  = 2*x_max - dq_SH/q_SH*q_CH/dq_CH ;
            
            }
            
          } while(x_m1 < x_m) ;
          
          for(j=0;j<4;j++) {
            q_CSH[j] = a_CSH[j]*pow(q_SH,y_CSH[j]) ;
          }
          
          {
            double y_m = 0,z_m = 0 ;
            double v_m = 0 ;
            double dgsdx = log(q_CH) ;
            double g = x_m*dgsdx + log(q_SH) ;
            for(j=0;j<4;j++) {
              y_m += q_CSH[j]*y_CSH[j] ;
              z_m += q_CSH[j]*z_CSH[j] ;
              v_m += q_CSH[j]*v_CSH[j] ;
            }
            sprintf(strchr(ltmp,'\n')," %e %e %e %e\n",dgsdx,z_m/y_m,v_m/y_m,g) ;
          }
        }
        fprintf(fic,"%s",ltmp) ;
      }
      /* n_courbes += 6 ; */
      n_courbes += 3 ;
      
    } else if(!strcmp(model,"Redlich-Kwong_CO2")) {
      double temperature ;
      
      sscanf(line,"{ %*s = %lf }",&temperature) ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double p,rho ;
          sscanf(ltmp,FCOLX(colx),&p) ;
          rho = MolarDensityOfCO2_RedlichKwong(p,temperature) ;
          sprintf(strchr(ltmp,'\n')," %e\n",rho) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else if(!strcmp(model,"Fenghour_CO2")) {
      double temperature ;
      
      sscanf(line,"{ %*s = %lf }",&temperature) ;

      while(fgets(ltmp,sizeof(ltmp),ftmp)) {
        if(ltmp[0] != '#') {
          double p,mu ;
          sscanf(ltmp,FCOLX(colx),&p) ;
          mu = ViscosityOfCO2_Fenghour(p,temperature) ;
          sprintf(strchr(ltmp,'\n')," %e\n",mu) ;
        }
        fprintf(fic,"%s",ltmp) ;
      }
      
    } else {
      fclose(ftmp) ;
      fclose(fic) ;
      arret("Curves_WriteCurves : fonction non connue") ;
      return(n_courbes) ;
    }
    fclose(ftmp) ;
    fclose(fic) ;
    line = strchr(line,'}') + 1 ;
  }
  return(n_courbes) ;
#undef MAX_COLX
#undef FCOLX
}

#endif
