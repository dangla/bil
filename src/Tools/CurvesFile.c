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
#include "Curves.h"
#include "CurvesFile.h"
#include "InternationalSystemOfUnits.h"


/* Shorthands of some units */
#define dm    (0.1*InternationalSystemOfUnits_OneMeter)
#define cm    (0.01*InternationalSystemOfUnits_OneMeter)
#define dm3   (dm*dm*dm)
#define cm3   (cm*cm*cm)



typedef double GenericFunction_t(double,va_list) ;

/*
static void  (CurvesFile_StoreFilePosition)(CurvesFile_t*) ;
*/
static void  (CurvesFile_PasteXaxisColumn)(CurvesFile_t*,const char*,double,double) ;
static void  (CurvesFile_PasteYaxisColumn)(CurvesFile_t*,const char*,const char*,int,GenericFunction_t*,...) ;



static GenericFunction_t Langmuir ;
static GenericFunction_t LangmuirN ;
static GenericFunction_t Freundlich ;
static GenericFunction_t MualemVanGenuchten ;
static GenericFunction_t VanGenuchten ;
static GenericFunction_t NavGenuchten ;
static GenericFunction_t Mualem_dry ;
static GenericFunction_t Mualem_gas ;
static GenericFunction_t Millington ;
static GenericFunction_t MonlouisBonnaire ;
static GenericFunction_t CSH3EndMembers ;
static GenericFunction_t CSHLangmuirN ;
static GenericFunction_t RedlichKwongCO2 ;
static GenericFunction_t FenghourCO2 ;
static GenericFunction_t FromCurve ; 
static GenericFunction_t Evaluate ; 
static GenericFunction_t Expressions ; 
static GenericFunction_t MolarDensityOfPerfectGas ;
static GenericFunction_t MolarDensityOfCO2 ;
static GenericFunction_t Affine ;
static GenericFunction_t VanGenuchten_gas ;


static double vangenuchten(double,double) ;
static double fraction_Silica(double*,double*,int,double) ;
static double MolarDensityOfCO2_RedlichKwong(double,double) ;
static double ViscosityOfCO2_Fenghour(double,double) ;
static double (langmuir)(double,double,double,double) ;



CurvesFile_t*   (CurvesFile_Create)(void)
{
  CurvesFile_t* curvesfile   = (CurvesFile_t*) malloc(sizeof(CurvesFile_t)) ;
  
  if(!curvesfile) arret("CurvesFile_Create(0)") ;
  

  /* Memory space for textfile */
  {
    TextFile_t* textfile = TextFile_Create(NULL) ;
    
    if(!textfile) {
      arret("CurvesFile_Create(1)") ;
    }
    
    CurvesFile_GetTextFile(curvesfile) = textfile ;
  }
  
  
  /* Memory space for the file position */
  /*
  {
    fpos_t* pos = (fpos_t*) malloc(sizeof(fpos_t)) ;
    
    if(!pos) {
      arret("CurvesFile_Create(3)") ;
    }
    CurvesFile_GetFilePositionStartingInputData(curvesfile) = pos ;
  }
  */
  
  
  /* Memory space for the text line to be read */
  {
    int n = CurvesFile_MaxLengthOfTextLine ;
    char *line = (char*) malloc(n*sizeof(char)) ;
    
    if(!line) {
      arret("CurvesFile_Create(4)") ;
    }
    
    CurvesFile_GetTextLine(curvesfile) = line ;
  }
  
  
  /* Initialization */
  CurvesFile_GetScaleType(curvesfile) = 'n' ;
  CurvesFile_GetNbOfCurves(curvesfile) = 0 ;
  CurvesFile_GetNbOfPoints(curvesfile) = 0 ;
  CurvesFile_GetCommandLine(curvesfile) = NULL ;
  CurvesFile_GetCurrentPositionInTheCommandLine(curvesfile) = NULL ;
  
    
  /* Space allocation for buffer */
  CurvesFile_GetBuffer(curvesfile) = Buffer_Create(CurvesFile_SizeOfBuffer) ;
  
  
  return(curvesfile) ;
}


void (CurvesFile_Delete)(CurvesFile_t** curvesfile)
{
  TextFile_Delete(&CurvesFile_GetTextFile(*curvesfile)) ;
  /* free(CurvesFile_GetFileName(*curvesfile)) ; */
  /* free(CurvesFile_GetFilePositionStartingInputData(*curvesfile)) ; */
  free(CurvesFile_GetTextLine(*curvesfile)) ;
  Buffer_Delete(&CurvesFile_GetBuffer(*curvesfile)) ;
  free(*curvesfile) ;
}


#ifdef NOTDEFINED
void (CurvesFile_MoveToFilePositionStartingInputData)(CurvesFile_t* curvesfile)
/** Set the file position of the stream to the beginning of the input data. */
{
  FILE *str  = CurvesFile_GetFileStream(curvesfile) ;
  fpos_t* pos = CurvesFile_GetFilePositionStartingInputData(curvesfile) ;
  
  /* Set the file position of the stream to the stored position */
  if(fsetpos(str,pos)) {
    arret("CurvesFile_MoveToFilePositionStartingInputData") ;
  }
}
#endif


int   (CurvesFile_Initialize)(CurvesFile_t* curvesfile,const char* cmdline)
/** Initialize the curvesfile from the command line "cmdline" 
 *  Return the nb of curves found in the file */
{
  /* The command line */
  CurvesFile_GetCommandLine(curvesfile) = cmdline ;
  
  /* Read the scale */
  {
    char  scale = 'n' ;
    const char* line = cmdline ;
    
    if(strstr(line,"_log")) scale = 'l' ;
    
    CurvesFile_GetScaleType(curvesfile) = scale ;
  }


  /* Read the file name */
  {
    char* filename = CurvesFile_GetFileName(curvesfile) ;

    {
      const char* line = strchr(cmdline,'=') + 1 ;
      sscanf(line," %s",filename) ;
    
      if(strlen(filename) > CurvesFile_MaxLengthOfFileName) {
        arret("CurvesFile_Initialize(1)") ;
      }
      
      line  = strstr(line,filename) + strlen(filename) ;
      CurvesFile_GetCurrentPositionInTheCommandLine(curvesfile) = line ;
    }
  }
  
  
  /* Does this file exist? If not return 0 */
  {
    if(CurvesFile_DoesNotExist(curvesfile)) {
      return(0) ;
    }
  }
  

  /* Position starting the input data */
  /*
  {
    char *line ;
    
    CurvesFile_OpenFile(curvesfile,"r") ;
    
    do {
      CurvesFile_StoreFilePosition(curvesfile) ;
      
      line = CurvesFile_ReadLineFromCurrentFilePosition(curvesfile) ;
      
    } while((line) && (line[0] == '#')) ;
      
    {
      fpos_t* pos = CurvesFile_GetFilePositionStartingInputData(curvesfile) ;

      *pos = *CurvesFile_GetFilePosition(curvesfile) ;
    }
    
    CurvesFile_CloseFile(curvesfile) ;
  }
  */
  
  

  /* Nb of curves */
  {
    CurvesFile_OpenFile(curvesfile,"r") ;
  
    {
      int n_curves = 0 ;
      char *line ;
      
      /* We skip the commented lines */
      do {
        line = CurvesFile_ReadLineFromCurrentFilePosition(curvesfile) ;
      } while((line) && (line[0] == '#')) ;
      
      /* We count the nb of curves in the first non-commented line */
      if(line) {
        char FS[] = CurvesFile_FieldDelimiters"\n" ;
        
        strtok(line,FS) ;
      
        while((char *) strtok(NULL,FS) != NULL) n_curves++ ;
      }
    
      CurvesFile_GetNbOfCurves(curvesfile) = n_curves ;
    }
    
    CurvesFile_CloseFile(curvesfile) ;
  }


  /* Nb of points */
  {
    CurvesFile_OpenFile(curvesfile,"r") ;
    
    {
      int n_points = 0 ;
      char *line ;
      
      do {
        line = CurvesFile_ReadLineFromCurrentFilePosition(curvesfile) ;
        
        /* We count the nb of non-commented lines */
        if(line != NULL && line[0] != '#') n_points++ ;
        
      } while(line) ;
    
      CurvesFile_GetNbOfPoints(curvesfile) = n_points ;
    }
    
    CurvesFile_CloseFile(curvesfile) ;
  }
  
  return(CurvesFile_GetNbOfCurves(curvesfile)) ;
}




int   (CurvesFile_WriteCurves)(CurvesFile_t* curvesfile)
/** Write curves in the file stored in curvesfile as discrete data
 *  Return the number of curves that has been writen */
{
  const char*  cmdline = CurvesFile_GetCommandLine(curvesfile) ;
  char   models[CurvesFile_MaxNbOfCurves+1][CurvesFile_MaxLengthOfKeyWord] ;
  char   labels[CurvesFile_MaxNbOfCurves+1][CurvesFile_MaxLengthOfKeyWord] ;
  const char*  line  = CurvesFile_GetCurrentPositionInTheCommandLine(curvesfile) ;
  char*  xmodel = models[0] ;
  char*  xlabel = labels[0] ;
  //char*  ymodel = models[1] ;
  //char*  ylabel = labels[1] ;
  int    acol ;
  int    ycol ;


  /* Write the first column x-axis */
  
  /* Read the range and the nb of points */
  if(sscanf(line," %s = %[^({]",xlabel,xmodel) == 2) {
      
    line = strchr(line,'{') ;
      
    if(!strcmp(xmodel,"Range")) {
      int    n_points ;
      double x_1,x_2 ;
        
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %d }",&x_1,&x_2,&n_points) ;
      
      CurvesFile_GetNbOfPoints(curvesfile) = n_points ;

      CurvesFile_PasteXaxisColumn(curvesfile,xlabel,x_1,x_2) ;
        
    } else {
      arret("CurvesFile_WriteCurves(1): The first key-word is not \"Range\"") ;
      return(0) ;
    }
      
    line = strchr(line,'}') + 1 ;
      
  } else {
    char* filename = CurvesFile_GetFileName(curvesfile) ;
    
    arret("CurvesFile_WriteCurves(2): no data to build \"%s\"",filename) ;
    return(0) ;
  }


  /* Write other columns y-axis */
  
#define YLABEL  (labels[ycol - 1])
#define YMODEL  (models[ycol - 1])
#define ALABEL  (labels[acol - 1])
    
  /* We initialize for the first next columns */
  ycol    = 2 ;
  //ylabel  = labels[ycol - 1] ;
  //ymodel  = models[ycol - 1] ;
  
  /* Read the key-word of the curve */
  while(sscanf(line," %s = %[^(](%d)",YLABEL,YMODEL,&acol) == 3) {
    //char* alabel = labels[acol - 1] ;
    
#define PasteColumn(a,...) \
        CurvesFile_PasteYaxisColumn(curvesfile,YLABEL,#a,acol,a,__VA_ARGS__)

    
    line = strchr(line,'{') ;

    /* Write a new column depending on the model name */
    if(!strcmp(YMODEL,"Freundlich")) {
      double alpha,beta ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&alpha,&beta) ;

      PasteColumn(Freundlich,alpha,beta) ;
      
    } else if(!strcmp(YMODEL,"Affine")) {
      double a,b ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&a,&b) ;

      PasteColumn(Affine,a,b) ;
      
    } else if(!strcmp(YMODEL,"Langmuir")) {
      double c0,ca_max ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&ca_max,&c0) ;

      PasteColumn(Langmuir,ca_max,c0) ;
      
    } else if(!strcmp(YMODEL,"LangmuirN")) {
      double n,x0,y0 ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf }",&y0,&x0,&n) ;

      PasteColumn(LangmuirN,y0,x0,n) ;
      
    } else if(!strcmp(YMODEL,"Mualem_wet") || !strcmp(YMODEL,"Mualem_liq")){
      double m ;
      
      sscanf(line,"{ %*s = %lf }",&m) ;
      
      PasteColumn(MualemVanGenuchten,m) ;
      
    } else if(!strcmp(YMODEL,"Mualem_dry")){
      double m_w,m_d,a_d,a_w ;
      FILE   *ftmp1 = tmpfile() ;
      double kh_max ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf , %*s = %lf }",&a_w,&m_w,&a_d,&m_d) ;

      /* Compute kh in the temporary file ftmp1 */
      {
#define MAX_COLX  5
#define FCOLX(i) (frmt + 4*(MAX_COLX + 1 - i))
        char   frmt[] = "%*s %*s %*s %*s %*s %lf" ;
        char   ltmp[CurvesFile_MaxLengthOfTextLine] ;
        double kh = 0. ;
        double p ;
        FILE *ftmp = CurvesFile_FileStreamCopy(curvesfile) ;
        
        do {
          fgets(ltmp,sizeof(ltmp),ftmp) ;
        } while(ltmp[0] == '#') ;
	
        sscanf(ltmp,FCOLX(acol),&p) ;
        fprintf(ftmp1,"%e %e\n",p,kh) ;
	
        while(fgets(ltmp,sizeof(ltmp),ftmp)) {
          if(ltmp[0] != '#') {
            double pn = p ;
            
            sscanf(ltmp,FCOLX(acol),&p) ;
            
            {
              double dp  = p - pn ;
              double s_d = vangenuchten(pn/a_d,m_d) ;
              double s_w = vangenuchten(pn/a_w,m_w) ;
              double hn = (s_d - s_w)/(1. - s_w) ;
              double h,dh ;
	      
              s_d = vangenuchten(p/a_d,m_d) ;
              s_w = vangenuchten(p/a_w,m_w) ;
              h   = (s_d - s_w)/(1. - s_w) ;
	      
              dh  = h - hn ;
              kh += dh/(p - 0.5*dp) ;
	      
              fprintf(ftmp1,"%e %e\n",p,kh) ;
            }
          }
        }
        
        kh_max = kh ;
        fclose(ftmp) ;
#undef FCOLX
#undef MAX_COLX
      }
      
      rewind(ftmp1) ;
      
      PasteColumn(Mualem_dry,a_w,m_w,a_d,m_d,kh_max,ftmp1) ;
      
      fclose(ftmp1) ;
      
    } else if(!strcmp(YMODEL,"Mualem_gas")){
      double m ;
      
      sscanf(line,"{ %*s = %lf }",&m) ;
      
      PasteColumn(Mualem_gas,m) ;
      
    } else if(!strcmp(YMODEL,"Van-Genuchten_gas")){
      double m,p,q ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf }",&m,&p,&q) ;
      
      PasteColumn(VanGenuchten_gas,m,p,q) ;
      
    } else if(!strcmp(YMODEL,"Van-Genuchten")){
      double a,m ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&a,&m) ;
      
      PasteColumn(VanGenuchten,a,m) ;
      
    } else if(!strcmp(YMODEL,"Nav-Genuchten")){
      double a,m ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&a,&m) ;
      
      PasteColumn(NavGenuchten,a,m) ;
      
    } else if(!strcmp(YMODEL,"Millington")){
      double b ;
      
      sscanf(line,"{ %*s = %lf }",&b) ;
      
      PasteColumn(Millington,b) ;
      
    } else if(!strcmp(YMODEL,"Monlouis-Bonnaire")){
      double m ;
      
      sscanf(line,"{ %*s = %lf }",&m) ;
      
      PasteColumn(MonlouisBonnaire,m) ;
      
    } else if(!strcmp(YMODEL,"Integral")){
      int    coly ;
      
      sscanf(line,"{ %*s = %d }",&coly) ;
      
      {
        int n_curves = CurvesFile_GetNbOfCurves(curvesfile) ;
        Curves_t* crvs = Curves_Create(n_curves) ;
        /* We read the existing curves in the now created file */
        int n_crvs = Curves_ReadCurves(crvs,cmdline) ;
        /* Curve to be integrated */
        int i = coly - 2 ;
        Curve_t* crvi = Curves_GetCurve(crvs) + i ;
        /* We create the integral */
        Curve_t* crvj = Curve_CreateIntegral(crvi) ;
        
        PasteColumn(FromCurve,crvj) ;
        
        /* Free memory */
        Curve_Delete(&crvj) ;
        Curves_Delete(&crvs) ;
      }
      
    } else if(!strcmp(YMODEL,"Evaluate")){
      char    expr[CurvesFile_MaxLengthOfTextLine] ;
      
      {
        const char* c = strchr(line,'{') ;
        
        strcpy(expr,c + 1) ;
      }
      
      {
        char* c = strchr(expr,'}') ;
        
        *c = '\0' ;
      }

      {
        Curves_t* curves = CurvesFile_GetCurves(curvesfile) ;
        
        PasteColumn(Evaluate,curves,expr,ALABEL) ;
      }
      
    } else if(!strcmp(YMODEL,"Expressions")){
      char    expr[CurvesFile_MaxLengthOfTextLine] ;
      
      {
        const char* c = strchr(line,'{') ;
        
        strcpy(expr,c + 1) ;
      }
      {
        char* c = strchr(expr,'}') ;
        
        *c = '\0' ;
      }
      
      {
        Curves_t* curves = CurvesFile_GetCurves(curvesfile) ;
        
        PasteColumn(Expressions,curves,expr,YLABEL,ALABEL) ;
      }
      
    } else if(!strcmp(YMODEL,"CSH3Poles")){
      double  y_Jen = 0.9 ;
      double  y_Tob = 1.8 ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf }",&y_Tob,&y_Jen) ;
      
      if(!strncmp(YLABEL,"X_CSH",5)) {
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"X_CSH") ;
      } else if(!strncmp(YLABEL,"Z_CSH",5)) {
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"Z_CSH") ;
      } else if(!strncmp(YLABEL,"V_CSH",5)) {
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"V_CSH") ;
      } else if(!strncmp(YLABEL,"S_SH",4)) {
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"S_SH") ;
      } else {
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"X_CSH") ;
        ycol += 1 ;
        strcpy(YLABEL,"Z_CSH") ;
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"Z_CSH") ;
        ycol += 1 ;
        strcpy(YLABEL,"V_CSH") ;
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"V_CSH") ;
        ycol += 1 ;
        strcpy(YLABEL,"S_SH") ;
        PasteColumn(CSH3EndMembers,y_Tob,y_Jen,"S_SH") ;
      }
      
    } else if(!strcmp(YMODEL,"CSHLangmuirN")){
      double  x1,s1,n1 ;
      double  x2,s2,n2 ;
      
      sscanf(line,"{ %*s = %lf , %*s = %lf , %*s = %lf \
                   , %*s = %lf , %*s = %lf , %*s = %lf }" \
                   , &x1,&s1,&n1,&x2,&s2,&n2) ;
      
      if(!strncmp(YLABEL,"X_CSH",1)) {
        PasteColumn(CSHLangmuirN,x1,s1,n1,x2,s2,n2,"X_CSH") ;
      } else if(!strncmp(YLABEL,"S_SH",1)) {
        PasteColumn(CSHLangmuirN,x1,s1,n1,x2,s2,n2,"S_SH") ;
      }
      
    } else if(!strcmp(YMODEL,"Redlich-Kwong_CO2")) {
      double temperature ;
      
      sscanf(line,"{ %*s = %lf }",&temperature) ;
      
      PasteColumn(RedlichKwongCO2,temperature) ;
      
    } else if(!strcmp(YMODEL,"MolarDensityOfCO2")) {
      double temperature ;
      
      sscanf(line,"{ %*s = %lf }",&temperature) ;
      
      PasteColumn(MolarDensityOfCO2,temperature) ;
      
    } else if(!strcmp(YMODEL,"MolarDensityOfPerfectGas")) {
      double temperature ;
      
      sscanf(line,"{ %*s = %lf }",&temperature) ;
      
      PasteColumn(MolarDensityOfPerfectGas,temperature) ;
      
    } else if(!strcmp(YMODEL,"Fenghour_CO2") || \
              !strcmp(YMODEL,"ViscosityOfCO2")) {
      double temperature ;
      
      sscanf(line,"{ %*s = %lf }",&temperature) ;
      
      PasteColumn(FenghourCO2,temperature) ;
      
    } else {
      arret("CurvesFile_WriteCurves : fonction non connue") ;
  
    }
    
    line = strchr(line,'}') + 1 ;
    
    /* Increment for the next column */
    ycol += 1 ;
    if(ycol > CurvesFile_MaxNbOfCurves) {
      arret("CurvesFile_WriteCurves: too many curves") ;
    }
    //ylabel  = labels[ycol - 1] ;
    //ymodel  = models[ycol - 1] ;
    
#undef PasteColumn
  }

#undef YLABEL
#undef YMODEL
#undef ALABEL
  
  return(CurvesFile_GetNbOfCurves(curvesfile)) ;
}



/* Intern functions */

void (CurvesFile_PasteXaxisColumn)(CurvesFile_t* curvesfile, const char *xlabel,double x_1,double x_2)
{
  FILE *fic = CurvesFile_OpenFile(curvesfile,"w") ;
  
  if(!fic) arret("CurvesFile_PasteXaxisColumn(1)") ;

  /* The first lines must be commented. At least two commented lines:
   * # Models: The name of models used to compute the columns 
   * # Labels: The labels of the columns */
  fprintf(fic,"# Models: X-axis(1)\n") ;
  fprintf(fic,"# Labels: %s(1)\n",xlabel) ;

  {
    char  scale = CurvesFile_GetScaleType(curvesfile) ;
    int n_points = CurvesFile_GetNbOfPoints(curvesfile) ;
    int i ;
    
    if(scale == 'l') {
      if(x_1 <= 0. || x_2 <= 0.) {
        arret("CurvesFile_PasteXaxisColumn(2)") ;
      }
    }
    
    /* Write the column x-axis */
    for(i = 0 ; i < n_points ; i++){
      int    n1 = n_points - 1 ;
      double n  = ((double) i)/n1 ;
      double x ;
      
      if(scale == 'l') {
        double ratio = x_2/x_1 ;
        
        x = x_1*pow(ratio,n) ;
        
      } else {
        double dx = (x_2 - x_1) ;
        
        x = x_1 + n*dx ;
      }
      
      fprintf(fic,"%e\n",x) ;   
    }
  }
  
  CurvesFile_CloseFile(curvesfile) ;
}


void (CurvesFile_PasteYaxisColumn)(CurvesFile_t* curvesfile,const char* ylabel,const char* model,int acol,GenericFunction_t* fct, ...)
/** Paste a new y-axis column in file "filename" by using the function "fct".
 *  acol   = is the column index to be used as the argument, ie y = fct(x) 
 *  model  = is the name of the model used as "fct"
 *  ylabel = is the label of the y-axis */
{
#define MAX_COLX  5
#define FCOLX(i) (frmt + 4*(MAX_COLX + 1 - i))
  char   frmt[] = "%*s %*s %*s %*s %*s %lf" ;
  int n_curves = CurvesFile_GetNbOfCurves(curvesfile) ;
    
  if(acol > MAX_COLX) {
    arret("CurvesFile_PasteYaxisColumn(1)") ;
  }

  {
    char   ltmp[CurvesFile_MaxLengthOfTextLine] ;
    FILE *ftmp = CurvesFile_FileStreamCopy(curvesfile) ;
    FILE *fic = CurvesFile_OpenFile(curvesfile,"w") ;
    
    if(!fic) {
      arret("CurvesFile_PasteYaxisColumn(2)") ;
    }
  
    while(fgets(ltmp,sizeof(ltmp),ftmp)) {
      
      if(ltmp[0] == '#') {
        /* Write the name of the model: "model" */
        if(!strncmp(ltmp,"# Models",8)) {
          char *c = strrchr(ltmp,'(') ;
          int n ;
          
          /* We read the number of columns */
          sscanf(c,"(%d)",&n) ;
          
          /* There is n-1 existing curves, so we increment */
          n_curves = n ;
    
          if(acol > n_curves) {
            arret("CurvesFile_PasteYaxisColumn(3)") ;
          }
          
          /* The index of this column is n + 1 */
          sprintf(strchr(ltmp,'\n')," %s(%d)\n",model,n + 1) ;
        
        /* Write the label of the y-axis: "ylabel" */
        } else if(!strncmp(ltmp,"# Labels",7)) {
          char *c = strrchr(ltmp,'(') ;
          int n ;
          
          sscanf(c,"(%d)",&n) ;
          
          n_curves = n ;
    
          if(acol > n_curves) {
            arret("CurvesFile_PasteYaxisColumn(4)") ;
          }
          
          sprintf(strchr(ltmp,'\n')," %s(%d)\n",ylabel,n + 1) ;
        }
        
      } else {
        double x,y ;
        va_list args ;

        va_start(args,fct) ;
      
        sscanf(ltmp,FCOLX(acol),&x) ;
        y = fct(x,args) ;
        sprintf(strchr(ltmp,'\n')," %e\n",y) ;
  
        va_end(args) ;
      }
      
      fprintf(fic,"%s",ltmp) ;
    }
    
    CurvesFile_CloseFile(curvesfile) ;
    fclose(ftmp) ;
  }
  
  CurvesFile_GetNbOfCurves(curvesfile) += 1 ;
  
  return ;

#undef FCOLX
#undef MAX_COLX
}






/* Predefined functions */
double (FromCurve)(double x,va_list args)
{
  Curve_t* curve = va_arg(args,Curve_t*) ;
  double y = Curve_ComputeValue(curve,x) ;
  
  return(y) ;
}


double (Evaluate)(double x,va_list args)
{
  Curves_t* curves = va_arg(args,Curves_t*) ;
  char*  expr = va_arg(args,char*) ;
  char*  xlabel = va_arg(args,char*) ;
  int NbOfCurves = Curves_GetNbOfCurves(curves) ;
  double y[10] = {0,0,0,0,0,0,0,0,0,0} ;
  char   line[CurvesFile_MaxLengthOfTextLine] ;
  char *c = expr ;
  
  /* Compute the curve values */
  for(c = expr ; *c ; c++) {
    if(*c == '$') {
      int i = *(c + 1) - '0' ;
      
      if(i > NbOfCurves || i > 10) {
        Message_RuntimeError("Evaluate: not valid curve index") ;
      }
      
      {
        Curve_t* curve = Curves_GetCurve(curves) + i - 1 ;
      
        y[i - 1] = Curve_ComputeValue(curve,x) ;
      }
    }
  }
  
  /* We copy expr in line by replacing 
   * xlabel with "($0)" and "$i" with "($i)" */
  {
    char *p = line ;
    
    c = expr ;
    
    while(*c) {
      
      if(*c == '$') {
        *p++ = '(' ;
        *p++ = '$' ;
        *p++ = *(c + 1) ;
        *p++ = ')' ;
        c += 2 ;
      } else if(!strncmp(c,xlabel,strlen(xlabel))) {
        *p++ = '(' ;
        *p++ = '$' ;
        *p++ = '0' ;
        *p++ = ')' ;
        c += strlen(xlabel) ;
      } else {
        *p++ = *c++ ;
      }
      
    }
    
    *p++ = '\0' ;
  }
  
  /* Substitute $i with the computed values */
  while((c = strchr(line,'$'))) {
    int i = *(c + 1) - '0' ;
      
    if(Curves_CannotAppendCurves(curves,i) || i > 10) {
      Message_RuntimeError("Evaluate") ;
    }
      
    *(c + 0) = '%' ;
    *(c + 1) = 'f' ;
      
    {
      char line1[CurvesFile_MaxLengthOfTextLine] ;
      
      if(i > 0) {
        sprintf(line1,line,y[i - 1]) ;
      } else if(i == 0) {
        sprintf(line1,line,x) ;
      }
      strcpy(line,line1) ;
    }
  }
  
  {
    double val = Math_EvaluateExpression(line) ;
    
    return(val) ;
  }
}


double (Expressions)(double x,va_list args)
{
  Curves_t* curves = va_arg(args,Curves_t*) ;
  char*  expr = va_arg(args,char*) ;
  char*  ylabel = va_arg(args,char*) ;
  char*  xlabel = va_arg(args,char*) ;
  int NbOfCurves = Curves_GetNbOfCurves(curves) ;
  double y[10] = {0,0,0,0,0,0,0,0,0,0} ;
  char   line[CurvesFile_MaxLengthOfTextLine] ;
  char *c = expr ;
  
  /* Compute the curve values */
  for(c = expr ; *c ; c++) {
    if(*c == '$') {
      int i = *(c + 1) - '0' ;
      
      if(i > NbOfCurves || i > 10) {
        Message_RuntimeError("Evaluate: not valid curve index") ;
      }
      
      {
        Curve_t* curve = Curves_GetCurve(curves) + i - 1 ;
      
        y[i - 1] = Curve_ComputeValue(curve,x) ;
      }
    }
  }
  
  /* We copy expr in line by 
   * 1. adding the expression "xlabel = x ;" 
   * 2. replacing "$i" with "($i)" 
   */
  {
    char *p = line ;
    
    c = expr ;
    
    /* 1. Add the expression "xlabel = x ;" */
    p += sprintf(p,"%s = %e ; ",xlabel,x) ;
    
    /* 2. Replace "$i" with "($i)" */
    while(*c) {
      
      if(*c == '$') {
        *p++ = '(' ;
        *p++ = '$' ;
        *p++ = *(c + 1) ;
        *p++ = ')' ;
        c += 2 ;
      } else {
        *p++ = *c++ ;
      }
      
    }
    
    *p++ = '\0' ;
  }
  
  /* Substitute $i with the computed values */
  while((c = strchr(line,'$'))) {
    int i = *(c + 1) - '0' ;
      
    if(Curves_CannotAppendCurves(curves,i) || i > 10) {
      Message_RuntimeError("Evaluate") ;
    }
      
    *(c + 0) = '%' ;
    *(c + 1) = 'e' ;
      
    {
      char line1[CurvesFile_MaxLengthOfTextLine] ;
      
      if(i > 0) {
        sprintf(line1,line,y[i - 1]) ;
      } else if(i == 0) {
        sprintf(line1,line,x) ;
      }
      strcpy(line,line1) ;
    }
  }
  
  {
    double val = Math_EvaluateExpressions(ylabel,line) ;
    
    return(val) ;
  }
}


double (Affine)(double x,va_list args)
/** y = a*x + b */
{
  double a = va_arg(args,double) ;
  double b = va_arg(args,double) ;
  double y  = a*x + b ;
  
  return(y) ; 
}


double (Langmuir)(double x,va_list args)
/** y = y0 * x / (x + x0) */
{
  double y0 = va_arg(args,double) ;
  double x0 = va_arg(args,double) ;
  double y  = (x > 0.) ? y0*x/(x0 + x) : 0. ;
  
  return(y) ; 
}


double (LangmuirN)(double x,va_list args)
/** y = y0 * (x/x0)**n / (1 + (x/x0)**n) */
{
  double y0  = va_arg(args,double) ;
  double x0  = va_arg(args,double) ;
  double n   = va_arg(args,double) ;
  double y   = langmuir(x,y0,x0,n) ;
  
  return(y) ; 
}


double (Freundlich)(double x,va_list args)
/** Freundlich: y = a*x**b */
{
  double a = va_arg(args,double) ;
  double b = va_arg(args,double) ;
  double y = (x > 0.) ? a*pow(x,b) : 0. ;
  
  return(y) ;
}


double (MualemVanGenuchten)(double s,va_list args)
/** Mualem-Van Genuchten: y = sqrt(x)*(1 - (1 - x**(1/m))**m)**2 */
{
  double m = va_arg(args,double) ;
  
  if(s <= 0) return(0) ;
  else if(s >= 1) return(1) ;
  else return(sqrt(s)*pow(1 - pow(1 - pow(s,1./m),m),2.)) ;
}


double (VanGenuchten)(double p,va_list args)
/** Van Genuchten: y = (1 + x**(1/(1-m)))**(-m) */
{
  double a = va_arg(args,double) ;
  double m = va_arg(args,double) ;
  
  if(p > 0) {
    double y = pow(1 + pow(p/a,1./(1-m)),-m) ;
    
    return(y) ;
  } else return(1) ;
}


double (NavGenuchten)(double s,va_list args)
/** Inverse of Van Genuchten: y = (x**(-1/m) - 1)**(1-m) */
{
  double a = va_arg(args,double) ;
  double m = va_arg(args,double) ;
  
  if(s >= 1) return(0.) ;
  else if(s <= 0) return(HUGE_VAL) ;
  else return(a*pow(pow(s,-1./m) - 1,1-m)) ;
}


double (Mualem_dry)(double p,va_list args)
{
  double a_w = va_arg(args,double) ;
  double m_w = va_arg(args,double) ;
  double a_d = va_arg(args,double) ;
  double m_d = va_arg(args,double) ;
  double kh_max = va_arg(args,double) ;
  FILE*  ftmp1 = va_arg(args,FILE*) ;
  double k_rd ;
  
  char   ltmp1[CurvesFile_MaxLengthOfTextLine] ;
  double kh ;
  double p1 ;
	  
      
  do {
    fgets(ltmp1,sizeof(ltmp1),ftmp1) ;
  } while(ltmp1[0] == '#') ;
  
  sscanf(ltmp1,"%lf %lf",&p1,&kh) ;
  
  if(p != p1) arret("Mualem_dry") ;
  
  kh = 1. - kh/kh_max ;
	  
  {
    double s_w  = vangenuchten(p/a_w,m_w) ;
    double kl   = 1 - pow(1 - pow(s_w,1./m_w),m_w) ;
	    
    double s_d  = vangenuchten(p/a_d,m_d) ;
	    
    k_rd = sqrt(s_d)*(kl + (1. - kl)*kh) ;
  }
          
  return(k_rd) ;
}


double (Mualem_gas)(double s,va_list args)
/** Mualem: y = sqrt(1-x)*(1 - x**(1/m))**(2*m) */
{
  double m = va_arg(args,double) ;
  
  if(s <= 0) return(1) ;
  else if(s >= 1) return(0) ;
  else return(sqrt(1 - s)*pow(1 - pow(s,1./m),2*m)) ;
}


double (VanGenuchten_gas)(double s,va_list args)
/** Mualem: y = ((1-x)**p)*(1 - x**(1/m))**(q*m) */
{
  double m = va_arg(args,double) ;
  double p = va_arg(args,double) ;
  double q = va_arg(args,double) ;
  
  if(s <= 0) return(1) ;
  else if(s >= 1) return(0) ;
  else return(pow(1 - s,p)*pow(1 - pow(s,1./m),q*m)) ;
}


double (Millington)(double s,va_list args)
/** Millington: y = (1 - x)**b */
{
  double b = va_arg(args,double) ;
  
  if(s <= 0) return(1) ;
  else if(s >= 1) return(0) ;
  else return(pow(1 - s,b)) ;
}

double (MonlouisBonnaire)(double s,va_list args)
/** Monlouis Bonnaire: y = (1 - x)**5.5*(1 - x**(1/m))**(2*m) */
{
  double m = va_arg(args,double) ;
  
  if(s <= 0) return(0) ;
  else if(s >= 1) return(1) ;
  else return(pow(1 - s,5.5)*pow(1 - pow(s,1./m),2*m)) ;
}


double (RedlichKwongCO2)(double Pa,va_list args)
{
  double T = va_arg(args,double) ;

  return(MolarDensityOfCO2_RedlichKwong(Pa,T)) ;
}


double (MolarDensityOfCO2)(double Pa,va_list args)
{
  double T = va_arg(args,double) ;
  double rho = MolarDensityOfCO2_RedlichKwong(Pa,T) ;

  return(1.e3*rho) ;
}


double (MolarDensityOfPerfectGas)(double Pa,va_list args)
/** Perfect gas law: y = Pa/RT */
{
  double T = va_arg(args,double) ;
  double R  = 8.3143 ;
  double RT = R*T ;
  
  return(Pa/RT) ;
}


double (FenghourCO2)(double Pa,va_list args)  
{
  double T = va_arg(args,double) ;
  
  return(ViscosityOfCO2_Fenghour(Pa,T)) ;
}

      
double (CSH3EndMembers)(double s_CH,va_list args)
{
  double k_CH  = 6.456e-6 ;
  /* Amorpheous Silica */
  double k_SH = 1.93642e-3 ;
  double x_SH = 0 ;
  double y_SH = 1 ;
  double z_SH = 2 ;
  double v_SH = (29. + 7.*z_SH)*cm3 ;
  /* Jennite (y_Jen = 1)         (y_Jen = 0.9)      */
  double k_Jen0 = 4.39e-18  ; /*  k_Jen = 2.39e-16  */
  double x_Jen0 = 1.66666   ; /*  x_Jen = 1.5       */
  double z_Jen0 = 2.66666   ; /*  z_Jen = 2.4       */
  double v_Jen0 = 81.8*cm3  ; /*  v_Jen = 73.6  (cm3/mol) */
  /* Tobermorite (y_Tob = 1)     (y_Tob = 1.8)      */
  double k_Tob0 = 1.684e-12 ; /*  k_Tob = 6.42e-22  */
  double x_Tob0 = 0.83333   ; /*  x_Tob = 1.5       */
  double z_Tob0 = 1.83333   ; /*  z_Tob = 3.3       */
  double v_Tob0 = 54.8*cm3  ; /*  v_Tob = 98.6  (cm3/mol) */
  
  double s_SH = 1. ;
  
  double y_Tob = va_arg(args,double) ;
  double y_Jen = va_arg(args,double) ;
  char *outputtype = va_arg(args,char*) ;
  
  /* Jennite */
  double x_Jen   = x_Jen0*y_Jen ;
  double z_Jen   = z_Jen0*y_Jen ;
  double k_Jen   = pow(k_Jen0,y_Jen) ;
  double v_Jen   = v_Jen0*y_Jen ;
  /* Tobermorite */
  double x_Tob   = x_Tob0*y_Tob ;
  double z_Tob   = z_Tob0*y_Tob ;
  double k_Tob   = pow(k_Tob0,y_Tob) ;
  double v_Tob   = v_Tob0*y_Tob ;
  
  double output ;

  {
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
          
          
    for(j = 0 ; j < 3 ; j++) {
      a_CSH[j] = pow(k_CH*s_CH,x_CSH[j])*pow(k_SH,y_CSH[j])/k_CSH[j] ;
    }
          
    s_SH = fraction_Silica(a_CSH,y_CSH,3,s_SH) ;
          
    for(j = 0 ; j < 3 ; j++) {
      q_CSH[j] = a_CSH[j]*pow(s_SH,y_CSH[j]) ;
    }
          
    {
      double x_m = 0,y_m = 0,z_m = 0 ;
      double v_m = 0 ;
      
      for(j = 0 ; j < 3 ; j++) {
        x_m += q_CSH[j]*x_CSH[j] ;
        y_m += q_CSH[j]*y_CSH[j] ;
        z_m += q_CSH[j]*z_CSH[j] ;
        v_m += q_CSH[j]*v_CSH[j] ;
      }
            
      x_m /= y_m ; z_m /= y_m ; v_m /= y_m ;
      
      if(!strcmp(outputtype,"S_SH")) {
        output = s_SH ;
      } else if(!strcmp(outputtype,"V_CSH")) {
        output = v_m ;
      } else if(!strcmp(outputtype,"Z_CSH")) {
        output = z_m ;
      } else if(!strcmp(outputtype,"X_CSH")) {
        output = x_m ;
      } else {
        arret("CSH3EndMembers") ;
      }
    }
  }
      
  return(output) ;
}

      
double (CSHLangmuirN)(double s,va_list args)
{
  double x1 = va_arg(args,double) ;
  double s1 = va_arg(args,double) ;
  double n1 = va_arg(args,double) ;
  double x2 = va_arg(args,double) ;
  double s2 = va_arg(args,double) ;
  double n2 = va_arg(args,double) ;
  char *outputtype = va_arg(args,char*) ;
  double output ;
  
  if(!strcmp(outputtype,"S_SH")) {
    double s_SH = pow(1 + pow(s/s1,n1),-x1/n1)*pow(1 + pow(s/s2,n2),-x2/n2) ;
    
    output = s_SH ;
  } else if(!strcmp(outputtype,"X_CSH")) {
    double x = langmuir(s,x1,s1,n1) + langmuir(s,x2,s2,n2) ;
    
    output = x ;
  } else {
    arret("CSHLangmuirN") ;
  }
      
  return(output) ;
}


/* Other intern functions */
double (vangenuchten)(double p,double m)
/** Van Genuchten: y = (1 + x**(1/(1-m)))**(-m) */
{
  if(p > 0) return(pow(1 + pow(p,1./(1-m)),-m)) ;
  else return(1) ;
}

double (fraction_Silica)(double *a,double *y,int n,double q0)
/**  */
{
  double err,tol = 1e-8 ;
  double q = q0 ;
  int    i = 0 ;
  
  do {
    double f  = - 1 ;
    double df = 0 ;
    double dq ;
    int    j ;
    
    for(j = 0 ; j < n ; j++) {
      double fj = a[j]*pow(q,y[j]) ;
      
      f  += fj ;
      df += y[j]*fj/q ;
    }
    
    dq = -f/df ;
    err = fabs(dq/q) ;
    q += dq ;
    
    if(i++ > 20) {
      printf("q0 = %e\n",q0) ;
      printf("q  = %e\n",q) ;
      arret("fraction_Silica : non convergence") ;
    }
  } while(err > tol) ;
  
  return(q) ;
}

double (MolarDensityOfCO2_RedlichKwong)(double Pa,double T)
/** Redlich Kwong EOS for scCO2 */
/* Implemented by J. Shen 
 * From Redlich-Kwong model of EOS for scCO2 (Spycher2003)
 * P = RT/(V - b) + a/(sqrt(T)*V*(V + b))
 * Input units: Pressure in Pascal, T in Kelvin 
 * Output units: mol/L */
{
  double R         = 83.14472 ;   /* cm3*bar/(K*mol) */
  double RT        = R*T ;
  double T05       = sqrt(T) ;
  double A_CO2     = 7.54e7 - 4.13e4*T ; /* bar*cm6*K0.5/mol2 */
  double B_CO2     = 27.80 ; /* cm3/mol */
  double AB_CO2    = A_CO2*B_CO2 ;
  double B2_CO2    = B_CO2*B_CO2 ;
  double P         = Pa/1.e5 ; /* unit of bar */
  double P0        = 1. ;
  
  double gc_co2 ;
  
  if(P > P0) {
    double RToP      = RT/P ;
    double PT05      = P*T05 ;
    double a         = 1. ;
    double b         = - RToP ;
    double c         = - (B_CO2*RToP - A_CO2/(PT05) + B2_CO2) ;
    double d         = - (AB_CO2/(PT05));
    double x[4]      ;
    
    x[0]  = a  ; x[1]  = b  ; x[2]  = c  ; x[3]  = d  ;
    Math_ComputePolynomialEquationRoots(x,3) ;
    
    {
      double V = x[0] ;    /* cm3/mol */
      gc_co2 = 1.e3/V ;    /* mol/dm3 */
    }
    
  } else {
    double RToP0     = RT/P0 ;
    double P0T05     = P0*T05 ;
    double a0        = 1. ;
    double b0        = - (RToP0) ;
    double c0        = - (B_CO2*RToP0 - A_CO2/(P0T05) + B2_CO2) ;
    double d0        = - (AB_CO2/(P0T05)) ;
    double x0[4]     ;
    
    x0[0] = a0 ; x0[1] = b0 ; x0[2] = c0 ; x0[3] = d0 ;
    Math_ComputePolynomialEquationRoots(x0,3) ;
    
    {
      double V0 = x0[0] ;       /* cm3/mol */
      gc_co2 = P/P0*1.e3/V0  ;  /* mol/dm3 */
    }
  }
  
  return(gc_co2) ;
}


double (ViscosityOfCO2_Fenghour)(double Pa,double T)  
/** Fenghour viscosity of scCO2: */
/* Implemented by J. Shen 
 * From	A. Fenghour (1997)
 * Input units: Pressure in Pa, Temperature in Kelvin
 * Output units: Pa.s 
 * Range of validity:
 *      P [0:300] MPa and T [200-1000] K
 *      P [0-30]  MPa and T [1000-1500] K 
*/
{
  double P		    = MAX(Pa,1.e5) ; /* mu_co2 cst for P < 1e5 Pa */
  double gc_co2 	= MolarDensityOfCO2_RedlichKwong(P,T) ;    /*mol/L */
  double rho    	= 44.* gc_co2 ;				/*kg/m3*/
  double d11  		= 0.4071119e-2 ;
  double d21		  = 0.7198037e-4 ;
  double d64 		  = 0.2411697e-16 ;
  double d81		  = 0.2971072e-22 ;
  double d82		  = -0.1627888e-22 ;
  double T_red		= T/251.196 ;
  double logT_red	= log(T_red) ;
  double a0			  = 0.235156 ;
  double a1			  = -0.491266 ;
  double a2			  = 5.211155e-2 ;
  double a3			  = 5.347906e-2 ;
  double a4			  = -1.537102e-2 ;
  double A_mu		  = a0 + a1*logT_red + a2*pow(logT_red,2.) + a3*pow(logT_red,3.)  + a4*pow(logT_red,4.) ;
  double g_mu   	= exp(A_mu) ;
  double mu_zero	= 1.00697*pow(T,0.5)/g_mu ; /*zero-density viscosity, uPa.s*/
  double rho2     = rho*rho ;
  double rho6     = rho2*rho2*rho2 ;
  double rho8     = rho6*rho2 ;
  double mu_excess 	= d11*rho+ d21*rho2 + d64*rho6/pow(T_red,3.) + d81*rho8 + d82*rho8/T_red;/* in uPa.s*/
  double mu 		  = mu_zero + mu_excess; /* uPa.s*/
  double mu_co2		= mu*1.e-6 ;           /* Pa.s*/
  return(mu_co2) ;
}



double (langmuir)(double x,double y0,double x0,double n)
/** y = y0 * (x/x0)**n / (1 + (x/x0)**n) */
{
  double a   = (x0 > 0.) ? x/x0 : 0. ;
  double p   = (a  > 0.) ? pow(a,n) : 0. ;
  double y   = (x0 > 0.) ? y0*p/(1 + p) : y0 ;
  
  return(y) ; 
}

