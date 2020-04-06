#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "OutputFiles.h"
#include "Message.h"
#include "Mesh.h"
#include "BilVersion.h"
#include "TextFile.h"
#include "Models.h"
#include "String.h"
#include "Mry.h"


//char   OutputFile_TypeOfCurrentFile ;


static void  (OutputFiles_PostProcessForGmshASCIIFileFormatVersion2_2)(OutputFiles_t*,DataSet_t*) ;
static void  (OutputFiles_PostProcessForGmshParsedFileFormatVersion2)(OutputFiles_t*,DataSet_t*) ;
static Views_t* (OutputFiles_CreateGlobalViews)(OutputFiles_t*,Models_t*,TextFile_t*) ;



OutputFiles_t*   (OutputFiles_Create)(char* filename,int n_dates,int n_points)
{
  OutputFiles_t* outputfiles = (OutputFiles_t*) Mry_New(OutputFiles_t) ;
  
  
  /* The file name */
  {
    int LengthOfName = strlen(filename) + 1 ;
    char* name = (char*) Mry_New(char[LengthOfName]) ;
    
    strcpy(name,filename) ;

    OutputFiles_GetDataFileName(outputfiles) = name ;
  }


  /* Date Files */
  OutputFiles_GetNbOfDateFiles(outputfiles) = n_dates ;

  /* Date Output Files */
  {
    int n = ceil(log10((double) n_dates+1)) ;
    int LengthOfName = strlen(filename) + 3 + n ;
    char* name = (char*) Mry_New(char[LengthOfName]) ;
      
    sprintf(name,"%s.t0",filename) ;
    
    {  
      OutputFile_t* outputfile = OutputFile_Create(name,n_dates) ;
    
      OutputFiles_GetDateOutputFile(outputfiles) = outputfile ;
    }
      
    free(name) ;
  }
  
  
  /* Point Files */
  OutputFiles_GetNbOfPointFiles(outputfiles) = n_points ;

  /* Point Output Files */
  {
    int n = ceil(log10((double) n_points+1)) ;
    int LengthOfName = strlen(filename) + 3 + n ;
    char* name = (char*) Mry_New(char[LengthOfName]) ;
      
    sprintf(name,"%s.p1",filename) ;
    
    {  
      OutputFile_t* outputfile = OutputFile_Create(name,n_points) ;
    
      OutputFiles_GetPointOutputFile(outputfiles) = outputfile ;
    }
      
    free(name) ;
  }
  
  
  /* Results */
  {
    Results_t* results = Results_Create(OutputFiles_MaxNbOfViews) ;
    
    OutputFiles_GetResults(outputfiles) = results ;
  }
  
  
  /* Text line */
  {
    char* line = (char*) Mry_New(char[OutputFiles_MaxLengthOfTextLine]) ;
    
    OutputFiles_GetTextLine(outputfiles) = line ;
  }
  
  
  return(outputfiles) ;
}



void   (OutputFiles_Delete)(void* self)
{
  OutputFiles_t** poutputfiles = (OutputFiles_t**) self ;
  OutputFiles_t*   outputfiles = *poutputfiles ;
  int n_dates = OutputFiles_GetNbOfDateFiles(outputfiles) ;
  int n_points = OutputFiles_GetNbOfPointFiles(outputfiles) ;
  
  free(OutputFiles_GetDataFileName(outputfiles)) ;
    
  OutputFile_Delete(&(OutputFiles_GetDateOutputFile(outputfiles)),n_dates) ;
  OutputFile_Delete(&(OutputFiles_GetPointOutputFile(outputfiles)),n_points) ;
  
  Results_Delete(&(OutputFiles_GetResults(outputfiles))) ;
  
  free(OutputFiles_GetTextLine(outputfiles)) ;
  
  free(outputfiles) ;
  *poutputfiles = NULL ;
}



void (OutputFiles_PostProcessForGmshASCIIFileFormat)(OutputFiles_t* outputfiles,DataSet_t* dataset)
{
  double version = 0 ;
  
  
  /* Read the version in the first line */
  {
    /* Assuming that there is at least one date file */
    OutputFile_t* outputfile = OutputFiles_GetDateOutputFile(outputfiles) ;
    TextFile_t* textfile = OutputFile_GetTextFile(outputfile) ;
    char* c ;
    
    TextFile_OpenFile(textfile,"r") ;
    c = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) ;
  
    //if((c = strstr(c,"Version") + strlen("Version"))) {
    if((c = String_FindAndSkipToken(c,"Version"))) {
      sscanf(c,"%lf",&version) ;
    }
  
    TextFile_CloseFile(textfile) ;
  }
  
  
  if(version < 2.2) {
    Materials_t* materials = DataSet_GetMaterials(dataset) ;
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
    
    if(n_usedmodels > 1) {
      arret("OutputFiles_PostProcessForGmshASCIIFileFormat") ;
    }
    
    OutputFiles_PostProcessForGmshASCIIFileFormatVersion2_2(outputfiles,dataset) ;
  } else {
    OutputFiles_PostProcessForGmshASCIIFileFormatVersion2_2(outputfiles,dataset) ;
  }
}



void (OutputFiles_PostProcessForGmshASCIIFileFormatVersion2_2)(OutputFiles_t* outputfiles,DataSet_t* dataset)
/* Create gmsh post-processing files in MSH ASCII file format version 2.2 */
{
  int    n_dates = OutputFiles_GetNbOfDateFiles(outputfiles) ;
  Mesh_t* mesh = DataSet_GetMesh(dataset) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  int    n_el = Mesh_GetNbOfElements(mesh) ;
    
  Materials_t* materials = DataSet_GetMaterials(dataset) ;
  Models_t* usedmodels = Materials_GetUsedModels(materials) ;

  OutputFile_t* outputfile = OutputFiles_GetDateOutputFile(outputfiles) ;
  TextFile_t* datefile = OutputFile_GetTextFile(outputfile) ;
  
  Views_t* globalviews = OutputFiles_CreateGlobalViews(outputfiles,usedmodels,datefile) ;
  
  int nbofglobalviews = Views_GetNbOfViews(globalviews) ;
  
  FILE*  fileofview[OutputFiles_MaxNbOfViews] ;
  
  int    nbofrecords[OutputFiles_MaxNbOfViews] ;
  
  
  if(nbofglobalviews > OutputFiles_MaxNbOfViews) {
    arret("OutputFiles_PostProcessForGmshASCIIFileFormatVersion2_2: too many views") ;
  }
   


  /* Nb of records */
  {
    int ie ;
    
    {
      int i ;
      
      for(i = 0 ; i < nbofglobalviews ; i++) {
        nbofrecords[i] = 0 ;
      }
    }
    
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t*  elt = el + ie ;
      Material_t* mat = Element_GetMaterial(elt) ;
      char* codename = Material_GetCodeNameOfModel(mat) ;
      int usedmodelindex = Models_FindModelIndex(usedmodels,codename) ;
      Model_t*  usedmodel = Models_GetModel(usedmodels) + usedmodelindex ;
      Views_t*  views = Model_GetViews(usedmodel) ;
      int nviews = Views_GetNbOfViews(views) ;
      View_t*   view  = Views_GetView(views) ;
      
      
      if(Element_IsSubmanifold(elt)) continue ;
      
      if(!mat) continue ;
      
      {
        int i ;
      
        for(i = 0 ; i < nviews ; i++) {
          int    jview = View_GetGlobalIndex(view + i) ;
      
          nbofrecords[jview] += 1 ;
        }
      }
    }
  }
  
  
  /* Open the date files for reading */
  {
    int i ;
    
    for(i = 0 ; i < n_dates ; i++) {
      TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i) ;
    
      TextFile_OpenFile(textfile,"r") ;
    }
  }


  /* Create and open the view files for writing */
  {
    char   *filename = OutputFiles_GetDataFileName(outputfiles) ;
    int n = ceil(log10((double) nbofglobalviews)) ;
    int i = strlen(filename) + 4 + n ;
    
    if(i > OutputFiles_MaxLengthOfFileName) {
      arret("OutputFiles_PostProcessForGmshASCIIFileFormat (1) : nom trop long") ;
    }
    
    for(i = 0 ; i < nbofglobalviews ; i++) {
      char   nom_pos[OutputFiles_MaxLengthOfFileName] ;
      FILE   *ficp ;

      sprintf(nom_pos,"%s.pos%d",filename,i + 1) ;
    
      ficp = fopen(nom_pos,"w") ;
    
      if(!ficp) {
        arret("erreur a l'ouverture du fichier %s",nom_pos) ;
      }
    
      fileofview[i] = ficp ;
    }
  }


  /* Record the results in the view files */
  {
    int i_temps ;
    
    /* First records */
    {
      int i ;
      
      for(i = 0 ; i < nbofglobalviews ; i++) {
        FILE*   ficp = fileofview[i] ;
        float   version_number = 2.2 ;
        int     file_type = 0 ;
        int     data_type = sizeof(double) ;
    
        fprintf(ficp,"$MeshFormat\n") ;
        fprintf(ficp,"%3.1f",version_number) ;
        fprintf(ficp," %d",file_type) ;
        fprintf(ficp," %d",data_type) ;
        fprintf(ficp,"\n") ;
        fprintf(ficp,"$EndMeshFormat\n") ;
      }
    }

    
    for(i_temps = 0 ; i_temps < n_dates ; i_temps++) {
      TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i_temps) ;
      double temps ;

      /* Read the time */
      {
        char*  pline ;
        
        while((pline = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile))) {
      
          if(pline[0] == '#') {
        
            /* Time */
            if(strstr(pline,"temps") || strstr(pline,"Time")) {
          
              pline = strchr(pline,'=') + 1 ;
          
              sscanf(pline,"%le",&temps) ;
            
              break ;
            }
        
          }
        }
      }
  
      /* Build a section $ElementNodeData/$EndElementNodeData for each time step */
      {
        int    ie ;
        

        /* First records of the section */
        {
          int    i ;
          
          for(i = 0 ; i < nbofglobalviews ; i++) {
            View_t* view = Views_GetView(globalviews) + i ;
            char*  nameofview = View_GetNameOfView(view) ;
            int    nbcompofview = View_GetNbOfComponents(view) ;
            FILE*  ficp = fileofview[i] ;
            int    nb_string_tags  = 1 ;
            int    nb_real_tags    = 1 ;
            int    nb_integer_tags = 4 ;
            int    partition = 0 ;

      
            fprintf(ficp,"$ElementNodeData\n") ;
            fprintf(ficp,"%d\n",nb_string_tags) ;
            fprintf(ficp,"\"%s\"\n",nameofview) ;
            fprintf(ficp,"%d\n",nb_real_tags) ;
            fprintf(ficp,"%e\n",temps) ;
            fprintf(ficp,"%d\n",nb_integer_tags) ;
            fprintf(ficp,"%d\n",i_temps) ;
            fprintf(ficp,"%d\n",nbcompofview) ;
            fprintf(ficp,"%d\n",nbofrecords[i]) ;
            fprintf(ficp,"%d\n",partition) ;
          }
        }

        
        /* Record the results per element */
        for(ie = 0 ; ie < n_el ; ie++) {
          Element_t*  elt = el + ie ;
          Material_t* mat = Element_GetMaterial(elt) ;
          char* codename = Material_GetCodeNameOfModel(mat) ;
          int usedmodelindex = Models_FindModelIndex(usedmodels,codename) ;
          Model_t*  usedmodel = Models_GetModel(usedmodels) + usedmodelindex ;
          Views_t*  views = Model_GetViews(usedmodel) ;
          int nviews = Views_GetNbOfViews(views) ;
          View_t*   view  = Views_GetView(views) ;
          int    nn = Element_GetNbOfNodes(elt) ;

          double val[OutputFiles_MaxNbOfViews][9*Element_MaxNbOfNodes] ;
          
          if(usedmodelindex < 0) {
            arret("OutputFiles_PostProcessForGmshASCIIFileFormat") ;
          }
          
          /* We skip elements of the boundaries 
           * (see OutputFiles_BackupSolutionAtTime) */
          if(Element_IsSubmanifold(elt)) continue ;

          if(!mat) continue ;

          
          /* Read the values */
          {
            int    in ;
            
            for(in = 0 ; in < nn ; in++) {
              char* pline = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) ;
              double x[3] ;
              int i ;
                  
              /* We skip the commented lines */
              while(pline[0] == '#') {
                pline = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) ;
              }

              /* Coordinates of nodes */
              for(i = 0 ; i < 3 ; i++) {
                sscanf(pline,"%le",x + i) ;
                pline  = strchr(pline,' ') ;
                pline += strspn(pline," ") ;
              }
            
              /* Values */
              for(i = 0 ; i < nviews ; i++) {
                int    nc = View_GetNbOfComponents(view + i) ;
                int    j ;
              
                for(j = 0 ; j < nc ; j++) {
                  sscanf(pline,"%le",val[i] + nc*in + j) ;
                  pline  = strchr(pline,' ') ;
                  pline += strspn(pline," ") ;
                }
              }
            }
          }

          /* Print these values in view files */
          {
            int i ;
            
            for(i = 0 ; i < nviews ; i++) {
              int    nc = View_GetNbOfComponents(view + i) ;
              int    jview = View_GetGlobalIndex(view + i) ;
              FILE*  ficp = fileofview[jview] ;
              int    j ;
            
              fprintf(ficp,"%d %d",ie + 1,nn) ;
            
              for(j = 0 ; j < nc*nn ; j++) {
                fprintf(ficp," %e",val[i][j]) ;
              }
            
              fprintf(ficp,"\n") ;
            }
          }
        }
      
      
        /* Last record of the views */
        {
          int i ;
        
          for(i = 0 ; i < nbofglobalviews ; i++) {
            FILE   *ficp = fileofview[i] ;
        
            fprintf(ficp,"$EndElementNodeData\n") ;
          }
        }
      }
    }
  }


  /* Close the view files */
  {
    int i ;
    
    for(i = 0 ; i < nbofglobalviews ; i++) {
      FILE   *ficp = fileofview[i] ;
    
      fclose(ficp) ;
    }
  }
  
  
  /* Close the date files */
  {
    int i ;
    
    for(i = 0 ; i < n_dates ; i++) {
      TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i) ;
    
      TextFile_CloseFile(textfile) ;
    }
  }

}



void (OutputFiles_PostProcessForGmshParsedFileFormat)(OutputFiles_t* outputfiles,DataSet_t* dataset)
{
  double version = 0 ;
  
  
  /* Read the version in the first line */
  {
    /* Assuming that there is at least one date file */
    OutputFile_t* outputfile = OutputFiles_GetDateOutputFile(outputfiles) ;
    TextFile_t* textfile = OutputFile_GetTextFile(outputfile) ;
    char* c ;
    
    TextFile_OpenFile(textfile,"r") ;
    c = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) ;
  
    //if((c = strstr(c,"Version") + strlen("Version"))) {
    if((c = String_FindAndSkipToken(c,"Version"))) {
      sscanf(c,"%lf",&version) ;
    }
  
    TextFile_CloseFile(textfile) ;
  }
  
  
  if(version < 2.2) {
    Materials_t* materials = DataSet_GetMaterials(dataset) ;
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
    
    if(n_usedmodels > 1) {
      arret("OutputFiles_PostProcessForGmshASCIIFileFormat") ;
    }
    
    OutputFiles_PostProcessForGmshParsedFileFormatVersion2(outputfiles,dataset) ;
  } else {
    OutputFiles_PostProcessForGmshParsedFileFormatVersion2(outputfiles,dataset) ;
  }
}



void (OutputFiles_PostProcessForGmshParsedFileFormatVersion2)(OutputFiles_t* outputfiles,DataSet_t* dataset)
/* Create gmsh post-processing files in parsed file format */
{
  int    n_dates = OutputFiles_GetNbOfDateFiles(outputfiles) ;
  
  Mesh_t* mesh = DataSet_GetMesh(dataset) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  int    n_el = Mesh_GetNbOfElements(mesh) ;
  int    dim = Mesh_GetDimension(mesh) ;
  
  Materials_t* materials = DataSet_GetMaterials(dataset) ;
  Models_t* usedmodels = Materials_GetUsedModels(materials) ;
  
  OutputFile_t* outputfile = OutputFiles_GetDateOutputFile(outputfiles) ;
  TextFile_t* datefile = OutputFile_GetTextFile(outputfile) ;
  
  Views_t* globalviews = OutputFiles_CreateGlobalViews(outputfiles,usedmodels,datefile) ;
  
  int nbofglobalviews = Views_GetNbOfViews(globalviews) ;

  FILE*  fileofview[OutputFiles_MaxNbOfViews] ;
  
  int    n_vertices[3][Element_MaxNbOfNodes] = {{1,2,2},{1,2,3,4,4,3,0,4},{1,2,3,4,4,0,0,8}} ;
  const char*  type[] = {"PL","PLTQ","PLTSYI7H"} ;
  char   svt[] = "S2V45678T" ;
  
  
  if(nbofglobalviews > OutputFiles_MaxNbOfViews) {
    arret("OutputFiles_PostProcessForGmshParsedFileFormatVersion2: too many views") ;
  }


  /* Open the date files for reading */
  {
    int i ;
    
    for(i = 0 ; i < n_dates ; i++) {
      TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i) ;
    
      TextFile_OpenFile(textfile,"r") ;
    }
  }


  /* Create and open the view files for writing */
  {
    char   *filename = OutputFiles_GetDataFileName(outputfiles) ;
    int n = ceil(log10((double) nbofglobalviews)) ;
    int i = strlen(filename) + 4 + n ;
    
    if(i > OutputFiles_MaxLengthOfFileName) {
      arret("OutputFiles_PostProcessForGmshParsedFileFormat (1) : nom trop long") ;
    }
    
    for(i = 0 ; i < nbofglobalviews ; i++) {
      char   nom_pos[OutputFiles_MaxLengthOfFileName] ;
      FILE   *ficp ;

      sprintf(nom_pos,"%s.pos%d",filename,i + 1) ;
    
      ficp = fopen(nom_pos,"w") ;
    
      if(!ficp) {
        arret("erreur a l'ouverture du fichier %s",nom_pos) ;
      }
    
      fileofview[i] = ficp ;
    }
  }
  
  
  /* Record the results in the view files */
  {
    int    ie ;
    
    
    /* The first record is the name of the view */
    {
      int i ;
      
      for(i = 0 ; i < nbofglobalviews ; i++) {
      View_t* view = Views_GetView(globalviews) + i ;
      char*   nameofview = View_GetNameOfView(view) ;
        FILE   *ficp = fileofview[i] ;
    
        fprintf(ficp,"View \"%s\" {\n",nameofview) ;
      }
    }


    /* Record the results per element */
    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t* elt = el + ie ;
      Material_t* mat = Element_GetMaterial(elt) ;
      char* codename = Material_GetCodeNameOfModel(mat) ;
      int usedmodelindex = Models_FindModelIndex(usedmodels,codename) ;
      Model_t*  usedmodel = Models_GetModel(usedmodels) + usedmodelindex ;
      Views_t*  views = Model_GetViews(usedmodel) ;
      int nviews = Views_GetNbOfViews(views) ;
      View_t*   view  = Views_GetView(views) ;
      int    nn = Element_GetNbOfNodes(elt) ;
      int    i_temps ;
      
      if(usedmodelindex < 0) {
        arret("OutputFiles_PostProcessForGmshParsedFileFormat") ;
      }
      
      /* We skip elements of the boundaries 
       * (see OutputFiles_BackupSolutionAtTime) */
      if(Element_IsSubmanifold(elt)) continue ;

      if(!mat) continue ;
    
    
      for(i_temps = 0 ; i_temps < n_dates ; i_temps++) {
        TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i_temps) ;
        double val[OutputFiles_MaxNbOfViews][9*Element_MaxNbOfNodes] ;
        double x_e[3*Element_MaxNbOfNodes] ;
        int    in ;
      
        for(in = 0 ; in < nn ; in++) {
          char* pline = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) ;
          int i ;
          
          /* We skip the commented lines */
          while(pline[0] == '#') {
            pline = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile) ;
          }
          
          for(i = 0 ; i < 3 ; i++) {
            sscanf(pline,"%le",x_e + 3*in + i) ;
            pline  = strchr(pline,' ') ;
            pline += strspn(pline," ") ;
          }
        
          for(i = 0 ; i < nviews ; i++) {
            int nc = View_GetNbOfComponents(view + i) ;
            int j ;
          
            for(j = 0 ; j < nc ; j++) {
              sscanf(pline,"%le",val[i] + nc*in + j) ;
              pline  = strchr(pline,' ') ;
              pline += strspn(pline," ") ;
            }
          }
        }
      
        if(i_temps == 0) {
          int i ;
          
          for(i = 0 ; i < nviews ; i++) {
            int    nvert = n_vertices[dim-1][nn-1] ;
            int    nc = View_GetNbOfComponents(view + i) ;
            int    jview = View_GetGlobalIndex(view + i) ;
            FILE*  ficp = fileofview[jview] ;
            int j ;

            if(nvert == 0) {
              arret("OutputFiles_PostProcessForGmshParsedFileFormat (2): type non traite") ;
            }

            /* Type d'enregistrement */
            fprintf(ficp,"%c",svt[nc - 1]) ;
            fprintf(ficp,"%c",type[dim-1][nvert-1]) ;
            if(nvert == nn/2) fprintf(ficp,"2") ; /* element P2 */

            /* Coordonnees */
            fprintf(ficp,"(%e",x_e[0]) ;
            for(j = 1 ; j < 3*nn ; j++) {
              fprintf(ficp,",%e",x_e[j]) ;
            }
            fprintf(ficp,")") ;

            /* Valeurs */
            fprintf(ficp,"{%e",val[i][0]) ;
            for(j = 1 ; j < nc*nn ; j++) {
              fprintf(ficp,",%e",val[i][j]) ;
            }
          }
        
        } else {
          int i ;
          
          for(i = 0 ; i < nviews ; i++) {
            int    nc = View_GetNbOfComponents(view + i) ;
            int    jview = View_GetGlobalIndex(view + i) ;
            FILE   *ficp = fileofview[jview] ;
            int j ;
          
            for(j = 0 ; j < nc*nn ; j++) {
              fprintf(ficp,",%e",val[i][j]) ;
            }
          }
        }
      }
      
      {
        int i ;
        
        for(i = 0 ; i < nviews ; i++) {
          int    jview = View_GetGlobalIndex(view + i) ;
          FILE   *ficp = fileofview[jview] ;
      
          fprintf(ficp,"};\n") ;
        }
      }
    }
    
    
    /* Last record of the views */
    {
      int i ;
      
      for(i = 0 ; i < nbofglobalviews ; i++) {
        FILE   *ficp = fileofview[i] ;
    
        fprintf(ficp,"};\n") ;
      }
    }
  
  }


  /* Close the view files */
  {
    int i ;
    
    for(i = 0 ; i < nbofglobalviews ; i++) {
      FILE   *ficp = fileofview[i] ;
    
      fclose(ficp) ;
    }
  }
  
  
  /* Close the date files */
  {
    int i ;
    
    for(i = 0 ; i < n_dates ; i++) {
      TextFile_t* textfile = OutputFile_GetTextFile(outputfile + i) ;
    
      TextFile_CloseFile(textfile) ;
    }
  }
}




void (OutputFiles_BackupSolutionAtTime_)(OutputFiles_t* outputfiles,DataSet_t* dataset,double t,int idate1)
/* Backup solutions at a given time in the appropriate output file */
{
  Mesh_t* mesh = DataSet_GetMesh(dataset) ;
  unsigned short int dim = Mesh_GetDimension(mesh) ;
  int n_el = Mesh_GetNbOfElements(mesh) ;
  Element_t* el = Mesh_GetElement(mesh) ;
  Materials_t* materials = DataSet_GetMaterials(dataset) ;
  Models_t* usedmodels = Materials_GetUsedModels(materials) ;
  int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
  
  OutputFile_t* outputfile = OutputFiles_GetDateOutputFile(outputfiles) ;
  
  Result_t* r_s = Results_GetResult(OutputFiles_GetResults(outputfiles)) ;
  
  /* Open the date file for writing */
  TextFile_t* textfile = OutputFile_GetTextFile(outputfile + idate1) ; 
  FILE *fict = TextFile_OpenFile(textfile,"w") ;
  
  OutputFile_TypeOfCurrentFile = 't' ;
  
  
  /* First lines: version and date */
  {
    time_t date ;
    time(&date) ;
    fprintf(fict,"# Version " BIL_VERSION ", %s",ctime(&date)) ;
    fprintf(fict,"# Time = %e\n",t) ;
  }
  
  
  {
    unsigned short int headings[10] = {0,0,0,0,0,0,0,0,0,0} ;
    int    ie ;
    
    if(n_usedmodels > 10) {
      arret("BackupSolutionAtTime (1): too many used models") ;
    }
  

    for(ie = 0 ; ie < n_el ; ie++) {
      Element_t* elt = el + ie ;
      Material_t* mat = Element_GetMaterial(elt) ;
      char*  codename = (mat) ? Material_GetCodeNameOfModel(mat) : NULL ;
      int usedmodelindex = Models_FindModelIndex(usedmodels,codename) ;
      unsigned short int entete = headings[usedmodelindex] ;
      
      /* if(!strcmp(codename,ucodename)) { */
      if(mat) {
        int  nn = Element_GetNbOfNodes(elt) ;
        int  in ;
        
        /* We skip elements of the boundaries */
        if(Element_IsSubmanifold(elt)) continue ;
        
        for(in = 0 ; in < nn ; in++) {
          double *x_s = Element_GetNodeCoordinate(elt,in) ;
          int    nso ;
          int    j,k ;
        
          Element_FreeBuffer(elt) ;
          nso = Element_ComputeOutputs(elt,t,x_s,r_s) ;
          
          if(nso > OutputFiles_MaxNbOfViews) {
            arret("BackupSolutionAtTime (1) : trop de valeurs") ;
          }

          if(nso == 0) continue ;
          
          /* Headings: Model and views */
          if(entete == 0) {
            entete = 1 ;
            headings[usedmodelindex] = 1 ;
            fprintf(fict,"# Model = %s\n",codename) ;
            fprintf(fict,"# Number of views = %d\n",nso) ;
            fprintf(fict,"# Numbers of components per view =") ;
            for(k = 0 ; k < nso ; k++) {
              int n = Result_GetNbOfValues(r_s + k) ;
              
              fprintf(fict," %d",n) ;
            }
            fprintf(fict,"\n") ;
            j = 1 ;
            fprintf(fict,"# Coordinates(%d)",j) ;
            j += 3 ;
            for(k = 0 ; k < nso ; k++) {
              char* name = Result_GetNameOfView(r_s + k) ;
              
              fprintf(fict," %s(%d)",name,j) ;
              j += Result_GetNbOfValues(r_s + k) ;
            }
            fprintf(fict,"\n") ;
          }
          
          /* Results per line */
          /* 1. Coordinates of node */
          for(j = 0 ; j < dim ; j++) {
            fprintf(fict,OutputFiles_RecordNumberFormat,x_s[j]) ;
          }
          for(j = dim ; j < 3 ; j++) {
            fprintf(fict,OutputFiles_RecordNumberFormat,0.) ;
          }
          
          /* 2. Components of views */
          for(k = 0 ; k < nso ; k++) {
            int n_r = Result_GetNbOfValues(r_s + k) ;
            
            for(j = 0 ; j < n_r ; j++) {
              fprintf(fict,OutputFiles_RecordNumberFormat,Result_GetValue(r_s + k)[j]) ;
            }
          }
          
          /* 3. End of line */
          fprintf(fict,"\n") ;
        }
      }
    }

  }
  
  
  /* Close the date file */
  TextFile_CloseFile(textfile) ;
}



void (OutputFiles_BackupSolutionAtPoint_)(OutputFiles_t* outputfiles,DataSet_t* dataset,double t,double t_0)
/* Backup solutions at given points in the approriate output files */
{
  Mesh_t* mesh = DataSet_GetMesh(dataset) ;
  Points_t* points = DataSet_GetPoints(dataset) ;
  unsigned short int dim = Mesh_GetDimension(mesh) ;
  int npt = Points_GetNbOfPoints(points) ;
  
  OutputFile_t* outputfile = OutputFiles_GetPointOutputFile(outputfiles) ;
  
  Result_t* r_s = Results_GetResult(OutputFiles_GetResults(outputfiles)) ;
  int    p ;
  
  OutputFile_TypeOfCurrentFile = 'p' ;
  
  
  /* Open point files for writing */
  if(t == t_0) {
    for(p = 0 ; p < npt ; p++) {
      TextFile_t* textfile = OutputFile_GetTextFile(outputfile + p) ;
      
      TextFile_OpenFile(textfile,"w") ;
    }
  }
  
  
  for(p = 0 ; p < npt ; p++) {
    TextFile_t* textfile = OutputFile_GetTextFile(outputfile + p) ;
    Point_t* point = Points_GetPoint(points) + p ;
    double *xp = Point_GetCoordinate(point) ;
    Element_t* elt = Point_GetEnclosingElement(point) ;
    FILE *ficp = TextFile_GetFileStream(textfile) ;
    
    /* First lines: version and point coordinates */
    if(t == t_0) {
      int    j ;
      time_t date ;
      time(&date) ;
    
      fprintf(ficp,"# Version " BIL_VERSION ", %s",ctime(&date)) ;
    
      fprintf(ficp,"# Point =") ;
      for(j = 0 ; j < dim ; j++) fprintf(ficp," %e",xp[j]) ;
      for(j = dim ; j < 3 ; j++) fprintf(ficp," 0") ;
      fprintf(ficp,"\n") ;
    }
    
    
    if(elt) {
      Material_t* mat = Element_GetMaterial(elt) ;
      
      if(mat) {
        char*  codename = Material_GetCodeNameOfModel(mat) ;
        int    nso ;
        int    i ;
        
        Element_FreeBuffer(elt) ;
        nso = Element_ComputeOutputs(elt,t,xp,r_s) ;
        
        if(nso > OutputFiles_MaxNbOfViews) {
          arret("BackupSolutionAtPoint: too much values") ;
        }
  
        /* Headings: Model and views */
        if(t == t_0) {
          int    j ;
          
          fprintf(ficp,"# Model = %s\n",codename) ;
          fprintf(ficp,"# Number of views = %d\n",nso) ;
          fprintf(ficp,"# Numbers of components per view =") ;
          for(i = 0 ; i < nso ; i++) {
            int n = Result_GetNbOfValues(r_s + i) ;
              
            fprintf(ficp," %d",n) ;
          }
          fprintf(ficp,"\n") ;
          j = 1 ;
          fprintf(ficp,"# Time(%d)",j) ;
          j += 1 ;
          for(i = 0 ; i < nso ; i++) {
            char* name = Result_GetNameOfView(r_s + i) ;
            
            fprintf(ficp," %s(%d)",name,j) ;
            j += Result_GetNbOfValues(r_s + i) ;
          }
          fprintf(ficp,"\n") ;
        }
  
        /* Results per line */
        /* 1. Time */
        fprintf(ficp,OutputFiles_RecordNumberFormat,t) ;
        
        /* 2. Components of views */
        for(i = 0 ; i < nso ; i++) {
          int n_r = Result_GetNbOfValues(r_s + i) ;
          int    j ;
          
          for(j = 0 ; j < n_r ; j++) {
            fprintf(ficp,OutputFiles_RecordNumberFormat,Result_GetValue(r_s + i)[j]) ;
          }
        }
        
        /* 3. End of line */
        fprintf(ficp,"\n") ;
      }
    }
    
    /* We clean the stream */
    TextFile_CleanTheStream(textfile) ;
  }
}



Views_t* (OutputFiles_CreateGlobalViews)(OutputFiles_t* outputfiles,Models_t* usedmodels,TextFile_t* textfile)
/* We use informations found in the date file textfile to build the global views and 
 * initialize the local views per used model including the index in the global views.
 * Return a pointer to Views_t.
 */
{
  Model_t*  usedmodel  = Models_GetModel(usedmodels) ;
  int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
  
  int    nbofviews ;
  char   nameofview[OutputFiles_MaxNbOfViews][OutputFiles_MaxLengthOfViewName] ;
  int    nbcompofview[OutputFiles_MaxNbOfViews] ;


  /* Open the date file for reading */
  TextFile_OpenFile(textfile,"r") ;
      
  {
    char*  pline ;
    int    imodel = -1 ;
    int    usedmodelindex ;
    
    nbofviews = 0 ;
      
    while((pline = OutputFiles_ReadLineFromCurrentFilePosition(outputfiles,textfile))) {
      
      if(pline[0] == '#') {
        
        /* Model */
        if(strstr(pline,"Model")) {
          char name[Model_MaxLengthOfKeyWord] ;
          
          pline = strchr(pline,'=') + 1 ;
          
          sscanf(pline,"%s",name) ;
          
          usedmodelindex = Models_FindModelIndex(usedmodels,name) ;
          
          imodel += 1 ;
            
        /* Nb of views */
        } else if(strstr(pline,"nombre de vues") || \
                  strstr(pline,"Number of views")) {
          int n ;
          
          pline = strchr(pline,'=') + 1 ;
          
          sscanf(pline,"%d",&n) ;
          
          /* We store the number of views in the current used model */
          {
            Model_t* model = usedmodel + usedmodelindex ;
            Views_t*  views = Model_GetViews(model) ;
            
            Views_GetNbOfViews(views) = n ;
          }
          
        /* Nb of components per view */
        } else if(strstr(pline,"nombres de composantes/vue") || \
                  strstr(pline,"Numbers of components per view")) {
          Model_t* model = usedmodel + usedmodelindex ;
          Views_t*  views = Model_GetViews(model) ;
          View_t*   view  = Views_GetView(views) ;
          int nviews = Views_GetNbOfViews(views) ;
          int i ;
          
          pline = strchr(pline,'=') + 1 ;
                    
          /* We store the number of components in the current used model */
          for(i = 0 ; i < nviews ; i++) {
            int n ;
            
            sscanf(pline,"%d",&n) ;
            
            View_GetNbOfComponents(view + i) = n ;
            
            pline = strpbrk(pline,"139") + 1 ;
          }
          
        /* Name of views */
        } else if(strstr(pline,"coordonnees") || \
                  strstr(pline,"Coordinates")) {
          Model_t* model = usedmodel + usedmodelindex ;
          Views_t*  views = Model_GetViews(model) ;
          View_t*   view  = Views_GetView(views) ;
          int nviews = Views_GetNbOfViews(views) ;
          int i ;
          
          pline = strstr(pline,"(1)") + strlen("(1)") ;
          
          /* We store the view names in the current used model */
          for(i = 0 ; i < nviews ; i++) {
            int    n = View_GetNbOfComponents(view + i) ;
            char*  name = View_GetNameOfView(view + i) ;
            char*  p ;
            int j ;
            
            sscanf(pline,"%s",name) ;
            
            if((p = strrchr(name,'('))) *p = '\0' ;
            
            /* Does this view have already been met? */
            for(j = 0 ; j < nbofviews ; j++) {
              if(!strcmp(nameofview[j],name) && (nbcompofview[j] == n)) break ;
            }
            
            /* If never met, we increment the nb of views */
            if(j == nbofviews) {
              strcpy(nameofview[j],name) ;
              nbcompofview[j] = n ;
              nbofviews += 1 ;
            }
            
            /* Here we store the index in the global views */
            View_GetGlobalIndex(view + i) = j ;
            
            pline = strchr(pline,')') + 1 ;
          }
          
          if(imodel == n_usedmodels - 1) break ;
        }
      }
    }
    
    if(imodel != n_usedmodels - 1) {
      Message_Direct("The nb of models doesn't fit that in the datafile\n") ;
      //arret("OutputFiles_PostProcessForGmshParsedFileFormat") ;
    }
    
  }
  
  TextFile_CloseFile(textfile) ;
  
  /* We build global views */
  {
    Views_t* views = Views_Create(nbofviews) ;
    View_t*  view  = Views_GetView(views) ;
    int i ;
    
    for(i = 0 ; i < nbofviews ; i++) {
      View_t* vi = view + i ;
      char* name = View_GetNameOfView(vi) ;
      
      strcpy(name,nameofview[i]) ;
      View_GetNbOfComponents(vi) = nbcompofview[i] ;
    }
    
    return(views) ;
  }
  
}
