#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "PosFilesForGMSH.h"
#include "Message.h"
#include "Mesh.h"
#include "BilVersion.h"
#include "TextFile.h"
#include "Models.h"
#include "String_.h"
#include "Mry.h"


static void (PosFilesForGMSH_ASCIIFileFormatVersion2_2)(PosFilesForGMSH_t*) ;
static void (PosFilesForGMSH_ParsedFileFormatVersion2)(PosFilesForGMSH_t*) ;



PosFilesForGMSH_t*   (PosFilesForGMSH_Create)(DataSet_t* dataset)
{
  PosFilesForGMSH_t* pf4gmsh = (PosFilesForGMSH_t*) Mry_New(PosFilesForGMSH_t) ;
  
  PosFilesForGMSH_GetDataSet(pf4gmsh) = dataset ;
  
  {
    DataFile_t* datafile = DataSet_GetDataFile(dataset) ;
    char* filename = DataFile_GetFileName(datafile) ;
    int n_dates = Dates_GetNbOfDates(DataSet_GetDates(dataset)) ;
    int n_points = Points_GetNbOfPoints(DataSet_GetPoints(dataset)) ;
    OutputFiles_t* outputfiles = OutputFiles_Create(filename,n_dates,n_points) ;
    
    PosFilesForGMSH_GetOutputFiles(pf4gmsh) = outputfiles ;
  
    PosFilesForGMSH_GetBilVersion(pf4gmsh) = OutputFiles_Version(outputfiles) ;
  }

  return(pf4gmsh) ;
}



void   (PosFilesForGMSH_Delete)(void* self)
{
  PosFilesForGMSH_t* pf4gmsh = (PosFilesForGMSH_t*) self ;
  
  {
    OutputFiles_t* outputfiles = PosFilesForGMSH_GetOutputFiles(pf4gmsh) ;
    
    if(outputfiles) {
      OutputFiles_Delete(outputfiles) ;
      free(outputfiles) ;
    }
  }
}



void (PosFilesForGMSH_ASCIIFileFormat)(PosFilesForGMSH_t* pf4gmsh)
{
  double version = PosFilesForGMSH_GetBilVersion(pf4gmsh) ;
  
  if(version < 2.2) {
    DataSet_t* dataset = PosFilesForGMSH_GetDataSet(pf4gmsh) ;
    Materials_t* materials = DataSet_GetMaterials(dataset) ;
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
    
    if(n_usedmodels > 1) {
      arret("PosFilesForGMSH_ASCIIFileFormat") ;
    }
    
    PosFilesForGMSH_ASCIIFileFormatVersion2_2(pf4gmsh) ;
  } else {
    PosFilesForGMSH_ASCIIFileFormatVersion2_2(pf4gmsh) ;
  }
}



void (PosFilesForGMSH_ASCIIFileFormatVersion2_2)(PosFilesForGMSH_t* pf4gmsh)
/* Create gmsh post-processing files in MSH ASCII file format version 2.2 */
{
  DataSet_t* dataset = PosFilesForGMSH_GetDataSet(pf4gmsh) ;
  OutputFiles_t* outputfiles = PosFilesForGMSH_GetOutputFiles(pf4gmsh) ;
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
    arret("PosFilesForGMSH_ASCIIFileFormatVersion2_2: too many views") ;
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
      arret("PosFilesForGMSH_ASCIIFileFormatVersion2_2: nom trop long") ;
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
            arret("PosFilesForGMSH_ASCIIFileFormatVersion2_2") ;
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



void (PosFilesForGMSH_ParsedFileFormat)(PosFilesForGMSH_t* pf4gmsh)
{
  double version = PosFilesForGMSH_GetBilVersion(pf4gmsh) ;
  
  if(version < 2.2) {
    DataSet_t* dataset = PosFilesForGMSH_GetDataSet(pf4gmsh) ;
    Materials_t* materials = DataSet_GetMaterials(dataset) ;
    Models_t* usedmodels = Materials_GetUsedModels(materials) ;
    int n_usedmodels = Models_GetNbOfModels(usedmodels) ;
    
    if(n_usedmodels > 1) {
      arret("PosFilesForGMSH_ParsedFileFormat") ;
    }
    
    PosFilesForGMSH_ParsedFileFormatVersion2(pf4gmsh) ;
  } else {
    PosFilesForGMSH_ParsedFileFormatVersion2(pf4gmsh) ;
  }
}




void (PosFilesForGMSH_ParsedFileFormatVersion2)(PosFilesForGMSH_t* pf4gmsh)
/* Create gmsh post-processing files in parsed file format */
{
  DataSet_t* dataset = PosFilesForGMSH_GetDataSet(pf4gmsh) ;
  OutputFiles_t* outputfiles = PosFilesForGMSH_GetOutputFiles(pf4gmsh) ;
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
    arret("PosFilesForGMSH_ParsedFileFormatVersion2: too many views") ;
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
    char* filename = OutputFiles_GetDataFileName(outputfiles) ;
    int n = ceil(log10((double) nbofglobalviews)) ;
    int i = strlen(filename) + 4 + n ;
    
    if(i > OutputFiles_MaxLengthOfFileName) {
      arret("PosFilesForGMSH_ParsedFileFormatVersion2: too lengthy name") ;
    }
    
    for(i = 0 ; i < nbofglobalviews ; i++) {
      char   nom_pos[OutputFiles_MaxLengthOfFileName] ;
      FILE* ficp ;

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
      int    nvert = n_vertices[dim-1][nn-1] ;
      int    i_temps ;
      
      if(nvert == 0) {
        arret("PosFilesForGMSH_ParsedFileFormatVersion2: unknown type") ;
      }
      
      if(usedmodelindex < 0) {
        arret("PosFilesForGMSH_ParsedFileFormatVersion2") ;
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
            int    nc = View_GetNbOfComponents(view + i) ;
            int    jview = View_GetGlobalIndex(view + i) ;
            FILE*  ficp = fileofview[jview] ;
            int j ;

            /* Type d'enregistrement */
            fprintf(ficp,"%c",svt[nc - 1]) ;
            fprintf(ficp,"%c",type[dim-1][nvert-1]) ;
            //if(nvert == nn/2) fprintf(ficp,"2") ; /* element P2 */

            /* Coordonnees */
            fprintf(ficp,"(") ;
            
            for(j = 0 ; j < 3*nvert ; j++) {
              if(j > 0) fprintf(ficp,",") ;
              fprintf(ficp,"%e",x_e[j]) ;
            }
            
            fprintf(ficp,")") ;

            /* Valeurs */
            fprintf(ficp,"{") ;
            
            for(j = 0 ; j < nc*nvert ; j++) {
              if(j > 0) fprintf(ficp,",") ;
              fprintf(ficp,"%e",val[i][j]) ;
            }
          }
        
        } else {
          int i ;
          
          for(i = 0 ; i < nviews ; i++) {
            int    nc = View_GetNbOfComponents(view + i) ;
            int    jview = View_GetGlobalIndex(view + i) ;
            FILE   *ficp = fileofview[jview] ;
            int j ;
          
            //for(j = 0 ; j < nc*nn ; j++) {
            for(j = 0 ; j < nc*nvert ; j++) {
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
    
    
    /* Times */
    {
      int i ;
      
      for(i = 0 ; i < nbofglobalviews ; i++) {
        FILE   *ficp = fileofview[i] ;
    
        fprintf(ficp,"TIME{\n") ;
        
        {
          Dates_t* dates = DataSet_GetDates(dataset) ;
          int nbofdates  = Dates_GetNbOfDates(dates) ;
          Date_t*    date   = Dates_GetDate(dates) ;
          int j ;
          
          for(j = 0 ; j < nbofdates ; j++) {
            double t = Date_GetTime(date + j) ;
            
            if(j > 0) fprintf(ficp,",") ;
            fprintf(ficp,"%e",t) ;
          }
        }
        
        fprintf(ficp,"};\n") ;
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
