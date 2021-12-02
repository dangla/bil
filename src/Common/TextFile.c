#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
//#include <stdarg.h>
#include <assert.h>

#include "Message.h"
#include "Tools/Math.h"
#include "Buffer.h"
#include "Mry.h"
#include "TextFile.h"


TextFile_t*   (TextFile_Create)(const char* filename)
{
  TextFile_t* textfile   = (TextFile_t*) Mry_New(TextFile_t) ;
  

  /* Initialization */
  {
    TextFile_GetFileStream(textfile) = NULL ;
  }
  

  /* Memory space for the file name */
  {
    char* name = (char*) Mry_New(char[TextFile_MaxLengthOfFileName]) ;
    
    if(filename) {
      if(strlen(filename) > TextFile_MaxLengthOfFileName) {
        arret("TextFile_Create(3)") ;
      }
      strcpy(name,filename) ;
    }
    
    TextFile_GetFileName(textfile) = name ;
  }
  
  
  /* Memory space for the file positions */
  {
    fpos_t* pos = (fpos_t*) Mry_New(fpos_t) ;
    
    TextFile_GetFilePosition(textfile) = pos ;
  }
  
  
  /* The pointer to the file content, intialized to NULL by default. */
  {
    TextFile_GetFileContent(textfile) = NULL ;
    TextFile_GetPreviousPositionInString(textfile) = 0 ;
    TextFile_GetCurrentPositionInString(textfile) = 0 ;
  }
  
  return(textfile) ;
}


void (TextFile_Delete)(void* self)
{
  TextFile_t* textfile = (TextFile_t*) self ;
  
  TextFile_CloseFile(textfile) ;
  
  {
    char* name = TextFile_GetFileName(textfile) ;
    
    if(name) {
      free(name) ;
    }
    
    TextFile_GetFileName(textfile) = NULL ;
  }
  
  {
    fpos_t* pos = TextFile_GetFilePosition(textfile) ;
    
    if(pos) {
      free(pos) ;
    }
    
    TextFile_GetFilePosition(textfile) = NULL ;
  }
  
  {
    char* c = TextFile_GetFileContent(textfile) ;
    
    if(c) {
      free(c) ;
    }
    
    TextFile_GetFileContent(textfile) = NULL ;
  }
}



int (TextFile_Exists)(TextFile_t* textfile)
{
  char* filename = TextFile_GetFileName(textfile) ;
  FILE* str ;
  
  /* After C, we can have multiple streams pointing to 
   * the same file open at the same time. So we open a
   * stream without worrying and close it afterwards. */
  str = fopen(filename,"r") ;
  
  if(!str) {
    return(0) ;
  }
  
  fclose(str) ;
  
  return(1) ;
}


FILE* (TextFile_OpenFile)(TextFile_t* textfile,const char* mode)
{
  char* filename = TextFile_GetFileName(textfile) ;
  FILE* str ;
  
  TextFile_CloseFile(textfile) ;
  
  str = fopen(filename,mode) ;
  
  TextFile_GetFileStream(textfile) = str ;
    
  if(!str) {
    Message_RuntimeError("TextFile_OpenFile: failed to open %s\n",filename) ;
  }
  
  return(str) ;
}


void (TextFile_CloseFile)(TextFile_t* textfile)
{
  FILE* str = TextFile_GetFileStream(textfile) ;
  
  if(!str) return ;
  
  fclose(str) ;
  
  TextFile_GetFileStream(textfile) = NULL ;
  TextFile_GetPreviousPositionInString(textfile) = 0 ;
  TextFile_GetCurrentPositionInString(textfile) = 0 ;
}


void (TextFile_CleanTheStream)(TextFile_t* textfile)
{
  FILE* str = TextFile_GetFileStream(textfile) ;
  
  if(!str) return ;
  
  fflush(str) ;
}


void (TextFile_StoreFilePosition)(TextFile_t* textfile)
/** Store the current file position of the stream for further reading. */
{
  FILE* str  = TextFile_GetFileStream(textfile) ;
  fpos_t* pos = TextFile_GetFilePosition(textfile) ;
  
  /* Store the current file position of the stream in pos */
  if(fgetpos(str,pos)) {
    arret("TextFile_StoreFilePosition") ;
  }
  
  TextFile_GetPreviousPositionInString(textfile) = TextFile_GetCurrentPositionInString(textfile) ;
  
}



void (TextFile_MoveToStoredFilePosition)(TextFile_t* textfile)
/** Set the file position of the stream to stored value. */
{
  FILE* str  = TextFile_GetFileStream(textfile) ;
  fpos_t* pos = TextFile_GetFilePosition(textfile) ;
  
  /* Set the file position of the stream to the stored position */
  if(fsetpos(str,pos)) {
    arret("TextFile_MoveToStoredFilePosition") ;
  }
  
  TextFile_GetCurrentPositionInString(textfile) = TextFile_GetPreviousPositionInString(textfile) ;
}



char* (TextFile_ReadLineFromCurrentFilePositionInString)(TextFile_t* textfile,char* line,int n)
/** Reads a line from the textfile at the current position of its string
 *  and stores it into the string pointed to by line. It stops when 
 *  either (n-1) characters are read, the newline character is read,
 *  or the end-of-file is reached, whichever comes first. 
 *  Return a pointer to the string line on success or NULL on failure. */
{
  char* c ;
  char* beg = TextFile_GetFileContent(textfile) ;
  char* end = (beg) ? String_FindEndOfString(beg) : NULL ;
  
  if(!beg) {
    arret("TextFile_ReadLineFromCurrentFilePositionInString") ;
  }
  
  
  /* Reads a non empty line from the current position of the string */
  do {
    char* cur = TextFile_GetCurrentPositionInFileContent(textfile) ;
    char* eol = String_FindEndOfLine(cur) ;
    char* nex = (eol) ? eol + 1 : end ;
    
    if(nex == cur) {
      return(NULL) ;
    }
    
    {
      int n0 = nex - cur ;
      int n1 = Math_Min(n0,n) ;
      char* cur1 = cur + n1 ;
      
      strncpy(line,cur,n1) ;
      c = line ;
      line[n1] = '\0' ;
    
      TextFile_GetCurrentPositionInString(textfile) = cur1 - beg ;
    }
    
    /* Eliminate the first blank characters */
    //if(*c == ' ') c += strspn(c," ") ;
    c = String_SkipBlankChars(c) ;
    
  //} while((*c == '\n')) ;
  } while(isspace(*c)) ;
  
  //String_FindEndOfLine(line)[0] = '\0' ;
  
  return(c) ;
}



char* (TextFile_ReadLineFromCurrentFilePosition)(TextFile_t* textfile,char* line,int n)
/** Reads a line from the stream of textfile at the current position 
 *  and stores it into the string pointed to by line. It stops when 
 *  either (n-1) characters are read, the newline character is read,
 *  or the end-of-file is reached, whichever comes first. 
 *  The file of textfile is assumed open for reading. 
 *  Return a pointer to the string line on success or NULL on failure. */
{
  char* c ;
  
  
  /* Reads a non empty line from the stream after the current position */
  {
    FILE* str  = TextFile_GetFileStream(textfile) ;
  
    do {
      if(feof(str)) return(NULL) ;
    
      c = fgets(line,n,str) ;
    
      if(!c) {
        return(NULL) ;
      }
    
      /* Eliminate the first blank characters */
      if(*c == ' ') c += strspn(c," ") ;
    
    /* } while((*c == '\n') || (*c == '#')) ; */
    //} while((*c == '\n') || (*c == '\r') || (*c == '\t') || (*c == '\v') || (*c == '\f')) ;
    } while(isspace(*c)) ;
    /* } while(0) ; */
  }
  
  return(c) ;
}



long int TextFile_CountNbOfEatenCharacters(TextFile_t* textfile)
{
  long int count = 0 ;
  
  if(TextFile_Exists(textfile)) {
    fpos_t* pos  = TextFile_GetFilePosition(textfile) ;
    fpos_t* cpos = NULL ;
    FILE* str = TextFile_OpenFile(textfile,"r") ;
    char c ;
  
    while((cpos != pos) && ((c = fgetc(str)) != EOF)) {
    
      /* Store the current file position of the stream in pos */
      if(fgetpos(str,cpos)) {
        arret("TextFile_CountNbOfEatenCharacters") ;
      }
    
      count++ ;
    
    }
  
    TextFile_CloseFile(textfile) ;
  }
  
  return(count) ;
}



long int TextFile_CountNbOfCharacters(TextFile_t* textfile)
/** Return the number of character of the file including end of file */
{
  long int count = 0 ;
  
  if(TextFile_Exists(textfile)) {
    FILE* str = TextFile_OpenFile(textfile,"r") ;
    char c ;
    
    count = 1 ;
  
    while((c = fgetc(str)) != EOF) {
    
      count++ ;
    
    }
  
    TextFile_CloseFile(textfile) ;
  }
  
  return(count) ;
}


int TextFile_CountTheMaxNbOfCharactersPerLine(TextFile_t* textfile)
/** Return the max number of character per line including end of line */
{
  int linelength = 0 ;
  
  if(TextFile_Exists(textfile)) {
    FILE* str = TextFile_OpenFile(textfile,"r") ;
    int ll = 0 ;
    char c ;
  
    while((c = fgetc(str)) != EOF) {
    
      if(c != '\n') {
        ll++ ;
      } else {
        if(ll > linelength) linelength = ll ;
        ll = 0 ;
      }
    }
  
    /* Include end of line character */
    linelength += 1 ;
  
    TextFile_CloseFile(textfile) ;
  }
  
  return(linelength) ;
}




char* (TextFile_FileCopy)(TextFile_t* textfile)
/** Copy the file in a temporary file
 *  Return its name. */
{
  char*   targetfile = tmpnam(NULL) ; /* temporary filename */
  FILE*   target  = fopen(targetfile,"w") ;
  //char    targetfile[] = "tmpXXXXXX" ;
  //int     descriptor   = mkstemp(targetfile) ;
  //FILE*   target  = (descriptor != -1) ? fopen(targetfile,"w") : NULL ;
  
  if(!target) {
    arret("TextFile_FileCopy(2)") ;
  }
  
  {
    char c ;
    FILE*   source  = TextFile_OpenFile(textfile,"r") ;
    
    while( (c = fgetc(source)) != EOF ) fputc(c,target) ;
  }
    
  TextFile_CloseFile(textfile) ;
  fclose(target) ;
  
  return(targetfile) ;
}


FILE* (TextFile_FileStreamCopy)(TextFile_t* textfile)
/** Copy the file in a temporary binary file
 *  Return the stream of the temporary file */
{
  FILE*   target = tmpfile() ; /* temporary file */
  FILE*   source = TextFile_OpenFile(textfile,"r") ;
  
  {
    char c ;
    
    while( (c = fgetc(source)) != EOF ) fputc(c,target) ;
  }
    
  TextFile_CloseFile(textfile) ;

  rewind(target) ;
  
  return(target) ;
}



char* (TextFile_StoreFileContent)(TextFile_t* textfile)
/** Allocate memory space for the file content in char type
 *  and copy the file content in it. 
 *  Return the pointer to the allocated space. */
{
  /* Free the space allocated for the content if any */
  {
    char* content = TextFile_GetFileContent(textfile) ;
    
    free(content) ;
    
    TextFile_GetFileContent(textfile) = NULL ;
  }
  
  /* Allocate the memory space for the content */
  {
    int n = TextFile_CountNbOfCharacters(textfile) ;
    
    if(n) {
      char* content = (char*) Mry_New(char[n]) ;
  
      TextFile_GetFileContent(textfile) = content ;
    }
  }
  
  /* Read the chars and fill in the allocated space */
  {
    char* c = TextFile_GetFileContent(textfile) ;
    
    if(c && TextFile_Exists(textfile)) {
      FILE* str  = TextFile_OpenFile(textfile,"r") ;
      int i ;
    
      //while( (*c = fgetc(str)) != EOF ) c++ ;
      while( (i = fgetc(str)) != EOF ) {
        c[0] = i ;
        c++ ;
      }
    
      c[0] = '\0' ;
    
      TextFile_CloseFile(textfile) ;
    }
  }

  return(TextFile_GetFileContent(textfile)) ;
}




/* From COS */
#if 0

static int rmlit = 0;

static void echo_1(FILE *in, FILE *out)
{
  fputc(fgetc(in), out);
}

static void echo_upto(FILE *in, FILE *out, int chr)
{
  int c;
 
  while ((c = fgetc(in)) != EOF) {
    fputc(c, out);
    if (c == '\\') { echo_1(in, out); continue; }
    if (c == chr ) return;
  }
}

static void skip_1(FILE *in)
{
  fgetc(in);
}

static void skip_upto(FILE *in, int chr)
{
  int c;
 
  while ((c = fgetc(in)) != EOF) {
    if (c == '\\') { skip_1(in); continue; }
    if (c == chr ) return;
  }
}

static void skip_line(FILE *in)
{
  int c;
 
  while ((c = fgetc(in)) != EOF && c != '\n') ;

  ungetc(c, in);
}

static void skip_cmt(FILE *in)
{
  int c, p;

  for(p = 0; (c = fgetc(in)) != EOF; p = c)
    if (c == '/' && p == '*')
      break;

  ungetc(' ', in);
}

static void remove_cmt(FILE *in , FILE *out)
{
  int c;
 
  while ((c = fgetc(in)) != EOF) {
    switch(c) {
      case '\'':
      case '"' : if (rmlit) skip_upto(in, c);
                 else { fputc(c, out); echo_upto(in, out, c); }
                 break;
      case '\\': fputc(c, out); echo_1(in, out); break;
      case '/' : c = fgetc(in);
                 if (c == '/') { skip_line(in); break; }
                 if (c == '*') { skip_cmt (in); break; }
                 fputc('/', out);
      default  : fputc(c  , out);
    }
  }
}

static void help(void)
{
  fprintf(stderr, "usage: coscmt [-h] [-l] [-o outfile] [infiles]\n");
  exit(EXIT_FAILURE);
}
 
int main(int argc, char *argv[])
{
  FILE *in  = stdin;
  FILE *out = stdout;
  int i, o = 0;
 
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (argv[i][1] == 0)
         continue;

      if (argv[i][1] == 'l')
        { argv[i][0] = 0; rmlit = 1; continue; }

      if (argv[i][1] == 'o')
        { argv[i][0] = 0; o = ++i; continue; }
         
      help();
    }
  }

  if (o) {
    if (argv[o][0] == '-')
      out = stdout;
      
    else {
      out = fopen(argv[o],"w");
      if (!out) {
        fprintf(stderr, "coscmt: unable to open file %s in write mode\n", argv[o]);
        exit(EXIT_FAILURE);
      }
    }

    argv[o][0] = 0;
  }

  for (i = 1; i < argc; i++) {
    if (argv[i][0] == 0)
      continue;

    if (argv[i][0] == '-')
      in = stdin;

    else {
      in = fopen(argv[i],"r");
      if (!in) {
        fprintf(stderr, "coscmt: unable to open file %s in read mode\n", argv[i]);
        exit(EXIT_FAILURE);
      }
    }

    remove_cmt(in, out);

    if (in != stdin)
      fclose(in);
  }

  if (out != stdout)
    fclose(out);

  return EXIT_SUCCESS;
}
#endif
