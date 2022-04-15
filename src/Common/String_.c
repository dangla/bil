#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "Message.h"
#include "String_.h"



static  char*  String_tokens[String_MaxNbOfKeyWords] ;
static  char   String_save[String_MaxLengthOfKeyWords] ;
static  char   String_line[String_MaxLengthOfLine] ;


#if 0
char* (String_Create)(const char* filename)
{
  //String_t* string = (String_t*) Mry_New(String_t) ;
  char* str ;
  int n = 0 ;
  
  /* Nb of characters in filename */
  {
    FILE* fp = fopen(filename,"r") ;
    char c ;
  
    while((c = fgetc(fp)) != EOF) {
      n++ ;
    }
    
    fclose(fp) ;
  }
  
  {
    str = (char*) malloc(n*sizeof(char)) ;
    assert(str) ;
  
    {
      FILE* fp = fopen(filename,"r") ;
      char* c = str ;
  
      while((*c = fgetc(fp)) != EOF) c++ ;
      
      c[0] = '\0' ;
    
      fclose(fp) ;
    }
    
  }
  
  return(str) ;
}





void (String_Delete)(void* self)
{
  char* str = (char*) self ;
  
  free(str) ;
}
#endif




char* String_FindToken3(const char* str,const char* cle,const char* del)
{
  char**    tok  = String_BreakIntoTokens(cle,del) ;
  short int ntok = String_NbOfTokens(tok) ;
  char* c = String_FindToken(str,tok[0]) ;

  if(ntok > 1) {
    int itok ;
    
    for(itok = 1 ; itok < ntok ; itok++) {
      char* c0 = String_FindToken(str,tok[itok]) ;

      if(!c) c = c0 ;

      if(c0 && (c0 < c)) c = c0 ;
    }
  }
  
  return(c) ;
}




char* String_FindAndSkipToken3(const char* str,const char* cle,const char* del)
{
  char**    tok  = String_BreakIntoTokens(cle,del) ;
  short int ntok = String_NbOfTokens(tok) ;
  char* c  = String_FindToken(str,tok[0]) ;
  char* c1 = String_FindAndSkipToken(c,tok[0]) ;

  if(ntok > 1) {
    int itok ;
    
    for(itok = 1 ; itok < ntok ; itok++) {
      char* c0 = String_FindToken(str,tok[itok]) ;

      if(!c) {
        c = c0 ;
        c1 = String_FindAndSkipToken(c,tok[itok]) ;
      }

      if(c0 && (c0 < c)) {
        c = c0 ;
        c1 = String_FindAndSkipToken(c,tok[itok]) ;
      }

      /* In case of same location select the longest token */
      if(c0 && (c0 == c)) {
        char* c2 = String_FindAndSkipToken(c,tok[itok]) ;
        
        if(c2 > c1) c1 = c2 ;
      }
    }
  }
  
  return(c1) ;
}




char* String_FindNthToken3(const char* str,const char* tok,const int n)
{
  char* c = String_FindToken(str,tok) ;
  int i = 1 ;
    
  while(c && i++ < n) {
    c = String_FindAndSkipToken(c,tok) ;
    if(c) c = String_FindToken(c,tok) ;
  }
  
  return(c) ;
}



#if 1
char* String_FindNthToken4(const char* str,const char* cle,const char* del,const int n)
{
  char**    tok  = String_BreakIntoTokens(cle,del) ;
  short int ntok = String_NbOfTokens(tok) ;
  char*  line = NULL ;
  
  {
    int itok = 0 ;
    
    while((itok < ntok) && !(line = String_FindNthToken(str,tok[itok],n))) itok++ ;
  }
  
  return(line) ;
}
#endif



#if 0
char* String_FindNthToken4(const char* str,const char* cle,const char* del,const int n)
{
  char* c = String_FindToken(str,cle,del) ;
  int i = 1 ;
    
  while(c && i++ < n) {
    c = String_FindAndSkipToken(c,cle,del) ;
    if(c) c = String_FindToken(c,cle,del) ;
  }
  
  return(c) ;
}
#endif




int String_CountTokens2(const char* str,const char* tok)
{
  const char* c = str ;
  int i = 0 ;

  while((c = String_FindAndSkipToken(c,tok))) i++ ;

  return(i) ;
}




int String_CountTokens3(const char* str,const char* cle,const char* del)
{
  char**    tok  = String_BreakIntoTokens(cle,del) ;
  short int ntok = String_NbOfTokens(tok) ;
  int     n = 0 ;
  
  /* Compute the nb of times we find the tokens in string */
  {
    int itok ;
    
    for(itok = 0 ; itok < ntok ; itok++) {
      n += String_CountTokens(str,tok[itok]) ;
    }
  }
  
  return(n) ;
}




int String_CountTokensAloneInOneLine2(const char* str,const char* tok)
{
  const char* c = str ;
  int i = 0 ;

  while((c = String_FindAndSkipToken(c,tok)) && (String_SkipBlankChars(c)[0] == '\n')) i++ ;

  return(i) ;
}




int String_CountTokensAloneInOneLine3(const char* str,const char* cle,const char* del)
{
  char**    tok  = String_BreakIntoTokens(cle,del) ;
  short int ntok = String_NbOfTokens(tok) ;
  int     n = 0 ;
  
  /* Compute the nb of times we find the tokens in string */
  {
    int itok ;
    
    for(itok = 0 ; itok < ntok ; itok++) {
      n += String_CountTokensAloneInOneLine(str,tok[itok]) ;
    }
  }
  
  return(n) ;
}





char** String_BreakIntoTokens(const char* str,const char* del)
{
  if(strlen(str) > String_MaxLengthOfKeyWords - 1) {
    arret("String_BreakIntoTokens: increase the length of String_save") ;
  }
  
  strcpy(String_save,str) ;
  
  if(!(String_tokens[0] = strtok(String_save,del))) return(NULL) ;
  
  {
    int ntok = 1 ;
  
    while((String_tokens[ntok] = strtok(NULL,del))) {
      ntok++ ;
  
      if(ntok > String_MaxNbOfKeyWords - 1) {
        arret("String_BreakIntoTokens: increase the length of String_tokens") ;
      }
    }
  }
  
  return(String_tokens) ;
}



int String_NbOfTokens(char** tok)
{    
  int ntok = 0 ;
  
  if(tok) {
    while(tok[ntok]) ntok++ ;
  }
  
  return(ntok) ;
}



char* String_CopyLine(const char* str)
{
  int   len = (str) ? strlen(str) : 0 ;
  char* eol = String_FindEndOfLine(str) ;
  int   n   = (eol) ? eol - str : len ;
  
  if(n > String_MaxLengthOfLine - 1) {
    arret("String_CopyLine: increase the length of String_line") ;
  }
  
  strncpy(String_line,str,n) ;
  String_line[n] = '\0' ;
  
  return(String_line) ;
}




const char* String_SkipRemainingComments(const char* str)
{
  const char* c = str ;
  
  while(*c != EOF) {
    if(c[-1] == '*' && c[0] == '/') return(c+1) ;
    c++ ;
  }

  return(str) ;
}



int String_NbOfUncommentedLines(const char* str,const char* cmt)
{
  int n = 0 ;
  
  {
    const char* c = str ;
    
    while(c && (*c != EOF)) {
      while(String_BeginsWithAnyChar(c,cmt)) c = String_SkipLine(c) ;
      c = String_SkipLine(c) ;
      n++ ;
    }
  }
  
  return(n) ;
}



int    (String_FindPositionIndex)(const char* str,const char* const* ss,const int n)
/** Return the position index in ss whose name is pointed to by str 
 *  or -1 if it fails. */
{
  int    i ;

  if(isdigit(str[0])) {
    
    i  = atoi(str) - 1 ;
    
  } else {
    
    for(i = 0 ; i < n ; i++) {
      if(String_Is(str,ss[i])) break ;
    }
    
    if(i == n) i = -1 ;
  }

  return(i) ;
}


char*   (String_RemoveComments)(char* src,char* dest)
/**  Remove the comments from src and copy it to dest.
 *   The pointer dest should point to an allocated space 
 *   large enough to accomodate the uncommented src. 
 *   We can give the same pointer, i.e. dest = src,
 *   in that case the string src will be modified.
 *   Return the pointer dest. */
{
  char* cin = src ;
  char* cou = dest ;
  
  if(src) {
    while(cin[0]) {
    
      if(String_BeginsWithSingleLineComment(cin)) {
      
        cin = String_FindEndOfLine(cin) ;
      
      } else if(String_BeginsWithMultiLineComment(cin)) {
        char* c = String_SkipMultiLineComment(cin) ;
      
        if(c) {
          cin = c ;
        } else {
          arret("String_RemoveComments: comment doesn't end") ;
        }
      }
    
      *cou++ = *cin++ ;
    }
  
    cou[0] = 0 ;
  }
  
  return dest ;
}



#if 0
#include "TextFile.h"

static int String_Test(int, char**) ;

int String_Test(int argc, char** argv)
{
  char* filename = argv[1] ;
  TextFile_t* textfile = TextFile_Create(filename) ;
  //char* str = String_Create(filename) ;
  char* str = TextFile_StoreFileContent(textfile) ;
  
  if(argc > 2) {
    char* tok = argv[2] ;
    int n ;
    char* c ;
    
    /* Test String_CountTokens */
    {
      n = String_CountTokens(str,tok) ;
    
      printf("nb of tokens = %d\n",n) ;
    }
  
    /* Test String_FindToken */
    {
      c = String_FindToken(str,tok) ;
    
      printf("token = %s\n%s\n",tok,c) ;
    }
  
    /* Test String_FindNthToken */
    {
      int j ;
      
      for(j = 0 ; j < n ; j++) {
        c = String_FindNthToken(str,tok,j+1) ;
    
        printf("token %d = %s\n%s\n",j,tok,c) ;
      }
    }
  }
    
  
  if(argc > 1) {
    /* Test String_RemoveComments */
    {
      printf("File content before removing the comments\n") ;
      printf("%s",str) ;
      
      String_RemoveComments(str,str) ;
      
      printf("File content after removing the comments\n") ;
      printf("%s",str) ;
    }
  }

  
  return(0) ;
}



int main(int argc, char** argv)
{
  return(String_Test(argc,argv)) ;
}
#endif








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
