#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <map>
#include <iostream>
#include <limits>
#include "core.h"


using namespace std;
typedef map<string,long int> KeyMap;







//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//


// Generic options
bool VERBOSE;
bool SCIENTIFIC;
float CUTOFF;
int DECIMALS;
bool UPPER;
bool GREATER;
bool EQUAL;
bool NAN_INIT;



// Global variables
const int BUFFER_SIZE = 100000000;








//-----Usage----------
//
void Usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: \n");
  fprintf(stderr, "  table OPERATION [OPTIONS] ENTRY_FILE\n"); 
  fprintf(stderr, "\n");
  fprintf(stderr, "OPERATIONS: \n");
  fprintf(stderr, "  * -c      create table from (key1,key2,value) entries\n");
  fprintf(stderr, "  * -x      extract (key1,key2,value) entries from table\n");
  fprintf(stderr, "  * -x0     extract (key1,key2,value) entries from tab-delimited table\n");
  fprintf(stderr, "  * -r      convert tab-delimited table to standard format\n");
  fprintf(stderr, "  * -rnorm  normalize rows\n");
  fprintf(stderr, "  * -cnorm  normalize columns\n");
  fprintf(stderr, "\n");
}






//-----InitCmdLine----------------
//
CmdLine *InitCmdLine(char *method)
{
  CmdLine *cmd_line = new CmdLine(); 
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-n", &DECIMALS, 2, "number of decimal points");
  cmd_line->AddOption("-s", &SCIENTIFIC, false, "scientific format");

  if (strcmp(method,"-c")==0) {
    cmd_line->AddOption("-nan", &NAN_INIT, false, "initialize to NaN");
  }
  else if (strcmp(method,"-x")==0) {
    cmd_line->AddOption("-c", &CUTOFF, 1, "cutoff value");
    cmd_line->AddOption("-g", &GREATER, false, "extract if greater that cutoff");
    cmd_line->AddOption("-e", &EQUAL, false, "extract if equal to cutoff");
    cmd_line->AddOption("-U", &UPPER, false, "extract upper triangle only");
  }
  else if (strcmp(method,"-x0")==0) {

  }
  else if (strcmp(method,"-r")==0) {

  }
  else {
    fprintf(stderr, "Unknown method '%s'!\n", method);
    exit(1);
  }
  return cmd_line;
}






//-----ReadToken----------------
//
char *ReadToken(char **pBUF, char DELIM)
{
  char *BUF = *pBUF;
  if (BUF==NULL) { fprintf(stderr, "Cannot read token from null buffer!\n"); exit(1); }
  while (BUF&&(BUF[0]==' ')) BUF++;                  // remove leading spaces
  int k = 0;
  for ( ; BUF[k]!=0; k++) if (BUF[k]==DELIM) break;          // read characters until delimiter or end of string is found
  if (BUF[k]==0) *pBUF = BUF+k;
  else { BUF[k] = 0; *pBUF = BUF+k+1; }
  return BUF;
}





//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  CmdLine *cmdLine = NULL;
  if (argc<2) { Usage(); exit(1); }
  char *method = argv[1];
  cmdLine = InitCmdLine(method);
  int next_arg = cmdLine->Read(argv+1,argc-1) + 1;

  // setup output format string
  char fmt[100];
  if (SCIENTIFIC) sprintf(fmt, "%%.%ie", DECIMALS); 
  else sprintf(fmt, "%%.%if", DECIMALS);

  
  //--------------------------------------------------------------------------------------//
  // OPTION -c: create table                                                              //
  //--------------------------------------------------------------------------------------//

  if (strcmp(method,"-c")==0) {
    if (argc<2) { cmdLine->Usage(argv[0],"-c [OPTIONS] ENTRY_FILE|<ENTRIES>"); exit(1); }
    MESSAGES(VERBOSE); 

    // initialize
    char *INPUT = next_arg>=argc ? NULL : argv[next_arg];
    FileBufferText *buffer;
    long int n_entries; 
    if (INPUT==NULL) buffer = new FileBufferText(LoadStdIn(&n_entries),BUFFER_SIZE);
    else { buffer = new FileBufferText(INPUT,BUFFER_SIZE); n_entries = buffer->CountLines(); }
    KeyMap key1_map, key2_map;
    long int n_key1 = 0, n_key2 = 0;
    
    // parse keys
    Progress PRG("Reading keys...",n_entries);
    for (long int n=0; n<n_entries; n++) {
      char *inp = buffer->Next();
      char *key1 = GetNextToken(&inp," \t");
      if (key1_map.find(key1)==key1_map.end()) key1_map[key1] = n_key1++;
      char *key2 = GetNextToken(&inp," \t");
      if (key2_map.find(key2)==key2_map.end()) key2_map[key2] = n_key2++;
      PRG.Check();
    }
    PRG.Done();
    if (VERBOSE) fprintf(stderr, "* Found %ld keys #1 and %ld keys #2.\n", n_key1, n_key2);

    // create table
    float **T;
    ALLOCATE1D(T,n_key1,float *);
    if (NAN_INIT==true) for (long int k1=0; k1<n_key1; k1++) { ALLOCATE1D_INIT(T[k1],n_key2,float,numeric_limits<float>::quiet_NaN()); }
    else for (long int k1=0; k1<n_key1; k1++) { ALLOCATE1D_INIT(T[k1],n_key2,float,0.0); }
    char **KEY1STR, **KEY2STR;
    ALLOCATE1D_INIT(KEY1STR,n_key1,char *,NULL);
    ALLOCATE1D_INIT(KEY2STR,n_key2,char *,NULL);
    buffer->Reset();
    Progress PRG2("Creating table...",n_entries);
    for (long int n=0; n<n_entries; n++) {
      char *inp = buffer->Next();
      char *key1 = GetNextToken(&inp," \t");
      char *key2 = GetNextToken(&inp," \t");
      float value = atof(GetNextToken(&inp,"\n"));
      long int x1 = key1_map[key1];
      long int x2 = key2_map[key2];
      if (KEY1STR[x1]==NULL) KEY1STR[x1] = StrCopy(key1);
      if (KEY2STR[x2]==NULL) KEY2STR[x2] = StrCopy(key2);
      if ((NAN_INIT==false)||(T[x1][x2]==T[x1][x2])) T[x1][x2] += value;
      else T[x1][x2] = value;
      //printf("[%ld][%ld] += %f => %f", key1_map[key1], key2_map[key2], value, T[key1_map[key1]][key2_map[key2]]); getchar();
      PRG2.Check();
    }
    PRG2.Done();
    
    // print table
    Progress PRG3("Printing...",n_key1);
    printf("\t"); for (long int k2=0; k2<n_key2; k2++) printf("%s ", KEY2STR[k2]); printf("\n");
    for (long int k1=0; k1<n_key1; k1++) {
      printf("%s\t", KEY1STR[k1]);
      for (long int k2=0; k2<n_key2; k2++) { printf(fmt, T[k1][k2]); printf(" "); }
      printf("\n");
      PRG3.Check();
    }
    PRG3.Done();

    
    // cleanup
    FREE2D(T,n_key1);
    FREE2D(KEY1STR,n_key1);
    FREE2D(KEY2STR,n_key2);
    delete buffer;
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -x: extract table entries                                                     //
  //--------------------------------------------------------------------------------------//

  else if (strcmp(method,"-x")==0) {
    if (argc<2) { cmdLine->Usage(argv[0],"-x [OPTIONS] TABLE_FILE|<TABLE>"); exit(1); }
    MESSAGES(VERBOSE); 

    // initialize
    char *INPUT = next_arg>=argc ? NULL : argv[next_arg];
    FileBufferText *buffer;
    long int n_lines; 
    if (INPUT==NULL) buffer = new FileBufferText(LoadStdIn(&n_lines),BUFFER_SIZE);
    else { buffer = new FileBufferText(INPUT,BUFFER_SIZE); n_lines = buffer->CountLines(); }
    if (n_lines>1) {
      // read header
      long int n_key1 = n_lines - 1;
      char *inp = buffer->Next();
      GetNextToken(&inp,'\t');
      long int n_key2 = CountTokens(inp,' ');
      if (VERBOSE) fprintf(stderr, "* Found %ld keys #1 and %ld keys #2.\n", n_key1, n_key2);
      char **KEY2STR;
      ALLOCATE1D(KEY2STR,n_key2,char *);
      for (long int k2=0; k2<n_key2; k2++) KEY2STR[k2] = StrCopy(GetNextToken(&inp," \n"));

      // extract entries
      Progress PRG("Extracting table entries...",n_key1);
      for (long int k1=0; k1<n_key1; k1++) {
        char *inp = buffer->Next();
        char *key1 = GetNextToken(&inp," \t");
        Vector V = Vec(inp);
	long int nn = (long int)V[0];
	if (nn!=n_key2) { fprintf(stderr, "Line %ld: vector length should be %ld!\n", k1+2, n_key2); exit(1); }
	long int start = UPPER==false ? 1 : k1+2;
        for (long int k=start; k<=nn; k++) 
          if ((EQUAL&&(V[k]==CUTOFF))||(GREATER&&(V[k]>CUTOFF))||((GREATER==false)&&(V[k]<CUTOFF))) {
            printf("%s\t%s\t", key1, KEY2STR[k-1]);
            printf(fmt, V[k]);
	    printf("\n");
	  }
	FREE1D(V);
	PRG.Check();
      }
      PRG.Done();
    
      // cleanup
      FREE2D(KEY2STR,n_key2);
      delete buffer;
    }
  }

  //--------------------------------------------------------------------------------------//
  // OPTION -x0: extract table entries (tab-delimited format)                             //
  //--------------------------------------------------------------------------------------//

  else if (strcmp(method,"-x0")==0) {
    if (argc<3) { cmdLine->Usage(argv[0],"-x0 [OPTIONS] ENTRY_FILE"); exit(1); }
    MESSAGES(VERBOSE); 

    // initialize
    char *INPUT = argv[argc-1];
    FileBufferText buffer(INPUT,BUFFER_SIZE);
    long int n_lines = buffer.CountLines();
    if (n_lines>1) {   
      // read header
      long int n_key1 = n_lines - 1;
      char *inp = buffer.Next();
      GetNextToken(&inp,'\t');
      long int n_key2 = CountTokens(inp,'\t');
      if (VERBOSE) fprintf(stderr, "* Found %ld keys #1 and %ld keys #2.\n", n_key1, n_key2);
      char **KEY2STR;
      ALLOCATE1D(KEY2STR,n_key2,char *);
      for (long int k2=0; k2<n_key2; k2++) KEY2STR[k2] = StrCopy(ReadToken(&inp,'\t'));

      // extract entries
      Progress PRG("Extracting table entries...",n_key1);
      for (long int k1=0; k1<n_key1; k1++) {
        char *inp = buffer.Next();
        char *key1 = ReadToken(&inp,'\t');
	bool err = inp[0]==0;
        for (long int k=0; k<n_key2; k++) {
          char *val = ReadToken(&inp,'\t');
          err = err&&(inp[0]==0);
	  if (err==true) { fprintf(stderr, "Line %ld: insufficient number of tokens!\n", k1+2); exit(1); }
          printf("%s\t%s\t%s\n", key1, KEY2STR[k], val);
	}
	PRG.Check();
      }
      PRG.Done();
   
      // cleanup
      FREE2D(KEY2STR,n_key2);
    }

  }


  //--------------------------------------------------------------------------------------//
  // OPTION -r: convert from tab-delimited to standard format                             //
  //--------------------------------------------------------------------------------------//

  else if (strcmp(method,"-r")==0) {
    if (argc<3) { cmdLine->Usage(argv[0],"-r [OPTIONS] ENTRY_FILE"); exit(1); }
    MESSAGES(VERBOSE); 

    // initialize
    char *INPUT = argv[argc-1];
    FileBufferText buffer(INPUT,BUFFER_SIZE);
    long int n_lines = buffer.CountLines();
    if (n_lines>1) {   
      // read header
      long int n_key1 = n_lines - 1;
      char *inp = buffer.Next();
      GetNextToken(&inp,'\t');
      long int n_key2 = CountTokens(inp,'\t');
      if (VERBOSE) fprintf(stderr, "* Found %ld keys #1 and %ld keys #2.\n", n_key1, n_key2);
      printf("\t"); for (long int k2=0; k2<n_key2; k2++) printf("%s ",ReadToken(&inp,'\t')); printf("\n");

      // extract entries
      Progress PRG("Extracting table entries...",n_key1);
      for (long int k1=0; k1<n_key1; k1++) {
        char *inp = buffer.Next();
        char *key1 = ReadToken(&inp,'\t');
        printf("%s\t", key1);
        bool err = inp[0]==0;
        for (long int k=0; k<n_key2; k++) {
          char *val = ReadToken(&inp,'\t');
	  UpperCase(val);
          err = err&&(inp[0]==0);
	  if (err==true) { fprintf(stderr, "Line %ld: insufficient number of tokens!\n", k1+2); exit(1); }
          if ((val[0]==0)||(strcmp(val,"NAN")==0)) printf("NaN");
	  else printf(fmt, atof(val));
	  printf(" ");
	}
	printf("\n");
	PRG.Check();
      }
      PRG.Done();
    }
  }


  else {
    fprintf(stderr, "Error: '%s' is an invalid method!\n", method);
    exit(1);
  }

  

  // clean up
  delete cmdLine;

  
  return 0;
}



