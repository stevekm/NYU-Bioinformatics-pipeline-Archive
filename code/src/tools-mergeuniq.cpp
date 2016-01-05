#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "core.h"




// options
bool VERBOSE;
char *SEPARATOR;
int DECIMALS;
char *VEC_FUNC;
char *FMT_STR;
bool FMT_INT;
bool AVERAGE;
const long int BUFFER_SIZE = 500000000;





//-------InitCmdLine-----------
//
CmdLine *InitCmdLine(char *method)
{
  CmdLine *cmd_line = new CmdLine(); 

  // common options
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");

  
  // Main options
  if (strcmp(method,"-add")==0) {

  }
  else if (strcmp(method,"-vec")==0) {
    cmd_line->AddOption("-func", &VEC_FUNC, "add", "function = {add,sub,prod,div}");
    cmd_line->AddOption("-f", &FMT_STR, "", "format string");
    cmd_line->AddOption("-i", &FMT_INT, false, "integer format");
    cmd_line->AddOption("-n", &DECIMALS, 2, "number of decimal points");
    cmd_line->AddOption("-a", &AVERAGE, false, "compute the average vector");
  }
  else if (strcmp(method,"-merge")==0) {
    cmd_line->AddOption("-t", &SEPARATOR, " ", "separator");
  }
  else {
    fprintf(stderr, "Unknown method '%s'!\n", method);
    exit(1);
  }
  
  return cmd_line;
}



//-----Usage----------
//
void Usage(char **argv)
{
  const int n_methods = 3;
  char methods[n_methods][20] = { "-add", "-vec", "-merge" };
  for (int m=0; m<n_methods; m++) {
    CmdLine *cmdLine = InitCmdLine(methods[m]);
    cmdLine->Read(argv,1);
    char msg[200];
    sprintf(msg, "%s %s [OPTIONS] FILE|<STDIN>", argv[0], methods[m]);
    cmdLine->Usage(msg);
    fprintf(stderr, "\n");
    delete cmdLine;
  }
}






// --------------------------------------------------------------------------------//
// MAIN	                                                                           //
// --------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // read command-line options 
  if (argc<2) { Usage(argv); exit(1); }
  char *METHOD = argv[1];
  CmdLine *cmdLine = InitCmdLine(METHOD);
  int next_arg = cmdLine->Read(argv+1,argc-1) + 1;
  MESSAGES(VERBOSE); 

  // open input file (or stdin)
  MESSAGES(VERBOSE); 
  char *INPUT = next_arg==argc ? NULL : argv[next_arg];
  FileBufferText buffer(INPUT,BUFFER_SIZE);
  long int n_lines = INPUT==NULL ? 1 : buffer.CountLines();
  
  

  //--------------------------------------------------------------------------------------//
  // OPTION -add:                                                                         //
  //--------------------------------------------------------------------------------------//

  if (strcmp(METHOD,"-add")==0) {
    unsigned long int nhits = 0;
    unsigned long int nseqs = 0;
    char *prev = NULL;
    while (buffer.Next()) {
      char *inp = buffer.Get();
      unsigned long int a = atol(GetNextToken(&inp,'\t'));
      unsigned long int b = atol(GetNextToken(&inp,'\t'));
      char *motif = GetNextToken(&inp,'\n');
      if (prev==NULL) prev = StrCopy(motif);
      if (strcmp(motif,prev)==0) { nhits += a; nseqs += b; }
      else { 
        printf("%lu\t%lu\t%s\n", nhits, nseqs, prev);
        nhits = a;
        nseqs = b;
        if (prev) free(prev);
        prev = StrCopy(motif);
      }
    }  

    printf("%lu\t%lu\t%s\n", nhits, nseqs, prev);

    // clean up
    free(prev);
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -vec:                                                                      //
  //--------------------------------------------------------------------------------------//

  else if (strcmp(METHOD,"-vec")==0) {    
    // setup output format
    char fmt[100];
    if (FMT_INT) strcpy(fmt,"%.0f");
    else if (strlen(FMT_STR)==0) sprintf(fmt, "%%.%if", DECIMALS); 
    else sprintf(fmt, "%s", FMT_STR);
    if (VERBOSE) fprintf(stderr, "* Using format '%s'\n", fmt);

    // process
    Progress PRG("Processing...",n_lines);
    char *key = NULL;
    long int cols = 0;
    Vector v = NULL;
    long int n_vec = 0;
    char *inp = buffer.Next();
    for (long int n=0; inp!=NULL; n++,inp=buffer.Next()) {
      // read/process input
      char *new_key = StrCopy(GetNextToken(&inp," \t"));
      Vector new_v = Vec(inp,' ');
      long int new_cols = (long int)new_v[0];
      
      // initialize first key
      if (key==NULL) { 
        key = new_key; 
        cols = new_cols; 
        v = new_v;
        n_vec = 1;
      }

      // if new key matches current key, add vectors
      else if (strcmp(key,new_key)==0) { 
        if (cols!=new_cols) { fprintf(stderr, "Error line %ld: column number mismatch!\n", n+1); exit(1); }
        if (strcmp(VEC_FUNC,"add")==0) for (long int k=1; k<=cols; k++) v[k] += new_v[k];
        else if (strcmp(VEC_FUNC,"sub")==0) for (long int k=1; k<=cols; k++) v[k] -= new_v[k];
        else if (strcmp(VEC_FUNC,"prod")==0) for (long int k=1; k<=cols; k++) v[k] *= new_v[k];
        else if (strcmp(VEC_FUNC,"div")==0) for (long int k=1; k<=cols; k++) v[k] /= new_v[k];
		else { fprintf(stderr, "Error: unknown vector function '%s'!\n", VEC_FUNC); exit(1); }
        n_vec++;
        FREE1D(new_key);
	    FREE1D(new_v);
      }

      // else print current key/vector, and initialize new key data
      else { 
        printf("%s\t", key); 
        if (AVERAGE) for (long int k=1; k<=cols; k++) { printf(fmt, v[k]/n_vec); printf(" "); }
        else for (long int k=1; k<=cols; k++) { printf(fmt, v[k]); printf(" "); }
        //printf("\t%ld", n_vec);
        printf("\n");
	
	FREE1D(key);
	FREE1D(v);
        key = new_key;
	cols = new_cols;
	v = new_v;
	n_vec = 1;
      }

      PRG.Check();
    }  
    PRG.Done();

    printf("%s\t", key); 
    if (AVERAGE) for (long int k=1; k<=cols; k++) { printf(fmt, v[k]/n_vec); printf(" "); }
    else for (long int k=1; k<=cols; k++) { printf(fmt, v[k]); printf(" "); }
    //printf("\t%ld", n_vec);
    printf("\n");

    // clean up
    FREE1D(key);
    FREE1D(v);
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -merge:                                                                       //
  //--------------------------------------------------------------------------------------//

  else if (strcmp(METHOD,"-merge")==0) {
    if (buffer.Next()!=NULL) {
      char *prev = NULL;
      for (long int n=0; ; n++) {
        char *inp = buffer.Get();
        char *key = GetNextToken(&inp," \t");
        if (prev==NULL) { prev = StrCopy(key); printf("%s\t%s", key, GetNextToken(&inp,'\n')); }
	else {
          if (strcmp(key,prev)==0) printf("%s%s", SEPARATOR, GetNextToken(&inp,'\n'));
          else { 
            printf("\n%s\t", key);
            printf("%s", GetNextToken(&inp,'\n'));
            if (prev) free(prev);
            prev = StrCopy(key);
          }
	}
	if (buffer.Next()==NULL) break;
      }  
      printf("\n");
      if (prev) free(prev);
    }

  }

  // clean up
  delete cmdLine;

  
  return 0;
}



