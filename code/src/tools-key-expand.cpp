#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_erf.h>
#include "core.h"



// options
bool VERBOSE;
bool HELP;



//-----InitCmdLine---------
//
CmdLine *InitCmdLine()
{
  CmdLine *cmd_line = new CmdLine(); 

  // options
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-h", &HELP, false, "help");
  
  return cmd_line;
}







// --------------------------------------------------------------------------------//
// MAIN	                                                                           //
// --------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  CmdLine *cmdLine = InitCmdLine();
  int next_arg = cmdLine->Read(argv,argc);
  if (HELP) { cmdLine->Usage("key-expand <STDIN>|FILE"); exit(1); }
  MESSAGES(VERBOSE); 
  char *INPUT = next_arg==argc ? NULL : argv[next_arg];

  // process
  FileBufferText buffer(INPUT);
  for (long int r=0; buffer.Next()!=NULL; r++) {
    char *inp = buffer.Get();
    char *key = GetNextToken(&inp," \t");
    while (inp[0]!=0) {
      char *val = GetNextToken(&inp," \t");
      printf("%s\t%s\n", key, val);
    }
  }

  // cleanup
  delete cmdLine;

  return 0;
}



