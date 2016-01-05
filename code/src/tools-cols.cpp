#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <string>
#include "core.h"
using namespace std;


// options
bool OFFSET;
char INPUT_DELIMITER;
bool INPUT_TAB;
bool APPEND;
char OUTPUT_DELIMITER;
bool OUTPUT_TAB;
char *NEW_DELIMITERS;
long int BUFFER_SIZE;
bool VERBOSE;




// ----- InitCmdLine ---------
//
CmdLine *InitCmdLine()
{
  CmdLine *cmd_line = new CmdLine(); 

  // main switches
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose");
  cmd_line->AddOption("-1", &OFFSET, false, "offset correction (subtract 1)");
  cmd_line->AddOption("-d", &INPUT_DELIMITER, ' ', "input delimiter");
  cmd_line->AddOption("-t", &INPUT_TAB, false, "use TAB as input delimiter");
  cmd_line->AddOption("-a", &APPEND, false, "append selected columns");
  cmd_line->AddOption("-D", &OUTPUT_DELIMITER, ' ', "output delimiter if append is applied");
  cmd_line->AddOption("-T", &OUTPUT_TAB, false, "use TAB as input delimiter");
  cmd_line->AddOption("-r", &NEW_DELIMITERS, "", "replace delimiter");
  cmd_line->AddOption("-b", &BUFFER_SIZE, 100000000, "input buffer size");
  
  // set default values
  cmd_line->Init();
  
  return cmd_line;
}





//-----ReadSequence---------
//
Sequence ReadSequence(char **s, long int n)
{
  Sequence Q = Seq(n);
  for (long int k=0; k<n; k++) { Q[k+1] = atol(s[k]); if (OFFSET) Q[k+1]--; }
  return Q;
}





// --------------------------------------------------------------------------------//
// MAIN	                                                                           //
// --------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  CmdLine *cmdLine = InitCmdLine();
  if (argc<2) { cmdLine->Usage("cols [OPTIONS] COLUMN_SEQUENCE <INPUT>"); exit(1); }
  int next_arg = cmdLine->Read(argv,argc);
  MESSAGES(VERBOSE);

  // read inputs
  if (strlen(NEW_DELIMITERS)==0) {

    if (INPUT_TAB==true) INPUT_DELIMITER = '\t';
    if (OUTPUT_TAB==true) OUTPUT_DELIMITER = '\t';
    if (VERBOSE) cerr << "input delimiter = " << '"' << INPUT_DELIMITER << '"' << '\n';
    if (VERBOSE) cerr << "output delimiter = " << '"' << OUTPUT_DELIMITER << '"' << '\n';
    Sequence cols = ReadSequence(&argv[next_arg],argc-next_arg);

    Progress PRG("Processing...",1);
    FileBufferText buffer((char*)NULL,BUFFER_SIZE);
    char *inp = buffer.Next();
    for (long int n=0; inp!=NULL; n++,inp=buffer.Next()) {
      if (APPEND) cout << inp << '\t';
      int n_cols = CountTokens(inp,INPUT_DELIMITER);
      string *S = new string[n_cols];
      for (int c=0; c<n_cols; c++) S[c] = GetNextToken(&inp,INPUT_DELIMITER);
      for (long int k=1; k<=cols[0]; k++) {
        if (cols[k]>=n_cols) { fprintf(stderr, "Line %ld: number of columns (%d) too small!\n", n+1, n_cols); exit(1); }
        cout << S[cols[k]];
        if (k<cols[0]) cout << (APPEND?OUTPUT_DELIMITER:INPUT_DELIMITER);
      }
      cout << '\n';
      delete [] S;
      PRG.Check();
    }
    PRG.Done();
    FREE1D(cols);
  }



  else {
    //
    // SOS: this is beta-version 
    //
    FileBufferText buffer((char*)NULL,BUFFER_SIZE);
    Progress PRG("Processing...",1);
    for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) {
      for (size_t c=0; c<strlen(NEW_DELIMITERS); c++) {
        char *s = GetNextToken(&inp,INPUT_DELIMITER);
	cout << s << NEW_DELIMITERS[c];
      }
      cout << inp << '\n';
      PRG.Check();
    }
    PRG.Done();

  }

  // clean up
  delete cmdLine;

  
  return 0;
}



