#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include "core.h"


#define BUFFER_SIZE 100000000


// options
int PERIOD;
bool MERGE;
bool OFFSET;
char *SEPARATOR;
bool EXCLUDE;
bool EMPTY;
bool VERBOSE;
long int N_SAMPLES;
char *ROW_FILE;
unsigned long int RND_SEED;
bool HEADER;
char *PREFIX = NULL;
long int PADDING;



// ----- InitCmdLine ---------
//
CmdLine *InitCmdLine(char *method=NULL)
{
  CmdLine *cmd_line = new CmdLine(); 

  cmd_line->AddOption("-v", &VERBOSE, false, "verbose");

  if (strcmp(method,"-permute")==0) {

  }
  if (strcmp(method,"-number")==0) {
    cmd_line->AddOption("-head", &HEADER, false, "input contains header");
    cmd_line->AddOption("-pref", &PREFIX, "", "prefix");
    cmd_line->AddOption("-pad", &PADDING, 10, "max number of leading zeroes");
  }
  else if (strcmp(method,"-resample")==0) {
    cmd_line->AddOption("-r", &N_SAMPLES, 0, "number of samples (default = number of input lines)");
    cmd_line->AddOption("--seed", &RND_SEED, 0, "seed for random number generator");
  }
  else {
    cmd_line->AddOption("-p", &PERIOD, 0, "period");
    cmd_line->AddOption("-m", &MERGE, false, "merge rows");
    cmd_line->AddOption("-t", &SEPARATOR, " ", "separator");
    cmd_line->AddOption("-1", &OFFSET, false, "offset correction (subtract 1)");
    cmd_line->AddOption("-x", &EXCLUDE, false, "exlude listed rows");
    cmd_line->AddOption("-e", &EMPTY, false, "print non-selected rows as empty lines");
    cmd_line->AddOption("-f", &ROW_FILE, "", "file containing rows to be selected");
  }
  
  // set default values
  cmd_line->Init();
  
  return cmd_line;
}



//-----Usage----------
//
void Usage()
{
  int n_methods = 4;
  char METHODS[][20] = { "-permute", "-resample", "-number", "" };
  for (int i=0; i<n_methods; i++) {
    CmdLine *cmd = InitCmdLine(METHODS[i]);
    char s[1000];
    sprintf(s, "rows %s [OPTIONS] <ROWS>", METHODS[i]);
    cmd->Usage(s);
    delete cmd;
  }
}





// ----- GetRows ---------
//
Sequence GetRows(char **s, long int n)
{
  Sequence Q = Seq(n);
  for (long int k=0; k<n; k++) { Q[k+1] = atol(s[k]); if (OFFSET) Q[k+1]--; }
  SeqUnique(Q);
  return Q;
}



//-----LoadRows-----
//
Sequence LoadRows(char *file)
{
  FileBufferText buffer(file);
  long int n = buffer.CountLines();
  Sequence Q = Seq(n);
  for (long int k=1; k<=n; k++) { Q[k] = atol(buffer.Next()); if (OFFSET) Q[k]--; }
  SeqUnique(Q);
  return Q;
}






// --------------------------------------------------------------------------------//
// MAIN	                                                                           //
// --------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  if (argc<2) { Usage(); exit(1); }
  char *method = argv[1];
  CmdLine *cmdLine = InitCmdLine(method);



  //--------------------------------------------------------------------------------------//
  // OPTION -permute: permute rows                                                        //
  //--------------------------------------------------------------------------------------//
  if (strcmp(method,"-permute")==0) { 
    // read options
    cmdLine->Read(argv+1,argc-1);
    MESSAGES(VERBOSE);

    // initialize random generator
    unsigned long int seed = RND_SEED+getpid()+time(NULL);
    gsl_rng *RANDOM_GENERATOR = InitRandomGenerator(seed);

    // read input
    long int n_lines;
    char **L;
    FILE *fp = LoadStdIn(&n_lines,BUFFER_SIZE);
    if (n_lines==0) return 0;
    ALLOCATE1D(L,n_lines,char *);
    FileBufferText *buffer = new FileBufferText(fp,BUFFER_SIZE);
    Progress PRG("Reading input rows...",n_lines);
    for (long int n=0; n<n_lines; n++) {
      L[n] = StrCopy(buffer->Next());
      PRG.Check();
    }
    PRG.Done();
    
    // permute
    long int *seq;
    ALLOCATE1D(seq,n_lines,long int);
    for (long int k=0; k<n_lines; k++) seq[k] = k;
    gsl_ran_shuffle(RANDOM_GENERATOR,seq,n_lines,sizeof(long int));

    // print output
    Progress PRG1("Printing output rows...",n_lines);
    for (long int n=0; n<n_lines; n++) {
      printf("%s\n", L[seq[n]]);
      PRG1.Check();
    }
    PRG1.Done();

    // cleanup
    FREE1D(seq);
    delete RANDOM_GENERATOR;
    delete buffer;
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -resample: resample rows                                                      //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-resample")==0) { 
    // read options
    cmdLine->Read(argv+1,argc-1);
    MESSAGES(VERBOSE);

    // initialize random generator
    unsigned long int seed = (RND_SEED==0)?(getpid()+time(NULL)):RND_SEED;
    gsl_rng *RANDOM_GENERATOR = InitRandomGenerator(seed);

    // read input
    long int n_lines;
    char **L;
    FILE *fp = LoadStdIn(&n_lines,BUFFER_SIZE);
    if (n_lines==0) return 0;
    ALLOCATE1D(L,n_lines,char *);
    FileBufferText *buffer = new FileBufferText(fp,BUFFER_SIZE);
    Progress PRG("Reading input rows...",n_lines);
    for (long int n=0; n<n_lines; n++) {
      L[n] = StrCopy(buffer->Next());
      PRG.Check();
    }
    PRG.Done();
    
    // print output
    if (N_SAMPLES<=0) N_SAMPLES = n_lines;
    Progress PRG1("Printing output rows...",N_SAMPLES);
    for (long int n=0; n<N_SAMPLES; n++) {
      long int k = gsl_rng_uniform_int(RANDOM_GENERATOR,n_lines);
      printf("%s\n", L[k]);
      PRG1.Check();
    }
    PRG1.Done();

    // cleanup
    delete RANDOM_GENERATOR;
    delete buffer;
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -number: do row numbering                                                     //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-number")==0) { 
    // read options
    cmdLine->Read(argv+1,argc-1);
    MESSAGES(VERBOSE);

    // read input
    FileBufferText buffer((char*)NULL,BUFFER_SIZE);
    Progress PRG("Reading input rows...",1);
    long int id = HEADER?0:1;
    for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) {
      if (id==0) printf("\t");
      else printf("%s%012ld\t", PREFIX, id);
      printf("%s\n", inp);
      id++;
      PRG.Check();
    }
    PRG.Done();    
  }


  //--------------------------------------------------------------------------------------//
  // OTHER OPTIONS                                                                        //
  //--------------------------------------------------------------------------------------//
  else {
    // read options
    int next_arg = cmdLine->Read(argv,argc);
    MESSAGES(VERBOSE);
  
    // allocate buffer
    char *BUFFER = (char *) malloc(BUFFER_SIZE*sizeof(char));
    if (BUFFER==NULL) { fprintf(stderr, "Out of memory!\n"); exit(1); }

    Sequence rows = strlen(ROW_FILE)==0 ? GetRows(&argv[next_arg],argc-next_arg) : LoadRows(ROW_FILE);

    // option 1: use period
    if (PERIOD>0) {
      int choose = MERGE ? -1 : atoi(argv[next_arg]);
      //printf("(PERIOD,NEXTARG,CHOOSE) = (%i,%i,%i)\n", PERIOD, next_arg, choose);
      for (int n=0; ; n++) {
        if (fgets(BUFFER,BUFFER_SIZE,stdin)==NULL) break;
        if (BUFFER[strlen(BUFFER)-1]!='\n') { 
          fprintf(stderr, "Error: line was only partially read (max chars=%i)!\n", BUFFER_SIZE); exit(1);
        }
        if (BUFFER[strlen(BUFFER)-1]=='\n') BUFFER[strlen(BUFFER)-1] = 0;
        if (MERGE) printf("%s%s", BUFFER, n%PERIOD==PERIOD-1?"\n":SEPARATOR); 
        else if (n%PERIOD==choose) printf("%s\n", BUFFER); 
      }
    }

    // option 2: choose rows
    else {
      if (EXCLUDE==true) {
        for (long int n=0,k=1; k<=rows[0]; n++) {
          fgets(BUFFER,BUFFER_SIZE,stdin);
          if (feof(stdin)) break;
          if (BUFFER[strlen(BUFFER)-1]!='\n') { fprintf(stderr, "Error: line was only partially read (max chars=%i)!\n", BUFFER_SIZE); exit(1); }
          if (n<rows[k]) printf("%s", BUFFER);
	  else { ++k; if (EMPTY) printf("\n"); }
        }
        while (true) {
          fgets(BUFFER,BUFFER_SIZE,stdin);
          if (feof(stdin)) break;
          if (BUFFER[strlen(BUFFER)-1]!='\n') { fprintf(stderr, "Error: line was only partially read (max chars=%i)!\n", BUFFER_SIZE); exit(1); }
          printf("%s", BUFFER);
        }
      }
      else {
        for (long int n=0,k=1; k<=rows[0]; n++) {
          fgets(BUFFER,BUFFER_SIZE,stdin);
          if (feof(stdin)) break;
          if (BUFFER[strlen(BUFFER)-1]!='\n') { fprintf(stderr, "Error: line was only partially read (max chars=%i)!\n", BUFFER_SIZE); exit(1); }
          if (n==rows[k]) { printf("%s", BUFFER); k++; }
	  else if (EMPTY) printf("\n"); 
        }
        if (EMPTY) {
          while (true) {
            fgets(BUFFER,BUFFER_SIZE,stdin);
            if (feof(stdin)) break;
            if (BUFFER[strlen(BUFFER)-1]!='\n') { fprintf(stderr, "Error: line was only partially read (max chars=%i)!\n", BUFFER_SIZE); exit(1); }
            printf("\n");
          }
        }
      }
    }

    // clean up
    FREE1D(rows);
    FREE1D(BUFFER);
  }


  // clean up
  delete cmdLine;

  
  return 0;
}



