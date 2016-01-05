//
// Copyright (c) 2014 Aristotelis Tsirigos
// All rights reserved. This program and the accompanying materials are made available under the terms of the GNU General Public License v2.0
// which accompanies this distribution, and is available at http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <map>
#include <string>
#include <list>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <thread>
#include "core.h"
#include "gtools-intervals.h"


using namespace std;


//---------------------------------------------------------------------------------//
// Global variables                                                                //
//---------------------------------------------------------------------------------//

const string PROGRAM = "gtools-threaded";
const long int BUFFER_SIZE = 10000;
gsl_rng *RANDOM_GENERATOR;







//---------------------------------------------------------------------------------//
// Command-line options                                                            //
//---------------------------------------------------------------------------------//

bool VERBOSE;
bool HELP;

char *OUTPUT_FILE;
char *RSCRIPT_INPUT_FILE_NAME;
bool REUSE;
char *OUT_PREFIX;
bool IGNORE_STRAND;
char *SHIFT;
char *LEGEND;
char *COLORS;
char *TITLE, *XLABEL, *YLABEL;
char *IMAGE_TYPE;
char *IMAGE_SIZE;
int IMAGE_RESOLUTION;
bool NORMALIZE_BY_REF_REGIONS;
bool NORMALIZE_BY_SIGNAL_REGIONS;
bool NORMALIZE_BY_BIN_SIZE;
bool NORMALIZE_RPKM;
char *PROFILE_TYPE;
int NBINS, NBINS_COMBINE;
double BIN_SIZE;
int NTHREADS;

char *GENOME_REG_FILE;
long int WIN_SIZE, WIN_DIST;
double PVALUE_CUTOFF;
float FDR;
float FOLD_CUTOFF;
long int FDR_BINS;
float OUTLIER_PROB;
char *SCALING;
char *NORMALIZATION;
double PSEUDOCOUNT;
char *LABELS;

bool SKIP_REF_GAPS;
bool NORM_REF_LEN;
double MAX_LABEL_VALUE_DOUBLE;
long int MAX_LABEL_VALUE_LINT;
char *VALUE_FORMAT; 
char *OVERLAP_OP;



//-------InitCmdLine-----------
//
CmdLineWithOperations *InitCmdLine(int argc, char *argv[], int *next_arg)
{
  // initialize
  CmdLineWithOperations *cmd_line = new CmdLineWithOperations();
  cmd_line->SetProgramName(PROGRAM,VERSION);

  cmd_line->AddOperation("matrix", "[OPTIONS] COMMA-SEPARATED-SIGNAL-REG-FILES REFERENCE-REGION-FILE", \
  "Computes overlaps with reference region set for multiple test region files; outputs values in matrix format.", \
  "* Input formats: REG, GFF, BED, SAM\n\
  * Operands: region, region-set\n\
  * Signal region requirements: single-interval\n\
  * Reference region requirements: chromosome/strand-compatible, sorted, non-overlapping\n" \
  );

  if (argc<2) {
    cmd_line->OperationSummary("OPERATION [OPTIONS] INPUT-FILES","Implements common genomics applications.");
    exit(1);
  }

  // current operation
  string op = argv[1];
  cmd_line->SetCurrentOperation(op);

  // common options
  cmd_line->AddOption("--help", &HELP, false, "help");
  cmd_line->AddOption("-h", &HELP, false, "help");
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");

  // Main options
  int n_args;
  bool image_output;
  if (op=="matrix") {
    n_args = 2;
    image_output = false;
    cmd_line->AddOption("-o", &OUTPUT_FILE, "", "output file name (default = standard output)");
    cmd_line->AddOption("-p", &NTHREADS, 1, "number of processes (threads)");
    cmd_line->AddOption("-i", &IGNORE_STRAND, false, "ignore signal region strand while finding overlaps");
    cmd_line->AddOption("-nbins", &NBINS, 10, "number of bins");
//    cmd_line->AddOption("--norm-by-total-reads", &NORMALIZE_BY_SIGNAL_REGIONS, false, "normalize by the total number of reads");
//    cmd_line->AddOption("--norm-by-bin-size", &NORMALIZE_BY_BIN_SIZE, false, "normalize by the bin size");
    cmd_line->AddOption("--overlap-op", &OVERLAP_OP, "hits", "determines how overlaps are counted = { hits, value }");
    cmd_line->AddOption("-rpkm", &NORMALIZE_RPKM, false, "converts read counts to RPKMs");
    cmd_line->AddOption("-profile", &PROFILE_TYPE, "none", "profile type {none, sum, mean}");
    cmd_line->AddOption("-format", &VALUE_FORMAT, "%.6e", "matrix value format string");
  }
  else {
    cerr << "Unknown operation '" << op << "'!\n";
    delete cmd_line;
    exit(1);
  }

  // process command line
  *next_arg = cmd_line->Read(argv+1,argc-1) + 1;
  if (HELP||(argc-*next_arg<n_args)) {
    cmd_line->OperationUsage();
    delete cmd_line;
    exit(1);
  }

  if (image_output&&(strcmp(IMAGE_TYPE,"pdf")!=0)&&(strcmp(IMAGE_TYPE,"tif")!=0)) {
    fprintf(stderr, "Error: unsupported image format '%s'!\n", IMAGE_TYPE);
    exit(1);
  }

  return cmd_line;
}








//-------Tokenize-----------
//
char **Tokenize(char *original_inp, char delim, int *n_tokens)
{
  if ((original_inp==NULL)||(strlen(original_inp)==0)) {
    *n_tokens = 0;
    return NULL;
  }
  char *inp = StrCopy(original_inp);
  *n_tokens = CountTokens(inp,delim);
  char **token = new char*[*n_tokens];
  char *p = inp;
  for (int k=0; k<*n_tokens; k++) token[k] = StrCopy(GetNextToken(&p,delim));
  free(inp);
  return token;
}





//-------PrintLogFile-----------
//
void PrintLogFile(string log_file_name)
{
  fprintf(stderr, "\n***********************************************\n");
  fprintf(stderr, "\n  R   S C R I P T   L O G   F I L E            \n");
  fprintf(stderr, "\n***********************************************\n");
  FileBufferText buffer(log_file_name.c_str());
  for (char *inp=buffer.Next(); inp!=NULL; inp=buffer.Next()) fprintf(stderr, "%s\n", inp);
}






//-----ThreadedOverlapMatrix----------
//
double ***ThreadedOverlapMatrix(
  int n_threads,
  char **signal_reg_file, 
  int n_signal_files, 
  char *ref_reg_file,
  long int n_bins,
  bool ignore_strand,
  char *overlap_op,
  bool rpkm)
{
  // initialize threads and data structures
  if (n_threads<1) n_threads = 1;
  if (n_threads>n_signal_files) n_threads = n_signal_files;
  int n = int(ceil((double)n_signal_files/n_threads));  // number of signal files to be processed per thread
  n_threads = int(ceil((double)n_signal_files/n));  // number of threads to be used needs to be adjusted
  double ****B = new double***[n_threads];  // pointers to the results of each thread (3-dim arrays)
  thread **overlap_thread = new thread*[n_threads];  // array of threads

  // process signal files in threads
  for (int t=0,s=0; t<n_threads; t++,s+=n) {
    int nn = t<n_threads-1?n:n_signal_files-s;
    overlap_thread[t] = new thread(CreateOverlapMatrixThread,&signal_reg_file[s],nn,ref_reg_file,n_bins,ignore_strand,overlap_op,rpkm,VERBOSE,&B[t]);
  }
  for (int t=0; t<n_threads; t++) overlap_thread[t]->join();

  // copy result pointers to final bins array
  double ***bins = new double**[n_signal_files];
  for (int t=0,s=0; t<n_threads; t++) for (int i=0; (i<n)&&(s<n_signal_files); i++,s++) bins[s] = B[t][i]; 
  delete B;
  if (VERBOSE) cerr << "All processes successfully completed!\n";
  
  // cleanup
  delete overlap_thread;
    
  return bins;
}





// ************************************************************************************
//
//   M  A  I  N
//
// ************************************************************************************
int main(int argc, char* argv[])
{
  // process command-line arguments
  int next_arg;
  CmdLineWithOperations *cmd_line = InitCmdLine(argc,argv,&next_arg);
  _MESSAGES_ = VERBOSE;

  // initilize random generator
  RANDOM_GENERATOR = InitRandomGenerator(getpid()+time(NULL));


  // ************************************************************************************
  //
  //     M  A  T  R  I  X
  //
  // ************************************************************************************

  if (cmd_line->current_cmd_operation=="matrix") {
    // input region files
    char *SIGNAL_REG_FILES = argv[next_arg];
    char *REF_REG_FILES = argv[next_arg+1];

    // process input parameters
    int n_signal_files;
    char **signal_reg_file = Tokenize(SIGNAL_REG_FILES,',',&n_signal_files);
    char *ref_reg_file = REF_REG_FILES;

    // create matrix
    double ***bins = ThreadedOverlapMatrix(NTHREADS,signal_reg_file,n_signal_files,ref_reg_file,NBINS,IGNORE_STRAND,OVERLAP_OP,NORMALIZE_RPKM);

    // post-process and print
    if (VERBOSE) fprintf(stderr, "Post-processing...\n"); 
    FILE *Fout = strlen(OUTPUT_FILE)==0?stdout:fopen(OUTPUT_FILE,"w");
    GenomicRegionSet RefRegSet(ref_reg_file,BUFFER_SIZE,false,true,true);
    long int n_ref = RefRegSet.n_regions;
    if (strcmp(PROFILE_TYPE,"none")==0) {
      // print header
      for (int s=0; s<n_signal_files; s++) for (int q=0; q<NBINS; q++) fprintf(Fout, "\t%s:bin=%d", signal_reg_file[s], q+1); fprintf(Fout, "\n");
      // print data
      for (int r=0; r<n_ref; r++) {
        fprintf(Fout, "%s\t", RefRegSet.R[r]->LABEL);
        for (int s=0; s<n_signal_files; s++) {
          for (int q=0; q<NBINS; q++) { fprintf(Fout, VALUE_FORMAT, bins[s][r][q]); fprintf(Fout, "%s", q!=NBINS-1?"\t":""); }
          fprintf(Fout, "%c", s!=n_signal_files-1?'\t':'\n');
        }
      }
    }
    else if ((strcmp(PROFILE_TYPE,"sum")==0)||(strcmp(PROFILE_TYPE,"mean")==0)) {
      for (int s=0; s<n_signal_files; s++) {
        fprintf(Fout, "%s\t", signal_reg_file[s]);
        for (int q=0; q<NBINS; q++) {
          double v = 0;
          for (int r=0; r<n_ref; r++) v += bins[s][r][q];
          if (strcmp(PROFILE_TYPE,"mean")==0) v /= n_ref;
          fprintf(Fout, "%.6e%c", v, q!=NBINS-1?'\t':'\n');
        }
      }
    }
    else { fprintf(stderr, "Error: unknown profile type!\n"); exit(1); }

    // cleanup
    for (int s=0; s<n_signal_files; s++) delete [] bins[s];
    delete bins;
    delete [] signal_reg_file;
    if (strlen(OUTPUT_FILE)>0) fclose(Fout);
  }

  // clean up
  delete cmd_line;
  delete RANDOM_GENERATOR;

  if (VERBOSE) fprintf(stderr, "Done.\n"); 

  return 0;
}








