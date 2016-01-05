#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_fit.h>
#include <string>
#include <iostream>
#include <limits>
#include "core.h"
using namespace std;



//---------------------------------------------------------------------------------//
// COMMAND-LINE OPTIONS                                                            //
//---------------------------------------------------------------------------------//


// Generic options
bool VERBOSE;
char *FMT_STR;
int DECIMALS;
float CORR_VAL;
long int N_TESTS;
double DIFF_PSEUDO;
double DIFF_FOLD;



//-----Usage----------
//
void Usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: \n");
  fprintf(stderr, "  matrix OPERATION [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2\n"); 
  fprintf(stderr, "\n");
  fprintf(stderr, "OPERATIONS: \n");
  fprintf(stderr, "  * -corr       correlations of all row vector pairs\n");
  fprintf(stderr, "  * -corrp      correlations and p-values of all row vector pairs \n");
  fprintf(stderr, "  * -diff       count pairwise differences based on fold-change\n");
  fprintf(stderr, "  * -div        divide\n");
  fprintf(stderr, "  * -prod       product vectors of all row vector pairs\n");
  fprintf(stderr, "  * -reg        compute all pairwise linear regressions\n");
  fprintf(stderr, "\n");
}



//-----Usage----------
//
void Usage(CmdLine *cmdLine, char *method)
{
  if (strcmp(method,"-diff")==0) cmdLine->Usage("matrix2 -diff [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2 => MATRIX"); 
  else if (strcmp(method,"-prod")==0) cmdLine->Usage("matrix2 -prod [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2 => MATRIX"); 
  else if (strcmp(method,"-div")==0) cmdLine->Usage("matrix2 -div [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2 => MATRIX"); 
  else if (strcmp(method,"-corr")==0) cmdLine->Usage("matrix2 -corr [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2 => MATRIX"); 
  else if (strcmp(method,"-corrp")==0) cmdLine->Usage("matrix2 -corrp [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2 => LIST"); 
  else if (strcmp(method,"-reg")==0) cmdLine->Usage("matrix2 -reg [OPTIONS] MATRIX_FILE_1 MATRIX_FILE_2 => LIST"); 
  else fprintf(stderr, "Unknown method '%s'!\n", method);
  exit(1); 
}



//--------InitCmdLine-----------
//
CmdLine *InitCmdLine(char *method)
{
  CmdLine *cmd_line = new CmdLine(); 

  // Processing
  cmd_line->AddOption("-v", &VERBOSE, false, "verbose mode");
  cmd_line->AddOption("-f", &FMT_STR, "", "format string");
  cmd_line->AddOption("-n", &DECIMALS, 2, "number of decimal points");

  if (strcmp(method,"-diff")==0) {
    cmd_line->AddOption("-pseudo", &DIFF_PSEUDO, 0.0, "add pseudo-count in fold-change computation");
    cmd_line->AddOption("-fold", &DIFF_FOLD, 2.0, "two-way fold-change cutoff");
  }

  else if (strcmp(method,"-corrp")==0) {
    cmd_line->AddOption("-c", &CORR_VAL, 0.5, "correlation value for significance test");
    cmd_line->AddOption("-t", &N_TESTS, 1000, "number of resampling tests");
  }

  return cmd_line;
}












//-----VectorProduct------
//
float *VectorProduct(float *V1, float *V2, long int n)
{
  float *V;
  ALLOCATE1D(V,n,float);
  for (long int k=0; k<n; k++) 
    if ((V1[k]==V1[k])&&(V2[k]==V2[k])) V[k] = 100*V1[k]*V2[k];
    else V[k] = numeric_limits<float>::quiet_NaN();    // 0.0/0.0;
  return V;
}


//-----VectorDiv------
//
float *VectorDiv(float *V1, float *V2, long int n)
{
  float *V;
  ALLOCATE1D(V,n,float);
  for (long int k=0; k<n; k++) 
    if ((V1[k]==V1[k])&&(V2[k]==V2[k])) V[k] = V1[k]/V2[k];
    else V[k] = numeric_limits<float>::quiet_NaN();
  return V;
}


//------DistL1------
//
float DistL1(float *A, float *B, unsigned long int M)
{
  float D = 0;
  int C = 0;
  for (unsigned long int i=0; i<M; i++) 
    if ((A[i]==A[i])&&(B[i]==B[i])) {
      D += fabs(A[i]-B[i]);
      C++;
    }
  return D/C;
}



//-------regression--------
//
float *regression(float *X, float *Y, int N, double *_chisq, double *_c0, double *_c1)
{
  // prepare the data
  double c0, c1, cov00, cov01, cov11, chisq;
  double *x, *y, *w;
  ALLOCATE1D(x,N,double);
  ALLOCATE1D(y,N,double);
  ALLOCATE1D(w,N,double);
  for (int k=0; k<N; k++) {
    if ((X[k]==X[k])&&(Y[k]==Y[k])) {
      x[k] = X[k];
      y[k] = Y[k];
      w[k] = 1;      
    }
    else x[k] = y[k] = w[k] = 0;
  }
  
  // run regression
  float *e;
  ALLOCATE1D(e,N,float);
  gsl_fit_wlinear(x, 1, w, 1, y, 1, N, &c0, &c1, &cov00, &cov01, &cov11, &chisq);
  chisq /= N;
  for (int k=0; k<N; k++) 
    if ((X[k]==X[k])&&(Y[k]==Y[k])) e[k] = y[k] - (c0 + c1*x[k]);
    else e[k] = numeric_limits<float>::quiet_NaN();
    
  // clean up memory
  FREE1D(x);
  FREE1D(y);
  FREE1D(w);

  // return values
  *_chisq = chisq;
  *_c0 = c0;
  *_c1 = c1;
  return e;
}




//-------Diff--------
//
long int Diff(float *X, float *Y, long int n, double pseudo, double fold)
{
  long int c = 0;
  for (long int k=0; k<n; k++) c += ((X[k]+pseudo)/(Y[k]+pseudo)>fold)||((Y[k]+pseudo)/(X[k]+pseudo)>fold);
  return c;
}




//---------------------------------------------------------------------------------//
// MAIN	                                                                           //
//---------------------------------------------------------------------------------//
int main(int argc, char* argv[]) 
{
  // process command-line arguments
  if (argc<2) { Usage(); exit(1); }
  char *method = argv[1];
  CmdLine *cmdLine = InitCmdLine(method);
  int next_arg = cmdLine->Read(argv+1,argc-1) + 1;
  if (next_arg>=argc) { Usage(cmdLine,method); exit(1); }
  MESSAGES(VERBOSE);

  // initialize random generator
  gsl_rng *RANDOM_GENERATOR = InitRandomGenerator(time(NULL));

  // load matrices
  MatrixTemplate<float> M1(argv[next_arg],VERBOSE);
  MatrixTemplate<float> M2(argv[next_arg+1],VERBOSE);

  // setup output format
  char fmt[1000];
  if (strlen(FMT_STR)==0) sprintf(fmt, "%%.%if", DECIMALS); 
  else sprintf(fmt, "%s", FMT_STR);
  if (VERBOSE) fprintf(stderr, "* Using format '%s'\n", fmt);



  //--------------------------------------------------------------------------------------//
  // OPTION -prod: pairwise row products                                                  //
  //--------------------------------------------------------------------------------------//
  if (strcmp(method,"-prod")==0) { 
    if ((M1.n_cols!=M2.n_cols)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    M1.PrintColLabels();
    for (int r1=0; r1<M1.n_rows; r1++) {
      for (int r2=0; r2<M2.n_rows; r2++) {
        if (M1.has_row_labels) cout << M1.row_labels[r1] << '|' << M2.row_labels[r2] << '\t';
        //float *V = VectorProduct(M1.val[r1],M2.val[r2],M1.n_cols);
        float *V = VectorDiv(M1.val[r1],M2.val[r2],M1.n_cols);
        VectorPrint(V,M1.n_cols,fmt);
        printf("\n");
	FREE1D(V);
        PRG.Check();
      }
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -div: divide                                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-div")==0) { 
    if ((M1.n_rows!=M2.n_rows)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    if (M1.has_col_labels) { 
      cout << '\t'; 
      for (int c1=0; c1<M1.n_rows; c1++) cout << M1.col_labels[c1] << ' '; 
      cout << '\n'; 
    }
    for (int r1=0; r1<M1.n_rows; r1++) {
      if (M1.has_row_labels) cout << M1.row_labels[r1] << '\t';
      for (int c1=0; c1<M1.n_cols; c1++) {
        printf(fmt, M1.val[r1][c1]/M2.val[r1][0]);
	printf(" ");
        PRG.Check();
      }
      cout << '\n';
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -corrp: pairwise row correlations and p-values                                //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-corrp")==0) { 
    if ((M1.n_cols!=M2.n_cols)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    for (int r1=0; r1<M1.n_rows; r1++) {
      for (int r2=0; r2<M2.n_rows; r2++) {
        if (M1.has_row_labels) cout << M1.row_labels[r1] << '\t';
        if (M2.has_row_labels) cout << M2.row_labels[r2] << '\t';
        double std;
        double c = Corr_resample(RANDOM_GENERATOR,M1.val[r1],M2.val[r2],M1.n_cols,N_TESTS,&std);
        printf(fmt, c); printf("\t");
        printf(fmt, std); printf("\t");
        double pval = gsl_cdf_gaussian_Q(c-CORR_VAL,std);
        printf("%.2e", pval);
        printf("\n");
        PRG.Check();
      }
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -corr: pairwise row correlations                                              //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-corr")==0) { 
    if ((M1.n_cols!=M2.n_cols)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    if (M2.has_row_labels) { 
      if (M1.has_row_labels) cout << '\t'; 
      for (int c=0; c<M2.n_rows; c++) cout << M2.row_labels[c] << ' '; 
      cout << '\n'; 
    }
    for (int r1=0; r1<M1.n_rows; r1++) {
      if (M1.has_row_labels) cout << M1.row_labels[r1] << '\t';
      for (int r2=0; r2<M2.n_rows; r2++) {
        printf(fmt, Corr(M1.val[r1],M2.val[r2],M1.n_cols));
	printf(" ");
        PRG.Check();
      }
      printf("\n");
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -diff: count pairwise differences based on fold-change                        //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-diff")==0) { 
    if ((M1.n_cols!=M2.n_cols)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    if (M2.has_row_labels) { 
      if (M1.has_row_labels) cout << '\t'; 
      for (int c=0; c<M2.n_rows; c++) cout << M2.row_labels[c] << ' '; 
      cout << '\n'; 
    }
    for (int r1=0; r1<M1.n_rows; r1++) {
      if (M1.has_row_labels) cout << M1.row_labels[r1] << '\t';
      for (int r2=0; r2<M2.n_rows; r2++) {
        printf("%ld", Diff(M1.val[r1],M2.val[r2],M1.n_cols,DIFF_PSEUDO,DIFF_FOLD));
	printf(" ");
        PRG.Check();
      }
      printf("\n");
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -reg: pairwise row linear regressions                                         //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-reg")==0) { 
    if ((M1.n_cols!=M2.n_cols)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    for (int r1=0; r1<M1.n_rows; r1++) {
      for (int r2=0; r2<M2.n_rows; r2++) {
        // print labels
        if (M1.has_row_labels) cout << M1.row_labels[r1] << '\t';
        if (M2.has_row_labels) cout << M2.row_labels[r2] << '\t';
	
	// compute
        double chi_sq, c0, c1; 
        float *e = regression(M1.val[r1],M2.val[r2],M1.n_cols,&chi_sq,&c0,&c1);
        
	// print results
	printf(fmt, chi_sq); printf("\t");
	printf(fmt, c0); printf("\t");
        printf(fmt, c1); printf("\t");
        for (int k=0; k<M1.n_cols; k++) { printf(fmt, e[k]); printf(" "); }
        printf("\n");

	// cleanup
	FREE1D(e);
        PRG.Check();
      }
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // OPTION -dist: pairwise L1 distances                                                  //
  //--------------------------------------------------------------------------------------//
  else if (strcmp(method,"-dist")==0) { 
    if ((M1.n_cols!=M2.n_cols)||(M1.has_row_labels!=M2.has_row_labels)||(M1.has_col_labels!=M2.has_col_labels)) { fprintf(stderr, "Incompatible matrices!\n"); exit(1); }
    Progress PRG("Processing row pairs...",M1.n_rows*M2.n_rows);
    if (M2.has_row_labels) { 
      if (M1.has_row_labels) cout << '\t'; 
      for (int c=0; c<M2.n_rows; c++) cout << M2.row_labels[c] << ' '; 
      cout << '\n'; 
    }
    for (int r1=0; r1<M1.n_rows; r1++) {
      if (M1.has_row_labels) cout << M1.row_labels[r1] << '\t';
      for (int r2=0; r2<M2.n_rows; r2++) {
        printf(fmt, DistL1(M1.val[r1],M2.val[r2],M1.n_cols));
	printf(" ");
        PRG.Check();
      }
      cout << '\n';
    }
    PRG.Done();   
  }


  //--------------------------------------------------------------------------------------//
  // INVALID OPTION                                                                       //
  //--------------------------------------------------------------------------------------//
  else {
    fprintf(stderr, "Error: '%s' is an invalid method!\n", method);
    exit(1);
  }


  // clean up
  gsl_rng_free(RANDOM_GENERATOR);
  delete cmdLine;

  return 0;
}



