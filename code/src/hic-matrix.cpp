#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// THIS IS HOW YOU COMPILE AND CALL FROM R:
//
// $ R CMD SHLIB hic_matrix.cpp
//
// > dyn.load("hic_matrix.so")
// > .C("hic_matrix",as.double(x),nrow(x),ncol(x),as.integer(dist),y=double(nnnnn))$y
//


//---- calc_df_explore -----
//
void calc_df_explore(double **x, int **y, int n, int m, int i, int j, double tolerance, int id) 
{
  y[i][j] = id;
  for (int di=-1; di<=1; di++) {
    for (int dj=-1; dj<=1; dj++) {
      if ((di==0)&&(dj==0)) continue;
      if ((i+di>=n)||(i+di<0)) continue;
      if ((j+dj>=m)||(j+dj<0)) continue;
      if ((y[i+di][j+dj]==-1)&&(fabs(x[i+di][j+dj]-x[i][j])<=tolerance)) calc_df_explore(x,y,n,m,i+di,j+dj,tolerance,id);
    }
  }
}




//---- calc_df -----
//
extern "C" void calc_df(double *v, int *n, int *m, int *max_dist, double *tolerance, int *df) 
{
  // convert input vector to matrix
  double **x = new double*[*n];
  if (x==NULL) { fprintf(stderr, "Error [CalcDF]: could not allocate memory!\n"); exit(1); }
  for (int i=0; i<*n; i++) { x[i] = new double[*m]; if (x[i]==NULL) { fprintf(stderr, "Error [CalcDF]: could not allocate memory!\n"); exit(1); } }
  for (int j=0,k=0; j<*m; j++) for (int i=0; i<*n; i++,k++) x[i][j] = v[k];
  //for (int i=0; i<*n; i++) { for (int j=0; j<*m; j++) printf("%.0f ", x[i][j]); printf("\n"); }

  // y will keep track of visited nodes (NOTE: values lower than tolerance and elements away from diagonal are marked as visited!)
  int **y = new int*[*n];
  if (y==NULL) { fprintf(stderr, "Error [CalcDF]: could not allocate memory!\n"); exit(1); }
  for (int i=0; i<*n; i++) { 
    y[i] = new int[*m];   
    if (y[i]==NULL) { fprintf(stderr, "Error [CalcDF]: could not allocate memory!\n"); exit(1); }
    for (int j=0; j<*m; j++) y[i][j] = ((abs(i-j)>*max_dist)||(x[i][j]<*tolerance))?0:-1;            // -1 means "not visited"; zeroes are not contributing to df
  }
  
  // initialize df to zero and run
  *df = 0;
  for (int i=0; i<*n; i++)
    for (int j=0; j<*m; j++) { 
      if (y[i][j]==-1) { 
        if (x[i][j]>0) *df = *df + 1;         // only positive values contribute to df
        calc_df_explore(x,y,*n,*m,i,j,*tolerance,x[i][j]>0?*df:0); 
      }
    }
  //for (int i=0; i<*n; i++) { for (int j=0; j<*m; j++) printf("%2d ", y[i][j]); printf("\n"); }

  // cleanup
  delete [] x;
  delete [] y;
}



//---- rotate45 -----
//
extern "C" void rotate45(double *v, int *n, int *m, int *dist, double *y) 
{
  // convert input vector to matrix
  double **x = new double*[*n];
  for (int i=0; i<*n; i++) x[i] = new double[*m];
  for (int j=0,k=0; j<*m; j++) for (int i=0; i<*n; i++,k++) x[i][j] = v[k];
  
  // rotate
  for (int q=0,d=0; d<=*dist; d++) {
    for (int r=1; r<=d/2; r++) y[q++] = 0;
    for (int k=0; k<*n-d; k++,q++) y[q] = x[k+d][k];
    for (int r=1; r<=(d+1)/2; r++) y[q++] = 0;
  }
//  for (int i=0; i<N; i++) { printf("%.0f ", y[i]); if ((i+1)%(*n)==0) printf("\n"); }

  // cleanup
  delete [] x;
}


//---- inverse_rotate45 -----
//
extern "C" void inverse_rotate45(double *v, int *n, int *m, double *y) 
{
  // convert input vector to matrix
  double **x = new double*[*n];
  for (int i=0; i<*n; i++) x[i] = new double[*m];
  for (int j=0,k=0; j<*m; j++) for (int i=0; i<*n; i++,k++) x[i][j] = v[k];
//  for (int i=0; i<*n; i++) { for (int j=0; j<*m; j++) printf("%.0f%c", x[i][j], j==*m-1?'\n':'\t'); }
  
  // initialize output matrix to zero
  int N = (*n)*(*n);
  double *ycurr=y, *ylast=y+N;
  while (ycurr!=ylast) *(ycurr++) = 0;

  // inverse-rotate  
  for (int d=0; d<*m; d++)
    for (int i=0; i<*n-d; i++) {
      int k = i+(*n)*(i+d);                     // TODO: this replaces upper triangle with lower triangle (does not matter if matrix is symmetric, as it should be)
      y[k] = x[i+d/2][d];
      //printf("x[%d][%d] => y[%d][%d] = %.0f\n", i+d/2, d, i, i+d, x[i+d/2][d]);
    }
//  for (int i=0; i<N; i++) printf("%.0f%c", y[i], (i+1)%(*n)==0?'\n':'\t');

  // cleanup
  delete [] x;
}



//---- calc_counts_per_distance -----
//
extern "C" void calc_counts_per_distance(double *v, int *n, double *y) 
{
  // convert input vector to matrix
  int m = *n;
  double **x = new double*[*n];
  for (int i=0; i<*n; i++) x[i] = new double[m];
  for (int j=0,k=0; j<m; j++) for (int i=0; i<*n; i++,k++) x[i][j] = v[k];
  
  // calculate mean value as a function of distance from diagonal
  for (int q=0,d=0; d<*n; d++) {
    y[d] = 0;
    for (int k=0; k<*n-d; k++,q++) y[d] += x[k][k+d] + x[k+d][k];
    y[d] /= 2*(*n-d);
  }

  // cleanup
  delete [] x;
}




//---- impute_rows -----
//
void impute_rows(double **mat, int n, double NaN)
{
  for (int r=0; r<n; r++) {
    int k = 0;
    while (k<n) {
      if (mat[r][k]==NaN) {
        int k0 = k;
        int k1 = k;
        for ( ; k1+1<n; k1++) if (mat[r][k1+1]!=NaN) break;
        if ((k0>0)||(k1<n)) { 
          double v0 = k0>0 ? mat[r][k0-1] : 0;
          double v1 = k1<n-1 ? mat[r][k1+1] : 0;
          double a = (v1-v0)/(k1-k0+2);
          if ((k0>0)||(k1<n-1)) {
            //fprintf(stderr, "Row %d: (k0,k1)=(%d,%d) (v0,v1)=(%.2f,%.2f) a=%.3f\n", r, k0, k1, v0, v1, a); 
            for (int z=k0; z<=k1; z++) {
              mat[r][z] = v0 + a*(z-k0+1);
              //fprintf(stderr, "-- imputing mat[%d][%d] <- %.2f\n", r, z, mat[r][z]); 
            }
          }
        }
        k = k1+1;
      } else {
        k = k+1;
      }
    }
  }
}


//---- impute -----
//
extern "C" void impute(double *v, int *n, double *NaN, double *y) 
{
  // NOTE: this operation is only valid for symmetric matrices
  double **x = new double*[*n];
  for (int i=0; i<*n; i++) x[i] = new double[*n];
  for (int j=0,k=0; j<*n; j++) for (int i=0; i<*n; i++,k++) x[i][j] = v[k];
  
  // impute
  impute_rows(x,*n,*NaN);
  
  // transpose
  for (int j=0; j<*n; j++) for (int i=j+1; i<*n; i++) { double t = x[i][j]; x[i][j] = x[j][i]; x[j][i] = t; }
  
  // impute again
  impute_rows(x,*n,*NaN);
  
  // convert to vector
  for (int j=0,k=0; j<*n; j++) for (int i=0; i<*n; i++,k++) y[k] = x[i][j];
    
  // cleanup
  delete [] x;
}






