// Time-stamp: <2006-06-30 15:47:31 zaykind> (written by Dmitri Zaykin)
//
// "Composite D-prime" code
//
// After: Zaykin DV 2004 Bounds and normalization of the composite
// linkage disequilibrium coefficient. Genet Epidemiol. 27(3):252-257.

#include "mcld.h"

// returns sample covariance estimated as Sxy/(n-1)
double Cov(int** x, int n) // x is n by 2
{ // this function is modified from pseudocode in
  // en.wikipedia.org/wiki/Correlation
  double ssq_xy = 0;
  double xbar = x[0][0];
  double ybar = x[0][1];
  for(int i=1; i<n; i++) {
    double dx = x[i][0] - xbar;
    double dy = x[i][1] - ybar;
    double ii = i+1;
    double sweep = (ii - 1.0) / ii;
    ssq_xy += dx * dy * sweep;
    xbar += dx / ii;
    ybar += dy / ii;
  }
  double Sxy = ssq_xy / (n-1.0);
  return Sxy;
}

int Sum(int **q, int n, int dim) {
  int s=0;
  for(int i=0; i<n; i++) s += q[i][dim];
  return s;
}

int Cnt(int **q, int n, int dim, int what) {
  int c=0;
  for(int i=0; i<n; i++) if(q[i][dim] == what) ++c;
  return c;
}

template <class T> const T& Min (const T& a, const T& b) {
  return a <= b ? a : b;
}

template <class T> const T& Max (const T& a, const T& b) {
  return a > b ? a : b;
}

// q is (n by 2) matrix of (-1;0;1) recoded values;
// two columns for two SNPs; n is the number of indiv
double DprIJ(int **q, int n) {
  double naa = Cnt(q,n,0,-1); // in R: length(q[q[,1] == -1,1])
  double nbb = Cnt(q,n,1,-1); // in R: length(q[q[,2] == -1,2])
  double nAa = Cnt(q,n,0, 0); // in R: length(q[q[,1] ==  0,1])
  double nBb = Cnt(q,n,1, 0); // in R: length(q[q[,2] ==  0,2])
  double nAA = Cnt(q,n,0, 1); // in R: length(q[q[,1] ==  1,1])
  double nBB = Cnt(q,n,1, 1); // in R: length(q[q[,2] ==  1,2])
  double pa = 0.5 - Sum(q,n,0)/(2.0*n);
  double pb = 0.5 - Sum(q,n,1)/(2.0*n);
  double delta = Cov(q, n);
  double d,s,mxd,ld = (n-1.0)/n * delta/2;
  if(delta > 0) {
    d = Min(nAA,nBB) + Min(naa,nbb);
    s = Max(n-d-nAa-nBb, 0.0);
    mxd = (d-s)/(2.0*n) - (1-2*pa)*(1-2*pb)/2;
  }
  else {
    d = Min(nAA,nbb) + Min(naa,nBB);
    s = Max(n-d-nAa-nBb, 0.0);
    mxd = (d-s)/(2.0*n) + (1-2*pa)*(1-2*pb)/2;
  }
  double dpr = ld/mxd;
  return dpr;
}
