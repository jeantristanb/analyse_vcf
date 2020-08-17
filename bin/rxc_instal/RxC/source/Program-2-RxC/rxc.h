// Time-stamp: <2008-08-21 16:26:24 zaykind>
// (written by Dmitri Zaykin)
//
// Note: "-lgsl -lgslcblas" and gsl includes below refer
// to "GNU Scientific Library",
// http://www.gnu.org/software/gsl/
// On Fedora Linux, install it as "yum install gsl gsl-devel"
// Compilation:
// g++ -o rxc.x rxc.cpp rxc1.cpp attic.cpp -Wno-deprecated -Wall -O3 -s -lgsl -lgslcblas

// GNU Scientific Library headers
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>

#include <time.h>
#if defined (__GNUG__)
   #include <unistd.h>
#else
   unsigned long getpid() { return 1; }
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "marsa.h"
#include "mymem.h"
#include "attic.h"
//#include "statstd.h"

struct RxC_labels { int r, c; };

double exact_kernel (int **dat, int rlen, int clen, bool zout);
void FillBag (RxC_labels *v2, int **x, int rlen, int clen);
int** Read(istream& Input_File, int &r, int &c, int& ssize);
void Perm (RxC_labels *v2, Rndm &rn, int n);
