// Time-stamp: <2010-01-12 14:18:17 zaykind> (Dmitri Zaykin)
//

#ifndef _ATTIC_H
#define _ATTIC_H

#include "mymem.h"
#include "marsa.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

double IntLog(int n);
void SpitMatrix (int **dat, int rlen, int clen, ostream &CoStr);
void SpitMatrix (int **dat, int *r, int *c, int rlen, int clen, ostream &);
int calc_margins(int **dat, int *r, int *c, int rlen, int clen);
int expected(int *r, int *c, int rlen, int clen, int SampleSize, double **e);
double chi(int **dat, double **e, int PutZeros, int rlen, int clen, int SampleSize);
void Spit(int **x, int rlen, int clen);
int** ReadTable(char *InputFile, int &r, int &c);
int** ReadTable(istream& InputFile, int &r, int &c);
bool Is_A_Number(const char *s);
bool Is_An_Integer_Gt_0(const char *s);
bool Is_An_Integer(const char *s);

#endif
