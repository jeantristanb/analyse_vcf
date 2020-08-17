// Time-stamp: <2008-08-21 16:26:55 zaykind>
// (written by Dmitri Zaykin)

#include "rxc.h"

template <class T> void Swap (T &a, T &b) { T temp = a; a = b; b = temp; }

void Perm (RxC_labels *v2, Rndm &rn, int n) {
  for(int i=0; i<n-1; i++) {
    Swap (v2[i].r, v2[rn(n-i)+i].r);
  }
}

void FillBag (RxC_labels *v2, int **x, int rlen, int clen)
{
  int i, j, v, w=0;
  for(i=0; i<rlen; i++) {
    for(j=0; j<clen; j++) {
      for(v=0; v<x[i][j]; v++) {
	v2[w].r = i;
	v2[w].c = j;
	++w;
      }
    }
  }
}

int Empty(int *Check, int n) {
  for(int i=0; i<n; i++) {
    if(Check[i] != 0) return 0;
  }
  return 1;
}

int** Read(istream& Input_File, int &r, int &c, int& ssize)
{
  Input_File >> r >> c; assert(Input_File.good());
  int i, j, **x = 0; x=Talloc(x,r,c);
  for(ssize=0, i=0; i<r; i++) {
    for(j=0; j<c; j++) {
      Input_File >> x[i][j]; assert(Input_File.good());
      ssize += x[i][j];
    }
  }

  int k, *Check = new int [r];
  for(j=0; j<c; j++) {
    for(k=0; k<r; k++) { Check[k] = 0; }
    for(i=0; i<r; i++) {
      Check[i] = x[i][j] ? 1:0;
    }
    if(Empty(Check,r)) {
      cout << "Empty column #" << (j+1) << " ...\n";
      Tfree (x, r);
      x = 0;
    }        
  }
  delete [] Check;
  return x;
}
