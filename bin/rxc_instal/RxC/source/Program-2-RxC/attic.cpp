// Time-stamp: <2008-08-21 16:24:55 zaykind> (Dmitri Zaykin)
//

#if defined (_Windows)
  #define _WinDog
#endif

#include "mymem.h"
#include "attic.h"

bool Is_A_Number(const char *s)
{
    while(*s && isspace(*s)) ++s;
    while(*s) {
        if(!isdigit(*s) && *s != '.' && *s != '-'
           && tolower(*s) != 'e') return false;
        ++s;
    }
    return true;   
}

bool Is_An_Integer(const char *s)
{
    while(*s && isspace(*s)) ++s;
    while(*s) {
        if(!isdigit(*s) && *s != '-') return false;
        ++s;
    }
    return true;   
}

bool Is_An_Integer_Gt_0(const char *s)
{
    while (*s && isspace(*s)) ++s;
    while (*s) {
        if(!isdigit(*s)) return false;
        ++s;
    }
    return true;
}

void SpitMatrix (int **dat, int rlen, int clen, ostream &CoStr)
{
    CoStr << "rows=" << rlen << "  ";
    CoStr << "cols=" << clen << "\n\nTable:";
    for(int i=0; i<rlen; i++) {
        CoStr << endl;
        for(int j=0; j<clen; j++) {
            CoStr << dat[i][j] << "  ";
        }
    }
    CoStr << endl;
}


void SpitMatrix (int **dat, int *r, int *c, int rlen, int clen, ostream &cost) 
{
    cost << "\nData has been read as:\n\n" << "rows=" << rlen << "  ";
    cost << "cols=" << clen << "\n\nTable:";
    int i, j;
    for(i=0; i<rlen; i++) {
        cost << endl;
        for(j=0; j<clen; j++) {
            cost << dat[i][j] << "  ";
        }
    }
    
    cost << "\n\nMargins are:" << endl;
    cost << "\nRows:\t";
    for(i=0; i<rlen; i++) { cost << r[i] << " "; }
    cost << "\nCols:\t";
    for(j=0; j<clen; j++) { cost << c[j] << " "; }
    cost << endl;
    
#if defined (DEBUG)
    cout << "\n\n" << "s=" << s;
    cout << " Prees a key...";
    cin.get();
#endif   
}

int calc_margins(int **dat, int *r, int *c, int rlen, int clen)
{
    int ssz,i,j;
    for(i=0; i<rlen; i++) { r[i]=0; }
    for(j=0; j<clen; j++) { c[j]=0; }
    for(i=0; i<rlen; i++)    // row totals
        for(j=0; j<clen; j++) { r[i] += dat[i][j]; }
    for(j=0; j<clen; j++)    // col totals
        for(i=0; i<rlen; i++) { c[j] += dat[i][j]; }
    for(ssz=0, i=0; i<clen; i++) ssz += r[i];
    return ssz;
}

int IsEmpty(int *Check, int n)
{
    for(int i=0; i<n; i++) {
        if(Check[i] != 0) return 0;  // not empty
    }
    return 1; // end reached, all zeros
}

int expected(int *r, int *c, int rlen, int clen, int SampleSize, double **e)
{
    int i, j;
    int retval=1;
    for(i=0; i<rlen; i++) {
        for(j=0; j<clen; j++) {
            e[i][j] = (double)r[i] * c[j] / SampleSize;
            if(e[i][j] == 0.0) retval=0;
        }}
    return retval;
}

double chi(int **dat, double **e, int PutZeros, int rlen, int clen, int SampleSize)
{
 double chi2=0.0;
 int i,j;
  for(i=0; i<rlen; i++) {
  for(j=0; j<clen; j++) {
    if(dat[i][j]) {
        chi2 += (double)dat[i][j]*dat[i][j] / e[i][j];
        if(PutZeros) dat[i][j]=0;
    }
  }}
  chi2 -= SampleSize;
  return chi2;
}


void Spit(int **x, int rlen, int clen)
{
    for(int i=0; i<rlen; i++) {
        for(int j=0; j<clen; j++) {
            cout << x[i][j] << "  ";
        }
        cout << endl;    
    }    
}

// read num of rows, num of cols, and then the matrix
int** ReadTable(char *InputFile, int &r, int &c)
{
    ifstream Input_File(InputFile); assert(!Input_File.fail());
    Input_File >> r >> c;
    assert(!Input_File.fail());
    int i, j, **x = 0; x = Talloc(x, r, c);
    for(i=0; i<r; i++) {
        for(j=0; j<c; j++) {
            Input_File >> x[i][j]; assert(!Input_File.fail());
        }
    }

    int k, *Check = new int [r];
    for(j=0; j<c; j++) {
        for(k=0; k<r; k++) { Check[k] = 0; }
        for(i=0; i<r; i++) {
            Check[i] = x[i][j] ? 1:0;
        }
        if(IsEmpty(Check,r)) {
            cout << "Empty column #" << (j+1) << ": exiting...\n";
            exit(1);
        }        
    }
    delete [] Check; 
    return x;
}

int** ReadTable(istream& Input, int &r, int &c)
{
    assert(Input.good());
    Input >> r >> c;
    assert(Input.good());
    int i, j, **x = 0; x = Talloc(x, r, c);
    for(i=0; i<r; i++) {
        for(j=0; j<c; j++) {
            Input >> x[i][j]; assert(Input.good());
        }
    }

    int k, *Check = new int [r];
    for(j=0; j<c; j++) {
        for(k=0; k<r; k++) { Check[k] = 0; }
        for(i=0; i<r; i++) {
            Check[i] = x[i][j] ? 1:0;
        }
        if(IsEmpty(Check,r)) {
            cout << "Empty column #" << (j+1) << ": exiting...\n";
            exit(1);
        }        
    }
    delete [] Check; 
    return x;
}
