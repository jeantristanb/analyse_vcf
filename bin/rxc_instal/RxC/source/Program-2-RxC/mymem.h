// Time-stamp: <2008-08-21 16:25:29 zaykind> (Dmitri Zaykin)
//

#ifndef DZ_MYMEM_H
#define DZ_MYMEM_H

#include "marsa.h"

#include <fstream>
#include <ostream>
#include <iomanip>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

template <class T>
int InRange(const T& mn, const T& mx, T& x) {
    return x >= mn && x <= mx;
}

template <class T>
void Zero2D (T **x, int rows, int cols)
{
    for (int i=0; i<rows; i++)
        memset(x[i], 0, cols*sizeof(T));
}

template <class T>
void Zero1D (T *x, int n)
{
    memset(x, 0, n*sizeof(T));
}

inline int IsOdd (int n)
{
    return n % 2;
}

template <class T> const T& Min (const T& a, const T& b)
{
    return a <= b ? a : b;
}

template <class T> const T& Max (const T& a, const T& b)
{
    return a > b ? a : b;
}

template <class T> const T Abs (const T& a)
{
    return a >= 0 ? a : -a;
}

template <class T>
void DzSwap (T &a, T &b)
{
    T temp = a;
    a = b;
    b = temp;
}

template<class T>
void TriFree (T ***x, int fir, int sec)
{
    int i;
    for(i=0; i < fir; i++) {
        for(int j=0; j < sec; j++) delete [] x[i][j];
    }
    for(i=0; i < fir; i++) delete [] x[i];
    delete [] x;
}

template <class T>
T*** TriAlloc (T*** d, int fir, int sec, int thi)
{
    d = new T** [fir];
    for(int i=0; i<fir; i++) {
        d[i] = new T* [sec];
        for(int j=0; j<sec; j++) {
            d[i][j] = new T [thi];
        }
    }
    return d;
}

template <class T>
T*** TriCopy(T*** dest, const T*** sour, int fir, int sec, int thi)
{
    for(int f=0; f<fir; f++) {
        for(int s=0; s<sec; s++) {
            for(int t=0; t<thi; t++) {
                dest[f][s][t] = sour[f][s][t];
            }
        }
    }
    return dest;
}

template<class T>
void FillVec (T *x, int n, T w)
{
    for(int i=0; i<n; i++) x[i]=w;
}

template<class T>
T *ArrDup (const T *s, int len)
{
    if(s==0) return 0;
    T *p = new T[len];
    for(int i=0; i<len; i++) p[i] = s[i];
    return p;
}

template<class T> 
T *ArrCpy (T *d, const T *s, int len)
{
    for(int i=0; i<len; i++) d[i] = s[i];
    return d;
}

template<class T>  // allocating 2D array
T **Talloc (T **d, int fi, int se) {
    d = new T*[fi];
    for(int i=0; i<fi; i++) { d[i] = new T[se]; }
    return d;
}

template<class T>  // allocating 2D array
T **TCalloc (T **d, int fi, int se, T with=0)
{
    d = new T*[fi];
    for(int i=0; i<fi; i++) {
        d[i] = new T[se];
        for(int j=0; j<se; j++) { d[i][j] = with; }
    }
    return d;
}

template<class T> // freeing 2D array
void Tfree (T **x, int fi) {
    for(int i=0; i < fi; i++) { delete [] x[i]; }
    delete [] x;
}

template<class T> 
T **DArrDup (T **s, int r, int c)
{
    if(s==0) return 0;
    T **p = 0;
    p = Talloc(p, r, c);
    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            p[i][j] = s[i][j];
        }
    }
    return p;
}

template<class T> 
T **DArrCpy (T **d, T **s, int r, int c)
{
    if(s==0) return 0;
    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            d[i][j] = s[i][j];
        }
    }
    return d;
}

enum { LineLen = 1000 };
enum { EOLINT = -9999 };

inline char *StrDup (const char *src)
{
    return strcpy(new char[strlen(src)+1], src);
}

template <class T>
struct Sarr {  // simple array
    T *w;
    int n;
    void Init (const T *s, int _n) {
        n = _n;
        if (w == 0) w = ArrDup (s, n);
        else { 
            delete [] w;
            w = ArrDup (s, n);
        } 
    }
    Sarr () : w(0), n(0) {}
    Sarr (int _n, bool init_0 = false) : n(_n) {
        w = new T [n];
        if(init_0 == true) Zero1D (w, n);
    }
    Sarr (const T *s, int _n) : n(_n) { w = 0; Init (s, n); }
    ~Sarr () { delete [] w; }
};

struct SimpleString {
    char *w;
    void Init (const char *s) { 
        if (w == 0) w = StrDup (s);
        else { 
            delete [] w;
            w = StrDup (s);
        } 
    }
    SimpleString () { w = 0; };
    SimpleString (int n) { Zero1D (w = new char [n], n); };
    SimpleString (const char *s) { w = 0; Init (s); }
    ~SimpleString () { delete [] w; }
};

template <class T>
void perm(T *t, int n, Rndm &rn)
{
   for(int i=0; i<n; i++) {
       DzSwap (t[i], t[i + rn(n-i)]);
   }
}

template <class T>
void BootStrap (T *to, T *from, int n, Rndm &rn)
{
    for (int i=0; i<n; i++) {
        to[i] = from[rn(n)];
    }
}

template <class T>
class twoD  // very simple 2D array that is 1D internally
{
    T *x;
    int r, c;
  public:
    twoD (int rw, int cl) : r(rw), c(cl) { x = new T [r*c]; }
    ~twoD () { delete [] x; }
    T& operator () (int i, int j) { return x[c*i+j]; }
};

template <class T> class triD  // very simple 3D array that is 1D internally
{
    T *x;
    int r, c, s;
  public:
    triD (int rw, int cl, int sl) : r(rw), c(cl), s(sl) {
        assert (double(c)*r*s <= INT_MAX); // for segmented dudes
        x = new T [r*c*s];
    }
    ~triD () { delete [] x; }
    T& operator () (int i, int j, int k) { return x[(i*c+j)*r+k]; }
};

double CumulSearch(double *dvect, double what, int len, int *position);
float CumulSearch(float *dvect, float what, int len, int *position);

int IsValidStr (const char *w, const char *set);
int IsEmpty(char *word);
int IsEmpty(int *Check, int n);
int CountChars(const char *w, char ch);
int GetNumWords(char *ww, const char *set = " \t");
double CumulSearch(double *dvect, double what, int len, int *position);
int* IntDup (const int *w, int * id_len = 0);
int IntLen (const int* x);
int IntCmp (const int* a, const int* b);

#endif
