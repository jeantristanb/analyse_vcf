// Time-stamp: <2006-06-30 15:59:54 zaykind>
// (written by Dmitri Zaykin)
//
// Compilation using GNU g++/gcc:
// g++ -c dcdflib.cpp -O3 -Wall
// g++ mcld1.cpp mcld2.cpp mcld3.cpp dcdflib.o -O3 -s -Wall -o mcld.x

#include "mcld.h"

template <class T>
void Swp (T &a, T &b) {
  T tmp = a; a = b; b = tmp;
}

template <class T>
void perm(T *t, int n, Rndm &rn) {
  for(int i=0; i<n; i++) {
    Swp(t[i], t[i + rn(n-i)]);
  }
}

template<class T>  // allocating 2D array
T **Talloc (T **d, int fi, int se) {
    d = new T*[fi];
    for(int i=0; i<fi; i++) { d[i] = new T[se]; }
    return d;
}

template<class T> // freeing 2D array
void Tfree (T **x, int fi) {
    for(int i=0; i < fi; i++) { delete [] x[i]; }
    delete [] x;
}

Vi* Read(const char *fnam, const char *miss, int *nvi)
{
  ifstream Cin(fnam);
  string mis = miss;
  if(!Cin.good()) {
    cerr << "couldn't open \"" << fnam << "\" for reading\n"; exit(1);
  }
  string s;
  Vi ncols;
  while(getline(Cin, s)) {
    istringstream is(s);
    string x;
    int k=0;
    while(is >> x) ++k;
    if(k) ncols.push_back(k);
  }
  Cin.close();

  int Errs=0;
  for(size_t i=1; i<ncols.size(); i++) {
    if(ncols[i] != ncols[i-1]) {
      ++Errs;
      cerr << "Error: unequal number of columns --\n";
      cerr << "\tRow " << i << " has " << ncols[i-1] << " columns\n";
      cerr << "\tRow " << (i+1) << " has " << ncols[i] << " columns\n";
    }
  }
  if(Errs) exit(1);
  int nc = ncols[0];
  if(nc % 2) { cerr << "Error: odd number of columns\n"; exit(1); }
  
  int nlo = nc/2;
  Ss *ss = new Ss [nlo];
  Vs *vs = new Vs [nlo*2];
  *nvi = nlo*2;
  Vi *vi = new Vi [nlo*2];
  Ms *ms = new Ms [nlo];

  //Cin.open(fnam); // <-- doesn't work on Cygwin !?!
  ifstream Cin1(fnam);
  while(Cin1) {
    for (int i=0; i < nlo*2; i++) {
      string s;
      if(!(Cin1 >> s)) break;
      ss[i/2].insert(s);
      vs[i].push_back(s);
    }
  }
  Cin1.close();

  for(int i=0; i<nlo; i++) {
    int v=0;
    for(Ss::iterator j=ss[i].begin(); j != ss[i].end(); j++) {
      if(*j != mis) ms[i][*j] = v++;
      else ms[i][*j] = -1;
    }
  }

  Vi tmp(vs[0].size());
  for (int i=0; i<nlo*2; i++) vi[i] = tmp;

  for(size_t j=0; j<vs[0].size(); j++) {
    for(int i=0; i<nlo*2; i+=2) {
      int a[2] = { 1+ms[i/2][vs[i][j]], 1+ms[i/2][vs[i+1][j]] };
      a[0] = a[0] == 0 ? -1 : a[0] - 1;
      a[1] = a[1] == 0 ? -1 : a[1] - 1;
      if (a[0] > a[1]) swap(a[0], a[1]);
      vi[i][j] = a[0];
      vi[i+1][j] = a[1];
    }
  }

  delete [] ss;
  delete [] vs;
  delete [] ms;

  return(vi);
}

double Cor(int** x, int n) // x is n by 2
{ // this function is modified from pseudocode in
  // en.wikipedia.org/wiki/Correlation
  double ssqx = 0, ssqy = 0, ssq_xy = 0;
  double xbar = x[0][0];
  double ybar = x[0][1];
  for(int i=1; i<n; i++) {
    double dx = x[i][0] - xbar;
    double dy = x[i][1] - ybar;
    double ii = i+1;
    double sweep = (ii - 1.0) / ii;
    ssqx += dx * dx * sweep;
    ssqy += dy * dy * sweep;
    ssq_xy += dx * dy * sweep;
    xbar += dx / ii;
    ybar += dy / ii;
  }
  double Sxx = sqrt(ssqx / n);
  double Syy = sqrt(ssqy / n);
  double Sxy = ssq_xy / n;
  double r = Sxy / (Sxx * Syy);
  return(r);
}

double ChisqPval(double x, double df) {
  double q, p, bound;
  int which=1, status;
  cdfchi(&which, &p, &q, &x, &df, &status, &bound);
  return 1.0 - p;
}

template <class T> const T Abs(const T& a) {return a >= 0 ? a : -a;}

double GetR(int **y, const int *ix, int **x, Vi& a1, Vi& a2,
	    int n, double *dpr, double *ApprPv)
{
  int n1 = a1.size(), n2 = a2.size();
  double r2 = 0;
  if(ApprPv != 0) *dpr=0;
  for(int i=0; i<n1; i++) {
    for(int j=0; j<n2; j++) {
      int ii = a1[i], jj = a2[j];
      for(int t=0; t<n; t++) {
	int k = ix[t];
	x[t][0] = ((y[t][0] == ii) + (y[t][1] == ii)) - 1;
	x[t][1] = ((y[k][2] == jj) + (y[k][3] == jj)) - 1;
      }
      double tr = Cor(x, n);
      r2 += tr*tr;
      if(ApprPv != 0) {
	double tdpr = DprIJ(x, n);
	*dpr += Abs(tdpr);
      }
    }
  }
  if(ApprPv != 0) {
    double df = (n1-1.0)*(n2-1.0);
    double ncells = double(n1)*n2;
    double x2 = (n*df/ncells) * r2;
    *ApprPv = ChisqPval(x2, df);
    *dpr /= ncells;
  }
  return (r2);
}

void permITa(int **y, int n, // yy is n by 4
	     Vi& a1, Vi& a2, int o, int seed, double PvThresh,
	     double *pr2, double *dpr, double *ppv, double *ApprPv)
{
  Rndm rn(seed);
  int *ix = new int [n];
  for(int i=0; i<n; i++) ix[i]=i;
  int **x = 0; x = Talloc(x, n, 2);
  int cn = 0;

  double t0 = GetR(y, ix, x, a1, a2, n, dpr, ApprPv);

  if(*ApprPv <= PvThresh) {
    for(int s=0; s<o; s++) {
      perm(ix, n, rn);
      double ti = GetR(y, ix, x, a1, a2, n, dpr, 0);
      if(ti >= t0) ++cn;
    }
    if(o) {
      double pv = double(cn)/o;
      *ppv = pv;
    }
  }
  int n1 = a1.size(), n2 = a2.size();
  *pr2 = t0 / (n1*n2);
  Tfree(x,n);
  delete [] ix;
}

int main(int ac, const char **av)
{
  string FileName, MissVal;
  int sim;
  double PvalThresh;
  unsigned long seed;
  ReadArgs(ac, av, FileName, MissVal, sim, PvalThresh, seed);

  int nvi;
  Vi *vi = Read(FileName.c_str(), MissVal.c_str(), &nvi);
  int nlo = nvi/2;
  int n = vi[0].size();
  int **y=0; y = Talloc(y, n, 4);

  cout << "Loci   Aver. |Delta-prime|   Aver. |Corr.|";
  cout << "     Approx. p-value   Permut. p-value\n";
  for(int i=0; i<nlo-1; i++) {
    for(int j=i+1; j<nlo; j++) {
      cout << (i+1) << "/" << (j+1) << " ";
      int c1=i*2,c2=c1+1,c3=j*2,c4=c3+1;
      int kk = 0;
      for(int k=0; k<n; k++) {
	if(vi[c1][k] != -1 && vi[c2][k] != -1 &&
	   vi[c3][k] != -1 && vi[c4][k] != -1)
	  { // skip over missings
	    y[kk][0] = vi[c1][k]; y[kk][1] = vi[c2][k];
	    y[kk][2] = vi[c3][k]; y[kk][3] = vi[c4][k];
	    ++kk;
	  }
      }

      // need to build allele sets here because they're
      // possibly reduced by skipping missing data
      Si s1, s2;
      for(int j=0; j<kk; j++) {
	s1.insert(y[j][0]); s1.insert(y[j][1]);
	s2.insert(y[j][2]); s2.insert(y[j][3]);
      }
      Vi a1, a2;
      for(Si::iterator s=s1.begin();s != s1.end();s++) {
	a1.push_back(*s);
      }
      for(Si::iterator s=s2.begin();s != s2.end();s++) {
	a2.push_back(*s);
      }

      double pr2=0, dpr=0, ApprPv=1, ppv = -1;
      if(a1.size() > 1 && a2.size() > 1) {
	permITa(y, kk, a1, a2, sim, seed, PvalThresh, 
		&pr2, &dpr, &ppv, &ApprPv);
      }
      cout << setw(17) << dpr << " ";
      cout << setw(17) << sqrt(pr2) << " ";
      cout << setw(17) << ApprPv << " ";
      if(ppv == -1) {
        cout << setw(12) << "> " << PvalThresh << endl;
      }
      else {
        cout << setw(17) << ppv << endl;
      }
    }
  }

  Tfree(y,n);
  delete [] vi;
}
