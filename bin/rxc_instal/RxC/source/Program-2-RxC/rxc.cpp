// Time-stamp: <2008-08-21 16:27:05 zaykind>
// (written by Dmitri Zaykin)

#include "rxc.h"

double chicdf(double df, double x2) { // chi-square p-value
    return (1.0 - gsl_cdf_chisq_P(x2, df));
}

double exact_kernel(int **dat, int rlen, int clen, bool zout)
{
  double prob=0.0;
  for(int i=0; i<rlen; i++) {
    for(int j=0; j<clen; j++) {
      if(!dat[i][j]) continue;
      prob += gsl_sf_lnfact(dat[i][j]);
      if(zout) dat[i][j]=0;
    }
  }
  return prob;
}

double LikRat(int **x, double **e, int r, int c)
{
  double x2=0;
  for(int i=0; i<r; i++) {
    for(int j=0; j<c; j++) {
      if(x[i][j]) x2 += x[i][j]*log(x[i][j]/e[i][j]);
    }
  }
  return 2.0*x2;
}

double avr2(int **x, int* rw, int* cl, int r, int c, int n)
{
  double r2=0, dn=n;
  for(int i=0; i<r; i++) {
    double p = rw[i]/dn;
    for(int j=0; j<c; j++) {
      double q = cl[j]/dn;
      double Pij = x[i][j]/dn;
      r2 += pow(Pij - p*q, 2) / (p*(1-p)*q*(1-q));
    }
  }
  r2 = r2 * dn * (r-1.0)/r *(c-1.0)/c;
  return r2;
}

double PowDiv(int **x, double **e, int r, int c, double dv=2/3.0)
{
  double x2=0;
  for(int i=0; i<r; i++) {
    for(int j=0; j<c; j++) {
      if(x[i][j]) {
	double d = x[i][j], f=e[i][j];
	x2 += d*(pow(d/f, dv) - 1.0);
      }
    }
  }
  double t = 2.0 / (dv*(dv+1.0));
  return t*x2;
}

int main(int ac, char *av[])
{
  if(ac < 2) {
    cerr << "run as:\n" << av[0] << " NumberOfPermutations [RandomSeed] < SomeFile" << "\n";
    cerr << "(RandomSeed is optional)" << "\n\n";
    cerr << "Input (SomeFile) should look like this:\n";
    cerr << "\t3 6\n";
    cerr << "\t4 1 0 2 1 0\n";
    cerr << "\t7 0 5 0 1 1\n";
    cerr << "\t1 1 1 5 0 4\n\n";
    cerr << "'3' refers to the number of rows,\n";
    cerr << "'6' -- to the number of columns.\n";
    cerr << "Then follows the 3x6 table.\n\n";
    cerr << "In the output, asymptotic and permutational tests are given\n";
    cerr << "for the following statistics:\n";
    cerr << "'avr2': average Correlation^2 over every cell vs. all other cells\n"; 
    cerr << "        (average R^2 over all collapsed 2x2 tables)\n";
    cerr << "'LikR': likelihood ratio statistic\n";
    cerr << "'PowD': Cressie & Read 1984 power divergence statistic with lambda=2/3\n";
    cerr << "'Pear': Pearson's chi-square\n";
    cerr << "'Fish': Fisher's exact test\n";
    return 1;
  }
  const int runs = atoi (av[1]);
  unsigned long s1 = ac==2 ? getpid()*time(0) : atol(av[2]);
  Rndm rn (s1);
  int **x, r, c, sz;
  if ((x = Read(cin, r, c, sz)) == 0) return 1;
  int *rw = new int [r];
  int *cl = new int [c];
  calc_margins(x, rw, cl, r, c);
  double **e=0; e=Talloc(e,r,c);
  expected(rw, cl, r, c, sz, e);

  RxC_labels *v2 = new RxC_labels [sz];
  FillBag (v2, x, r, c);
  double pv, df = (r-1.0)*(c-1.0);

  double w0r = avr2(x, rw, cl, r, c, sz);
  pv = chicdf(df, w0r);
  cout << "avr2 asym = " << pv << "\n";

  double w0g = LikRat(x, e, r, c);
  pv = chicdf(df, w0g);
  cout << "LikR asym = " << pv << "\n";

  double w0d = PowDiv(x, e, r, c);
  pv = chicdf(df, w0d);
  cout << "PowD asym = " << pv << "\n";

  double w0c = chi(x, e, false, r, c, sz);  
  pv = chicdf(df, w0c);
  cout << "Pear asym = " << pv << "\n";

  double w0e = exact_kernel (x, r, c, true);

  int i, cnr=0, cng=0, cnd=0, cnc=0, cne=0;
  for(i=0; i<runs; i++) {
    Perm (v2, rn, sz);
    for(int j=0; j<sz; j++) ++x[v2[j].r][v2[j].c];
    double wir = avr2(x, rw, cl, r, c, sz);
    double wig = LikRat(x, e, r, c);
    double wid = PowDiv(x, e, r, c);
    double wic = chi(x, e, false, r, c, sz);
    double wie = exact_kernel (x, r, c, true);
    if (wir >= w0r) ++cnr;
    if (wig >= w0g) ++cng;
    if (wid >= w0d) ++cnd;
    if (wic >= w0c) ++cnc;
    if (wie >= w0e) ++cne;
  }
  double dr = runs;
  cout << "avr2 perm = " << cnr/dr << "\n";
  cout << "LikR perm = " << cng/dr << "\n";
  cout << "PowD perm = " << cnd/dr << "\n";
  cout << "Pear perm = " << cnc/dr << "\n";
  cout << "Fish perm = " << cne/dr << "\n";
  delete [] v2;
  delete [] rw;
  delete [] cl;
  Tfree(x, r);
  Tfree(e, r);
  return 0;
}
