// Time-stamp: <2008-04-01 19:41:50 zaykind>
// (written by Dmitri Zaykin)
//
// Command line processing code

#include "mcld.h"

int InString(const string& t, const string& stuff) {
  return t.find_first_of(stuff, 0) < t.length();   
}

void split_string (const string& tok, const string& sep, Vs& w)
{
  w.resize(0);
  int n = tok.length();
  int start=tok.find_first_not_of(sep), stop;
  while (start >= 0 && start < n) {
    stop = tok.find_first_of(sep, start);
    if (stop < 0 || stop > n) stop = n;
    w.push_back(tok.substr(start, stop-start));
    start = tok.find_first_not_of(sep, stop+1);
  }
}

struct Par {
  string s, v;
  Par() {s=v="";}
  Par(const string s_, const string v_) : s(s_), v(v_) {}
  Par(const string s_, unsigned long v_) : s(s_) {
    ostringstream os;
    os << v_;
    v = os.str();
  }
  int operator < (const Par& w) const { return s < w.s; }
  int operator != (const Par& w) const { return s != w.s; }
  int operator == (const Par& w) const { return s == w.s; }
};

typedef map <string, Par, less<string> > Ps;

template <class T> void stri2t(const string& s, T& x)
{
  istringstream i(s);
  if (!(i >> x)) {
    string err = "Cannot convert \"" + s + "\"";
    cerr << err << "\n";
    exit(1);
  }
}

void Instruct(const char *av0, Ss &ss, Ps &ps) 
{
  cerr << "Need command line arguments entered in any order, e.g.\n\n";

  Vs wav;
  split_string(av0, "\\/", wav);
  cerr << wav[wav.size()-1] << " ";

  ldiv_t qt = ldiv((long)time(0), (long)100000); // just making "seed" shorter,
  Par tmp("-seed", (unsigned long)qt.rem);       // for the screen output
  ps["-seed"] = tmp;

  for(Ss::iterator j=ss.begin(); j != ss.end(); j++) {
    cerr << ps[*j].s << "=" << ps[*j].v << " ";
  }
  cerr << "> Results.txt\n\n";
  cerr << "At a minimum, you have to supply -file=YourDataFile.\n";
  cerr << "Other parameters will take default values as shown.\n";
  cerr << "Input file format is a text file with rows as individuals,\n";
  cerr << "and two columns per locus, holding allele labels.\n";
  cerr << "Allele labels and missing data label can be any strings, or numbers.\n";

  cerr << "\n";
  cerr << "-minp= denotes the p-value threshold below which\n";
  cerr << "\tthe permutational test is performed.\n";

  cerr << "\n";
  cerr << "-miss= denotes the missing data value\n";

  cerr << "\n";
  cerr << "-perm= denotes the number of permutaions\n";
  cerr << "\tfor the permutational test.\n";

  cerr << "\n";
  cerr << "-seed= is an integer seed for the random generator;\n";
  cerr << "\tif omitted, it is set to the computer time.\n";

  cin.get();
  exit(1);
}

void ReadArgs(int ac, const char** av, 
	      string & FileName,
	      string & MissVal,
	      int &sim,
	      double &PvalThresh,
	      unsigned long &seed)
{
  const int NumPar=5;
  Par par[NumPar] = {
    Par("-file", "***"),
    Par("-miss", "?"),
    Par("-perm", "19999"),
    Par("-minp", "0.05"),
    Par("-seed", time(0))
  };

  Ps ps;
  Ss sp;

  for(int i=0; i<NumPar; i++) {
    ps[par[i].s] = par[i];
    sp.insert(par[i].s);
  }

  if(ac==1) Instruct(av[0], sp, ps);

  for(int i=1; i<ac; i++) {
    string s(av[i]);
    if(!InString(s, "=")) Instruct(av[0], sp, ps);
    Vs w;
    split_string(s, "=", w);
    if(sp.find(w[0]) == sp.end()) Instruct(av[0], sp, ps);
    else {
      Par tmp(w[0], w[1]);
      ps[tmp.s] = tmp;
    }
  }

  stri2t(ps["-file"].v, FileName);
  if(FileName == "***") {
    cerr << "Need an input file parameter set as:\n";
    cerr <<  av[0] << " -file=***" << "\n";
    exit(1);
  }
  stri2t(ps["-miss"].v, MissVal);
  stri2t(ps["-perm"].v, sim);
  stri2t(ps["-minp"].v, PvalThresh);
  stri2t(ps["-seed"].v, seed);
}
