The programs here are to accompany

Zaykin, DV, Pudovkin AI, Weir BS 2008. Correlation-based inference for
linkage disequilibrium with multiple alleles. Genetics (in press).

They are compiled for the Windows's "Command Line". The C++ source
code is also provided that was developed and should compile under
Linux, or any other system with the GNU compiler (g++)

dz-rxc.exe -- computes correlation-based tests for 
RxC contingency tables (as well as several other tests)

mcld.exe -- computes correlation-based tests for linkage
disequilibrium. A sample input file is "loci.txt" in this directory
(this file is modelled after LD found in STRs in Rosenberg NA et al.,
Science 2002; 298:2981-2985)

To try out the programs -

1) Start the "Command Prompt" window.
(In Windows XP it is in Programs -> Accessories -> Command Prompt)

2) Change to the directory where you stored the programs, e.g. type
"cd C:\MyStuff"

3) Run one of the programs to get help screen on its usage, e.g.
type at the prompt "dz-rxc.exe" or "mcld.exe"

4) To run sample files downloaded from this directory, type
"mcld.exe -file=loci.txt" or "dz-rxc.exe 10000 < rxcdat.txt"

5) you can direct output to a file as e.g.
"dz-rxc.exe 10000 < rxcdat.txt > Results.txt"
