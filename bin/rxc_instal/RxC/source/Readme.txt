(1) Correlation-based tests for linkage disequilibrium

Compilation under Linux
-----------------------

At the terminal window, cd to the directory with the files and type

 g++ mcld1.cpp mcld2.cpp mcld3.cpp dcdflib.cpp -O3 -s -Wall -o mcld.x

That will produce the executable, "mcld.x"


Compilation under Windows
-------------------------

Install Cygwin (http://www.cygwin.com/). During the setup, add g++ to
the packages to be installed. Start the Cygwin bash shell, cd to the
directory with the files, and type

 g++ mcld1.cpp mcld2.cpp mcld3.cpp dcdflib.cpp -O3 -s -Wall -o mcld.exe -mno-cygwin

That will produce the executable, "mcld.exe"


(2) Correlation-based tests for RxC contingency tables

This program uses two functions from "GNU Scientific Library"
(http://www.gnu.org/software/gsl/): log of factorial and the chisquare
CDF. This library needs to be installed. 

Compilation under Linux
-----------------------

To install the GNU Scientific Library on Fedora Linux, just type 
"yum install gsl gsl-devel". There should be a similar to do this
via a package manager with any major Linux distribution.

Then the compilation command is

  g++ -o rxc.x rxc.cpp rxc1.cpp attic.cpp -Wno-deprecated -Wall -O3 -s -lgsl -lgslcblas

"-lgsl -lgslcblas" is to link the GNU Scientific Library functions.
That will produce the executable, "rxc.x"

Compilation under Windows
-------------------------

Install Cygwin (http://www.cygwin.com/). During the setup, add g++ and
GSL to the packages to be installed. Start the Cygwin bash shell, cd
to the directory with the files, and type

  g++ -o rxc.exe rxc.cpp rxc1.cpp attic.cpp -Wno-deprecated -Wall -O3 -s -lgsl -lgslcblas -mno-cygwin

"-lgsl -lgslcblas" is to link the GNU Scientific Library functions.
That will produce the executable, "rxc.exe".
