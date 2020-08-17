wget  https://www.niehs.nih.gov/research/atniehs/labs/assets/docs/q_z/rxc.zip
unzip rxc.zip
cd RxC/source/Program-1-LD
g++ mcld1.cpp mcld2.cpp mcld3.cpp dcdflib.cpp -O3 -s -Wall -o mcld
mv mcld ../../../../
cd ../../../
cd RxC/source/Program-2-RxC/
g++ -o rxc rxc.cpp rxc1.cpp attic.cpp -Wno-deprecated -Wall -O3 -s -lgsl -lgslcblas -L /home/jeantristan/bin/gsl/lib -I /home/jeantristan/bin/gsl/include/ 
mv ./rxc  ../../../../
#cd ../../../
#rm -rf rxc.zip RxC


