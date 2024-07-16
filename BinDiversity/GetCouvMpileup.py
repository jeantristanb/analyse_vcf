#!/usr/bin/python
import fileinput
import sys


if len(sys.argv)!=2 :
   print "Exe SortieCouv SortieResume"
   sys.exit()


try :
    EcrireSortieCouv=open(sys.argv[1], 'w')
except :
   print "file "+ sys.argv[1] + " Can't open"



cmt=1
BaliseEcrire=False
SommeProfTot=0
NbLocusCouvTot=0
SommeLastBase=0

SommeProfChro=0
NbLocusCouvChro=0
LastBase=0
Chro=""
for line in sys.stdin:
    linesplit=line.split()
    if len(linesplit)==6 :
       if Chro != linesplit[0] :
           if NbLocusCouvChro>0 :
               EcrireSortieCouv.write(Chro + "\t" + str(SommeProfChro) + "\t" +str(NbLocusCouvChro) +"\t"+ str(LastBase) +"\n")
               SommeProfTot+=SommeProfChro
               NbLocusCouvTot+=NbLocusCouvChro
               SommeLastBase+=int(LastBase)
           SommeCouvChro=0
           SommeProfChro=0
           NbLocusCouvChro=0
           LastBase=0
           Chro=linesplit[0]
       SommeProfChro+=float(linesplit[3])/1000.0
       NbLocusCouvChro+=1 
       LastBase=linesplit[1]

EcrireSortieCouv.write(Chro + "\t" + str(SommeProfChro) + "\t" +str(NbLocusCouvChro) +"\t"+ str(LastBase) +"\n")
EcrireSortieCouv.write("Total\t" + str(SommeProfTot) + "\t" +str(NbLocusCouvTot) +"\t"+ str(SommeLastBase) +"\n")
EcrireSortieCouv.close()
    
