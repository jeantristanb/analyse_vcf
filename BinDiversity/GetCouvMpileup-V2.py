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

ListeCouvMin=[1,3,6]
SommeProf=[0]*(max(ListeCouvMin)+1)
SommeCouv=[0]*(max(ListeCouvMin)+1)
LastBase=0
Chro=""
EcrireSortieCouv.write("Chro\tLong\tCouvMin\tProfSum\tCouvSum\n")
for line in sys.stdin:
    linesplit=line.split()
    if len(linesplit)==6 :
       Couv=float(linesplit[3])
       if Chro != linesplit[0] :
           if SommeProf[ListeCouvMin[0]]>0 :
               #EcrireSortieCouv.write(Chro + "\t" +str(LastBase))
               for CouvMin in ListeCouvMin :
                   EcrireSortieCouv.write(Chro + "\t" +str(LastBase)+"\t"+str(CouvMin)+"\t"+str(SommeProf[CouvMin]) +"\t" + str(SommeCouv[CouvMin])+"\n")
               SommeLastBase+=int(LastBase)
           for CouvMin in ListeCouvMin :
               SommeProf[CouvMin]=0
               SommeCouv[CouvMin]=0
           LastBase=0
           Chro=linesplit[0]
       for CouvMin in ListeCouvMin :
           if Couv >= CouvMin:
               SommeProf[CouvMin]+=Couv/1000.0
               SommeCouv[CouvMin]+=1 
       LastBase=linesplit[1]

for CouvMin in ListeCouvMin :
    EcrireSortieCouv.write(Chro + "\t" +str(LastBase)+"\t"+str(CouvMin)+"\t"+str(SommeProf[CouvMin]) +"\t" + str(SommeCouv[CouvMin])+"\n")


