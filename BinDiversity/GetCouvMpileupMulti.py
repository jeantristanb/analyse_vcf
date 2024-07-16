#!/usr/bin/python
import fileinput
import sys

def GetMatrice(Nrow, Ncol):
    M = [[0 for j in range(0,Ncol)] for i in range(0,Nrow)]
    return M

if len(sys.argv)!=2 :
   print "Exe SortieCouv"
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

ListeCouvMin=[1,3,6,15, 30,50,100]
#SommeProf=[0]*(max(ListeCouvMin)+1)
#SommeCouv=[0]*(max(ListeCouvMin)+1)
LastBase=0
Chro=""
#EcrireSortieCouv.write("Chro\tLong\tCouvMin\tProfSum\tCouvSum\n")
Cmt=0
for line in sys.stdin:
    linesplit=line.split('\t')
    if Cmt==0 :
       NbInd=int((len(linesplit)-3)/3.0)
       SommeProf=GetMatrice(NbInd,max(ListeCouvMin)+1)  
       SommeCouv=GetMatrice(NbInd,max(ListeCouvMin)+1)  
       Chro=linesplit[0]
    if Chro != linesplit[0] :
       for CouvMin in ListeCouvMin :
           EcrireSortieCouv.write(Chro + "\t" +str(LastBase)+"\t"+str(CouvMin))
           for CmtInd in range(0,NbInd) :
               EcrireSortieCouv.write("\t"+str(SommeProf[CmtInd][CouvMin]) +"\t" + str(SommeCouv[CmtInd][CouvMin]))
               SommeProf[CmtInd][CouvMin]=0
               SommeCouv[CmtInd][CouvMin]=0
           EcrireSortieCouv.write("\n")
           SommeLastBase+=int(LastBase)
       LastBase=0
       Chro=linesplit[0]
    for CmtInd in range(0,NbInd): 
       try :
          Couv=float(linesplit[(CmtInd+1)*3])
       except :
          print CmtInd, Cmt, linesplit[(CmtInd+1)*3],(CmtInd+1)*3
          sys.exit()
       for CouvMin in ListeCouvMin :
           if Couv >= CouvMin:
               SommeProf[CmtInd][CouvMin]+=Couv/1000.0
               SommeCouv[CmtInd][CouvMin]+=1 
    LastBase=linesplit[1]
    Cmt+=1

for CouvMin in ListeCouvMin :
    EcrireSortieCouv.write(Chro + "\t" +str(LastBase)+"\t"+str(CouvMin))
    for CmtInd in range(0,NbInd) :
        EcrireSortieCouv.write("\t"+str(SommeProf[CmtInd][CouvMin]) +"\t" + str(SommeCouv[CmtInd][CouvMin]))
    EcrireSortieCouv.write("\n")


