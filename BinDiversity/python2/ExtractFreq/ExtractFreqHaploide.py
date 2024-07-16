#-*- coding: utf-8 -*-
import sys
import os 
sys.path.append(os.environ['HOME']+"/bin/python2")
import ComputePiHaploide
from ExtractInfoFile import GetPop, GetNomVCF, GetGenoPositionPopHaploide

def GetMatrice0(NbLigne, NbColonne) :
    return [ [ 0 for i in range(NbColonne) ] for j in range(NbLigne) ] 

def GetMatriceVec(NbLigne, NbColonne) :
    return [ [ 0 for i in range(NbColonne) ] for j in range(NbLigne) ] 

## P sVCFIndPop contient la position de chacun des individu pour chacune des populations [[Ind1,ind3, ind5], [ind2, ind4]]
def ExtractStatPiD(ListeNomFile, Entete, PosVcfIndPop, MinFreqInd, NomPop, PosAncetre=None):
   NbPop=len(PosVcfIndPop)
   ListeCmtPop=range(len(PosVcfIndPop))
   CmtFenetre=0
   ListeFreq=[[0]*(len(PosVcfIndPop[x])+1) for x in range(NbPop)]
   ListeEcrire={}
   ListeFreq2D=GetMatriceVec(NbPop, NbPop)
   for CmtPop in range(0,NbPop) :
       for Cmt2Pop in range(0,NbPop) :
           ListeFreq2D[CmtPop][Cmt2Pop]=GetMatrice0((len(PosVcfIndPop[CmtPop])+1), (len(PosVcfIndPop[Cmt2Pop])+1))
   for CmtPop in range(0,len(NomPop)) :
      ListeEcrire[CmtPop]=open(Entete+"."+NomPop[CmtPop]+".out",'w') 
   for File in ListeNomFile :
       Lire=open(File)
       for Ligne in Lire :
           if Ligne[0]!="#" and len(Ligne)>5 :
              SplitLigne=Ligne.split()
              if len(SplitLigne[3])==1 and len(SplitLigne[4])==1:
                 ListeNbInd=[]
                 for CmtPop in ListeCmtPop :
                     TempPop=GetGenoPositionPopHaploide(SplitLigne, PosVcfIndPop[CmtPop], PosAncetre)
                     NbInd=sum(TempPop)
                     NbIndBon=sum(TempPop[1::])
                     NbAllele1=TempPop[2]
                     if NbInd>0 and NbIndBon/float(NbInd)>MinFreqInd :
                        ListeFreq[CmtPop][NbAllele1]+=1          
                        ListeNbInd.append(NbAllele1)
                     else :
                        ListeNbInd.append(-1)
                 for CmtPop in range(0,NbPop-1) :
                     for Cmt2Pop in range(CmtPop+1,NbPop) :
                         if ListeNbInd[CmtPop]!=-1 and ListeNbInd[Cmt2Pop]!=-1 :
                            try :
                               ListeFreq2D[CmtPop][Cmt2Pop][ListeNbInd[CmtPop]][ListeNbInd[Cmt2Pop]]+=1 
                            except :
                               print CmtPop, Cmt2Pop
                               print ListeFreq2D[CmtPop][Cmt2Pop]
                               print len(ListeFreq2D),  len(ListeFreq2D[0]), "Freq2"
                               print CmtPop, Cmt2Pop, ListeNbInd[CmtPop], ListeNbInd[Cmt2Pop], len(ListeFreq2D[CmtPop][Cmt2Pop]), len(ListeFreq2D[CmtPop][Cmt2Pop][0])
                               sys.exit()
       Lire.close()
   for CmtPop in range(0,NbPop) :
        ListeEcrire[CmtPop].write("\n".join([str(x) for x in ListeFreq[CmtPop]]))
 
   for CmtPop in range(0,NbPop-1) :
       for Cmt2Pop in range(CmtPop+1,NbPop) :
           Ecrire=open(Entete+"2d."+NomPop[CmtPop]+"."+NomPop[Cmt2Pop]+".out",'w')       
           ##
           for X in range(0,len(ListeFreq2D[CmtPop][Cmt2Pop])):
               Ecrire.write(" ".join([str(x) for x in ListeFreq2D[CmtPop][Cmt2Pop][X]])+"\n")
           Ecrire.close()
#SF=Spectrum_mod.Spectrum([10,1,2,3])
#SF.Tajima_D()
#SF.S()
#SF.
Cmd="Exe.py NomIndResume Sortie MinFreq FileVCF"
if(len(sys.argv)<5):
  sys.exit(Cmd)


MinFreqInd=float(sys.argv[3])
FilePop=sys.argv[1]
### on recupere les populations
DicInfoPop=GetPop(FilePop)
ListeNomVCF=sys.argv[4::]
print ListeNomVCF
LireVCF=open(ListeNomVCF[0])
### on recupere dans les entetes les positions des individus
[InfoNomPosPop,PosOutGroup]=GetNomVCF(LireVCF, DicInfoPop)
## On verifie si il y a un groupe Externe
if PosOutGroup==None and DicInfoPop['OutGroup']!=None :
   sys.exit('Ind Outgroup : '+ DicInfoPop['OutGroup'] +' not find in the vcf : ' + ListeNomVCF[0])
if PosOutGroup!=None :
   print "OutGroup defined"

VecteurNbPop=range(len(InfoNomPosPop))
NomPop=DicInfoPop['PopIndic']
LireVCF.close()


ExtractStatPiD(ListeNomVCF,sys.argv[2],InfoNomPosPop, MinFreqInd, NomPop, PosOutGroup)

