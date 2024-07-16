# -*- coding: utf-8 -*-
import sys
import os 
sys.path.append(os.environ['HOME']+"/Travail/bin/python")
import ComputePiHaploide
from ExtractInfoFile import GetPop, GetNomVCF
import Spectrum_mod
import gzip

def ValueGeno(a):
    if a=='.':
       return 0
    elif a=='0':
       return 1
    elif a=='1':
       return 2
    else :
       return 0

## Permet pour une position donnee, a une pop donne de recuperer les genotypes
def GetGenoPositionPopHaploide(LigneVCFSplit, Geno,ParcoursInd) :
    Geno[0]=0; Geno[1]=0;Geno[2]=0
    for Cmt in ParcoursInd :
        A1=ValueGeno(LigneVCFSplit[Cmt][0])
        A2=ValueGeno(LigneVCFSplit[Cmt][2])
        ### Cas ou on enregistre que les homozygotes
        if A1==A2 :
           Geno[A1]+=1
        else :
           Geno[0]+=1


def Ouvrir(FileVCF, IsGZ=False) :
    try :
       if  IsGZ :
         Lire=gzip.open(FileVCF)
       else :
         Lire=open(FileVCF)
    except :
        sys.exit("File : "+FileVCF+"not found")
    return Lire

### http://www.cell.com/cms/attachment/599136/4711220/mmc1.pdf
def ComputeDiv(ListeVariant1, ListeVariant2 , MinFreqBon, NInd, PreComputea1, PreComputea2):
    Taj=ComputePiHaploide.TajimaDWithMissingData(ListeVariant1,ListeVariant2,PreComputea1,PreComputea2,NInd, MinFreqBon, 4)
    return "\t".join([str(x) for x in Taj])
      

def RecupereDiv(ListeVCF,PosVcfIndPop, TailleFenetre,MinFreqBon,IsGZ, EnteteSortie, NomEnteteDtaj):
    ListeCmtPop=range(len(PosVcfIndPop))
    NbPop=len(PosVcfIndPop)
    BaliseBegin=True
    MaxNbInd=max([len(x) for x in PosVcfIndPop])
    PreComputea1=ComputePiHaploide.PrecomputingDa1(MaxNbInd)
    PreComputea2=ComputePiHaploide.PrecomputingDa2(MaxNbInd)
    ListeAll1=[[] for x in ListeCmtPop]
    ListeAll2=[[] for x in ListeCmtPop]
    EcrireDtajFen=open(EnteteSortie+".win.dtaj",'w')
    Geno=[0,0,0]
    EcrireDtajFen.write("Chro\tPosDeb\tPosFin"+NomEnteteDtaj+"\n")
    for FileVCF in ListeVCF :
       Lire=Ouvrir(FileVCF, IsGZ)
       CmtLigne=0
       FenetreGeno=[]
       PosEnregistre=[]
       for Ligne in Lire :
          SplitLigne=Ligne.split()
          if Ligne[0]!="#" and len(SplitLigne)>8 and len(SplitLigne[3])==1 and len(SplitLigne[4])==1:
             Position=int(SplitLigne[1])
             Chro=SplitLigne[0]
             if BaliseBegin :
                PositionDebFen=0#int(Position/TailleFenetre)*TailleFenetre
                PositionFinFen=PositionDebFen+TailleFenetre#PositionDebFen+TailleFenetre
                ChroEnCours=Chro
                BaliseBegin=False
                ListeVariant2=[[] for x in ListeCmtPop]
                ListeVariant1=[[] for x in ListeCmtPop]
                #ListeVariant0=[[] for x in ListeCmtPop]
             elif Position > PositionDebFen or Chro != ChroEnCours :
                Chaine=Chro+"\t"+str(PositionDebFen)+"\t"+str(PositionFinFen)
                for Pop in ListeCmtPop :
                    Chaine+="\t"+ComputeDiv(ListeVariant1[Pop], ListeVariant2[Pop] , MinFreqBon, len(PosVcfIndPop[Pop]), PreComputea1, PreComputea2)
                EcrireDtajFen.write(Chaine+"\n")
                ChroEnCours=Chro    
                PositionDebFen=PositionFinFen
                PositionFinFen=PositionDebFen+SizeWindows
                ListeVariant2=[[] for x in ListeCmtPop]
                ListeVariant1=[[] for x in ListeCmtPop]
             for CmtPop in ListeCmtPop :
                 GetGenoPositionPopHaploide(SplitLigne, Geno,PosVcfIndPop[CmtPop]) 
                 ListeVariant2[CmtPop].append(Geno[2])
                 ListeVariant1[CmtPop].append(Geno[1])
                 ListeAll1[CmtPop].append(Geno[2])
                 ListeAll2[CmtPop].append(Geno[1])
    Chaine=Chro+"\t"+str(PositionDebFen)+"\t"+str(PositionFinFen)
    for Pop in ListeCmtPop :
        Chaine+="\t"+ComputeDiv(ListeVariant1[Pop], ListeVariant2[Pop] , MinFreqBon, len(PosVcfIndPop[Pop]), PreComputea1, PreComputea2)
    EcrireDtajFen.write(Chaine+"\n")
    EcrireDtajFen.close()
    EcrireDtaAll=open(EnteteSortie+".all.dtaj",'w')
    EcrireDtaAll.write(NomEnteteDtaj[1::]+"\n")
    Chaine=""
    for Pop in ListeCmtPop :
        Chaine+="\t"+ComputeDiv(ListeAll1[Pop], ListeAll2[Pop] , MinFreqBon, len(PosVcfIndPop[Pop]), PreComputea1, PreComputea2)
    EcrireDtaAll.write(Chaine+"\n")
    EcrireDtaAll.close()
     

ListeFileVCF=["/home/jeantristan/Travail/Data/DataTrip/SnpNoGeneWithTripAA.vcf.gz"]
#ListeFileVCF=["/home/jeantristan/Travail/Data/AllDataWithGene/SnpAll.vcf.gzip"]
FilePop="../NomInd-Conc"
EstGZ=True
SizeWindows=100000
MinFreqBon=0.7
EnteteSortie="TripWindow"


### Lecture des Ind et de leur pop
DicInfoPop=GetPop(FilePop)


LireVCF=Ouvrir(ListeFileVCF[0], EstGZ)
### on recupere les positions des noms des VCF en fonctio de ceux des individus
(InfoNomPosPop,PosExterne)=GetNomVCF(LireVCF, DicInfoPop)
###
LireVCF.close()

NomPop=DicInfoPop['PopIndic']
ChaineEntete=""
NomEnteteDtaj=['Dtj', 'S', 'pi', 'TetWat', 'n', 'NbPos', 'NbPoly']
for Pop in range(0,len(NomPop)):
    ##[(pihat - theta)/C, S, pihat, theta, n]
    for Entete in NomEnteteDtaj :
       ChaineEntete+= "\t"+Entete+"_"+NomPop[Pop]


RecupereDiv(ListeFileVCF,InfoNomPosPop,SizeWindows,MinFreqBon,EstGZ, EnteteSortie, ChaineEntete)









