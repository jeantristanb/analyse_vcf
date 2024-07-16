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
      
def ComputeFstReichFen(ListeVariant1Pop1, ListeVariant2Pop1, N1,ListeVariant1Pop2, ListeVariant2Pop2, N2, MinFreqBon) :
    Nf=0.0
    Df=0.0
    N1=float(N1)
    N2=float(N2)
    CmtLocus=0
    for CmtVariant in range(0,len(ListeVariant1Pop1)):
       n1=float(ListeVariant1Pop1[CmtVariant]+ListeVariant2Pop1[CmtVariant])
       n2=float(ListeVariant1Pop2[CmtVariant]+ListeVariant2Pop2[CmtVariant])
       a1=float(ListeVariant1Pop1[CmtVariant])
       a2=float(ListeVariant1Pop2[CmtVariant])
       if n1/N1 > MinFreqBon and n2/N2 > MinFreqBon and a1+a2!=0 and a1+a2!=n1+n2:
          (N,D)=ComputePiHaploide.ComputeFstReich(a1,a2, n1, n2)
          Nf+=N 
          Df+=D 
          CmtLocus+=1
    #print Nf, Df, a1,a2,n1, n2
    if CmtLocus==0 :
       return ["NA", 0, 0, CmtLocus]
    return [Nf/Df, Nf, Df, CmtLocus]
    

def RecupereDiv(ListeVCF,PosVcfIndPop, TailleFenetre,MinFreqBon,IsGZ, EnteteSortie, NomEnteteDtaj, ChaineEnteteFst):
    ListeCmtPop=range(len(PosVcfIndPop))
    NbPop=len(PosVcfIndPop)
    BaliseBegin=True
    MaxNbInd=max([len(x) for x in PosVcfIndPop])
    PreComputea1=ComputePiHaploide.PrecomputingDa1(MaxNbInd)
    PreComputea2=ComputePiHaploide.PrecomputingDa2(MaxNbInd)
    ListeAll1=[[] for x in ListeCmtPop]
    ListeAll2=[[] for x in ListeCmtPop]
    EcrireDtajFen=open(EnteteSortie+".win.dtaj",'w')
    EcrireDtajFen.write("Chro\tPosDeb\tPosFin"+NomEnteteDtaj+"\n")
    EcrireFstFen=open(EnteteSortie+".win.fst",'w')
    EcrireFstFen.write("Chro\tPosDeb\tPosFin"+ChaineEnteteFst+"\n")
    Geno=[0,0,0]
    AllFst=[[] for x in ListeCmtPop]
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
                Chaine=Chro+"\t"+str(PositionDebFen)+"\t"+str(PositionFinFen)
                CmtPop=0
                while CmtPop<NbPop-1 :
                      CmtPop2=CmtPop+1
                      while CmtPop2<NbPop :
                            Tmp=ComputeFstReichFen(ListeVariant1[CmtPop], ListeVariant2[CmtPop], len(PosVcfIndPop[CmtPop]),ListeVariant1[CmtPop2], ListeVariant2[CmtPop2], len(PosVcfIndPop[CmtPop2]), MinFreqBon)
                            Chaine+="\t"+"\t".join([str(x) for x in Tmp])
                            CmtPop2+=1
                      CmtPop+=1
                EcrireFstFen.write(Chaine+"\n")
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
    Chaine=Chro+"\t"+str(PositionDebFen)+"\t"+str(PositionFinFen)
    CmtPop=0
    while CmtPop<NbPop-1 :
          CmtPop2=CmtPop+1
          while CmtPop2<NbPop :
                Tmp=ComputeFstReichFen(ListeVariant1[CmtPop], ListeVariant2[CmtPop], len(PosVcfIndPop[CmtPop]),ListeVariant1[CmtPop2], ListeVariant2[CmtPop2], len(PosVcfIndPop[CmtPop2]), MinFreqBon)
                Chaine+="\t"+"\t".join([str(x) for x in Tmp])
                CmtPop2+=1
          CmtPop+=1
    EcrireFstFen.write(Chaine+"\n")
    EcrireFstFen.close()
    EcrireDtaAll=open(EnteteSortie+".all.dtaj",'w')
    EcrireDtaAll.write("Type"+NomEnteteDtaj+"\n")
    Chaine="All"
    for Pop in ListeCmtPop :
        Chaine+="\t"+ComputeDiv(ListeAll1[Pop], ListeAll2[Pop] , MinFreqBon, len(PosVcfIndPop[Pop]), PreComputea1, PreComputea2)
    EcrireDtaAll.write(Chaine+"\n")
    EcrireDtaAll.close()
    EcrireFstAll=open(EnteteSortie+".all.fst",'w')
    EcrireFstAll.write("Type"+ChaineEnteteFst+"\n")
    Chaine="All"
    CmtPop=0
    while CmtPop<NbPop-1 :
          CmtPop2=CmtPop+1
          while CmtPop2<NbPop :
                Tmp=ComputeFstReichFen(ListeAll1[CmtPop], ListeAll2[CmtPop], len(PosVcfIndPop[CmtPop]),ListeAll1[CmtPop2], ListeAll2[CmtPop2], len(PosVcfIndPop[CmtPop2]), MinFreqBon)
                Chaine+="\t"+"\t".join([str(x) for x in Tmp])
                CmtPop2+=1
          CmtPop+=1
    EcrireFstAll.write(Chaine)

    

     

#ListeFileVCF=["/home/jeantristan/Travail/Data/DataTrip/SnpNoGeneWithTripAA.vcf.gz"]
#ListeFileVCF=["/home/jeantristan/Travail/Data/AllDataWithGene/SnpAll.vcf.gzip"]
#FilePop="../NomInd-Conc"
#EstGZ=True
#SizeWindows=100000
#MinFreqBon=0.7
#EnteteSortie="TripWindow"
if len(sys.argv)<7 :
   sys.exit("exe FilePop EstGZ[T/F] SizeWindows MinFreqAvailable Entete ListeFileVCF\n"+"\t".join(sys.argv))
FilePop=sys.argv[1]
EstGZ=sys.argv[2]
if EstGZ=="F" :
   EstGZ=False
elif EstGZ=="T":
   EstGZ=True
else: 
   sys.exit("EstGZ : " + EstGZ + "not T or F")

SizeWindows=int(sys.argv[3])

MinFreqBon=float(sys.argv[4])
EnteteSortie=sys.argv[5]
ListeFileVCF=sys.argv[6::]


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

NbPop=len(NomPop)
ChaineEnteteFst=""
NomEnteteDtaj=['Fst', 'N', 'D', 'NbPoly']
CmtPop=0
while CmtPop<NbPop-1 :
      CmtPop2=CmtPop+1
      while CmtPop2<NbPop :
         for Entete in NomEnteteDtaj :
             ChaineEnteteFst+="\t"+Entete+"_"+NomPop[CmtPop]+"_"+NomPop[CmtPop2]
         CmtPop2+=1
      CmtPop+=1

RecupereDiv(ListeFileVCF,InfoNomPosPop,SizeWindows,MinFreqBon,EstGZ, EnteteSortie, ChaineEntete, ChaineEnteteFst)









