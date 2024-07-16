# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys
import math
from lib_diversity import *
#sys.path.append("/home/gear/jeantristan/bin/python")
#import ComputePiHaploide


def nCr(n,r):
    fact = math.factorial
    return fact(n) / fact(r) / fact(n-r)

def GetPi(n1,n2):
    ### n1 => allele 1 
    ### n2 => allele 2
    ### formule (1 - (2 parmis n1 + 2 parmis n2))/ (2 parmis n)
    ### Population Genomics of Parallel Adaptation in Threespine Stickleback using Sequenced RAD Tags
    if n1 < 2 or n2 < 2 :
       return 0
    return 1- (nCr(n1,2)+nCr(n2,2))/float(nCr(n1+n2,2) )

### nj : nombre d'allele
def GetFst_Cresko(NAllelePop1, NAllelePop2, PiPop1, PiPop2, PiAll):
    ### Population Genomics of Parallel Adaptation inThreespine Stickleback using Sequenced RAD Tags
    ### formul :
    ###   Fst = 1- (sommepop( 2 parmis nj * pii_j)/pi_a * somme( 2 parmis nj)
    ### nj is the number of alleles sampled in population j,
    FactPop1=nCr(NAllelePop1, 2)
    FactPop2=nCr(NAllelePop2, 2)
    return 1- ( FactPop1*PiPop1 + FactPop2*PiPop2)/float((PiAll*(FactPop1+ FactPop2)))


## NomInd NomPop
### on recupere les informations par population
def GetPop(FilePop):
   LirePop=open(FilePop) 
   DicInd={}
   DicPop={}
   DicPopByIndic={}
   CmtPop=0
   for Ligne in LirePop :
       SplitLigne=Ligne.split()
       Pop=SplitLigne[1].replace('\n','')
       if Pop not in DicPop :
          DicPop[Pop]=CmtPop
          DicPopByIndic[CmtPop]=Pop
          CmtPop+=1
       #DicInd[SplitLigne[0].upper()]=Pop
       DicInd[SplitLigne[0]]=Pop
   return {'Pop':DicPop, 'Ind':DicInd, 'PopIndic':DicPopByIndic}

### on recupere les noms dans le VCF
def GetNomVCF(LireVCF,DicInfoPop ) :
    Balise=True
    InfoInd=DicInfoPop['Ind'] 
    InfoPop=DicInfoPop['Pop']
    IndPop=[[] for x in range(len(InfoPop))]
    while(Balise):
       Ligne=LireVCF.readline() 
       if Ligne[0:6]=="#CHROM":
          SplitLigne=Ligne.split()
          Cmt =9
          for Ind in SplitLigne[9::]:
              Ind=Ind.replace('\n','')
              #IndPop[InfoPop[InfoInd[Ind.upper()]]].append(Cmt)
              if Ind in InfoInd:
                 IndPop[InfoPop[InfoInd[Ind]]].append(Cmt)
              Cmt+=1 
          Balise=False
    return IndPop

def GetGeno(Char):
    A1=Char[0]
    A2=Char[2]
    if A1 == '.':
       return 0
    if A1==A2:
       if A1=='0':
          return 1
       return 2
    return 3
def NettoyerVec2D(Vecteur) :
    for Vec in Vecteur :
        for Cmt in range(len(Vec)):
           Vec[Cmt]=0

def NettoyerVec(Vec):
    Cmt=0
    NbElem=len(Vec)
    while Cmt<NbElem:
        Vec[Cmt]=0
        Cmt+=1

def GetP(Tab, SumInd):
   return (Tab[1]*2 + Tab[3])/(float(SumInd*2))


def GetP(HAA,HET,SumInd):
   if SumInd==0 :
      return -1
   return (HAA*2 + HET)/(float(SumInd*2))

def GetHetTheo(p) :
    if p < 0 :
       return -1 
    return 2*p*(1-p)

def GetDPInd(InfoInd):
    TmpSplit=InfoInd.split(':')
    if len(TmpSplit)<14 :
       return -1
    return int(TmpSplit[4])+int(TmpSplit[5])

def GetMeanDPByPop(InfoInd, InfoNomPosPop, VecteurNbPop):
    MeanDP=[0]*len(InfoNomPosPop)
    MeanDPHet=[0]*len(InfoNomPosPop)
    MeanDPHom=[0]*len(InfoNomPosPop)
    SumIndDP=[0]*len(InfoNomPosPop)
    SumIndDPHet=[0]*len(InfoNomPosPop)
    for Pop in VecteurNbPop:
        for Ind in InfoNomPosPop[Pop] :
            DPInd=GetDPInd(InfoInd[Ind])
            if DPInd>0:
               MeanDP[Pop]+=DPInd
               SumIndDP[Pop]+=1
               Geno=GetGeno(InfoInd[Ind])
               if Geno==3 :
                  MeanDPHet[Pop]+=DPInd
                  SumIndDPHet[Pop]+=1
               elif Geno>0:
                  MeanDPHom[Pop]+=DPInd
        if SumIndDP[Pop]>0:
           MeanDP[Pop]/=float(SumIndDP[Pop])
           if SumIndDPHet[Pop]>0 :
              MeanDPHet[Pop]/=float(SumIndDPHet[Pop])
           else :
              MeanDPHet[Pop]=-1
           if SumIndDP[Pop] - SumIndDPHet[Pop]> 0 :
              MeanDPHom[Pop]/=float(SumIndDP[Pop]-SumIndDPHet[Pop])
           else :
              MeanDPHom[Pop]=-1
        else :
          MeanDP[Pop]=-1
          MeanDPHom[Pop]=-1
          MeanDPHet[Pop]=-1
    return [MeanDP, MeanDPHet, MeanDPHom]

        
### Extrait de SnpStat
def ComputeFst(FreqPop1, FreqPop2, FreqAll, NIndPop1, NIndPop2, NAll):  #(FreqPop,NbIndByPop,Freq2Pop, NTot, NumPop1, NumPop2):
   Ys=FreqAll*(1-FreqAll)*(NAll/(NAll-1))
   NNMoin1Pop1=(NIndPop1/float(NIndPop1-1))
   NNMoin1Pop2=(NIndPop2/float(NIndPop2-1))
   Xs1=FreqPop1*(1-FreqPop1)*NNMoin1Pop1
   Xs2=FreqPop2*(1-FreqPop2)*NNMoin1Pop2
   Poid1=NNMoin1Pop1/(NNMoin1Pop1+NNMoin1Pop2)
   Poid2=NNMoin1Pop2/(NNMoin1Pop1+NNMoin1Pop2)
   Xs=Poid1*Xs1 + Poid2*Xs2
   Fs=1-(Xs/Ys)
   return (Fs,Ys)


### Lecture du VCF
def ExtractStat(LireVCF,InfoNomPosPop, NomPop, FileSortie):
    ## Genotype (0 : None, 1 : HomRef, 2 : HomAlt, 3 HEt for each populaion 
    StatGeno=[[0]*4 for x in range(len(InfoNomPosPop))]
    FstPop=[[0]*len(InfoNomPosPop) for x in range(len(InfoNomPosPop))]
    FstPopCresko=[[0]*len(InfoNomPosPop) for x in range(len(InfoNomPosPop))]
    PoidFst=[[0]*len(InfoNomPosPop) for x in range(len(InfoNomPosPop))]
    ## Genotype (0 : None, 1 : HomRef, 2 : HomAlt, 3 HEt for all individu
    StatAll=[0]*4
    EcrireSortie=open(FileSortie, 'w')
    ## p by population
    FreqPop=[0]*len(InfoNomPosPop)
    NbIndByPop=[0]*len(InfoNomPosPop)
    HoPop=[0]*len(InfoNomPosPop)
    HsPop=[0]*len(InfoNomPosPop)
    FisPop=[0]*len(InfoNomPosPop)
    NbIndNonCouvert=[0]*len(InfoNomPosPop)
    PiPop=[0]*len(InfoNomPosPop)
    ## Nb Ind for all population
    ## Vector kept individual where number individual as information
    NbPop=len(InfoNomPosPop)
    EcrireInfoInd={}
    Freq=[0]*NbPop
    VecteurNbPop=range(NbPop)
    #for Pop in VecteurNbPop:
    #    EcrireInfoInd[Pop]=open(FileSortie+"_"+NomPop[Pop]+".out", 'w')
    ChaineSortie="Chro\tPos\tNbPop\tNbIndTot\tPt\tHoT\tHT\tHSMean\tHoMean\tNbHet\tNbHom\t"
    Cmt=1
    CharInfoInd=[[""]*len(x) for x in InfoNomPosPop]
    for Pop in VecteurNbPop:
        ChaineSortie+="\t"+"NbInd_"+NomPop[Pop]+"\tP_"+NomPop[Pop]+"\tHo_"+NomPop[Pop]+"\tHs_"+NomPop[Pop]+"\tDP_"+NomPop[Pop] + "\tDPHet_"+NomPop[Pop]+"\tDPHom_"+NomPop[Pop]+"\tPi_"+NomPop[Pop]
        #ChaineSortie+="\t"+"NbInd_"+NomPop[Pop]+"\tP_"+NomPop[Pop]+"\tHo_"+NomPop[Pop]+"\tHs_"+NomPop[Pop]#+"\tDP_"+NomPop[Pop] + "\tDPHet_"+NomPop[Pop]+"\tDPHom_"+NomPop[Pop]
    for Pop in VecteurNbPop :
           for Pop2 in VecteurNbPop[(Pop+1)::]:
               ChaineSortie+="\t"+"Fst_"+NomPop[Pop]+"_"+NomPop[Pop2]+"\t"+"FstW_"+NomPop[Pop]+"_"+NomPop[Pop2]
    for Pop in VecteurNbPop :
           for Pop2 in VecteurNbPop[(Pop+1)::]:
               ChaineSortie+="\t"+"FstCres_"+NomPop[Pop]+"_"+NomPop[Pop2]
    ChaineSortie+="\n"
    for Ligne in LireVCF:
       SplitLigne=Ligne.split()
       NbPopWithInd=0
       HoSum=0.0
       HsSum=0.0
       NbPopWithPop=0
       for Pop in VecteurNbPop:
           ### on recupere les Geno pour chaque populations
           for PosInd in InfoNomPosPop[Pop]:
               Geno=GetGeno(SplitLigne[PosInd])
               StatGeno[Pop][Geno]+=1
           ## Nombre d'individu par pop
           NbIndByPop[Pop]=sum(StatGeno[Pop][1:4])
           NbIndNonCouvert[Pop]=StatGeno[Pop][0]
           #if NbIndByPop[Pop]> 0:
           if NbIndNonCouvert[Pop]==0 :
              PiPop[Pop]=GetPi(2*StatGeno[Pop][1]+StatGeno[Pop][3],2*StatGeno[Pop][2]+StatGeno[Pop][3] )
              FreqPop[Pop]=GetP(StatGeno[Pop][1], StatGeno[Pop][3], NbIndByPop[Pop]) 
              HoPop[Pop]=StatGeno[Pop][3]/float(NbIndByPop[Pop])
              HsPop[Pop]=GetHetTheo(FreqPop[Pop])
              if HsPop[Pop]==0 or HsPop[Pop]==1 :
                  FisPop[Pop]=2
              else :
                 FisPop[Pop]=(HsPop[Pop]-HoPop[Pop])/HsPop[Pop]
              HsSum+=HsPop[Pop]
              HoSum+=HoPop[Pop]
              NbPopWithPop+=1
              #if FisPop[Pop] < 0.2 :
              #   for PosInd in InfoNomPosPop[Pop]:
              #       EcrireInfoInd[Pop].write(SplitLigne[PosInd]+"\t")
              #   EcrireInfoInd[Pop].write("\n")
              for CmtStat in range(4):
                  StatAll[CmtStat]+= StatGeno[Pop][CmtStat]
           else :
              FreqPop[Pop]=-1
              HoPop[Pop]=-1
              HsPop[Pop]=-1
              PiPop[Pop]=-1
       ## On calcule la diversite : Fst 
       LimCalcul=-10**3
       for Pop in VecteurNbPop :
           for Pop2 in VecteurNbPop[(Pop+1)::]:
               NbInd0=min(StatGeno[Pop][1]+StatGeno[Pop2][1], StatGeno[Pop][2]+StatGeno[Pop2][2])
               if NbIndNonCouvert[Pop]==0 and NbIndNonCouvert[Pop2]==0:
                     ### Calcul des frequences de 2 Pop
                     #print Freq2Pop
                     #if StatGeno[Pop2][3]==0 and StatGeno[Pop][3]==0 and NbInd0>1 :
                     if NbInd0>1 :
                        Freq2Pop=GetP(StatGeno[Pop][1]+StatGeno[Pop2][1], StatGeno[Pop][3]+StatGeno[Pop2][3], NbIndByPop[Pop]+NbIndByPop[Pop2]) 
                        Pi2Pop=GetPi(StatGeno[Pop][1]*2+StatGeno[Pop2][1]*2 + StatGeno[Pop][3]+StatGeno[Pop2][3],  StatGeno[Pop][2]*2+StatGeno[Pop2][2]*2 + StatGeno[Pop][3]+StatGeno[Pop2][3])
                        HT2Pops=GetHetTheo(Freq2Pop)
                        ### New metho to compute Fst
                        #def ComputeFst(FreqPop1, FreqPop2, FreqAll, NIndPop1, NIndPop2, NAll)  
                        (FstPop[Pop][Pop2],PoidFst[Pop][Pop2])=ComputeFst(FreqPop[Pop],FreqPop[Pop2],Freq2Pop,NbIndByPop[Pop], NbIndByPop[Pop2], NbIndByPop[Pop]+NbIndByPop[Pop2])
                        #def GetFst_Cresko(NAllelePop1, NAllelePop2, PiPop1, PiPop2, PiAll):
                        if PiPop[Pop]!=-1 and PiPop[Pop2]!=-1:
                           FstPopCresko[Pop][Pop2]=GetFst_Cresko(StatGeno[Pop][1]*2 + StatGeno[Pop][3] + StatGeno[Pop][2]*2,StatGeno[Pop2][1]*2 + StatGeno[Pop2][3] + StatGeno[Pop2][2]*2,PiPop[Pop],PiPop[Pop2] ,Pi2Pop)
                        else :
                           FstPopCresko[Pop][Pop2]=-1
                     else :
                        FstPop[Pop][Pop2]=-1
                        PoidFst[Pop][Pop2]=-1
                        FstPopCresko[Pop][Pop2]=-1
               else :
                   PoidFst[Pop][Pop2]=-1
                   FstPop[Pop][Pop2]=-1
                   FstPopCresko[Pop][Pop2]=-1
                
       MeanDPPop=GetMeanDPByPop(SplitLigne, InfoNomPosPop, VecteurNbPop)
       ### On calcule les stat
       ## Heterozygotie tot pop
       NbIndTot=sum(StatAll[1:4])
       PT=GetP(StatAll[1],StatAll[3], NbIndTot)
       HT=GetHetTheo(PT)
       #print NbIndTot
       if NbIndTot>0:
          HoT=StatAll[3]/float(NbIndTot)
       else :
          HoT=-1
       if NbPopWithPop==0 :
          MeanHs=-1
          MeanHo=-1
       else :
          MeanHs=HsSum/NbPopWithPop
          MeanHo=HoSum/NbPopWithPop
       ## Stat pour toutes les populations
       ChaineSortie+=SplitLigne[0]+"\t"+SplitLigne[1]+"\t"+str(NbPopWithPop)+"\t"+str(NbIndTot)+"\t"+str(PT)+"\t"+str(HoT)+"\t"+str(HT)+"\t"+str(MeanHs)+"\t"+str(MeanHo) +"\t"+str(StatAll[3])+"\t"+str(StatAll[2])
       for Pop in VecteurNbPop:
           ChaineSortie+="\t"+str(NbIndByPop[Pop])+"\t"+str(FreqPop[Pop])+"\t"+str(HoPop[Pop])+"\t"+str(HsPop[Pop]) +"\t"+str(MeanDPPop[0][Pop])+"\t"+str(MeanDPPop[1][Pop])+"\t"+str(MeanDPPop[2][Pop])+"\t"+str(PiPop[Pop])
       for Pop in VecteurNbPop :
           for Pop2 in VecteurNbPop[(Pop+1)::]:
               ChaineSortie+="\t"+str(FstPop[Pop][Pop2])+"\t"+str(PoidFst[Pop][Pop2])
       for Pop in VecteurNbPop :
           for Pop2 in VecteurNbPop[(Pop+1)::]:
               ChaineSortie+="\t"+str(FstPopCresko[Pop][Pop2])
       ChaineSortie+="\n"
       ## To clean vector (no new allocation)
       NettoyerVec(StatAll)
       NettoyerVec(FisPop)
       NettoyerVec(FreqPop)
       NettoyerVec(NbIndByPop)
       NettoyerVec(HoPop)
       NettoyerVec(HsPop)
       NettoyerVec(NbIndNonCouvert)
       NettoyerVec(PiPop)
       ## 
       NettoyerVec2D(StatGeno)
       NettoyerVec2D(FstPop)
       NettoyerVec2D(FstPopCresko)
       NettoyerVec2D(PoidFst)
       if Cmt%10000==0:
          EcrireSortie.write(ChaineSortie)
          ChaineSortie="\n"
       Cmt+=1
    EcrireSortie.write(ChaineSortie)

FilePop="NomIndPopMin.out"
#Snp_Het_AvecFiltreV2-Binom-50000-NewDefHet.vcf
#FileVCF="Snp_Het_AvecFiltreV2-Binom-50000-NewDefHet.vcf"
#FileSortie="Sortie-NewDefHet.out"
if len(sys.argv )!=4 :
   sys.exit("exe FileVCF FilePop FileSortie")

FileVCF=sys.argv[1]
FilePop=sys.argv[2]
FileSortie=sys.argv[3]



DicInfoPop=GetPop(FilePop)
LireVCF=open(FileVCF)
InfoNomPosPop=GetNomVCF(LireVCF, DicInfoPop)
ExtractStat(LireVCF,InfoNomPosPop,  DicInfoPop['PopIndic'], FileSortie)

