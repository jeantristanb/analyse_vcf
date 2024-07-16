# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys
import math
import Spectrum_mod

def _tajimas_d(num_sequences, avg_num_pairwise_differences, num_segregating_sites):

    ### VERIFICATION ###
    ###
    ### Given: num_sequences = 10, num_pairwise_differences = 3.888889, num_segregating_sites = 16
    ###  i.e.: tajimas_d(10, 3.888889, 16)  == -1.44617198561
    ###  Then:    a1 == 2.82896825397
    ###           a2 == 1.53976773117
    ###           b1 == 0.407407407407
    ###           b2 == 0.279012345679
    ###           c1 == 0.0539216450284
    ###           c2 == 0.0472267720013
    ###           e1 == 0.0190605338016
    ###           e2 == 0.0049489277699
    ###           D ==  -1.44617198561

    a1 = sum([1.0/i for i in range(1, num_sequences)])
    a2 = sum([1.0/(i**2) for i in range(1, num_sequences)])
    b1 = float(num_sequences+1)/(3*(num_sequences-1))
    b2 = float(2 * ( (num_sequences**2) + num_sequences + 3 )) / (9*num_sequences*(num_sequences-1))
    c1 = b1 - 1.0/a1
    c2 = b2 - float(num_sequences+2)/(a1 * num_sequences) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    D = (
        float(avg_num_pairwise_differences - (float(num_segregating_sites)/a1))
        / math.sqrt(
            (e1 * num_segregating_sites )
          + ((e2 * num_segregating_sites) * (num_segregating_sites - 1) ))
        )
    return D

def ComputeDa1(N) :
    if N<2 :
       return None
    return sum([1.0/i for i in range(1, N)])


def ComputeDa2(N) :
    if N<2 :
       return None
    return sum([1.0/(i**2) for i in range(1, N)])

## 
def PrecomputingDa1(NMax) :
    return [ComputeDa1(x) for x in range(0,(NMax+1))]

def PrecomputingDa2(NMax) :
    return [ComputeDa2(x) for x in range(0,(NMax+1))]
    
### http://www.cell.com/cms/attachment/599136/4711220/mmc1.pdf
def TajimaDWithMissingData(Liste1,Liste2,PreComputea1,PreComputea2,MaxN, MinFreqNa, MinFreqS, SumThetaPi=None, SumThetaW=None, S=None, NbPosition=None) :
    if SumThetaW==None and SumThetaPi==None and S==None:
       SumDifTheta=0
       SumThetaPi=0
       SumThetaW=0
       S=0.0
       NbPosition=0
       for Cmt in range(0, len(Liste1)):
           N=Liste1[Cmt]+Liste2[Cmt]
           FreqNa=N/float(MaxN)
           if FreqNa > MinFreqNa :
              Freq=Liste1[Cmt]/float(N)
              if Freq>0 and Freq<1.0 :
                 Pi=N/(float(N)-1.0)*2.0*(Freq*(1-Freq))
                 Wat=1/PreComputea1[int(N)]
                 SumThetaPi+=Pi
                 SumThetaW+=Wat
                 SumDifTheta+=Pi-Wat
                 S+=1.0
              NbPosition+=1
    else :
       SumDifTheta=SumThetaPi-SumThetaW
    if S< MinFreqS :
       return ["NA",S, SumThetaPi, SumThetaW,  MaxN, NbPosition, S]
    a1 = PreComputea1[MaxN]
    a2 = PreComputea2[MaxN]
    b1 = float(MaxN+1)/(3*(MaxN-1))
    b2 = float(2 * ( (MaxN**2) + MaxN + 3 )) / (9*MaxN*(MaxN-1))
    c1 = b1 - 1.0/a1
    c2 = b2 - float(MaxN+2)/(a1 * MaxN) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    Var=math.sqrt((e1 * S)+((e2 * S) * (S- 1) ))
    if Var==0 :
       return ["NA",S, SumThetaPi, SumThetaW,  MaxN, NbPosition, S]
    return [SumDifTheta/Var,S, SumThetaPi, SumThetaW, MaxN, NbPosition, S]

### http://www.cell.com/cms/attachment/599136/4711220/mmc1.pdf
def TajimaDWithMissingDataV2(Listea,Listen,PreComputea1,PreComputea2,MaxN, MinFreqNa, MinFreqS) :
    SumDifTheta=0
    SumThetaPi=0
    SumThetaW=0
    #MaxN2=float(MaxN)
    print MaxN, MinFreqNa
    S=0.0
    NbPosition=0
    for Cmt in range(0, len(Listea)):
        N=Listen[Cmt]#+Liste2[Cmt]
        FreqNa=N/float(MaxN)
        if FreqNa > MinFreqNa :
           Freq=Listea[Cmt]/float(N)
           if Freq>0 and Freq<1.0 :
              Pi=N/(float(N)-1.0)*2.0*(Freq*(1-Freq))
              Wat=1/PreComputea1[int(N)]
              SumThetaPi+=Pi
              SumThetaW+=Wat
              SumDifTheta+=Pi-Wat
              S+=1.0
           NbPosition+=1
    if S< MinFreqS :
       return ["NA",str(S), "NA", "NA",  MaxN, NbPosition, S]
    a1 = PreComputea1[MaxN]
    a2 = PreComputea2[MaxN]
    b1 = float(MaxN+1)/(3*(MaxN-1))
    b2 = float(2 * ( (MaxN**2) + MaxN + 3 )) / (9*MaxN*(MaxN-1))
    c1 = b1 - 1.0/a1
    c2 = b2 - float(MaxN+2)/(a1 * MaxN) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    Var=math.sqrt((e1 * S)+((e2 * S) * (S- 1) ))
    if Var==0 :
       return ["NA",str(S), "NA", "NA",  MaxN, NbPosition, S]
    return [SumDifTheta/Var,S, SumThetaPi, SumThetaW, MaxN, NbPosition, S]





def nCr(n,r):
    fact = math.factorial
    return fact(n) / fact(r) / fact(n-r)

def GetPiCresko(n1,n2):
    ### n1 => allele 1 
    ### n2 => allele 2
    ### formule (1 - (2 parmis n1 + 2 parmis n2))/ (2 parmis n)
    ### Population Genomics of Parallel Adaptation in Threespine Stickleback using Sequenced RAD Tags
    if n1<2  or n2 < 2 :
       return 0
    return 1- (nCr(n1,2)+nCr(n2,2))/float(nCr(n1+n2,2) )

def ComputeFstGreskoAll(a1,a2,n1,n2) :
    a1=a1*2
    a2=a2*2
    n1=n1*2
    n2=n2*2
    PiPop1=GetPiCresko(a1, n1-a1)
    PiPop2=GetPiCresko(a2, n2-a2)
    PiAll=GetPiCresko(a2+a1, (n1-a1) +(n2-a2))
    FactPop1=nCr(n1, 2)
    FactPop2=nCr(n2, 2)
    return ( 1- ( FactPop1*PiPop1 + FactPop2*PiPop2)/float((PiAll*(FactPop1+ FactPop2))), PiPop1, PiPop2) 

def GetFreq(Tab, SumInd):
   return (Tab[1]*2 + Tab[3])/(float(SumInd*2))

### Extrait de SnpStat Bioconductor
#def ComputeFstBioCond(FreqPop1, FreqPop2, FreqAll, NIndPop1, NIndPop2, NAll):  #(FreqPop,NbIndByPop,Freq2Pop, NTot, NumPop1, NumPop2):
#http://www.bioconductor.org/packages/release/bioc/vignettes/snpStats/inst/doc/Fst-vignette.pdf
def ComputeFstBioConduct(a1,a2,n1,n2) :
   NAll=n1+n2
   F1=a1/n1
   F2=a2/n2
   FAll=(a1+a2)/NAll
   Ys=FAll*(1-FAll)*(NAll/(NAll-1))
   NNMoin1Pop1=(n1/float(n1-1))
   NNMoin1Pop2=(n2/float(n2-1))
   Xs1=F1*(1-F1)*NNMoin1Pop1
   Xs2=F2*(1-F2)*NNMoin1Pop2
   Poid1=n1*NNMoin1Pop1/(n1*NNMoin1Pop1+n2*NNMoin1Pop2)
   Poid2=n2*NNMoin1Pop2/(n1*NNMoin1Pop1+n2*NNMoin1Pop2)
   Xs=Poid1*Xs1 + Poid2*Xs2
   Fs=1-(Xs/Ys)
   return Fs


#### Calcul du Fst selon reich, nature 2009
def ComputeFstReich(a1,a2, n1, n2) :
    if (a1+a2)==(n1+n2) or a1+a2==0 :
       sys.exit("(a1+a2)==(n1+n2) a1+a2==0") 
    h1=(a1*(n1-a1))/(n1*(n1-1))
    h2=(a2*(n2-a2))/(n2*(n2-1))
    N=(a1/n1 - a2/n2)**2 -h1/n1 -h2/n2
    D=N+h1+h2
    return (N,D)

        
def ComputeDivDtajima(Spectre) :
    SF=Spectrum_mod.Spectrum(Spectre)
    return SF.Tajima_D()

###    DataMap$FreqEF<-apply(DataEF,1, function(x)length(which(x==0))/length(which(x!=9)))
#    DataMap$HeEF<-DataMap$FreqEF*(1-DataMap$FreqEF)*2
#    DataMap$NEF<-apply(DataEF,1, function(x)length(which(x!=9)))
#    DataMap$FreqNF<-apply(DataNF,1, function(x)length(which(x==0))/length(which(x!=9)))
#    DataMap$HeNF<-DataMap$FreqNF*(1-DataMap$FreqNF)*2
#    DataMap$NNF<-apply(DataNF,1, function(x)length(which(x!=9)))
#    DataMap$Hs<- (DataMap$HeEF*DataMap$NEF + DataMap$HeNF*DataMap$NNF)/(DataMap$NNF+DataMap$NEF)
#    DataMap$FreqAll<-apply(cbind(DataNF,DataEF),1, function(x)length(which(x==0))/length(which(x!=9)))
#    DataMap$Ht<-2*DataMap$FreqAll*(1-DataMap$FreqAll)
#    DataMap$Fst<-(DataMap$Ht - DataMap$Hs)/(DataMap$Ht)
#    DataMap
#
#
def ComputeFstTrad(a1,a2, n1,n2) :
     F1=a1/n1     
     F2=a2/n2     
     He1=F1*(1-F1)*2
     He2=F2*(1-F2)*2
     FAll=(a1+a2)/(n1+n2)
     Hs=(He1*n1 + He2*n2)/(n1+n2)
     Ht=2*FAll*(1-FAll)
     Fst=(Ht - Hs)/Ht
     return (Fst, He1, He2)

#def ComputeDivOneLocus(Listea1,Listen1) :
    #for Cmt in range(0,len(Listen1)) :
    #   if a1<n1 and a1> 0 :
    #      F1=a1/float(n1)
    #      He1=F1*(1-F1)*2
    #      PiGre=GetPiCresko(a, n-a)
     

def ComputeAllFst(Listea1, Listea2, Listen1, Listen2, N1, N2, MinNbInd1, MinNbInd2, MinPosition) :
    SumN=0.0
    SumD=0.0
    CmtPosition=0.0
    FstTrad=0.0
    FstGe=0.0
    FstBC=0.0
    for Cmt in range(len(Listea1)) :
        if Listea1[Cmt]+Listea2[Cmt]>0 and Listen1[Cmt]>=MinNbInd1 and Listen2[Cmt]>=MinNbInd2 and (Listen1[Cmt]+Listen2[Cmt])-(Listea1[Cmt]+Listea2[Cmt])>0:
            ## Calcul du Fst de Reich
            (N,D)=ComputeFstReich(Listea1[Cmt],Listea2[Cmt], Listen1[Cmt], Listen2[Cmt])
            SumN+=N
            SumD+=D
            ### Compute le Fst Trad
            (FstTradTm,He1, He2)=ComputeFstTrad(Listea1[Cmt],Listea2[Cmt], Listen1[Cmt], Listen2[Cmt])
            FstTrad+=FstTradTm
            (FstGreTm, PiGre1Tm, PiGre2Tm)=ComputeFstGreskoAll(Listea1[Cmt],Listea2[Cmt], Listen1[Cmt], Listen2[Cmt])
            FstGe+=FstGreTm
            (FstBCTm)=ComputeFstBioConduct(Listea1[Cmt],Listea2[Cmt], Listen1[Cmt], Listen2[Cmt])
            FstBC+=FstBCTm
            CmtPosition+=1.0
    if CmtPosition > MinPosition :
       FstReich= SumN/SumD
       FstTrad=FstTrad/CmtPosition
       FstGe=FstGe/CmtPosition
       FstBC=FstBC/CmtPosition
    else :
       FstReich="NA"
       FstTrad="NA"
       FstGe="NA"
       FstBC="NA"
       SumN="NA"
       SumD="NA"
    return [CmtPosition,FstReich,SumN,SumD, FstTrad, FstGe, FstBC]

def ComputeAllFstOneLocus(a1, a2, n1, n2, N1, N2, MinNbInd1,MinNbInd2, MinPosition=None) :
    if a1+a2>0 and n1>=MinNbInd1 and n2>=MinNbInd2 and (n1+n2)-(a1+a2)>0:
       ## Calcul du Fst de Reich
       (N,D)=ComputeFstReich(a1,a2, n1, n2)
       FstReich=N/D
       ### Compute le Fst Trad
       (FstTrad,He1, He2)=ComputeFstTrad(a1,a2, n1, n2)
       (FstGre, PiGre1, PiGre2)=ComputeFstGreskoAll(a1,a2, n1, n2)
       (FstBC)=ComputeFstBioConduct(a1,a2, n1, n2)
       MinPosition[0]=1
    else :
       FstReich="NA"
       FstTrad="NA"
       FstGre="NA"
       FstBC="NA"
       N="NA"
       D="NA"
       Cmt=0
       MinPosition[0]=0
    return [FstReich,N,D, FstTrad, FstGre, FstBC]


def ComputeAllPi(Listea, Listen, N, MinNbInd, MinNbPosition) :
   SpectFreq=[0]*(N+1)
   PiGre=0.0
   CmtPosition=0
   CmtPosition2=0
   for Cmt in range(len(Listea)) :
        if Listen[Cmt]>=MinNbInd :
           PiGre+=GetPiCresko(Listea[Cmt], Listen[Cmt]-Listea[Cmt])
           CmtPosition+=1
        if Listen[Cmt]==N :
           SpectFreq[int(Listea[Cmt])] +=1
           CmtPosition2+=1
   if CmtPosition> MinNbPosition :
      PiGre=PiGre/CmtPosition
   else :
      PiGre="NA"
   if CmtPosition2> MinNbPosition :
      PiDadi=ComputeDivDtajima(SpectFreq)
   else :
      PiDadi=["NA","NA","NA","NA", "NA"]
   return (PiGre, CmtPosition, PiDadi[2], PiDadi[0] , CmtPosition2)

def ComputeAllPiOneLocus(a, n, N, MinNbInd, MinNbPosition=None) :
   SpectFreq=[0]*(N+1)
   PiGre=0.0
   if n>=MinNbInd :
      PiGre=GetPiCresko(a, n-a)
   else :
      PiGre="NA"
   return (PiGre, a, n)
