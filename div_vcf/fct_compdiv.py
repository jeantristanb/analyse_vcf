#### *** Authors : Jean-Tristan Brandenburg ***
# -*- coding: utf-8 -*-
#!/usr/bin/python
import fct_div


def ComputeAllFst(Listea1, Listea2, Listen1, Listen2, N1, N2, MinNbInd1, MinNbInd2, MinS) :
    '''
    input :
      Listea1 : list conatains for each position count of allele 1 for pop 1
      Listea2 : list contains for each position count of allele 1 for pop 2
      Listen1 : list conatains for each position count of all genotype for pop 1
      Listen2 : list conatains for each position count of all genotype for pop 2
      N1 : inital number of allele of pop 1
      N2 : inital number of allele of pop 1
      MinNbInd1 : Min nb ind for pop 1 otherwise position discarded due too much missing data
      MinNbInd2 : Min nb ind for pop 2  otherwise position discarded due too much missing data
      MinS : min position polymorphics to compute FST otherwise return NA
    output :
         vector contains values :  [S, FstReich, NReich, DReich, FstTrad, FstGeeco, FstBioconducto]
         S : position count used to computed Fst 
         Fst Reich : Fst Reich
         NReich : Value of N computed by reich
         DReich : Value of D computed by Reich
         FstTrad : Fst used 
         FstGeeco : Fst of Geeco
         FstBioconducto : Fst computed in bio conductor
    what done :  
          
    '''
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
    if CmtPosition > MinS :
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
      #PiDadi=ComputeDivDtajima(SpectFreq)
      PiDadi=["NA","NA","NA","NA", "NA"]
   else :
      PiDadi=["NA","NA","NA","NA", "NA"]
   return (PiGre, CmtPosition, PiDadi[2], PiDadi[0] , CmtPosition2)





