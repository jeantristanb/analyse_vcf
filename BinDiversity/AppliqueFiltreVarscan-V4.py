#!/usr/bin/python
# -*- coding: utf8 -*-
## 25 fevrier : developpement du test de KS
import sys
from scipy import stats
#from fisher import pvalue
from scipy.stats import kstest
import scipy
from Binom import Binom





## On peut s'en sortir avec la moyenne
#global GlobalFreqA=0.265538268593245## 0
#global GlobalFreqT=0.265537395432029## 1
#global GlobalFreqC=0.23456418978987## 2
#global GlobalFreqG=0.234360146184855## 3
#global GlobalVecFreq=[0.265538268593245, 0.265537395432029, 0.23456418978987, 0.234360146184855, 0]

def GetMatrice(Nrow, Ncol):
    M = [[0 for j in range(0,Ncol)] for i in range(0,Nrow)]
    return M

def GetStatKS(Vector, lim=20) :
   VectorA=scipy.array([x for x in Vector if x>=0])
   if len(VectorA) < lim :
       return True
   SumVec=sum(VectorA)
   if SumVec==0:
      return True
   MeanVec=SumVec/float(len(VectorA))
   test_stat= kstest(VectorA, 'poisson', args=(MeanVec,), alternative = 'greater')  
   return test_stat[1]
    



###Ouverture d'un fichier
def OuvrirFichier(Fichier, Type='r'):
    try :
       if Type =='r' :
           LireFich=open(Fichier)
       else :
           LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " open in mode " + Type +" can't open\n")
    return LireFich

### Ces calculs ne sont pas genials puisque l'erreur a 1/4 chance d'etre PEAlt
def GetPvalueBinomHom(PE, PEI, BHom):
  return BHom.Test(int(PEI-PE),int(PEI))

def GetPvalueBinomHet(PERef,PEAlt, PEI, BHet):
   return  BHet.Test(int(PEI-max(PERef,PEAlt)), int(PEI)) 

#def GetPvalueBinomHomLynchModif(NbPEI, PERef,PEAlt, p=0.015):
#   return  stats.binom_test(NbPEI-max(PERef,PEAlt), NbPEI, p=p) 

#def GetPvalueBinomHomLynchModif(NbPEI, PERef,PEAlt, p=0.015):
#   SommePEAltRef=(PERef+PEAlt)
#   PError=stats.binom_test(NbPEI-(PERef+PEAlt) , NbPEI, p=2.0/3.0*p)
   #return  PError*stats.binom_test( PERef ,SommePEAlt, p=0.5) +  PError*stats.binom_test(PEAlt,

def GetPvalueFisherExactRef(PERef,PEAlt, ProbRef=1, ProbAlt=0):
   n=PERef+PEAlt
   PeTheoAlt=float(ProbAlt*n)
   PeTheoRef=n-PeTheoAlt
   return pvalue(PeTheoRef, PeTheoAlt, PERef, PEAlt).right_tail

def GetPvalueFisherExactHet(PERef,PEAlt, ProbRef=0.5, ProbAlt=0.5) :
   n=PERef+PEAlt
   PeTheoAlt=float(ProbAlt*n)
   PeTheoRef=n-PeTheoAlt
   return pvalue(PeTheoRef, PeTheoAlt, PERef, PEAlt).two_tail

def GetPvalueFisherExactRef2(PERef,PEAlt, ProbRef=1, ProbAlt=0):
   n=PERef+PEAlt
   PeTheoAlt=float(ProbAlt*n)
   PeTheoRef=n-PeTheoAlt
   return pvalue(PeTheoRef, PeTheoAlt, PERef, PEAlt).two_tail

def GetPvalueFisherExactAlt(PERef,PEAlt, ProbRef=0, ProbAlt=1) :
   n=PERef+PEAlt
   PeTheoAlt=float(ProbAlt*n)
   PeTheoRef=n-PeTheoAlt
   return pvalue(PeTheoRef, PeTheoAlt, PERef, PEAlt).left_tail

def GetPvalueFisherExactAlt2(PERef,PEAlt, ProbRef=0, ProbAlt=1) :
   n=PERef+PEAlt
   PeTheoAlt=float(ProbAlt*n)
   PeTheoRef=n-PeTheoAlt
   return pvalue(PeTheoRef, PeTheoAlt, PERef, PEAlt).two_tail


### Calu de l'erreur en fonctio de l'individu=> Pour les stats
## En pourcentage
### Erreur definit comme NbPEI - NbPE avec :
###  NbPE : nombre de PE qui differe du 
def GetErrorRatePerc(NbAllele, NbPEI, NbAlt,NbRef):
    ErrorRate=-1.0
    if NbAllele==2:
        ErrorRate=float(NbPEI-(NbRef+NbAlt))/float(NbPEI)
    elif NbAllele==1:
        ErrorRate=float(NbPEI-(NbAlt))/float(NbPEI)
    elif NbAllele==3:
        ErrorRate=float(NbPEI-(NbRef))/float(NbPEI)
    return ErrorRate
    
def GetErrorRate(NbAllele, NbPEI, NbAlt,NbRef):
    ErrorRate=-1.0
    if NbAllele==2:
        ErrorRate=float(NbPEI-(NbRef+NbAlt))
    elif NbAllele==1:
        ErrorRate=float(NbPEI-(NbAlt))
    elif NbAllele==3:
        ErrorRate=float(NbPEI-(NbRef))
    return ErrorRate

## Controle de la 
## Return 1 quand cok
## 2 si inferieur a couverture
## 3 si redondance trop fore
def IsCorrectRed(NbRef, NbAlt,NbPEI, NbPeRed, NbPeIRed, MinCouv, MaxCouv, PercRed):
     if NbRef+NbAlt<MinCouv or (MaxCouv>0 and NbRef+NbAlt>MaxCouv ):
        return 2
     if float(NbPeRed)/float(NbPeIRed)<PercRed:
        return 1
     else :
        return 3
     return 1


def IsCorrectNonRed(NbRef, NbAlt,NbPEI, NbPeRed, NbPeIRed, MinCouv, MaxCouv, PercRed):
     if NbRef+NbAlt<MinCouv or (MaxCouv>0 and NbRef+NbAlt>MaxCouv ):
        return 2
     return 1

#return 0 si Homozygote Ref
##      1 Si Homozygote Alt
##      2 Si Het Ref/Alt
#def IsHomozygoteFiltreX(NbRef, NbAlt, NbPEI, NbPeRed)
 
## V1 :
  ## Homozygote si allele >0.75
  ## sinon homozygote
def IsHomozygoteFiltreV1(NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
    if PercAlt<= 0.25:
       return 0
    elif PercAlt>=0.75 :
       return 1
    else :
       return 2

### Filtre V2 et V2B sur des tests de Ficher :
##  * on definit qu'un locus est hétérozygote si critere filtre 1 pour V2 pas pour V2B
##  * Si la distribution suit un tableau de fisher avec comme hypothese l'homozygotie
def IsHomozygoteFiltreV2 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
   if PercAlt<=0.25 or GetPvalueFisherExactRef(NbRef,NbAlt)>PValue :
      return 0
   elif PercAlt >= 0.75 or GetPvalueFisherExactAlt(NbRef, NbAlt)>PValue :
      return 1
   else :
       return 2
     
def IsHomozygoteFiltreV2B (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
   if GetPvalueFisherExactRef(NbRef,NbAlt)>PValue :
      return 0
   elif GetPvalueFisherExactAlt(NbRef, NbAlt)>PValue :
      return 1
   else :
      return 2

### Filtre V3 :
## Test our verifier si la repartition et statistqiuement different d'une distribution 0.5  
def  IsHomozygoteFiltreV3 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet): 
   if GetPvalueFisherExactHet(NbRef,NbAlt)<PValue/2.0:
      if PercAlt<=0.5:
         return 0
      else : 
         return 1
   else :
     return 2

## 3 tests de Fishers, Pvalue la plus importante => Gagne
## Pas sur que cpomparer les P values sont adequates
def IsHomozygoteFiltreV4 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
    PvalueRef=GetPvalueFisherExactRef2(NbRef,NbAlt)
    PvalueHom=GetPvalueFisherExactAlt2(NbRef,NbAlt)
    ### Cas ou on est dans des faibles effectifs
    if NbPEI<25 :
       PvalueHet=GetPvalueFisherExactHet(NbRef,NbAlt)
       if PvalueRef>PvalueHom and PvalueRef>PvalueHet:
          return 0
       if PvalueHom>PvalueHet:
          return 1
       return 2
    else :
       if PvalueRef>PValue :
          return 0
       if PvalueHom>PValue:
          return 1
       return 2

def IsHomozygoteFiltreV5 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
    PvalueRef=GetPvalueBinomHom(NbRef,NbRef+NbAlt, BHom)
    PvalueHom=GetPvalueBinomHom(NbAlt, NbRef+NbAlt, BHom)
    if NbPEI<25 :
       PvalueHet=GetPvalueBinomHet(NbRef,NbAlt, NbRef+NbAlt, BHet)
       if PvalueRef>= PvalueHom and PvalueRef>=PvalueHet:
          return 0
       if PvalueHom>PvalueHet:
          return 1
       return 2
    else :
       if PvalueRef>PValue:
          return 0
       if PvalueHom>PValue:
          return 1
       return 2

##Permet de definir l'etat de la position d'un individu
## Avec une loi binomiale
## Test H0 : Est il Homozygote Reference
## Sinon Test H0 : Est il Homozygote Alternatif
## Si Test H1 : Est il heterozygote
## Sinon Pas d'assigantion
def IsHomozygoteFiltreV6 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
    PvalueRef=GetPvalueBinomHom(NbRef,NbRef+NbAlt, BHom)
    if PvalueRef>PValue:
       return 0
    PvalueHom=GetPvalueBinomHom(NbAlt, NbRef+NbAlt,BHom)
    if PvalueHom>PValue :
       return 1
    PvalueHet=GetPvalueBinomHet(NbRef,NbAlt, NbRef+NbAlt, BHet)
    if PvalueHet>PValue:
       return 2
    return -1

##Test de comparaison de vraisemblance
##

def IsHomozygoteFiltreV7 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
    NbRef=int(NbRef)
    NbAlt=int(NbAlt)
    ProbLogRef=BHom.Test(NbAlt,NbRef+NbAlt)
    ProbLogAlt=BHom.Test(NbRef,NbRef+NbAlt)
    ProbLogHet=BHet.Test(NbAlt,NbRef+NbAlt)
    PHom=max(ProbLogAlt, ProbLogRef)
    X2Score=-2*(min(ProbLogHet,PHom) - max(ProbLogHet,PHom))
    if X2Score<PValue :
       #print "Dist NC ", NbRef, NbAlt, X2Score
       return -1
    MaxProbLog=max([ProbLogRef, ProbLogAlt, ProbLogHet])
    if MaxProbLog==ProbLogRef :
       #print "Dist HomRef ", NbRef, NbAlt, X2Score
       return 0
    elif MaxProbLog==ProbLogAlt :
       #print "Dist AltRef ", NbRef, NbAlt, X2Score
       return 1
    elif MaxProbLog==ProbLogHet:
       #print "Dist HetRef ", NbRef, NbAlt, X2Score
       return 2
    else :
       sys.exit("tu t'es plante JT")
    return -1


def GetInfoIndCSVFiltre(StrInfoInd, MinCouv, MaxCouv, PercPEError, PercRed,PValue ,IsHomozygoteFiltreX, BHom, BHet, BaliseHetRed1Loc):
    A1=0
    A2=0
    InfoInd=StrInfoInd.split(':')
    TempAllele=InfoInd[0].split('/')
    A1Str=TempAllele[0]
    A2Str=TempAllele[1]
    ## Balise si Locus bon
    LocusBon=1
    ### Cas ou il n'y a pas de couverture
    if len(InfoInd)<14 :
       LocusBon=2
       NbPEI=0
       NbRef=0
       NbAlt=0
       NbAllele=0
    else :
        NbPEI=float(InfoInd[3])
        NbRef=float(InfoInd[4])
        NbAlt=float(InfoInd[5])
        if NbRef+NbAlt==0 or A1Str=='.' or A2Str=='.':
           LocusBon=2
        else :
           A1=int(A1Str)+1
           A2=int(A2Str)+1
           if PercRed==1:
             NbPeRed=0
             NbPeIRed=0
             LocusBon=(IsCorrectNonRed(NbRef, NbAlt,NbPEI, NbPeRed, NbPeIRed, MinCouv,MaxCouv,  PercRed))
           else :
             NbPeRed=float(InfoInd[15])
             NbPeIRed=float(InfoInd[14])
             LocusBon=(IsCorrectRed(NbRef, NbAlt,NbPEI, NbPeRed, NbPeIRed, MinCouv,MaxCouv ,PercRed))
    ### Cas ou la couverture est pas assez important
    if LocusBon==2 or NbRef+NbAlt==0:
    #if LocusBon!=1 or NbRef+NbAlt==0:
        A1=0
        A2=0
        NbAllele=0
        PercAlt=-1 
    elif BaliseHetRed1Loc==False and LocusBon==3 : 
          A1=0
          A2=0
          NbAllele=0
          PercAlt=-1
    else :
        if NbRef+NbAlt==0:
           sys.exit("Erreur NbRef+NbAlt==0 \n" + StrInfoInd + "sys exit\n" )
        PercAlt=NbAlt/(NbRef+NbAlt)
        ## Les test ne sont fait que si le pourcentage de l'alternatif > 0 ou < 1
        if PercAlt==0 :
           AlleleStatut=0
        elif PercAlt==1:
           AlleleStatut=1
        else :
###def IsHomozygoteFiltreV6 (NbRef, NbAlt, NbPEI, PercAlt, PValue,ErrorRate, BHom, BHet):
           AlleleStatut=IsHomozygoteFiltreX(NbRef, NbAlt, NbPEI, PercAlt, PValue, PercPEError, BHom, BHet) 
        if AlleleStatut==-1 :
           A1=0
           A2=0
           NbAllele=0
           PercAlt=-1
        elif AlleleStatut==0:
           NbAllele=3
           A1=1
           A2=1
        elif AlleleStatut==1:
           NbAllele=1
           A1=A2
        else :
           if A1==A2 :
          ## ypothes que 2 alleles pas toujours vrai...
              A1=1
              A2=2 
           NbAllele=2
    ### Traitement de la redondance
    ### Dans ce scenario on supprime juste les heterozygotes redondant (plus les homozygotes)
    #if BaliseHetRed1Loc==False :
    #   if LocusBon==3 and A1!=A2 :
    #      A1=0
    #      A2=0
    #      NbAllele=0
    #      PercAlt=-1 
    return (NbPEI, NbRef, NbAlt, PercAlt,NbAllele, A1, A2, ":"+":".join(InfoInd[1::]), GetErrorRate(NbAllele, NbPEI, NbAlt,NbRef), LocusBon)


def GetAlleleAltOld(FreqAllele,ListeAllele,NbAlt) :
    #if NbAlt==1 : 
    #   ListeAlt=ListeAllele[0]
    #   NbBase=2
    #else :
    if True :
       ListeAlleleNew=[]
       Cmt=0
       for FreqAlt in FreqAllele[2::] :
           if FreqAlt>0 :
              ListeAlleleNew.append(ListeAllele[Cmt])
           Cmt+=1
       NbAlt=len(ListeAlleleNew)
       if NbAlt==0:
          return [".",1]
       ListeAlt=ListeAllele[0]
       NbBase=len(ListeAlleleNew)+1
       for Alt in ListeAlleleNew[1::]:
          ListeAlt+=","+Alt
    return [ListeAlt, NbBase]
### A
def GetAlleleAlt(FreqAllele, ListeAllele,NbAlt, ListeAlleleInd, ListeInfoIndStr):
    ListeAlleleNew=[]
    CorrespondancePosition=["-1"]*(NbAlt+2)
    CorrespondancePosition[0]="." 
    CorrespondancePosition[1]="0" 
    Cmt=0
    CmtNewAllele=1
    for FreqAlt in FreqAllele[2::]:
       if FreqAlt>0 :
          ListeAlleleNew.append(ListeAllele[Cmt])
          CorrespondancePosition[Cmt+2]=str(CmtNewAllele)
          CmtNewAllele+=1
       Cmt+=1 
    NbAlt=len(ListeAlleleNew)
    if NbAlt==0:
       ListeAlt="."
       NbBase=1
    else :
        ListeAlt=ListeAllele[0]
        NbBase=len(ListeAlleleNew)+1
    for Alt in ListeAlleleNew[1::]:
        ListeAlt+=","+Alt
    Chaine=""
    for CmtAlleleInd in range(0,len(ListeAlleleInd)):
       Chaine+="\t"+CorrespondancePosition[ListeAlleleInd[CmtAlleleInd][0]]+"/"+CorrespondancePosition[ListeAlleleInd[CmtAlleleInd][1]]+ListeInfoIndStr[CmtAlleleInd]
    return (ListeAlt, NbBase, Chaine, ListeAlleleNew)

def GetNewStatStringPop(OldString, FreqInd, CouvMean):
    Temp=OldString.split(";")
    return "ADP="+str(int(CouvMean))+";WT="+str(FreqInd[3])+";HET="+str(FreqInd[2])+";HOM="+str(FreqInd[1])+";NC="+str(FreqInd[0])
## Donnee  
##   FreqInd : Frequence Het / Hom / Het
##   MeanCouvLocus : Moyenne de la couverture pour le locus
##   NbPErrorHom : Nombre de PE en Error pour les Hom alternatif et Ref
##   NbPeHom : Nombre de PE total par les homozygotes alt et ref
##   NbAlt : Nombre d'alternatif
## Param pour choix
##   MaxCouvMean : moyenne maximal souhate pour le locus
##   MinNbHomAlt : Nombre minimum d'homozyogte Alternatif souhaite
##   PercPEError => Pourcentage attendu pour l'erreur sur un locus
##   
def IsBonLocus(FreqInd, MeanCouvLocus,NbPerrorHom, NbPeHom, ListeAlt, Ref,MaxCouvMean ,MinNbHomAlt, PercPEError,NbLocusMauvaisDupIDD, NbAlt, PValue, TypePoly, PercNaMax, DistributionErreur, BaliseKs, BaliseHetRed1Loc, BalisePrintError=False) :
    ## Si pas d'hom alt ou de Het => Le locus	ne nous interesse pas 
#   if  FreqInd[1]+FreqInd[2] ==0 :
#       return [False]
    ## Si le locus est inferieur a MinNbHomAlt
   NbNC=FreqInd[0]
   #BalisePrint=True
   ## si presence des deux alleles
   if (FreqInd[1] <MinNbHomAlt or  FreqInd[3]<MinNbHomAlt) :
        if BalisePrintError :
            sys.stderr.write('No Diversity')
        return [False]
   PercNaObs=float(FreqInd[0])/float(sum(FreqInd))
   if PercNaObs>PercNaMax :
       if BalisePrintError :
            sys.stderr.write('Na Obs')
       return [False]
   ## Cas d'une couverture trop importante  
   if MaxCouvMean!=-1 and MeanCouvLocus>=MaxCouvMean :
       if BalisePrintError :
            sys.stderr.write('Coverage too important')
       return [False]
   ## suppression des locus ou il y a de la redondance et des individus heterozygotes
   if BaliseHetRed1Loc and NbLocusMauvaisDupIDD>0 and FreqInd[2]>0 :
      if BalisePrintError :
            sys.stderr.write('Depth Too important')
      return [False]
   NbAllele=len(ListeAlt)
   if  FreqInd[3]>0 :
       NbAllele+=1
   if NbAlt!=-1 and  NbAllele>NbAlt  :
       if BalisePrintError :
            sys.stderr.write('Number Alt not good')
       return [False]
   if TypePoly=='S':
      if len(Ref)>1 :
         if BalisePrintError :
            sys.stderr.write('Poly != S')
         return [False]
      for Alt in ListeAlt:
          if len(Alt)>1:
             if BalisePrintError :
                sys.stderr.write('More one 1 Alt')
             return [False] 
   if TypePoly=='I':
      BaliseSRef=False
      if len(Ref)>1 :
         BaliseSRef=True
      for Alt in ListeAlt:
          if len(Alt)>1:
             BaliseSRef=True
      if BaliseSRef==False :
         if BalisePrintError :
            sys.stderr.write('Poly != I')
         return [False]
   if NbPeHom>0 and PercPEError!=-1:
      PercPErrorObs=float(NbPerrorHom)/float(NbPeHom)
      if PercPErrorObs> PercPEError  and stats.binom_test(NbPerrorHom,NbPeHom, p= PercPEError)<PValue:
         if BalisePrintError :
            sys.stderr.write('Perc error false')
         return [False]
   if BaliseKs :
      StatKs=GetStatKS(DistributionErreur)
      if StatKs!=None and StatKs <0.05 :
        if BalisePrintError :
            sys.stderr.write('error doesnt Ks')
        return [False]
   return [True]

def IsBonLocusStat(FreqInd, MeanCouvLocus,NbPerrorHom, NbPeHom, ListeAlt, Ref,MaxCouvMean ,MinNbHomAlt, PercPEError, NbLocusMauvaisDupIDD, NbAlt, PValue, TypePoly, PercNaMax, DistributionErreur, BaliseKs, BaliseHetRed1Loc) :
    ## Si pas d'hom alt ou de Het => Le locus	ne nous interesse pas 
   Balise=True
   VecTest=["-"]*9
   #if  FreqInd[1]+FreqInd[2] ==0 :
   #     Balise=[False]
   #     VecTest[0]="F"
   #else: 
   #     VecTest[0]="T"
    ## Si le locus est inferieur a MinNbHomAlt
    ## Si le locus est inferieur a MinNbHomAlt
   NbNC=FreqInd[0]
   NbC=FreqInd[1]+FreqInd[2]+FreqInd[3]
   if NbC==0:
        Balise=False
        VecTest[1]="F"
   else :
        VecTest[1]="T"
   ## Cas plusieur indidivu + Cas ou j'ai un individu 
   if ((FreqInd[1] <MinNbHomAlt or  FreqInd[3]<MinNbHomAlt) and NbNC+NbC>1) or (NbC+NbNC==1 and FreqInd[1]<MinNbHomAlt) :
        Balise=False
        VecTest[1]="F"
   else: 
        VecTest[1]="T"
   PercNaObs=float(FreqInd[0])/float(sum(FreqInd))
   if PercNaObs>PercNaMax :
        Balise=False
        VecTest[2]="F"
   else: 
        VecTest[2]="T"
   ## Cas d'une couverture trop importante  
   if MaxCouvMean!=-1 and MeanCouvLocus>=MaxCouvMean :
        Balise=False
        VecTest[3]="F"
   else: 
        VecTest[3]="T"
   NbAllele=len(ListeAlt)
   if  FreqInd[3]>0 :
       NbAllele+=1
   if NbAlt!=-1 and  NbAllele>NbAlt  :
        Balise=False
        VecTest[4]="F"
   else: 
        VecTest[4]="T"
   if BaliseHetRed1Loc and NbLocusMauvaisDupIDD>0 and FreqInd[2]>0 :
      Balise=False
      VecTest[8]="F"
   else :
      VecTest[8]="T"
   if TypePoly=='S':
      if len(Ref)>1 :
          Balise=False
          VecTest[5]="F"
      else: 
          VecTest[5]="T"
      for Alt in ListeAlt:
          if len(Alt)>1:
             Balise=False
             VecTest[5]="F"

   if PercPEError!=-1 :
      if NbPeHom>0 :
         PercPErrorObs=float(NbPerrorHom)/float(NbPeHom)
         if PercPErrorObs> PercPEError  and stats.binom_test(NbPerrorHom,NbPeHom, p= PercPEError)<PValue:
            Balise=False
            VecTest[6]="F"
         else: 
            VecTest[6]="T"
      else :
         VecTest[6]="T"
   if BaliseKs :
      StatKs=GetStatKS(DistributionErreur)
      if StatKs <0.05 :
        Balise=False
        VecTest[7]="F"
      else: 
        VecTest[7]="T"
   return [Balise, VecTest] 

def GetInfoAllIndividu(SplitLigne, MinCouv, MaxCouv,PercPEError,MaxPercRed,Pvalue ,FoncFiltre, PositionAvecIDHet, BaliseHetID, BHom, BHet, BaliseHetRed1Loc) :
    Chro=SplitLigne[0]
    Pos=SplitLigne[1]
    Ref=SplitLigne[3]
    Alt=SplitLigne[4]
    CmtInd=0
    ListeAllele2=['0',Ref]
    ListeAllele2.extend(Alt.split(','))
    ListeInfoInd=[]
    for Cmt in range(9, len(SplitLigne)) :
        InfoInd=GetInfoIndCSVFiltre(SplitLigne[Cmt], MinCouv, MaxCouv,PercPEError,MaxPercRed,Pvalue ,FoncFiltre, BHom, BHet, BaliseHetRed1Loc)
        ListeInfoInd.append(InfoInd)
        if BaliseHetID :
           if InfoInd[4]==2 and (len(ListeAllele2[InfoInd[5]])>1 or len(ListeAllele2[InfoInd[6]])>1) :
              PositionAvecIDHet[CmtInd].append(int(Pos))
        CmtInd+=1
    ListeAllele=Alt.split(',')
    return [[Chro, Pos, Ref, ListeAllele], ListeInfoInd, SplitLigne]

def CheckInsertionDeletion(ListePos,PosCourante, Dist) :
     Cmt=0
     #for Cmt in range(len(ListePos)) :
     for Pos in ListePos:
         DiffPos=PosCourante-Pos
         if abs(DiffPos)<Dist :
            return False
         if DiffPos>Dist :
            ListePos.remove(Pos)
         elif DiffPos<Dist :
            return True
     return True


def GetTransforIndividu(ListeInfoInd, PositionAvecIDHet, EcrireStat, Ecrire, PosInt, FonctionFiltreBonLocus, MaxCouvMean ,MinNbHomAlt, PercPEError, NbAltParam,  Pvalue, TypePoly, PercNaMax,  BaliseKs, StatParInd, PrintAllLocus, BaliseHetID, BaliseHetRed1Loc):
           ## Recupere le chro / pos la valeur des alt et ref
           ## On intialise les vecteurs qui COntinnent
           ## Frequence de l'allele
           SplitLigne=ListeInfoInd[2]
           Chro=ListeInfoInd[0][0]
           Pos=ListeInfoInd[0][1]
           Ref=ListeInfoInd[0][2]
           ListeAllele=ListeInfoInd[0][3]
           NbAlt=len(ListeAllele)
           FreqAllele=[0]*(NbAlt+2)
           ## Si le locus est bon sinon pk redondance / couverture
           DistTypeFiltre=[0]*5
           ## Frequence des ind : Het NA Alt Ref
           FreqInd=[0]*4
           ChaineTemp=""
           ## Nombre de PE en erreur
           NbPEErreur=0
           ## Nombre de Pe tot
           NbPEITot=0
           ## Nb individu Couvert
           CmtCouv=0.0
           ## Nombe de Pe Erreur pour les Homozygotes
           NbPEErreurHom=0
           #Nombre dePe tot calcule sur les Hom alt et ref
           NbPEITotHom=0
           ## Nb Individu sur lequel cela a ete calcule
           CmtCouvHom=0
           ## Pour chaque individu
           ListeInfoIndStr=[]
           ListeInfoAllele=[]
           ## Pour chaque individu du vcf
           DistributionErreur=[]
           StatParLocuParInd=[]
           ### Nombre de locus mauvais
           NbLocusMauvaisDupID=0
           for Cmt in range(len(ListeInfoInd[1])):
               InfoInd=ListeInfoInd[1][Cmt]
               if BaliseHetID and int(InfoInd[4])==2:
                  BaliseID=CheckInsertionDeletion(PositionAvecIDHet[Cmt], PosInt, 500)
                  if BaliseID==False :
                     #NbLocusMauvaisDupID+=1
                     InfoInd=(InfoInd[0],InfoInd[1], InfoInd[2],-1, 0, 0,0,InfoInd[7],-1,4 )
               ### Calcul du nombre d'individu ou il y a de la redondance
               if InfoInd[9]==3 :
                  NbLocusMauvaisDupID+=1
               # on recupere les information et applique les filtres
               ## On enregistre les frequences des alleles
               DistTypeFiltre[InfoInd[9]]+=1
               FreqAllele[int(InfoInd[5])]+=1
               FreqAllele[int(InfoInd[6])]+=1
               FreqInd[int(InfoInd[4])]+=1
               StatParLocuParInd.append(int(InfoInd[4]))
               ## On recupere l'information des individus
               ListeInfoIndStr.append(InfoInd[7])
               ListeInfoAllele.append([InfoInd[5],InfoInd[6]])
               DistributionErreur.append(InfoInd[8]) 
               ## On enregistre le taux d erreur pour les homozygotes
               if InfoInd[8]!=-1 and InfoInd[4]!=2:
                  CmtCouvHom+=1
                  NbPEErreurHom+=InfoInd[8]
                  NbPEITotHom+=InfoInd[0]
               if InfoInd[8]!=-1 :
                  CmtCouv+=1
                  NbPEErreur+=InfoInd[8]
                  NbPEITot+=InfoInd[0]
           ListeAlt=GetAlleleAlt(FreqAllele,ListeAllele,NbAlt, ListeInfoAllele, ListeInfoIndStr)
           if CmtCouv>0 :
               MeanCouvLocus=float(NbPEITot)/float(CmtCouv)
           else :
               MeanCouvLocus=0
           ## Definit si le locus est bon => 
           BalisePrintError=False
           BonLocus=FonctionFiltreBonLocus(FreqInd,MeanCouvLocus,NbPEErreurHom,  NbPEITotHom, ListeAlt[3], Ref,MaxCouvMean ,MinNbHomAlt, PercPEError,NbLocusMauvaisDupID ,NbAltParam,  Pvalue, TypePoly, PercNaMax, DistributionErreur, BaliseKs, BaliseHetRed1Loc, BalisePrintError)
           if BonLocus[0]==False and BalisePrintError :
              sys.stderr.write(SplitLigne[0]+"\t"+SplitLigne[1]+"\t"+SplitLigne[2]+"\t"+ SplitLigne[3]+"\n")
           ## On recalcu
           if BonLocus[0] :
               FinalChaine=SplitLigne[0]+"\t"+SplitLigne[1]+"\t"+SplitLigne[2]+"\t"+ SplitLigne[3]+ "\t" + ListeAlt[0]+"\t"+SplitLigne[5]+"\t"+SplitLigne[6]+"\t"+GetNewStatStringPop(SplitLigne[7], FreqInd, MeanCouvLocus)+"\t"+SplitLigne[8]+ListeAlt[2]
               Ecrire.write(FinalChaine+"\n")
               for CmtInd in range(0,len(ListeInfoInd[1])) :
                   StatParInd[CmtInd][StatParLocuParInd[CmtInd]]+=1
           elif PrintAllLocus==True :
                FinalChaine=SplitLigne[0]+"\t"+SplitLigne[1]+"\t"+SplitLigne[2]+"\t"+ SplitLigne[3]+ "\t" + ListeAlt[0]+"\t"+SplitLigne[5] +"\t"+"NOPASS"+"\t"+GetNewStatStringPop(SplitLigne[7], FreqInd, MeanCouvLocus)+"\t"+SplitLigne[8]+ListeAlt[2]
                Ecrire.write(FinalChaine+"\n")
           if EcrireStat!=None :
             EcrireStat.write(SplitLigne[0]+"_"+SplitLigne[1]+"\t"+str((ListeAlt[1]))+"\t" +str(FreqInd[3])+"\t"+str(FreqInd[2])+"\t"+str(FreqInd[1])+"\t"+str(FreqInd[0]) +"\t"+str(NbPEErreur)+"\t"+str(NbPEITot)+"\t"+str(CmtCouv)+"\t"+str(NbPEErreurHom)+"\t"+str(NbPEITotHom)+"\t"+str(CmtCouvHom)+"\t"+str(DistTypeFiltre[2])+"\t"+str(DistTypeFiltre[3])+"\t"+str(BonLocus[0]))
             EcrireStat.write("\t"+"\t".join(BonLocus[1]))
             EcrireStat.write("\n")
     
    

         

## chromosome:AGPv2:10:1:150189435:1       6419    .       T       G,C     .       PASS  ADP=60;WT=6;HET=30;HOM=8;NC=0   GT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR    0/1:32:16:15:3:8:53,33%:5,16E-4:39:35:3:0:3:5
### Filtre sur chacun des locus
### Modifie un fichier de sortie par Varscan
### 
#MaxCouvMean ,MinNbHomAlt, PercPEError, MaxNbAltBaliseBalise
def TransformMultiCSV(Fichier, FichierSortie, FichierSortieStat, FichierSortieStatParInd,FoncFiltre,MinCouv, MaxCouv, MaxCouvMean, NbAltParam,MinNbHomAlt,MaxPercRed,PercPEError ,ListeInd, Pvalue, PrintAllLocus, TypePoly,NomModel, PercNaMax, BaliseKs,BaliseHetID, BaliseHetRed1Loc):
    ### Lecture du vcf 
    Lire=OuvrirFichier(Fichier, 'r')
    ## Sortie des nouveaux VCF
    Ecrire=OuvrirFichier(FichierSortie,'w')
    EcrireStat=None
    ## Si la fonction est V6 on initialise les stat
    if FoncFiltre==IsHomozygoteFiltreV6 or FoncFiltre==IsHomozygoteFiltreV5:
       BHom=Binom(200,PercPEError)
       BHet=Binom(200,0.5)
    elif FoncFiltre==IsHomozygoteFiltreV7:
       BHom=Binom(200,PercPEError, 'probLog')
       BHet=Binom(200,0.5, 'probLog')
       Pvalue=stats.chi2.ppf(1-Pvalue,1)
    else: 
       BHom=None
       BHet=None
    ## On ecrit les stats que si on a un fichier 
    if FichierSortieStat!=None :
      EcrireStat=OuvrirFichier(FichierSortieStat,'w')
      FonctionFiltreBonLocus=IsBonLocusStat
    else :
      FonctionFiltreBonLocus=IsBonLocus
    ### Ecriture des entetes
    if FichierSortieStat!=None :
      EcrireStat.write("Pos\tNbAlt\tNbIndRef\tNbIndHet\tNbIndHom\tNbNonBon\tNbPEErreurHom\tNbPEIHom\tNbIndNonHet\tNbPEErreur\tNbPEI\tNbInd\tNbNC\tNbRed\tPass")
      EcrireStat.write("\tBdiv\tBNbHom\tBNA\tBMaxCouv\tBnbAlt\tBIns\tNErr\tBks\tBDupl")
      EcrireStat.write("\n")
    ### Pour chaque ligne
    ListeLocus=[]
    CmtLocus1=0
    for Ligne in Lire :
       ## On verifie qu'il ne s'agit pas d'une entete
        if Ligne[0]!="#" and len(Ligne)>5 : 
           ## Decoupage lignie
           SplitLigne=Ligne.split()
           ## On sauvegarde les informations sur les individu
           ListeInfo1Locus=GetInfoAllIndividu(SplitLigne, MinCouv, MaxCouv,PercPEError,MaxPercRed,Pvalue ,FoncFiltre, PositionAvecIDHet, BaliseHetID, BHom,BHet, BaliseHetRed1Loc)
           ListeLocus.append(ListeInfo1Locus)
           Position=int(ListeInfo1Locus[0][1])
           if CmtLocus1==0 :
              MinPosition=Position
              ChroI=ListeInfo1Locus[0][0]
           if Position-MinPosition>5000 :
              Chro=ListeLocus[0]
              while Position-MinPosition>500 and ChroI==Chro:
                   MinPosition=int(ListeLocus[0][0][1])
                   GetTransforIndividu(ListeLocus[0], PositionAvecIDHet, EcrireStat, Ecrire, MinPosition, FonctionFiltreBonLocus,MaxCouvMean ,MinNbHomAlt, PercPEError, NbAltParam,  Pvalue, TypePoly, PercNaMax,  BaliseKs, StatParInd,PrintAllLocus, BaliseHetID, BaliseHetRed1Loc) 
                   del(ListeLocus[0])
                   if len(ListeLocus) >0 : 
                      Chro=ListeLocus[0]
              ChroI=Chro
           CmtLocus1+=1
        else :
          #MinCouv, MaxCouv, MaxCouvMean, NbAltParam,MinNbHomAlt,MaxPercRed,PercPEError ,ListeInd, Pvalue, PrintAllLocus, TypePoly,NomModel
          if Ligne[0:6]=="#CHROM" :
             #Ecrire.write("#Param Model:"+NomModel+";Min Couv by individu :"+str(MinCouv)+";Max Couv by individu :"+str(MaxCouv)+";Max Couv mean by locus :"+str(MaxCouvMean)+";Min Couv by individu :"+str(MinCouv)+";\n")
             #Ecrire.write("#Number Allele by locus"+str(NbAltParam)+";Percentage max of redondance"+str(MaxPercRed)+";Percentage of error"+str(PercPEError)+";Alpha risk"+str(Pvalue)+"; Print All locus:"+str(PrintAllLocus)+";Type Of polymorphims : "+str(TypePoly)+";Perc Na Max "+str(PercNaMax)+";\n")
             if ListeInd!=None:
                Ligne="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
                for Ind in ListeInd:
                    Ligne+="\t"+Ind
                Ligne+="\n"
                ListeIndFin=ListeInd
             else :
                ListeIndFin=Ligne.split()[9::]
             NbInd=len(ListeIndFin)
             StatParInd=GetMatrice(NbInd, 4)
             if BaliseHetID:
                PositionAvecIDHet=[[] for x in range(NbInd)]
             else :
                PositionAvecIDHet=None
          Ecrire.write(Ligne)
    for Locus in ListeLocus :
        GetTransforIndividu(Locus, PositionAvecIDHet, EcrireStat, Ecrire, MinPosition, FonctionFiltreBonLocus,MaxCouvMean ,MinNbHomAlt, PercPEError, NbAltParam,  Pvalue, TypePoly, PercNaMax,  BaliseKs, StatParInd,PrintAllLocus, BaliseHetID, BaliseHetRed1Loc)
    Lire.close()
    Ecrire.close()
    if FichierSortieStat!=None :
       EcrireStat.close()
    if FichierSortieStatParInd!=None :
        EcrireStatInd=open(FichierSortieStatParInd,'w')
    else :
        EcrireStatInd=sys.stdout
    EcrireStatInd.write("Ind\tNbNA\tNbAlt\tNbHet\tNbRef\n")
    for CmtInd in range(0,NbInd) :
        EcrireStatInd.write(ListeIndFin[CmtInd])
        for Stat in StatParInd[CmtInd] :
             EcrireStatInd.write("\t"+str(Stat))
        EcrireStatInd.write("\n")
    if FichierSortieStatParInd!=None :
        EcrireStatInd.close()
    #return Dico

##
def GetUnArgument(ListeArgv, opt, nMin, nMax, Cmd):
    BaliseSave=False
    Argc=len(ListeArgv)
    cmt=0
    ListeArg=[]
    while cmt<Argc:
        if opt==ListeArgv[cmt].lower() and BaliseSave==False:
           BaliseSave=True
        elif BaliseSave and ListeArgv[cmt][0] == '-' :
           return ListeArg
        elif BaliseSave :
           ListeArg.append(ListeArgv[cmt])
        cmt+=1
    if len(ListeArg)>nMax or len(ListeArg)<nMin:
       sys.exit("Error with arg"+opt+"\n Cmd :"+Cmd)
    if len(ListeArg)==0:
       return [None]
    return ListeArg


def Launch():
    ListeArg=["-f", "-i","-o","-os", "-e", "-snf", "-cmin", "-prm", "-pv", "-pa", "-cmax", "-cmaxmean", "-nballmax", "-nbindalt", "-percnamax", "-osi", "-tks", "-hetdlid", "-hetinddup", '-t']
    ### Parametre a integrer dans la ligne de commande
    # Balise pour supprimer les heterozygotes lorsqu'il y a une I/D heterozygote a + - 500
    for cmt in range(1,len(sys.argv),2):
       if sys.argv[cmt].lower() not in ListeArg:
            sys.exit("Arg "+ sys.argv[cmt] + " not known \n"+Cmd)
    TypeFiltre=GetUnArgument(sys.argv, "-f", 1, 1, Cmd)[0]
    FileIn=GetUnArgument(sys.argv, "-i", 1,1, Cmd)[0]
    FileOut=GetUnArgument(sys.argv, "-o", 1,1, Cmd)[0]
    FileOutStat=GetUnArgument(sys.argv, "-os",0,1, Cmd)[0]
    PercPEError=GetUnArgument(sys.argv, "-e", 0,1,Cmd)[0]
    if PercPEError==None:
       PercPEError=0.01
    else :
       PercPEError=float(PercPEError)
    FileOutStatInd=GetUnArgument(sys.argv, "-osi",0,1, Cmd)[0]
    NomEchantillon= GetUnArgument(sys.argv, "-snf", 0, 1,Cmd)[0]
    ListeInd=None
    if NomEchantillon!=None:
       ListeInd=[]
       LireInd=open(NomEchantillon)
       for Ligne in LireInd:
          ListeInd.append(Ligne.replace('\n','')) 
       LireInd.close()
    if TypeFiltre=="V1":
        FoncFiltre=IsHomozygoteFiltreV1
    elif TypeFiltre=="V2":
        FoncFiltre=IsHomozygoteFiltreV2
    elif TypeFiltre=="V3":
        FoncFiltre=IsHomozygoteFiltreV3
    elif TypeFiltre=="V4":
        FoncFiltre=IsHomozygoteFiltreV4
    elif TypeFiltre=="V5":
        FoncFiltre=IsHomozygoteFiltreV5
    elif TypeFiltre=="V2B":
        FoncFiltre=IsHomozygoteFiltreV2B
    elif TypeFiltre=="V6":
        FoncFiltre=IsHomozygoteFiltreV6
    elif TypeFiltre=="V7":
        FoncFiltre=IsHomozygoteFiltreV7
    else :
        sys.exit("Modele no known")
    ## Couverture Minimum pour un individu
    MinCouv=GetUnArgument(sys.argv, "-cmin", 0,1,Cmd)[0]
    if MinCouv==None:
       MinCouv=0
    MinCouv=float(MinCouv)
    MaxPercRed=GetUnArgument(sys.argv, "-prm", 0,1, Cmd)[0]
    if MaxPercRed==None:
       MaxPercRed=1.0
    MaxPercRed=float(MaxPercRed)
    Pvalue=GetUnArgument(sys.argv, "-pv", 0,1, Cmd)[0]
    if Pvalue==None:
       Pvalue=0.05
    Pvalue=float(Pvalue)
    ### Couverture Maximum 
    MaxCouv=GetUnArgument(sys.argv, "-cmax", 0,1, Cmd)[0]
    if MaxCouv==None:
       MaxCouv=-1
    MaxCouv=float(MaxCouv)
    MaxCouvMean=GetUnArgument(sys.argv, "-cmaxmean", 0,1, Cmd)[0]
    if MaxCouvMean==None:
       MaxCouvMean=-1
    MaxCouvMean=float(MaxCouvMean)
    PrintAllLocus=GetUnArgument(sys.argv, "-pa", 0,1, Cmd)[0]
    if PrintAllLocus==None:
       PrintAllLocus=False
    else:
       PrintAllLocus=PrintAllLocus.lower()
       if PrintAllLocus=='f' or PrintAllLocus=='false':
          PrintAllLocus=False
       else :
          PrintAllLocus=True
    ## NbAllMax=> Nombre d'alternatif maximum => A mettre a -1
    NbAllMax=GetUnArgument(sys.argv, "-nballmax", 0,1, Cmd)[0]
    if NbAllMax==None:
       NbAllMax=-1
    NbAllMax=float(NbAllMax)
    if NbAllMax==1 :
       print "warning : Nombre d'allele au locus demande est de 1"
    
    MinNbHomAlt=GetUnArgument(sys.argv, "-nbindalt", 0,1, Cmd)[0]
    if MinNbHomAlt==None :
       MinNbHomAlt=0
    MinNbHomAlt=float(MinNbHomAlt)
    
    TypePoly=GetUnArgument(sys.argv, "-t", 0,1,Cmd)[0]
    if TypePoly==None :
       TypePoly="S"
    if TypePoly not in ['S','I', 'A'] :
       sys.exit("TypePoly not known : "+TempTypePoly)


    PercNaMax=GetUnArgument(sys.argv, "-percnamax", 0,1, Cmd)[0]
    if PercNaMax==None:
       PercNaMax=1.0
    PercNaMax=float(PercNaMax)
    
    BaliseKs=GetUnArgument(sys.argv, "-tks", 0,1, Cmd)[0]
    if BaliseKs==None:
       BaliseKs=False
    else:
       BaliseKs=BaliseKs.lower()
       if  BaliseKs=='f' or  BaliseKs=='false':
           BaliseKs=False
       elif BaliseKs=='true' or BaliseKs=='t':
          print "Tks test perfomed"
          BaliseKs=True
       else :
          sys.exit('arg tks : t/f arg given unknown: '+BaliseKs )

    BaliseHetID=GetUnArgument(sys.argv, "-hetdlid",0,1, Cmd)[0]
    if BaliseHetID==None or BaliseHetID=='F' or BaliseHetID=='f':
       BaliseHetID=False
    elif BaliseHetID=='T' or BaliseHetID=='t':
       print "Delete heterozygote in LD with I/D heterozygous"
       BaliseHetID=True
    else :
       sys.exit("Param -hetdlId false T/F")

    BaliseHetRed1Loc=GetUnArgument(sys.argv, "-hetinddup", 0,1, Cmd)[0]
    if BaliseHetRed1Loc==None:
       BaliseHetRed1Loc=False
    else:
       BaliseHetRed1Loc=BaliseHetRed1Loc.lower()
       if  BaliseHetRed1Loc=='f' or  BaliseHetRed1Loc=='false':
           BaliseHetRed1Loc=False
       elif BaliseHetRed1Loc=='true' or BaliseHetRed1Loc=='t':
          print "Het Red 1 Loc test perfomed"
          BaliseHetRed1Loc=True
       else :
          sys.exit('arg hetinddup : t/f arg given unknown: '+BaliseHetRed1Loc )
    ### MinNbHomAlt => Nombre minimum d'homozygote pour un locus )=> defaut 0
    DataVCF=TransformMultiCSV(FileIn, FileOut,FileOutStat,FileOutStatInd,FoncFiltre,MinCouv, MaxCouv,MaxCouvMean,NbAllMax, MinNbHomAlt,MaxPercRed,PercPEError,ListeInd, Pvalue, PrintAllLocus, TypePoly, TypeFiltre, PercNaMax, BaliseKs, BaliseHetID, BaliseHetRed1Loc)



if __name__ == '__main__':
    Cmd="exe -f FiltreType -i VCF File -o Out CSV File [-os Out Stat] [-e Error Rate] [-cmin coverage min for one individual default : 0] [-cmax coveraga max for one individual default : NA] [-cmaxmean mean coverage max for one locus  ] [-nballmax Nombre d'allele maximal by locus : default NA] [-nbindalt : individual number with alternatif allele default : 0]  [-prm percentage of maximal redondance] [-pv value] [-pa print All locus false] [-snf nom des echantillons optionnel] [-hetdlId : Supression of position by individu where Het is in LD with Heterozygous I/D Default : F]  -hetinddup"
    Launch()

