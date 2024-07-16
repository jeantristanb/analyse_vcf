#### *** Authors : Jean-Tristan Brandenburg ***
# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys
import math

### compute Pi of gesco
def nCr(n,r):
    fact = math.factorial
    return fact(n) / fact(r) / fact(n-r)

def GetPiCresko(a1,a2):
    '''
    Pi developpend by Cresko for one position
      *  formule (1 - (2 parmis n1 + 2 parmis n2))/ (2 parmis n)
      see : Population Genomics of Parallel Adaptation in Threespine Stickleback using Sequenced RAD Tags
    input :
      n1  nb allele 1
      n2 nb allele 2
    output :
       pi value 
    '''
    if a1<2  or a2 < 2 :
       return 0
    return 1- (nCr(a1,2)+nCr(a2,2))/float(nCr(a1+a2,2) )



def TajimasD(n, avg_pair_diff, s):
    '''
    compute tajima_d 
    input :
      n : individu number
      avg_pair_diff : mean of difference two by two
      s : site number segregating 
      extract of : to foun
    output :
      D values
    '''

    a1 = sum([1.0/i for i in range(1, n)])
    a2 = sum([1.0/(i**2) for i in range(1, n)])
    b1 = float(n+1)/(3*(n-1))
    b2 = float(2 * ( (n**2) + n + 3 )) / (9*n*(n-1))
    c1 = b1 - 1.0/a1
    c2 = b2 - float(n+2)/(a1 * n) + float(a2)/(a1 ** 2)
    e1 = float(c1) / a1
    e2 = float(c2) / ( (a1**2) + a2 )
    D = (
        float(avg_pair_diff - (float(s)/a1))
        / math.sqrt(
            (e1 * s )
          + ((e2 * s) * (s - 1) ))
        )
    return D


##
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

def TajimaDWithMissingData(ListeA1,ListeA2,PreComputea1,PreComputea2,MaxN, MinFreqNa, MinFreqS, SumThetaPi=None, SumThetaW=None, S=None, NbPosition=None) :
    '''
    compute d tajima with missing values
    see : http://www.cell.com/cms/attachment/599136/4711220/mmc1.pdf
    input : 
      * ListeA1 : list contains for each position nb allele of A1
      * ListeA2 : list contains for each position nb allele of A2
      * PreComputea1 : value of a1 precomputed 
      * PreComputea2 : value of a2 precomputed 
      * MaxN : number of allele ( 2 * n individu)
      * MinFreqNa : cut off for missing values for each position
      * MinFreqS : cut off for position polymorphics
      * another way fo computed tajima D :
          if SumThetaPi, SumThetaW, S, NbPosition is allowed, values will be used to computed tajima D
    output :
       D as SumDifTheta/Var 
       S as position polymorphic 
       SumThetaPi :
       SumThetaW : 
       MaxN : value N give in entrance3 
       NbPosition : NbPositon analyse
       S : site polymorphic 

    '''
    if SumThetaW==None and SumThetaPi==None and S==None:
       SumDifTheta=0
       SumThetaPi=0
       SumThetaW=0
       S=0.0
       NbPosition=0
       for Cmt in range(0, len(ListeA1)):
           N=ListeA1[Cmt]+ListeA2[Cmt]
           FreqNa=N/float(MaxN)
           if FreqNa > MinFreqNa :
              Freq=ListeA1[Cmt]/float(N)
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

def ComputeDivDtajima(Spectre) :
    SF=Spectrum_mod.Spectrum(Spectre)
    return SF.Tajima_D()

def FstGreskoAll(a1,a2,n1,n2) :
    '''
    compute gresco FST for one position
    input :
      * a1 : nb genotype of allele1 for pop 1 
      * a2 : nb genotype of allele1 for pop 2 
      * n1 : nb genotype total for pop 1 
      * n2 : nb genotype total pop 2 
    output :
      * return fst values for locus
    '''
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



def FstBioConduct(a1,a2,n1,n2) :
    '''
    compute FST for one position
    input :
      * a1 : nb genotype of allele1 for pop 1 
      * a2 : nb genotype of allele1 for pop 2 
      * n1 : nb genotype total for pop 1 
      * n2 : nb genotype total pop 2 
    output :
      * return fst values for locus
    '''
   NAll=n1+n2
   F1=a1/n1
   F2=a2/n2
   FAll=(a1+a2)/NAll
   Ys=FAll*(1-FAll)*(NAll/(NAll-1))
   NNMoin1Pop1=(n1/float(n1-1))
   NNMoin1Pop2=(n2/float(n2-1))
   Xs1=F1*(1-F1)*NNMoin1Pop1
   Xs2=F2*(1-F2)*NNMoin1Pop2
   Weigth1=n1*NNMoin1Pop1/(n1*NNMoin1Pop1+n2*NNMoin1Pop2)
   Weight2=n2*NNMoin1Pop2/(n1*NNMoin1Pop1+n2*NNMoin1Pop2)
   Xs=Weigth1*Xs1 + Weigth2*Xs2
   Fs=1-(Xs/Ys)
   return Fs




#### Calcul du Fst selon reich, nature 2009
def ComputeFstReich(a1,a2, n1, n2) :
    '''
    compute FST for one position with reich , Fst is obtained by division of N/D
    to do added publication
    input :
      * a1 : nb genotype of allele1 for pop 1 
      * a2 : nb genotype of allele1 for pop 2 
      * n1 : nb genotype total for pop 1 
      * n2 : nb genotype total pop 2 
    output :
      * return N value 
      * return D values 
    '''
    if (a1+a2)==(n1+n2) or a1+a2==0 :
       sys.exit("(a1+a2)==(n1+n2) a1+a2==0")
    h1=(a1*(n1-a1))/(n1*(n1-1))
    h2=(a2*(n2-a2))/(n2*(n2-1))
    N=(a1/n1 - a2/n2)**2 -h1/n1 -h2/n2
    D=N+h1+h2
    return (N,D)


