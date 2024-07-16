#-*- coding: utf-8 -*-

##
##pa1 : frequence de l'allele 1 la position a
##pb2 : frequence de l'allele 1 la position b
## n11 : nombre d'haplotype a1 b1 
## n12 : nombre d'haplotype a1 b2 
## n22 : nombre d'haplotype a2 b2 
## n21 : nombre d'haplotype a2 b1 
#def ComputeLD(pa1,pb1,p11, n12, n22, n21) :
def ComputeLDPos(p11,p22,p12,p21, pa1,pb1):
     #raw difference in frequency between the observed number of AB pairs and the expected number:
     D = p11*p22-p12*p21
     pa2=1-pa1
     pb2=1-pb1
     if pa1==0 or  pa2==0 or pb1==0 or  pb2==0 :
        return (1, 1, pa1,pb1)
     if D> 0 :
        ## Dmax = min( p(A)p(b), p(a)p(B) ) 
        Dmax = min(pa1*pb2, pa2*pb1)
     else :
        #Dmax = max( -p(A)p(B), -p(a)p(b) ) 
        Dmax = min(-pa1*pb1, -pa2*pb2)
     Dprime=D/Dmax
     Rcar=D/((pa1*pa2*pb1*pb2)**0.5)
     Rcar=Rcar**2
     return (Dprime,Rcar,pa1,pb1)

## prevu pour des ge
## Snp1 : distribution du snp 1 pour la pop : 0 donnee manquante, 1 Allele1 2 allele2
def ComputeLDDeuxPos(Snp1,Snp2) :
    NbInd=float(len(Snp1))
    Cmt=0
    Dist=[[0,0,0],[0,0,0],[0,0,0]]
    while Cmt<NbInd:
         Dist[Snp1[Cmt]][Snp2[Cmt]]+=1
         Cmt+=1
    pa1=(Dist[1][1]+Dist[1][2])/NbInd
    pb1=(Dist[1][1]+Dist[2][1])/NbInd
    p11=Dist[1][1]/NbInd
    p12=Dist[1][2]/NbInd
    p22=Dist[2][2]/NbInd
    p21=Dist[2][1]/NbInd
    return ComputeLDPos(p11,p22,p12,p21 ,pa1,pb1)


### format [Snp1, Snp2, Snp3], Snp1=[1,2,2,1]
def ComputeLD1Segment(ListeSnp, MinFreq): 
    NbSnp=len(ListeSnp)
    ListeDprime=[]
    ListeRcarre=[]
    ListeDprimeMinFreq=[]
    ListeRcarreMinFreq=[]
    NbInd=len(ListeSnp[0])
    for CmtSnp in range(0,len(ListeSnp)-1) : 
        for CmtSnp1 in range(CmtSnp+1,len(ListeSnp)) : 
            (DPrime, Rcar,pa1,pb1)=ComputeLDDeuxPos(ListeSnp[CmtSnp],ListeSnp[CmtSnp1])
            if pa1>MinFreq and pb1>MinFreq and (1-pa1)>MinFreq and (1-pb1)>MinFreq:
               ListeDprimeMinFreq.append(DPrime)
               ListeRcarreMinFreq.append(Rcar)
            ListeDprime.append(DPrime)
            ListeRcarre.append(Rcar)
    NbCmp=float(len(ListeRcarre))
    if NbCmp>0:
       MeanDprime=sum(ListeDprime)/NbCmp
       MeanRcarre=sum(ListeRcarre)/NbCmp
    else :
       MeanDprime=0
       MeanRcarre=0

    NbCmpMinFreq=float(len(ListeRcarreMinFreq))
    if NbCmpMinFreq>0:
       MeanDprimeMinFreq=sum(ListeDprimeMinFreq)/NbCmpMinFreq
       MeanRcarreMinFreq=sum(ListeRcarreMinFreq)/NbCmpMinFreq
    else :
       MeanDprimeMinFreq=0
       MeanRcarreMinFreq=0
    NbSignif=0
    for RCarre in ListeRcarre:
        if RCarre*NbInd> 3.84 :
           NbSignif+=1
    NbSignifMinFreq=0
    for RCarre in ListeRcarreMinFreq:
        if RCarre*NbInd> 3.84 :
           NbSignifMinFreq+=1
    return (MeanDprime,MeanRcarre, NbSignif,NbCmp,MeanDprimeMinFreq,MeanRcarreMinFreq, NbSignifMinFreq,NbCmpMinFreq)



def CompareLD2Segment(ListeSnp1, ListeSnp2, MinFreq) :
    NbSnp1=len(ListeSnp1)
    NbSnp2=len(ListeSnp2)
    ListeDprime=[]
    ListeRcarre=[]
    ListeDprimeMinFreq=[]
    ListeRcarreMinFreq=[]
    NbInd=len(ListeSnp1[0])
    for CmtSnp1 in range(0,NbSnp1):
        for CmtSnp2 in range(0,NbSnp2):
            (DPrime, Rcar,pa1,pb1)=ComputeLDDeuxPos(ListeSnp1[CmtSnp1],ListeSnp2[CmtSnp2])
            if pa1>MinFreq and pb1>MinFreq and (1-pa1)>MinFreq and (1-pb1)>MinFreq:
               ListeDprimeMinFreq.append(DPrime)
               ListeRcarreMinFreq.append(Rcar)
            ListeDprime.append(DPrime)
            ListeRcarre.append(Rcar)
    NbCmp=float(len(ListeRcarre))
    if NbCmp>0:
       MeanDprime=sum(ListeDprime)/NbCmp
       MeanRcarre=sum(ListeRcarre)/NbCmp
    else :
       MeanDprime=0
       MeanRcarre=0

    NbCmpMinFreq=float(len(ListeRcarreMinFreq))
    if NbCmpMinFreq>0:
       MeanDprimeMinFreq=sum(ListeDprimeMinFreq)/NbCmpMinFreq
       MeanRcarreMinFreq=sum(ListeRcarreMinFreq)/NbCmpMinFreq
    else :
       MeanDprimeMinFreq=0
       MeanRcarreMinFreq=0
    NbSignif=0
    for RCarre in ListeRcarre:
        if RCarre*NbInd> 3.84 :
           NbSignif+=1
    NbSignifMinFreq=0
    for RCarre in ListeRcarreMinFreq:
        if RCarre*NbInd> 3.84 :
           NbSignifMinFreq+=1
    return (MeanDprime,MeanRcarre, NbSignif,NbCmp,MeanDprimeMinFreq,MeanRcarreMinFreq, NbSignifMinFreq,NbCmpMinFreq)
