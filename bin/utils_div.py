# coding=utf-8

def compute_ld_pos_pos(p11,p22,p12,p21, pa1,pb1, n, Na="NA"):
     #raw difference in frequency between the observed number of AB pairs and the expected number:
     D = p11*p22-p12*p21
     pa2=1-pa1
     pb2=1-pb1
     if pa1==0 or  pa2==0 or pb1==0 or  pb2==0 :
        return (Na, Na, pa1,pb1)
     if D> 0 :
        ## Dmax = min( p(A)p(b), p(a)p(B) ) 
        Dmax = min(pa1*pb2, pa2*pb1)
     else :
        #Dmax = max( -p(A)p(B), -p(a)p(b) ) 
        Dmax = min(-pa1*pb1, -pa2*pb2)
     Dprime=D/Dmax
     Rcar=D/((pa1*pa2*pb1*pb2)**0.5)
     Rcar=Rcar**2
     return (Dprime,Rcar,pa1,pb1, n)

## prevu pour des ge
## Snp1 : distribution du snp 1 pour la pop : 0 donnee manquante, 1 Allele1 2 allele2
def compute_geno_bin(Snp1,Snp2, MinFreqBon=0, Na="NA") :
    NbIndI=float(len(Snp1))
    Cmt=0
    Dist=[[0,0,0],[0,0,0],[0,0,0]]
    ##
    for All in [2,3]:
      Cmt=1
      while Cmt<NbIndI:
         Dist[Snp1[Cmt][All]][Snp2[Cmt][All]]+=1
         Cmt+=1
    NBon=Dist[1][1] + Dist[1][2] + Dist[2][1] +Dist[2][2]
    PBon=NBon/NbIndI
    if PBon>=MinFreqBon :
       NbInd=float(NBon)
       pa1=(Dist[1][1]+Dist[1][2])/NbInd
       pb1=(Dist[1][1]+Dist[2][1])/NbInd
       if pa1==0 or pb1==0 or pa1==1 or pb1==1:
          return (Na, Na, pa1, pb1, NBon)
       p11=Dist[1][1]/NbInd
       p12=Dist[1][2]/NbInd
       p22=Dist[2][2]/NbInd
       p21=Dist[2][1]/NbInd
       return compute_ld_pos_pos(p11,p22,p12,p21 ,pa1,pb1, NBon, Na=Na)
    return (Na, Na, Na, Na, NBon)

