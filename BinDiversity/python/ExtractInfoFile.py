import sys

def ValueGeno(a):
    if a=='.':
       return 0
    elif a=='0':
       return 1
    elif a=='1':
       return 2
    else :
       return None

## Permet pour une position donnee, a une pop donne de recuperer les genotypes
def GetGenoPositionPopHaploide(LigneVCFSplit, ParcoursInd, PosAncetre=None) :
    Geno=[0]*3
    if PosAncetre==None :
        for Cmt in ParcoursInd :
            A1=ValueGeno(LigneVCFSplit[Cmt][0])
            A2=ValueGeno(LigneVCFSplit[Cmt][2])
            ### Cas ou on enregistre que les homozygotes
            try :
               if A1==A2 :
                  Geno[A1]+=1
               else :
                  Geno[0]+=1
            except :
               sys.exit("Problem with genotyp : "+ str(Cmt)+" "+LigneVCFSplit[Cmt]+"\n"+"\t".join(LigneVCFSplit))
    else :
        AlleleAncetre=ValueGeno(LigneVCFSplit[PosAncetre][0])
        AlleleAncetre2=ValueGeno(LigneVCFSplit[PosAncetre][2])
        if AlleleAncetre!=AlleleAncetre2 or AlleleAncetre==None or AlleleAncetre==0:
           for Cmt in ParcoursInd :
               Geno[0]+=1
           ## Allele ancestrale => 1 
        else :
           for Cmt in ParcoursInd :
               A1=ValueGeno(LigneVCFSplit[Cmt][0])
               A2=ValueGeno(LigneVCFSplit[Cmt][2])
               if A1==A2 :
                  if A1==0 :
                     A1=0
                  elif A1==AlleleAncetre:
                     A1=1
                  else :
                     A1=2
               else :
                     A1=0
               Geno[A1]+=1
    return Geno

## NomInd NomPop
### on recupere les informations par population
def GetPop(FilePop):
   LirePop=open(FilePop)
   DicInd={}
   DicPop={}
   DicPopByIndic={}
   CmtPop=0
   DicPopPourInd={}
   NomExt=None
   for Ligne in LirePop :
       SplitLigne=Ligne.split()
       Pop=SplitLigne[1].replace('\n','')
       if Pop.upper() == "OUT":
          NomExt=SplitLigne[0].upper()
       else :      
          if Pop not in DicPop :
             DicPop[Pop]=CmtPop
             DicPopByIndic[CmtPop]=Pop
             DicPopPourInd[Pop]=[]
             CmtPop+=1
          DicPopPourInd[Pop].append(SplitLigne[0].upper())
   return {'Pop':DicPop, 'Ind':DicInd, 'PopIndic':DicPopByIndic, 'ListeIndParPop':DicPopPourInd, 'OutGroup':NomExt}

### on recupere les noms dans le VCF
def GetNomVCF(LireVCF,DicInfoPop,GetLigneVCF=False ) :
    Balise=True
    InfoPop=DicInfoPop['Pop']
    InfoPopParInd=DicInfoPop['ListeIndParPop']
    IndPop=[[] for x in range(len(InfoPop))]
    NomExt=DicInfoPop['OutGroup']
    PosExterne=None
    while(Balise):
       Ligne=LireVCF.readline()
       if Ligne[0:6]=="#CHROM":
          SplitLigne=Ligne.split()
          Cmt =9
          for Ind in SplitLigne[9::]:
              Ind=Ind.replace('\n','')
              if Ind.upper() == NomExt :
                 PosExterne=Cmt
              else :
                  for Pop in InfoPopParInd.keys() :
                     if Ind.upper() in InfoPopParInd[Pop] :
                        IndPop[InfoPop[Pop]].append(Cmt)
              Cmt+=1
          Balise=False
    if GetLigneVCF : 
       return (IndPop, PosExterne, SplitLigne)
    else :
        return [IndPop, PosExterne]



