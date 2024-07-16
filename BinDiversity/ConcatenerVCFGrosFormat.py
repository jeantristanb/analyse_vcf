#!/usr/bin/python
# -*- coding: utf-8 -*-
## pour fusionner deux VCF meme chro, meme pos en 1=> avec individu different
## Pour varscan :)
import sys
import gzip

def intersect(b1,b2):
    return list(set(b1).intersection( set(b2) ))
def GetEntete(Lire):
    CmtLigne=0
    Entete=""
    for Ligne in Lire :
        if Ligne[0:6]=="#CHROM" :
           return (Entete,Ligne)
        if Ligne[0]!='#' :
           sys.exit('Pattern #CHROM not found')
        Entete+=Ligne
        CmtLigne+=1
    sys.exit("FileWithout Information")


def LireVCFBrut(File):
    Lire=open(File)
    (Entete,NomInd)=GetEntete(Lire)
    Dico={'Entete':Entete,'NomInd':NomInd.split()[9::],'InfoLocus':[]}
    InfoAll={}
    for Ligne in Lire :
        SplitLigne=Ligne.split()
        if SplitLigne[0] not in InfoAll :
           InfoAll[SplitLigne[0]]={}
        NbPoly=SplitLigne[7].split(';')
        ## on recupere le nombre d'ADP, d'het
        #Info=[int(x.split('=')[1]) for x in NbPoly]
        Info=None
         
        InfoAll[SplitLigne[0]][int(SplitLigne[1])]=[SplitLigne[2:9],Info,"\t".join(SplitLigne[9::])]
    Dico['InfoLocus']=InfoAll
    return Dico

def GetAltInfo(Alt1,Alt2):
    ListeAlt1=Alt1.split(',') 
    ListeAlt2=Alt2.split(',') 
    ListeAlt1.extend(ListeAlt2)
    ListeF=list(set(ListeAlt1))
    if '.' in ListeF :
      ListeF.remove('.')
    if len(ListeF)==0 :
       ListeF=['.']
    return ListeF

def CheckRef(Ref1,Ref2) :
    if Ref1!=Ref2:
       return False
    return True
### INFO pour vascan
##ADP=345;WT=1;HET=0;HOM=0;NC=0
def GetInfo(Info1, Info2): 
    SumNbInd1=float(sum(Info1[1::]))
    SumNbInd2=float(sum(Info2[1::]))
    ADP="ADP="+str(int((Info1[0]*SumNbInd1+ Info2[0]*SumNbInd2)/(SumNbInd1+SumNbInd2)))
    #ADP=345;WT=1;HET=0;HOM=0;NC=0
    WT="WT="+str(Info1[1]+Info2[1])
    HET="WT="+str(Info1[2]+Info2[2])
    HOM="HOM="+str(Info1[3]+Info2[3])
    NC="NC="+str(Info1[4]+Info2[4])
    return ADP+";"+WT+";"+HET+";"+HOM+";"+NC
      
def ConcatenerVCF(Dico1,Dico2, Ecrire=None, MaxAlt=1, BaliseEntete=True):
    EcrireErreur=sys.stderr
    if Ecrire==None :
       Ecrire=sys.stdout
    if BaliseEntete :
       Ecrire.write(Dico1['Entete'] )
       Ecrire.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	"+"\t".join(Dico1['NomInd'])+"\t"+"\t".join(Dico2['NomInd'])+"\n")
    ChroCommun=intersect(Dico1['InfoLocus'].keys(), Dico2['InfoLocus'].keys())
    for Chro in ChroCommun :
        ## On recupere les informations des locus
        Chro1=Dico1['InfoLocus'][Chro]
        Chro2=Dico2['InfoLocus'][Chro]
        NaValueChr1="\t".join(["./."]*len(Chro1[Chro1.keys()[0]][2].split()))
        NaValueChr2="\t".join(["./."]*len(Chro2[Chro2.keys()[0]][2].split()))
        ## Positions communnes
        #PosCommun=intersect(Chro1.keys(), Chro2.keys())
        PosAll=Chro1.keys()
        PosAll.extend(Chro2.keys())
        PosAll=list(set(PosAll))
        PosAll.sort()
        for Pos in PosAll :
           if Pos in Chro1 and Pos in Chro2 :
              PtrPos1=Chro1[Pos]
              PtrPos2=Chro2[Pos]
              ListeAlt=GetAltInfo(PtrPos1[0][2],PtrPos2[0][2])
              INFO="." #GetInfo(PtrPos1[1], PtrPos2[1])
              ### on verifie que la position avec fusion n'a pas plus de de Max alternatif et que les references sont les memes
              if len(ListeAlt)<=MaxAlt and CheckRef(PtrPos1[0][1],PtrPos2[0][1]):
                 if len(ListeAlt)==1:
                    Alt=ListeAlt[0]
                 else :
                    Alt=",".join(ListeAlt)
                 Ecrire.write(Chro + "\t" + str(Pos)+"\t"+".\t"+PtrPos1[0][1]+"\t"+",".join(ListeAlt)+"\t.\tPASS\t"+INFO +"\t"+PtrPos1[0][6]+"\t"+PtrPos1[2]+ "\t"+PtrPos2[2] +"\n")
              else : 
                  EcrireErreur.write('Locus False : Alt1 : '+PtrPos1[0][2]+"\t"+PtrPos2[0][2]+'\tRef 1: '+PtrPos1[0][1]+'\t Ref 2 '+PtrPos2[0][1]+"\n")
           elif Pos in Chro1 :
              PtrPos1=Chro1[Pos]
              Alt=PtrPos1[0][2]
              INFO="." 
              Ecrire.write(Chro + "\t" + str(Pos)+"\t"+".\t"+PtrPos1[0][1]+"\t"+Alt+"\t.\tPASS\t"+INFO +"\t"+PtrPos1[0][6]+"\t"+PtrPos1[2]+ "\t"+NaValueChr2+"\n")
           elif Pos in Chro2 :
              PtrPos2=Chro2[Pos]
              Alt=PtrPos2[0][2]
              INFO="." 
              Ecrire.write(Chro + "\t" + str(Pos)+"\t"+".\t"+PtrPos2[0][1]+"\t"+Alt+"\t.\tPASS\t"+INFO +"\t"+PtrPos2[0][6]+"\t"+"\t"+NaValueChr1+"\t"+PtrPos2[2]+"\n")

def GetInfoVCF(Lire) :
    Balise=True 
    while Balise :
      Ligne=Lire.readline()
      if Ligne=="" :
        Balise=False
      elif Ligne[0]!="#":
          Balise=False
    try :
       if Ligne!="" :
         SplitLigne=Ligne.split()
         Chro=SplitLigne[0]
         Pos=int(SplitLigne[1])
         Ref=SplitLigne[3]
         Alt=SplitLigne[4]
         Info=SplitLigne[9::]
         return (True, Chro, Pos, (Ref, Alt, Info))
       else :
         return (False, None, None, (None, None,None))
    except :
       sys.exit(Ligne)

def TransformValue(Ref, InfoAlt,InfoInd):
    AllAltern=[".",Ref]
    AllAltern.extend(InfoAlt.split(","))
    NewListe1=[]
    NewListe2=[]
    for Ind in InfoInd :
       All1=Ind[0]
       All2=Ind[2]
       if All1 == '.':
         NewListe1.append(All1) 
         NewListe2.append(All2) 
       else :
         NewListe1.append(AllAltern[int(All1)+1]) 
         NewListe2.append(AllAltern[int(All2)+1]) 
    return (AllAltern, NewListe1, NewListe2)

def ConcatenerLigneVCF(Chro,Pos, Info1, Info2) :
   ## On check les references 
   if  CheckRef(Info1[0],Info2[0]) ==False :
        print "warning : "+Chro+" "+ str(Pos) + " reference varie "+ Info1[0] +" "+ Info2[0]
        return (False,"") 
   ## On regarde les alternatifs 
   ListeAlt=[Info1[0]]
   (AllAltern1, NewListePop1A1, NewListePop1A2)=TransformValue(Info1[0], Info1[1],Info1[2])
   (AllAltern2, NewListePop2A1, NewListePop2A2)=TransformValue(Info2[0], Info2[1],Info2[2])
   ListeAlt=AllAltern1[2::] 
   ListeAlt.extend(AllAltern2[2::] )
   NewListeA1=NewListePop1A1
   NewListeA1.extend(NewListePop2A1)
   NewListeA2=NewListePop1A2
   NewListeA2.extend(NewListePop2A2)
   DicAlt={'.':'.'}
   DicAlt[Info1[0]]='0'
   Cmt=1
   NewListeAlt=[]
   for Alt in ListeAlt :
      if Alt not in DicAlt :
         DicAlt[Alt]=str(Cmt)
         Cmt+=1
         NewListeAlt.append(Alt)
   NewListeAlt=",".join(NewListeAlt)
   FinalListe=[] 
   for Cmt in range(0, len(NewListeA1)):
      FinalListe.append(DicAlt[NewListeA1[Cmt]]+"/"+DicAlt[NewListeA2[Cmt]]) 
   FinalLines=Chro+"\t"+str(Pos)+  "\t.\t"+Info1[0]+"\t"+NewListeAlt+"\t.\tPASS\t.\tGT\t"+"\t".join(FinalListe)
   return (True,FinalLines)
 

def ConcatenerDeuxVCF(File1, File2, FileOut, TGZ) :
    print File1,File2, FileOut
    if TGZ==False:
       Lire1=open(File1, "r")
       Lire2=open(File2, "r")
    else :
      Lire1=gzip.open(File1)
      Lire2=gzip.open(File2)
    (Entete1,CharNomInd1)=GetEntete(Lire1)
    (Entete2,CharNomInd2)=GetEntete(Lire2)
    Lire1.close()
    Lire2.close()
    #print Entete1 
    #print "Entete1"
    #print Entete2
    if TGZ==False:
       Lire1=open(File1, "r")
       Lire2=open(File2, "r")
    else :
      Lire1=gzip.open(File1)
      Lire2=gzip.open(File2)
    NomInd1=CharNomInd1.split()
    NomInd2=CharNomInd2.split()
    EnteteF="\t".join(NomInd1[0:9])+"\t"
    EnteteF+="\t".join(NomInd1[9:])+"\t"
    EnteteF+="\t".join(NomInd2[9:])
    Ecrire=open(FileOut, 'w')
    Ecrire.write(EnteteF+"\n")
    Balise=True
    Cmt=1
    while Balise :
      if Cmt ==1 :
         (Balise1, Chro1, Pos1, Info1)=GetInfoVCF(Lire1)
         (Balise2, Chro2, Pos2, Info2)=GetInfoVCF(Lire2)
      elif Balise1==False or Balise2==False :
         Balise=False
         print "end of File "
      elif Chro1!=Chro2 :
         sys.exit("Chro1!=Chro2")
      elif Pos1==Pos2 :
          (BaliseVCFPrint, LigneVCF)=ConcatenerLigneVCF(Chro1,Pos1, Info1, Info2)
          (Balise1, Chro1, Pos1, Info1)=GetInfoVCF(Lire1)
          (Balise2, Chro2, Pos2, Info2)=GetInfoVCF(Lire2)
          if BaliseVCFPrint :
            Ecrire.write(LigneVCF+"\n")
      elif Pos1>Pos2 :
         (Balise2, Chro2, Pos2, Info2)=GetInfoVCF(Lire2)
      else :
         (Balise1, Chro1, Pos1, Info1)=GetInfoVCF(Lire1)
      Cmt+=1
         
       


File1=sys.argv[1]
File2=sys.argv[2]
FileSortie=sys.argv[3]
TGZ=True

ConcatenerDeuxVCF(File1,File2,FileSortie, TGZ)
