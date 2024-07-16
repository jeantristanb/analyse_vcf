#!/usr/bin/python
# -*- coding: utf-8 -*-
## pour fusionner deux VCF meme chro, meme pos en 1=> avec individu different
## Pour varscan :)
import sys

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

def ConcatenerDeuxVCF(File1, File2) :
    Lire1=open(File1)
    Lire2=open(File2)
    (Entete1,NomInd1)=GetEntete(Lire1)
    (Entete2,NomInd2)=GetEntete(Lire1)


File1=sys.argv[1]
File2=sys.argv[2]
TmpBaliseEntete=sys.argv[3]
if TmpBaliseEntete=='T':
   BaliseEntete=True
else :
   BaliseEntete=False
if len(sys.argv)==5:
   FileSortie=sys.argv[4]
   Ecrire=open(FileSortie,'w')
else :
   Ecrire=None

Dico1=LireVCFBrut(File1)
Dico2=LireVCFBrut(File2)
ConcatenerVCF(Dico1,Dico2, Ecrire, MaxAlt=1)

