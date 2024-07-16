#-*- coding: utf-8 -*-
import sys
import os
sys.path.append(os.environ['HOME']+"/bin/python2")
from ExtractInfoFile import GetPop, GetNomVCF

def ValueGeno(a):
    if a=='.':
       return 0
    elif a=='0':
       return 1
    elif a=='1':
       return 2
    else :
       return 0

def GetGeno(Ind,AncestralAll=None) :
    Geno1=ValueGeno(Ind[0])
    Geno2=ValueGeno(Ind[2])
    if Geno1!=Geno2 :
       Geno=0
    else :
       Geno=Geno1
    if AncestralAll==None or Geno==0 or  AncestralAll==0 :
      return Geno
    if AncestralAll==Geno :
      return 1
    else :
      return 2

def GetGFF(FileGFF) :
    LireGFF=open(FileGFF) 
    DataGFF={}
    for Ligne in LireGFF :
        SplitLigne=Ligne.split()
        Chro=SplitLigne[0]
        PosI=int(SplitLigne[3])
        PosF=int(SplitLigne[4])
        if Chro not in DataGFF :
           DataGFF[Chro]=[[] for x in range(6)]
        DataGFF[Chro][0].append(PosI)
        DataGFF[Chro][1].append(PosF)
        DataGFF[Chro][3].append([])
        DataGFF[Chro][4].append("")
        DataGFF[Chro][5].append([])
    return DataGFF

def GetCorrMs(Geno) :
    if Geno==1:
       return '0'
    if Geno==2 :
        return '1'
    return 'N'

def EstDansGFF(DicGFF, Chro, Pos):
    if Chro not in DicGFF :
       return None
    GFFChroDeb=DicGFF[Chro][0]
    GFFChroFin=DicGFF[Chro][1]
    NbElem=len(GFFChroDeb)
    Cmt=0
    while Cmt <NbElem :
      if GFFChroDeb[Cmt]<=Pos and Pos<=GFFChroFin[Cmt] :
         return Cmt
      Cmt+=1
    return None

def PrintMS(DicGFF, FileVCF, FileMS, FileIndMS, InfoNomPosPop, SplitLigneHeader, DicInfoPop, PosAnc=None) : 
    EcrireMS=open(FileMS, 'w')
    EcrireVCF=open(FileVCF, 'w')
    EcrireFileIndMS=open(FileIndMS, 'w')
    ChaineEntete="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    ChaineIndMS=""
    CmtPop=0
    for Pop in InfoNomPosPop:
        for CmtInd in Pop :
            ChaineEntete+="\t"+SplitLigneHeader[CmtInd]
            ChaineIndMS+=SplitLigneHeader[CmtInd]+"\t"+DicInfoPop['PopIndic'][CmtPop]+"\n" 
        CmtPop+=1
    if PosAnc!=None :
       ChaineEntete+="\t"+SplitLigneHeader[PosAnc]
    ChaineEntete+="\n"
    EcrireFileIndMS.write(ChaineIndMS)
    EcrireVCF.write(ChaineEntete)
    for Chro in DicGFF :
        for Cmt in range(len(DicGFF[Chro][1])) :
             EcrireVCF.write("".join(DicGFF[Chro][5][Cmt]))
             ChaineMS="\n//\n"
             if len(DicGFF[Chro][3][Cmt])==0 :
                ChaineMS+="segsites: 0\n"
                ChaineMS+="\n\n"
             else :
                ChaineMS+="segsites: "+str(len(DicGFF[Chro][3][Cmt][0]))+"\n"
                ChaineMS+="positions: "+DicGFF[Chro][4][Cmt]+"\n"
                ChaineMS+="\n".join(DicGFF[Chro][3][Cmt])+"\n"
             EcrireMS.write(ChaineMS)
              
            

def AppendValueInGFF(DicGFF, Chro, Pos, PosInGFF, InfoNomPosPop, SplitVCF, MinFreqBon, PosAnc=None) :
    NbInd=sum([len(x) for x in InfoNomPosPop])
    if len(DicGFF[Chro][3][PosInGFF])==0 :
       DicGFF[Chro][3][PosInGFF]=["" for x in range(NbInd)]
    ListeMsValue=DicGFF[Chro][3][PosInGFF]
    CmtInd=0
    Freq=[0,0,0]
    TempInfo=[]
    AncestralAll=None
    if PosAnc!=None :
       AncestralAll=GetGeno(SplitVCF[PosAnc])
       if AncestralAll==0 :
          return 
    CharVCF="\t".join(SplitVCF[:9])
    for Pop in InfoNomPosPop :
        for PosInd in Pop :
            Ind=SplitVCF[PosInd]
            CharVCF+="\t"+Ind
            Geno=GetGeno(Ind,AncestralAll)
            Freq[Geno]+=1
            TempInfo.append(GetCorrMs(Geno))
    if PosAnc!=None :
       CharVCF+="\t"+SplitVCF[PosAnc]+"\n"
    PBon=sum(Freq[1::])/float(sum(Freq))
    if PBon>= MinFreqBon :
       Pa= Freq[1]/float(sum(Freq[1::]))
       if Pa>0 and Pa<1:
          DicGFF[Chro][5][PosInGFF].append(CharVCF)
          CmtInd=0
          DicGFF[Chro][4][PosInGFF]+=str(Pos)+"\t"
          while CmtInd< NbInd :
               DicGFF[Chro][3][PosInGFF][CmtInd]+=TempInfo[CmtInd]
               CmtInd+=1
       
      

def PrintInMSGFF(DicGFF, FileVCF, InfoNomPosPop, MinFreqBon, Entete, SplitLigneHeader, DicInfoPop, PosAnc=None):
    Lire=open(FileVCF)
    for Ligne in Lire :
        SplitLigne=Ligne.split()
        if Ligne[0]!="#" and len(SplitLigne)>5 and len(SplitLigne[3])==1 and len(SplitLigne[4])==1:
           Chro=SplitLigne[0]
           Pos=int(SplitLigne[1])
           PosInGFF=EstDansGFF(DicGFF, SplitLigne[0], Pos)
           if PosInGFF!=None :
              AppendValueInGFF(DicGFF, Chro, Pos, PosInGFF, InfoNomPosPop, SplitLigne, MinFreqBon, PosAnc)
    PrintMS(DicGFF, Entete+".vcf", Entete+".ms", Entete+".ind", InfoNomPosPop, SplitLigneHeader, DicInfoPop, PosAnc)
              
#python ../ExtractInVCFMs.py $FileInd $FileGFF  $Entete $MinFreqBon FileVCF

Cmd="../ExtractInVCFMs.py $FileInd $FileGFF  $Entete $MinFreqBon ListeFileVCF"
print sys.argv
if(len(sys.argv)<5):
  sys.exit(Cmd)

FilePop=sys.argv[1]
FileGFF=sys.argv[2]
Entete=sys.argv[3]
MinFreqBon=float(sys.argv[4])
ListeNomVCF=sys.argv[5::]
### on recupere les populations
DicInfoPop=GetPop(FilePop)
### on recupere dans les entetes les positions des individus
LireVCF=open(ListeNomVCF[0])
[InfoNomPosPop,PosAncetre, SplitLigneHeader]=GetNomVCF(LireVCF, DicInfoPop, True)
print InfoNomPosPop
print DicInfoPop
LireVCF.close()
print DicInfoPop['PopIndic']
print DicInfoPop['PopIndic'][0]
print PosAncetre


 
DicGFF=GetGFF(FileGFF)
PrintInMSGFF(DicGFF, ListeNomVCF[0], InfoNomPosPop, MinFreqBon, Entete, SplitLigneHeader, DicInfoPop, PosAncetre)




