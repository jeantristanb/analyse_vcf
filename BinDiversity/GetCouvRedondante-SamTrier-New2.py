#!/usr/bin/python
# -*- coding: utf8 -*-
import sys
import string

def GetUnArgument(ListeArgv, opt):
    BaliseSave=False
    Argc=len(ListeArgv)
    cmt=0
    ListeArg=[]
    while cmt<Argc:
        if opt==ListeArgv[cmt] and BaliseSave==False:
           BaliseSave=True
        elif BaliseSave and ListeArgv[cmt][0] == '-' :
           return ListeArg
        elif BaliseSave :
           ListeArg.append(ListeArgv[cmt])
        cmt+=1
    return ListeArg

def OuvrirFichier(Fichier, OptionOuverture='r'):
   try :
      if OptionOuverture=='r':
          LE=open(Fichier)
      else :
          LE=open(Fichier,OptionOuverture)
   except :
      sys.exit("File "+Fichier+" with " + OptionOuverture + "Can't open")
   return LE



ListeFlag=set(['I', 'M', 'D', 'N'])
def GetInfoFlag(FlagSam):
   Chaine=""
   Longueur=0
   for x in FlagSam :
      if x in ListeFlag :
          if x!='I':
              Longueur+=int(Chaine)  
          Chaine=""
      else :
          Chaine+=x
   return Longueur

def GetNbBase(ChaineCigar) :
    if ChaineCigar=='101M' or ChaineCigar=='=':
       return 101
    else :
       return GetInfoFlag(ChaineCigar) 

def ChroCorrInt(Chro):
   if Chro=='chromosome:AGPv2:10:1:150189435:1':
       return 10
   if Chro=='chromosome:AGPv2:1:1:301354135:1':
       return 1
   if Chro=='chromosome:AGPv2:2:1:237068873:1':
       return 2
   if Chro=='chromosome:AGPv2:3:1:232140174:1':
       return 3
   if Chro=='chromosome:AGPv2:4:1:241473504:1':
       return 4
   if Chro=='chromosome:AGPv2:5:1:217872852:1':
       return 5
   if Chro=='chromosome:AGPv2:6:1:169174353:1':
       return 6
   if Chro=='chromosome:AGPv2:7:1:176764762:1':
       return 7
   if Chro=='chromosome:AGPv2:8:1:175793759:1':
       return 8
   if Chro=='chromosome:AGPv2:9:1:156750706:1':
       return 9
   if Chro=='chromosome:AGPv2:mitochondrion:1:569630:1':
       return 11
   if Chro=='chromosome:AGPv2:chloroplast:1:140384:1':
       return 12
   if Chro=='chromosome:AGPv2:UNKNOWN:1:7140151:1':
       return 13


def GetTaille(Chro):
   if Chro==10:
       return 150189435
   if Chro==1:
       return 301354135
   if Chro==2:
       return 237068873
   if Chro==3:
       return 232140174
   if Chro==4:
       return 241473504
   if Chro==5:
       return 217872852
   if Chro==6:
       return 169174353
   if Chro==7:
       return 176764762
   if Chro==8:
       return 175793759
   if Chro==9:
       return 156750706
   if Chro==11:
       return 569630
   if Chro==12:
       return 140384
   if Chro==13:
       return 7140151


    

def GetValueXA(splitligne):
    try :
        for Info in splitligne[11::]:
           if Info[0:3]=='XS:': 
              return True
    except :
      sys.exit(splitligne +Info)
    
def IsRedondantScoreStampy(splitligne1, splitligne2, LimScore):
    Score1=int(splitligne1[4])
    Score2=int(splitligne2[4])
    if Score1<LimScore and Score2<LimScore :
        return True
    return False


def IsRedondant2(splitligne1, splitligne2, ScoreBowtie=7, ScoreStampy=11):
    XS1=100
    XS2=100
    AS1=100
    AS2=100
    try :
        for Info in splitligne1[11::]:
           if Info[0:3]=='XS:': 
              XS1=int(Info[5::])
           if Info[0:3]=='AS:': 
              AS1=int(Info[5::])
        if AS1==100 :
           return IsRedondantScoreStampy(splitligne1, splitligne2, ScoreStampy)       
        if XS1==100 :
             return False
        for Info in splitligne2[11::]:
           if Info[0:3]=='XS:': 
              XS2=int(Info[5::])
           if Info[0:3]=='AS:': 
              AS2=int(Info[5::])
        if XS2==100 :
             return False
        if AS1==100 or AS2==100 :
            sys.exit("AS not find")
        if (AS1-XS1) + (AS2-XS2)< ScoreBowtie :
            return True
        else :
            return False
    except :
      sys.exit(splitligne +Info)
    return False


def IsRedondant(splitligne):
    try :
        for Info in splitligne[11::]:
           if Info[0:3]=='XS:': 
              return True
    except :
      sys.exit(splitligne +Info)
    return False

#AppendRedondantPosFrwd(splitligne, ListeNonRed1, ListeNonRed2,PosI1,PosI2)
def AppendRedondantPosFrwd(splitligne, ListeNonRed1, ListeNonRed2, PosI1, PosI2):
    Chro=ChroCorrInt(splitligne[2])
    Pos=int(splitligne[3])
    NbBase=GetNbBase(splitligne[5])
    MaxPos=Pos+NbBase
    for CmtPos in range(Pos-PosI1,min(PosI2,MaxPos)-PosI1):
        ListeNonRed1[CmtPos]+=1
    for CmtPos in range(max( PosI2, Pos)-PosI2, max(PosI2, MaxPos)-PosI2):
        ListeNonRed2[CmtPos]+=1

def AppendCouvFrwd(splitligne):
    NbBase=GetNbBase(splitligne[5])
    Chro=ChroCorrInt(splitligne[2])
    Pos=int(splitligne[3])
    MaxPos=Pos+NbBase
    for CmtPos in range(Pos,MaxPos):
        InfoBase[Chro][CmtPos]+=1


def AppendRedondantPosRev(splitligne):
    NbBase=GetNbBase(splitligne[5])
    Chro=ChroCorrInt(splitligne[2])
    Pos=int(splitligne[3])
    MinPos=Pos-NbBase+1
    for CmtPos in range(MinPos,Pos+1):
        InfoBase[Chro][CmtPos]+=1
        InfoBaseDup[Chro][CmtPos]+=1

def AppendCouvRev(splitligne):
    NbBase=GetNbBase(splitligne[5])
    Chro=ChroCorrInt(splitligne[2])
    Pos=int(splitligne[3])
    MinPos=Pos-NbBase+1
    for CmtPos in range(MinPos,Pos+1):
        InfoBase[Chro][CmtPos]+=1

#AfficherDataAll(ListeNonRed1, PosI1, Chro)
def AfficherDataAll(ListeRed,ListeNonRed, PosI, Chro, Ecrire):
    Chaine=""
    Cmt=0
    for Cmt in range(0, len(ListeNonRed)):
        Chaine+=(Chro+" "+str(Cmt+PosI)+" "+str(ListeNonRed[Cmt]+ListeRed[Cmt]) + " "+ str(ListeRed[Cmt])+"\n")
    Ecrire.write(Chaine)

def NettoyerVect(Liste):
    for Cmt in range(0, len(Liste)):
        Liste[Cmt]=0
    return Liste

def ReadInfoSnp(Fichier):
    Lire = open(Fichier)
    Dico={}
    for Ligne in Lire :
       if Ligne[0]!="#" :
          splitligne=Ligne.split()
          if splitligne[0] not in Dico:
             Dico[splitligne[0]]=[]
          Dico[splitligne[0]].append(int(splitligne[1]))
    for Chro in Dico.keys():
        Dico[Chro].sort()
    return Dico

def ReadInfoSnpPosDebFin(Fichier) :
    Lire = open(Fichier)
    Dico={}
    for Ligne in Lire :
       if Ligne[0]!="#" :
          splitligne=Ligne.split()
          if splitligne[0] not in Dico:
             Dico[splitligne[0]]=[]
          Dico[splitligne[0]].append([int(splitligne[1]), int(splitligne[2])])
    for Chro in Dico.keys():
        Dico[Chro][0].sort()
        Dico[Chro][1].sort()
    return Dico




def EcrireListeSnp(ListeNonRed, ListeRed, PosI, Pas, Chro, DicoListe, Ecrire) :
    ListeSnp=DicoListe[Chro]
    MinPos=PosI
    MaxPos=PosI+Pas
    ListeSnp2=list(ListeSnp)
    for Snp in ListeSnp2:
        if Snp >= MaxPos :
           return 0
        if Snp>= MinPos :
           PosListe=Snp-PosI
           Ecrire.write(Chro+" "+str(Snp)+" "+str(ListeNonRed[PosListe]+ListeRed[PosListe])+" "+ str(ListeRed[PosListe])+"\n")
           DicoListe[Chro].remove(Snp)
       
def EcrireListeSnpPosDebutFin(ListeNonRed, ListeRed, PosI, Pas, Chro, DicoListe, Ecrire) :
    MinPos=PosI
    ## Maximum postion of data analysed
    MaxPos=PosI+Pas
    ListeWrite=set()
    ListePortionGeneToRemove=[]
    for PortionGene in range(len(DicoListe[Chro])) :
        Debut=DicoListe[Chro][PortionGene][0]
        Fin=DicoListe[Chro][PortionGene][1]
        for Snp in range(Debut, Fin+1) :
           if Snp >= MaxPos :
              for Remove in ListePortionGeneToRemove :
                  if Remove in DicoListe[Chro] :
                     DicoListe[Chro].remove(Remove)
              return 0
           if Snp>= MinPos :
              PosListe=Snp-PosI
              if PosListe not in ListeWrite :
                 Ecrire.write(Chro+" "+str(Snp)+" "+str(ListeNonRed[PosListe]+ListeRed[PosListe])+" "+ str(ListeRed[PosListe])+"\n")
                 ListeWrite.add(PosListe)
        if MaxPos>Fin :
           ListePortionGeneToRemove.append([Debut,Fin])
    for Remove in ListePortionGeneToRemove :
        if Remove in DicoListe[Chro] :
           DicoListe[Chro].remove(Remove)

FichierListePos=GetUnArgument(sys.argv, '-i')
ScoreBowtie=GetUnArgument(sys.argv, '-sb')
ScoreStampy=GetUnArgument(sys.argv, '-st')
FichierOut=GetUnArgument(sys.argv, '-o')
OptBaliseHeader=GetUnArgument(sys.argv, '-h')
if len(FichierListePos)>1 or len(FichierOut)!=1 or len(ScoreBowtie)!=1 or len(ScoreStampy)!=1:
    sys.exit("Sam Info  | Exec -i ListePos FILE -o OutFile -st [int] -sb [int] -h T/F")

if len(OptBaliseHeader)==0:
    BalHeader=True

TypeFichierPos="Pos"


if len(FichierListePos)==0:
   BalPrintAll=True
   InfoPos=None
else :
   FichierListePos=FichierListePos[0]
   BalPrintAll=False
   if TypeFichierPos=="Pos" :
      InfoPos=ReadInfoSnp(FichierListePos)  
      FoncEcrireSnp=EcrireListeSnp
   else:
      InfoPos=ReadInfoSnpPosDebFin(FichierListePos)
      FoncEcrireSnp=EcrireListeSnpPosDebutFin

if len(OptBaliseHeader)==0:
    BalHeader=True
else :
    if OptBaliseHeader[0]=='T':
        BalHeader=True
    elif  OptBaliseHeader[0]=='F':
        BalHeader=False
    else :
        sys.exit(" Balise Header non reconnue T or F\nSam Info  | Exec -i ListePos FILE -o OutFile -st [int] -sb [int] -h T/F")
       
ScoreBowtie=int(ScoreBowtie[0])
ScoreStampy=int(ScoreStampy[0])
#ScoreBowtie=6
#ScoreStampy=11
FichierOut=FichierOut[0]
BaliseEcrire=False



Cmt=0
if BalHeader==True :
    for line in sys.stdin:
       if line[0]=='@':
           break
       Cmt+=1
       if Cmt > 150 :
          sys.exit("Erreur nb header ne peut depasser 150")

PasPourSave=500000

ListeRed1=[0]*PasPourSave
ListeNonRed1=[0]*PasPourSave
ListeRed2=[0]*PasPourSave
ListeNonRed2=[0]*PasPourSave
PosI1=1
PosI2=PosI1+PasPourSave

print "Fin Allocation memory et debut lecture "
EcrireSortie=OuvrirFichier(FichierOut ,'w')
#Ecrire2=open('Test.sort.byname.sam', "w")

cmtTot=0
CmtBon=0
ListePE=set([])
Dico={}
CmtRed=0
for line in sys.stdin:
    if line[0]!='@' :
       splitligne=line.split()
       if cmtTot%10000000==0:
          print cmtTot
       if len(splitligne)>5 and(splitligne[1]=='99' or splitligne[1]=='83' or splitligne[1]=='163' or splitligne[1]=='147'):
           if CmtBon==0 :
              PosI1=int(splitligne[7]) 
              PosI2=PosI1+PasPourSave
           if splitligne[0] not in ListePE :
              ListePE.add(splitligne[0])
              Dico[splitligne[0]]=splitligne
              if CmtBon==0:
                 Chro = splitligne[2]
           else :
               ListePE.remove(splitligne[0])
               splitligne2=Dico[splitligne[0]]
               del Dico[splitligne[0]]
               Pos1=int(splitligne[7])
               Pos2=int(splitligne2[7])
               MaxPos=max(Pos1, Pos2)
               MinPos=min(Pos1, Pos2)
               if (MaxPos >= PosI2+PasPourSave and MinPos >= PosI2+PasPourSave) or Chro!=splitligne[2]:
                 if BalPrintAll :
                    AfficherDataAll(ListeRed1,ListeNonRed1, PosI1, Chro, EcrireSortie)
                    AfficherDataAll(ListeRed2,ListeNonRed2, PosI2, Chro, EcrireSortie)
                 else :
                    FoncEcrireSnp(ListeNonRed1, ListeRed1, PosI1, PasPourSave, Chro, InfoPos, EcrireSortie)
                    FoncEcrireSnp(ListeNonRed2, ListeRed2, PosI2, PasPourSave, Chro, InfoPos, EcrireSortie)
                 PosI1=MinPos
                 PosI2=PosI1+PasPourSave
                 Chro=splitligne[2]
                 ListeNonRed1=NettoyerVect(ListeNonRed1)
                 ListeNonRed2=NettoyerVect(ListeNonRed2)
                 ListeRed1=NettoyerVect(ListeRed1)
                 ListeRed2=NettoyerVect(ListeRed2)
               elif MinPos>PosI2+5000 :
                 if BalPrintAll :
                     AfficherDataAll(ListeRed1,ListeNonRed1, PosI1, Chro, EcrireSortie) ## uniquement pour 1
                 else :
                      FoncEcrireSnp(ListeNonRed1, ListeRed1, PosI1, PasPourSave, Chro, InfoPos, EcrireSortie)
                 Temp=ListeNonRed1
                 ListeNonRed1=ListeNonRed2
                 ListeNonRed2=NettoyerVect(Temp)
                 Temp=ListeRed1
                 ListeRed1=ListeRed2
                 ListeRed2=NettoyerVect(Temp)
                 PosI1=PosI2
                 PosI2=PosI1+PasPourSave
               if IsRedondant2(splitligne,splitligne2, ScoreBowtie, ScoreStampy):
                   CmtRed+=1
                   AppendRedondantPosFrwd(splitligne, ListeRed1, ListeRed2,PosI1,PosI2)
                   AppendRedondantPosFrwd(splitligne2, ListeRed1, ListeRed2, PosI1, PosI2)
               else :
                   AppendRedondantPosFrwd(splitligne, ListeNonRed1, ListeNonRed2,PosI1,PosI2)
                   AppendRedondantPosFrwd(splitligne2, ListeNonRed1, ListeNonRed2, PosI1, PosI2)
           CmtBon+=1
       cmtTot+=1


if BalPrintAll :
  AfficherDataAll(ListeRed1,ListeNonRed1, PosI1, Chro, EcrireSortie)
  AfficherDataAll(ListeRed2,ListeNonRed2, PosI2, Chro, EcrireSortie)
else :
  FoncEcrireSnp(ListeNonRed1, ListeRed1, PosI1, PasPourSave, Chro, InfoPos, EcrireSortie)
  FoncEcrireSnp(ListeNonRed2, ListeRed2, PosI2, PasPourSave, Chro, InfoPos, EcrireSortie)


#ecrire=OuvrirFichier(FichierOut, 'w')
#Chaine=""
#Cmt=1
#for ligne in LireListePos:
#    if ligne[0]!='#':
#       splitligne=ligne.split()
#       Chro=ChroCorrInt(splitligne[0])
#       Pos=int(splitligne[1])
#       Chaine+=(splitligne[0]+" "+str(Pos)+" "+str(InfoBase[Chro][Pos])+" "+str(InfoBaseDup[Chro][Pos])+"\n")
#       if Cmt%10000==0 :
#          ecrire.write(Chaine)   
#          Chaine=""
#       Cmt+=1
#ecrire.write(Chaine)
#ecrire.close()

#LireListePos.close()
