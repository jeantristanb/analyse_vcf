#!/usr/bin/python
# coding=utf-8
import argparse
import re
import os
import sys
import glob


def Str(Number, Lim=2):
    if type(Number)==float :
       return str(round(Number,Lim))
    return str(Number)

def Concatener(Vecteur, sep="\t"):
    Chaine=""
    for El in Vecteur :
       Chaine+=str(El)+sep
    return Chaine


def GetMatrix(Col, Row):
    return [[0 for x in xrange(Col)] for x in xrange(Row)]
def intersect(liste1, liste2):
    return list(set(liste1) & set(liste2))

def diff(liste1, liste2):
    return list(set(liste1) - set(liste2))


def GetInfoVcf(lines, postosave):
    def getgeno(x,l):
      if x[0]=='.'
        return [None,None]
     # if '/' in x :
     #   p='/'
     # else :
     #   p='|'
      splx=re.split([':|/'])
      return [splx[int(x[0])], splx[int(x[1])]]
    TmpSplit=lines.split()
    Chro=TmpSplit[0]
    Pos=TmpSplit[1]
    NomPos=Chro+"_"+Pos
    Ref=TmpSplit[3]
    ListeAlt=TmpSplit[4].split(',')
    TmpPos=[Ref] 
    TmpPos+=ListAlt
    return (Chro, Pos, Ref, ListAlt,[getgeno(TmpSplitp[x], ListAlt)  for x in postosave ])

def GetHeaderVcf(File):
    VcfRead=OpenFile(File)
    Entete=[]
    for line in Read:
       if line[0]=="#" :
         Entete.append(line.lower().replace("\n",""))
       else :
         Read.close()
         return Entete
    return Entete


def OpenFile(Fichier, Type='r'):
    try :
       LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " OuvrirFichier in mode " + Type + "can't OuvrirFichier")
    return LireFich

## 
def parseArguments():
    parser = argparse.ArgumentParser(description='fill in missing bim values')
    parser.add_argument('--vcf1',type=str,required=True)
    parser.add_argument('--vcf2',type=str,required=True)
    args = parser.parse_args()
    return args


def PosNotFound(Pos, Chro,Vcf):
    return 0

def ComparePos(Allele1,Allele2,Vcf1,Vcf2):
  Allele1+=Allele2
  Auni=list(set(Allele1))
  freqVcf1={} 
  GenoVcf1={}
  freqVcf2={} 
  GenoVcf2={}
  for Al in Auni :
    freqVcf1[Al]=0
    GenoVcf1[Al]={}
    freqVcf2[Al]=0
    GenoVcf2[Al]={}
    for Al2 in Auni :
     GenoVcf1[Al][Al2]=0
     GenoVcf2[Al][Al2]=0
  def IsDiff(Ia1,Ia2, Ib1, Ib2): 
   ## homozyogous case
    if Ia1==Ia2 and Ia1==Ib1 :
      ResumePos[]
      return 1
    if Ia1!=Ia2 and Ib1!=Ib2 :
        if (Ia1==Ib1 and Ia2==Ib2) or  (Ia1==Ib2 and Ia2==Ib1):	
           return 1
    return 0
  if len(Vcf1)!=len(Vcf2):
    sys.exit('len vcf1!=vcf2') 
  Head=["NbVcf1","NbVcf2","Vcf1A1","Vcf2A2","Vcf1NbA1","Vcf2NbA2","Vcf1G1","Vcf2G2","Vcf1NbG1","Vcf2NbG2","NbHom1Hom2Well","NbHom1Hom2False","NbHom1Het2","NbHet1Hom2","NbId"]
  ResumePos=[0]*len(Head)
  PosNbVcf1=0
  PosNbVcf2=1
  PosVcf1A1=2
  PosVcf2A2=3
  PosVcf1NbA1=4
  PosVcf2NbA2=5
  PosVcf1G1=6
  PosVcf2G2=7
  PosVcf1NbG1=8
  PosVcf2NbG2=9
  PosNbHom1Hom2Well=10
  PosNbHom1Hom2False=11
  PosNbHom1Het2=12
  PosNbHet1Hom2=13
  PosNbId=14
  Cmt=0
  for Geno1 in Vcf1
     Geno2=Vcf2[Cmt]
     if Geno1[0]!=None :
       Res[PosNbVcf1]+=1
     if Geno2[0]!=None :
       Res[PosNbVcf2]+=1
     if Geno1[0]!=None and Geno2[0]!=None:
       GenoVcf1[Geno1[0]][Geno1[1]]+=1
       GenoVcf2[Geno2[0]][Geno2[1]]+=1
       freqVcf1[Geno1[0]]+=1
       freqVcf1[Geno1[1]]+=1
       freqVcf2[Geno2[0]]+=1
       freqVcf2[Geno2[1]]+=1
       IsDiff=IsDiff(Geno1[0],Geno1[1],Geno2[0],Geno2[1])
            
          
          
     Cmt+=1
args = parseArguments()
HeadVcf1=GetHeaderVcf(args.vcf1)
HeadVcf2=GetHeaderVcf(args.vcf2)

ListInd1=HeadVcf1[-1].split()
ListInd2=HeadVcf2[-1].split()
lind=intersect(ListInd1[9::], ListInd2[9::])
if len(lind)==0:
   print('common ind between file null')
   print(" ".join(ListInd1[9::])
   print(" ".join(ListInd2[9::])
   sys.exit(2)
posind1=[]
posind2=[]
for ind in posind1:
  posind1.append(ListInd1.index(ind))
  posind2.append(ListInd2.index(ind))

##
lire1=OpenFile(args.vcf1,'r')
for line1 in lire1:
  if line1[0]!='#':
     break
     
lire2=OpenFile(args.vcf2,'r')
for line2 in lire2:
  if line2[0]!='#':
     break
balise=True

(Chro1, Pos1, Ref1, ListAlt1)=GetInfoVcf(lines1, posind1)
(Chro2, Pos2, Ref2, ListAlt2)=GetInfoVcf(lines2, posind2)
while balise:
  if Pos1==Pos2:
    ComparePos(ListAlt1, ListAlt2, StatByInd)
    line1=lire1.readline()  
    line2=lire2.readline()
    if not line1 or not line2 :
     break
    (Chro1,Pos1, Ref1, ListAlt1)=GetInfoVcf(lines1, posind1)
    (Chro2,Pos2, Ref2, ListAlt2)=GetInfoVcf(lines2, posind2)
  elif Pos1>Pos2:
    line2=lire2.readline()
    if not line1 or not line2 :
       break
    (Chro2,Pos2, Ref2, ListAlt2)=GetInfoVcf(lines2, posind2)
  elif Pos1<Pos2:
    line2=lire2.readline()
    if not line1 or not line2 :
       break
