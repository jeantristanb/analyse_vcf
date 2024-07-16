#!/usr/bin/python
## Version 1.2
## Modifier de la version 1.1 de ~/bin/ConvertSamPE-BUOU.py 
## Ajout d un  nombre de PE ou le fichier d ecriture est enleve
import fileinput
import sys


###Lecture de VCF a plusieurs Individu
def OuvrirFichier(Fichier, Type='r'):
    try :
       if Type =='r' :
           LireFich=open(Fichier)
       else :
           LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " in mode " + Type + "can't OuvrirFichier")
    return LireFich




def Rev(x) :
   if x=='A':
      return 'T'
   if x=='T':
      return 'A'
   if x=='C':
      return 'G'
   if x=='G':
      return 'C'
   return x


def RevSeq(Sequence) :
    Seq=""
    for x in Sequence :
      Seq=Rev(x)+Seq
    return Seq


def PrintFormatIllumina(EcrireSortie,splitline, FlagSam):
   #splitline=line.split()
   EcrireSortie.write("@"+splitline[0]+"\n")
   if FlagSam == '89' or FlagSam=='153':
       EcrireSortie.write(RevSeq(splitline[9])+"\n")
   else :
       EcrireSortie.write(splitline[9]+"\n")
   EcrireSortie.write("+\n")
   EcrireSortie.write(splitline[10]+"\n")


if len(sys.argv)!=5 :
   print "Exe SortieR1 SortieR2 SAMFILE NbLigne"
   sys.exit()

### Nombre de sequence qui sera contenue dans chaque fichier
NbSequence=int(sys.argv[4])
SortieR1=sys.argv[1]
SortieR2=sys.argv[2]


ecrireSam=OuvrirFichier(sys.argv[3], 'w')


Flag1=['73','89','121','69','77','101','117']
Flag2=['133','165','181','137','141', '153','185']
#Flag1=['73','89','121','69','77','101','117', '133','165','181','137','141', '153','185']
cmt=1
BaliseEcrire=False
BalHeader=True
Cmt=0
if BalHeader==True :
    for line in sys.stdin:
       if line[0]=='@':
           break
       Cmt+=1
print Cmt

if BalHeader==True :
    ecrireSam.write(line)

CmtPEEnregistre=0
CmtFichier=0
EcrireR1=OuvrirFichier(SortieR1+"_"+str(CmtFichier)+"_PE.txt", 'w')
EcrireR2=OuvrirFichier(SortieR2+"_"+str(CmtFichier)+"_PE.txt", 'w')
CmtFichier+=1

for line in sys.stdin:
    if line[0]!='@' or BaliseEcrire:
       if cmt%2==0: 
          if splitline1[1] in Flag1:
               splitline2=line.split()
               if splitline1[0] != splitline2[0] :
                   sys.exit("Nom PE different entre 1 et 2 : "+splitline1[0]+"\t"+splitline2[0])
               PrintFormatIllumina(EcrireR2,splitline2, splitline2[1])
               PrintFormatIllumina(EcrireR1,splitline1, splitline1[1])
               CmtPEEnregistre+=1
          elif splitline1[1] in Flag2 :
               splitline2=line.split()
               if splitline1[0] != splitline2[0] :
                   sys.exit("Nom PE different entre 1 et 2 : "+splitline1[0]+"\t"+splitline2[0])
               PrintFormatIllumina(EcrireR1,splitline2, splitline2[1])
               PrintFormatIllumina(EcrireR2,splitline1, splitline1[1])
               CmtPEEnregistre+=1
          else :
               ecrireSam.write(line)
               ecrireSam.write(line1)
          if CmtPEEnregistre==NbSequence:
              EcrireR1.close()
              EcrireR2.close()
              EcrireR1=OuvrirFichier(SortieR1+"_"+str(CmtFichier)+"_PE.txt", 'w')
              EcrireR2=OuvrirFichier(SortieR2+"_"+str(CmtFichier)+"_PE.txt", 'w')
              CmtPEEnregistre=0
              CmtFichier+=1
       else :
          line1=line
          splitline1=line.split()
       BaliseEcrire=True
       cmt+=1
    else :
       ecrireSam.write(line)

EcrireR1.close() 
EcrireR2.close() 
ecrireSam.close()
