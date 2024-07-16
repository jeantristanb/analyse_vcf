import sys
import glob

def GetFileCor(Fichier):
    LireFichier=open(Fichier)
    DicCorr={}
    for Ligne in LireFichier:
        if Ligne[0]!="#" :
           TmpLigne=Ligne.replace('\n',"").split()
           DicCorr[TmpLigne[0].upper()]=TmpLigne[1]
    return DicCorr

#DicCorChro=GetFileCor(sys.argv[1])
DicCorNom=GetFileCor(sys.argv[1])
LireVCF=open(sys.argv[2])
FichierSortie=sys.argv[3]
print glob.glob(FichierSortie)
if len(glob.glob(FichierSortie))>0 :
   sys.exit("File "+ FichierSortie+" Exist")
Ecrire=open(FichierSortie, 'w')


for Ligne in LireVCF :
    if Ligne[0]!="#" :
       #SplitLigne=Ligne.replace("\n","").split()
       #SplitLigne[0]=DicCorChro[SplitLigne[0]]
       Chaine=Ligne
    else :
      if Ligne[0:6]=="#CHROM" :
         SplitLigne=Ligne.replace("\n","").split()
         Chaine="\t".join(SplitLigne[0:9])
         for OldInd in SplitLigne[9::] :
             Chaine+="\t"+DicCorNom[OldInd.upper()]
         Chaine+="\n"
      else :
         Chaine=Ligne
    Ecrire.write(Chaine)

