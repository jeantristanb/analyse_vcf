import glob
import sys



def GetPosition(FichierPosition):
    LireFichierPosition=open(FichierPosition)
    ListePoly=[]
    tem=LireFichierPosition.readline()
    for Ligne in LireFichierPosition :
        if Ligne[0]!="#" :
           TmpSplit=Ligne.split() 
           ListePoly.append(TmpSplit[0]+"_"+TmpSplit[1])
    LireFichierPosition.close()
    ListePoly=set(ListePoly)
    return ListePoly

def GetIndividu(FichierInd) :
    Lire=open(FichierInd)
    ListeInd=[]
    for Ligne in Lire :
        if Ligne[0]!="#":
          ListeInd.append(Ligne.replace("\n",""))
    return ListeInd

def GetPosInd(Lire,ListeInd) :
    Entete=""
    for Ligne in Lire :
        if Ligne[0:6]=="#CHROM" :
           SplitLigne=Ligne.split()
           ListePos=[]
           CmtPos=9
           ListeIndFind=[]
           for Ind in SplitLigne[9::] :
               if Ind in ListeInd :
                  ListePos.append(CmtPos) 
                  ListeIndFind.append(Ind)
               CmtPos+=1
           break
        else :
           Entete+=Ligne
    if len(ListePos)!=len(ListeInd) :
       print ListePos
       print ListeInd
       print ListeIndFind
       sys.exit("all Ind not find")
    return (set(ListePos), Entete+"\t".join(SplitLigne[0:9])+"\t"+"\t".join(ListeIndFind)+"\n")

FichierPosition=sys.argv[1]
FichierInd=sys.argv[2]
FichierSortie=sys.argv[3]
ListeFichier=sys.argv[4::]
if len(glob.glob(FichierSortie))>0 :
   sys.exist("File Sortie exist : " + FichierSortie)

ListePosition=GetPosition(FichierPosition)
ListeInd=GetIndividu(FichierInd)


ListeSortie=[]
CmtFichier=0
Chaine=""
for Fichier in ListeFichier :
    LireFichier=open(Fichier)
    print CmtFichier
    if CmtFichier==0 :
        (ListePosVCFInd, Entete)=GetPosInd(LireFichier,ListeInd)
    for Ligne in LireFichier :
        if Ligne[0]!="#":
           TmpSplit=Ligne.split()
           if TmpSplit[0]+"_"+TmpSplit[1] in ListePosition :
              Chaine+="\t".join(TmpSplit[0:9])
              for CmtInd in ListePosVCFInd :
                  Chaine+="\t"+TmpSplit[CmtInd]
              Chaine+="\n"
    LireFichier.close()
    CmtFichier+=1 

Ecrire=open(FichierSortie, "w")
Ecrire.write(Entete)
Ecrire.write(Chaine)
Ecrire.close()
