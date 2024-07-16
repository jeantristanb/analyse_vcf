#!/usr/bin/python
import sys
import gzip

def OuvrirFichier(Fichier, Type='r'):
    try :
       if Type =='r' :
           LireFich=open(Fichier)
       else :
           LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " OuvrirFichier in mode " + Type + "can't OuvrirFichier")
    return LireFich

def OuvrirFichierGzip(Fichier, Type='r'):
    try :
       if Type =='r' :
           LireFich=gzip.open(Fichier)
       else :
           LireFich=gzip.open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " OuvrirFichier in mode " + Type + "can't OuvrirFichier")
    return LireFich


#def GetCouvRed(Fichier):
#   Lire=open(Fichier)
#   Ligne=Lire.readline()
#   ListeInd=[]
#   Dico={}
#   SplitLigne=Ligne.split()
#   for Ind in SplitLigne:
#      ListeInd.append(Ind) 
#   for Ligne in Lire :
#      SplitLigne=Ligne.split()
#      NomPos=SplitLigne[1]
#      Dico[NomPos]={}
#      for cmt in range(2,len(ListeInd)):
#          Dico[NomPos][ListeInd[cmt]]=SplitLigne[cmt]
#   return Dico


def Concatener(Vecteur):
    Chaine=""
    for El in Vecteur :
       Chaine+=El+'\t'
    return Chaine

def GetLigneRed(LireRed,Pos):
    Balise=True
    while Balise: 
       PosFile=LireRed.tell()
       Ligne=LireRed.readline()
       SplitLigne=Ligne.split()
       if len(SplitLigne)<2 :
          return -1
       try :
          if SplitLigne[1]==Pos :
            return SplitLigne
          elif int(Pos)<int(SplitLigne[1]):
             ### on revient a l'ancienne pos dans le fichier
             LireRed.seek(PosFile)
       except :
            print Pos
            sys.exit()
       return None

def GetEnteteRed(LireRed):
   Ligne=LireRed.readline()
   SplitLigne=Ligne.split()
   ListeInd={}
   Cmt=0
   for Ind in SplitLigne:
      ListeInd[Ind]=Cmt
      Cmt+=1
   return ListeInd 
          

def ModifVCF(Fichier, FichierSortie, FichierSortieErreur, FichierRed, ListeInd, FichierEntete=None, BaliseGZVCF=False, BaliseGZRed=False):
    if BaliseGZVCF==False :
       Lire=OuvrirFichier(Fichier, 'r')
    else :
       Lire=OuvrirFichierGzip(Fichier, 'r')
    Ecrire=OuvrirFichier(FichierSortie,'w')
    EcrireErreur=OuvrirFichier(FichierSortieErreur,'w')
    if BaliseGZRed==False :
       LireRed=OuvrirFichier(FichierRed)
    else :
       LireRed=OuvrirFichierGzip(FichierRed)
    if FichierEntete==None :
       DicoEntete=GetEnteteRed(LireRed)
    else :
       LireTemp=open(FichierEntete)
       DicoEntete=GetEnteteRed(LireTemp)
       LireTemp.close()
    CmtLigne=0
    Chaine=""
    NbCase=len(ListeInd)+9
    for Ligne in Lire :
       if Ligne[0]!='#' and len(Ligne)>2:
           SplitLigne=Ligne.replace('\n','').split()
           Chro=SplitLigne[0]
           Pos=SplitLigne[1]
           SplitLigneRed=GetLigneRed(LireRed, Pos)
           if SplitLigneRed==-1 :
              Ecrire.write(Chaine)
              EcrireErreur.write(Chaine)
              print "Not Pos more after "+str(Pos)
              return -1
           if SplitLigneRed==None or len(SplitLigne)!=NbCase:
              if len(SplitLigne)!=NbCase :
                 print "Pos "+ Pos +" "+ Chro+" Wrong formatage" + str(len(SplitLigne))
                 EcrireErreur.write(Ligne)
              else :
                 print "warning Pos"+Pos +"not found "
           else :
              SplitLigne[8]+=":PR:CR"
              Cmt=9
              for Ind in ListeInd :
                  try :
                     SplitLigne[Cmt]+=":"+SplitLigneRed[DicoEntete[Ind]]
                  except :
                     print DicoEntete.keys()
                     sys.exit()
                  Cmt+=1
              Chaine+=Concatener(SplitLigne)+"\n"
       else :
          if Ligne[0:6]=="#CHROM" and len(ListeInd)>0:
             Ligne="#CHROM      POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  "
             Ligne+=Concatener(ListeInd)+"\n"
          Chaine+=Ligne
       if CmtLigne==100000 :
             Ecrire.write (Chaine)
             Chaine=""
             CmtLigne=0
       CmtLigne+=1
    Ecrire.write(Chaine)
    EcrireErreur.write(Chaine)

       
def GetListeInd(Fichier) :
   ListeInd=[]
   LireInd=open(Fichier)
   for Ligne in LireInd:
      ListeInd.append(Ligne.replace('\n',''))
   LireInd.close()
   return ListeInd


#BaliseGZVCF=True
#BaliseGZRed=True

BaliseGZVCF=False
BaliseGZRed=False



Cmd="exe FichierVCF FichierRed FichierVCFSortie FIchierVCFErreur FichierInd"
if len(sys.argv)not in [6,7] :
   sys.exit(Cmd)

FichierVCF=sys.argv[1]
FichierRed=sys.argv[2]
FichierVCFSortie=sys.argv[3]
FichierVCFSortieErreur=sys.argv[4]
FichierInd=sys.argv[5]
if len(sys.argv)==7:
   FichierEntete=sys.argv[6]
else: 
   FichierEntete=None
ListeInd=GetListeInd(FichierInd)
#def ModifVCF(Fichier, FichierSortie, FichierSortieErreur, FichierRed, ListeInd, FichierEntete=None, BaliseGZVCF=False, BaliseGZRed=False):
ModifVCF(FichierVCF, FichierVCFSortie, FichierVCFSortieErreur, FichierRed, ListeInd, FichierEntete, BaliseGZVCF=BaliseGZVCF, BaliseGZRed=BaliseGZRed)



