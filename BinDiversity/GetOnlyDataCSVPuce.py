#!/usr/bin/python
import sys


def GetListeSnpPuceV1(Fichier) :
    PosChro=1
    PosPos=2
    ListePos=[]
    LireFich=OuvrirFichier(Fichier, 'r')
    for ligne in LireFich:
        if ligne[0]!='#':
            splitligne=ligne.split()
            NomPos=splitligne[PosChro]+"-"+splitligne[PosPos]
            ListePos.append(NomPos)
    print ListePos[0]
    LireFich.close()
    return(list(ListePos))

def GetListeSnpPuceV2(Fichier) :
    PosChro=1
    PosPos=2
    ListePos=[]
    Dic={}
    LireFich=OuvrirFichier(Fichier, 'r')
    for ligne in LireFich:
        if ligne[0]!='#':
            splitligne=ligne.split()
            #NomPos=splitligne[PosChro]+"-"+splitligne[PosPos]
            if splitligne[PosChro] not in Dic.keys():
               Dic[splitligne[PosChro]]=[]
            Dic[splitligne[PosChro]].append(int(splitligne[PosPos]))
    for key in Dic.keys():
       Dic[key].sort()
       Dic[key]=set( Dic[key])
    LireFich.close()
    return(Dic)


def OuvrirFichier(Fichier, OptionOuverture) :
   try :
      LE=open(Fichier,OptionOuverture)
   except :
      sys.exit("File "+Fichier+" with " + OptionOuverture + "Can't open")
   return LE

def NouveauCSVAvecPuce(FichierCSVIn, FichierCSVOut, ListePos) :
    LireFich=OuvrirFichier(FichierCSVIn, 'r')
    Ecrire=OuvrirFichier(FichierCSVOut, 'w')
    Cmt=0
    for ligne in LireFich:
       if Cmt%1000==0:
          print Cmt
       if ligne[0]!='#':
           splitligne=ligne.split()
           Chro=splitligne[0]
           Pos=int(splitligne[1])
           if Chro in ListePos.keys() :
              if Pos in ListePos[Chro]:
                 Ecrire.write(ligne)
       else : 
           Ecrire.write(ligne)
       Cmt+=1
    Ecrire.close()
    LireFich.close()


##
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


def LaunchCommande():
     Commande="exe -ci [CSV] -co [CSVOut] -p [50K Puce]"
     FichierOut=GetUnArgument(sys.argv, '-co')
     if len(FichierOut) != 1 :
        sys.exit(Commande)
     FichierOut=FichierOut[0]
     FichierEntree=GetUnArgument(sys.argv, '-ci')
     if len(FichierEntree) != 1 :
        sys.exit(Commande)
     FichierEntree=FichierEntree[0]
     if FichierEntree == FichierOut :
        sys.exit("File Entree same that file out"+FichierOut+"\n"+Commande)
     FichierPuce=GetUnArgument(sys.argv, '-p') 
     if len(FichierPuce) != 1 :
        sys.exit(Commande)
     FichierPuce=FichierPuce[0]
     ListePos=GetListeSnpPuceV2(FichierPuce) 
     NouveauCSVAvecPuce(FichierEntree, FichierOut, ListePos)

LaunchCommande()
