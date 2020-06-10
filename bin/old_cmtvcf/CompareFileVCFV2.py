#!/usr/bin/python
# coding=utf-8
import os
import sys
import glob


def Str(Number, Lim=2):
    if type(Number)==float :
       return str(round(Number,Lim))
    return str(Number)



def OpenFile(Fichier, Type):
    try :
       LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " OpenFile in mode " + Type + "can't OpenFile")
    return LireFich


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



def CreateDirectory(Dir):
    try:
       os.mkdir(Dir)
    except OSError:
       pass

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


def PrintStatBrutInd(DataStatInd, Fichier) :
    Ecrire=open(Fichier, 'w')
    Cmt=0
    for Ind in DataStatInd :
        if Cmt ==0:
           Ecrire.write('Ind')
           for NomEntet in DataStatInd[Ind]: 
               Ecrire.write("\t"+NomEntet)
           Ecrire.write('\n')
        Ecrire.write(Ind)
        for Stat in DataStatInd[Ind]:
            Ecrire.write("\t"+str(DataStatInd[Ind][Stat]))
        Ecrire.write("\n")
        Cmt+=1 



def GetInfoVCF(Fichier, FichierStatInd,FichierStatLocus ,Balisebcftools=False,ListePosCom=[], Entete=None, SaveInfo=True, BaliseAll=False) :
    Lire=open(Fichier)
    if (SaveInfo) :
       Dico={}
       DicoLocus={}
    else :
       Dico=None
       DicoLocus=None
    DicoInfoInd={}
    EcrireLocus=open(FichierStatLocus, 'w')
    EcrireLocus.write("Chro\tPos\tNbAl\tNbIndHom\tNbIndHet\tNbIndRef\tNbIndNC\n")
    for Ligne in Lire :
       if Ligne[0]!='#':
          TmpSplit=Ligne.split()
          Chro=TmpSplit[0]
          Pos=TmpSplit[1]
          NomPos=Chro+"_"+Pos
          Ref=TmpSplit[3]
          ListeAlt=TmpSplit[4].split(',')
          NbAltLocus=len(ListeAlt)
          NbIndHom=0
          NbIndHet=0
          NbIndRef=0
          NbIndNC=0
          if len(TmpSplit)!= len(ListeInd)+9:
             print "Warning nb ind not good "+Chro +" " + Pos + " " + str( len(TmpSplit)) + " " + str(len(ListeInd)+9)
             print TmpSplit, ListeInd, Entete
          elif (len(ListePosCom)==0 or NomPos in ListePosCom) and (TmpSplit[6]=='PASS'  or BaliseAll):
              if SaveInfo:
                 DicoLocus[NomPos]={'Ref':Ref, 'Alt':ListeAlt, 'Chro':Chro, 'Pos':Pos, 'Freq':[0]*(len(ListeAlt)+1), 'FreqHaplotype':GetMatrix(len(ListeAlt)+1,len(ListeAlt)+1)}
              for CmtInd1 in range(0,len(ListeInd)) :
                  CmtInd=CmtInd1+9
                  SumProb=10
                  if Balisebcftools:
                     TmpProb=TmpSplit[CmtInd].split(':')[1].split(',')
                     SumProb=sum([int(x) for x in TmpProb])
                  if TmpSplit[CmtInd][0]=='.' or SumProb==0:
                     A1S='.'
                     A2S='.'
                     DicoInfoInd[ListeInd[CmtInd1]]['NbNC']+=1
                     NbIndNC+=1
                  else :
                     try :
                       TempSplit2=TmpSplit[CmtInd].split(':')[0].split('/')
                       A1=int(TempSplit2[0])
                       A2=int(TempSplit2[1])
                     except :
                       print TmpSplit[CmtInd]
                       print TmpSplit
                       sys.exit()
                     if SaveInfo :
                        DicoLocus[NomPos]['Freq'][A1]+=1
                        DicoLocus[NomPos]['Freq'][A2]+=1
                        DicoLocus[NomPos]['FreqHaplotype'][A1][A2]+=1
                        DicoLocus[NomPos]['FreqHaplotype'][A2][A1]+=1
                     if A1 == 0:
                        A1S=Ref
                     else :
                        A1S=ListeAlt[A1-1]  
                     if A2 == 0:
                        A2S=Ref
                     else :
                        A2S=ListeAlt[A2-1]  
                     if A1!=A2 :
                         DicoInfoInd[ListeInd[CmtInd1]]['NbHet']+=1
                         NbIndHet+=1
                     elif A1==0:
                         DicoInfoInd[ListeInd[CmtInd1]]['NbHomRef']+=1
                         NbIndRef+=1
                     else :
                         DicoInfoInd[ListeInd[CmtInd1]]['NbHomAlt']+=1
                         NbIndHom+=1
                  if SaveInfo :
                     Dico[ListeInd[CmtInd1]][NomPos]=[A1S,A2S]#, TmpSplit[CmtInd]]
              EcrireLocus.write(str(Chro)+"\t"+str(Pos)+"\t"+str(NbAltLocus+1)+"\t"+str(NbIndHom)+"\t"+str(NbIndHet)+"\t"+str(NbIndRef)+"\t"+str(NbIndNC)+"\n")
       else :
          if Ligne[0:6]=="#CHROM":
             ListeInd=[]
             TmpSplit=Ligne.split()
             if Entete==None :
                TempListeInd=TmpSplit[9::]
             else :
                 TempListeInd=Entete
                 if len(TempListeInd)!=len(Entete) :
                    sys.exit("Nombre d'entete donnee au prog different de celui du fichier")
             for Ind in TempListeInd:
                 ListeInd.append(Ind.lower())
                 if SaveInfo :
                    Dico[Ind.lower()]={}
                 DicoInfoInd[Ind.lower()]={'NbHet':0, 'NbHomAlt':0, 'NbNC':0, 'NbHomRef':0}
    EcrireLocus.close()
    PrintStatBrutInd(DicoInfoInd, FichierStatInd)
    return {'DataAll':Dico, 'InfoSnp':DicoLocus, 'InfoInd':DicoInfoInd}

def CmpLocusPuce(InfoLocus1, InfoLocusPuce):
    NomLocus=intersect(InfoLocus1.keys(), InfoLocusPuce.keys()) 
    Chaine="" 
    NbLocusPuce=len(InfoLocusPuce.keys())
    NbLocusNGS=len( NomLocus)
    Chaine+="===Resume===\n  * Nombre de Locus de la Puce : "+str( NbLocusPuce) + "\n  * Nombre de locus NGS : "+ str(NbLocusNGS)+"\n"
    ### Pour chaque locus
    NbLocusAltPucePresent=0
    NbLocusAltNGSNonPresent=0
    NbLocusAltNGSNonPresentAltPucePresent=0
    NbLocusAltIdentique=0
    NbLocusDifferent=0
    for Locus in NomLocus :
       Info1Locus=InfoLocus1[Locus]
       Info1LocusPuce=InfoLocusPuce[Locus]
       NewAltLocus=[]
       ### On verifie les frequences
       for cmt in range(1,len(Info1Locus['Freq'])) :
           if Info1Locus['Freq'][cmt] > 0 :
               NewAltLocus.append(InfoLocus1[Locus]['Alt'][cmt-1])
       Diff1=diff(NewAltLocus, InfoLocusPuce[Locus]['Alt'])  
       DiffPuce=diff(InfoLocusPuce[Locus]['Alt'], NewAltLocus)
       Info1Locus['NewAlt']=NewAltLocus
       if len(DiffPuce)==0 :
           NbLocusAltPucePresent+=1
       if len(Diff1) > 0 :
           NbLocusAltNGSNonPresent+=1
       if len(Diff1) > 0 and len(DiffPuce)==0:
           NbLocusAltNGSNonPresentAltPucePresent+=1
       if len(Diff1) == 0 and len(DiffPuce)==0:
           NbLocusAltIdentique+=1
       if len(Diff1) != 0 and len(DiffPuce)!=0:
           NbLocusDifferent+=1
    Chaine+="=== Comparaison des Bases ALternatives ===\n  * Nb Locus identique :"+str( NbLocusAltIdentique) +"  \n  * Nombre de locus ayant le meme Alt  : " + str(NbLocusAltPucePresent) + "\n  * Nombre de Locus ayant un alt different " + str(NbLocusAltNGSNonPresent)+"\n  *  Nombre de locus ayant le meme Alt et de NGS ayant un autre Alt "+ str(NbLocusAltNGSNonPresentAltPucePresent) +"\n  * Nb Locus Different : "+str(NbLocusDifferent)+ "\n"
    Stat={'Stat':[NbLocusPuce, NbLocusNGS, NbLocusAltPucePresent, NbLocusAltNGSNonPresent, NbLocusAltNGSNonPresentAltPucePresent, NbLocusAltIdentique, NbLocusDifferent] ,'ChaineWiki':Chaine}
    return Stat


def NbDifference(A1,B1,A2,B2):
   NbDiff=0
   if A1!=A2 :
     NbDiff+=1
   if A1!=B2:
     NbDiff+=1
   if B1!=A2:
     NbDiff+=1
   if B1!=B2:
     NbDiff+=1
   return NbDiff

def GetScoreLocus(GenLocus1,GenLocusPuce, Ref, Alt):
     ## Cas ou les deux sont Non Couverts :
    if GenLocus1[0]=='.' and GenLocusPuce[0]=='.':
       return 0
    ## Cas ou le NGS est NC
    if GenLocus1[0]=='.':
       return 1
    ## Cas ou la puce est NC
    if GenLocusPuce[0]=='.':
        return 2
    NbDiff=NbDifference(GenLocus1[0],GenLocus1[1],GenLocusPuce[0],GenLocusPuce[1]) 
    ## Cas ou les deux sont homogotes
    if GenLocus1[0]==GenLocus1[1] and GenLocusPuce[0]==GenLocusPuce[1]:
       #Cas ou ce n'est ni la reference ni l'alternatif
       if (GenLocus1[0]!=Alt and GenLocus1[0]!=Ref) or (GenLocus1[1]!=Alt and GenLocus1[1]!=Ref):
          return 15
       ## Cas ou c'est identique
       if NbDiff==0 :
          if Ref==GenLocusPuce[0]:
             return 3
          else :
             return 4 
       ## Cas ou cela est different
       elif NbDiff==4:
          if Ref==GenLocusPuce[0]:
             return 5
          else :
             return 6 
    if GenLocus1[0]!=GenLocus1[1] and  GenLocusPuce[0]!=GenLocusPuce[1]  :
       if NbDiff==2 :
          return 7
       else :
          return 8
    ## Cas ou la puce est heterozygote
    if GenLocusPuce[0]!=GenLocusPuce[1]:
       if NbDiff==2:
          if GenLocus1[0]==Ref:
             return 9
          else :
             return 10
       else :
          return 11
    if GenLocusPuce[0]==GenLocusPuce[1]:
       if NbDiff==2:
          if GenLocusPuce[0]==Ref:
             return 12
          else :
             return 13
       else :
          return 14


def CmpIndLocus(Dico1, DicoPuce, DicInfoCorPuce) :
    ListeLocus=intersect(Dico1['InfoSnp'].keys(), DicoPuce['InfoSnp'].keys())
    if DicInfoCorPuce==None :
       ListeInd=intersect(Dico1['DataAll'].keys(),DicoPuce['DataAll'].keys())
    else :
       ListeCorr=DicInfoCorPuce.keys()
       for Corr in ListeCorr :
          if DicInfoCorPuce[Corr] not in DicoPuce['DataAll'].keys() : 
             print Corr, DicInfoCorPuce[Corr]
             del DicInfoCorPuce[Corr]
       ListeInd=intersect(Dico1['DataAll'].keys(),DicInfoCorPuce.keys())
       ListeDNA=intersect(DicInfoCorPuce.values(), DicoPuce['DataAll'].keys())
       if len(ListeDNA)==0: 
           sys.exit("Pas de dna commun d'entete entre puce et vcf")
    if len(ListeInd)==0: 
       sys.exit("Pas de nom commun d'entete entre puce et vcf")
    if DicInfoCorPuce==None :
       print diff(DicoPuce['DataAll'].keys(), Dico1['DataAll'].keys())
       print diff(Dico1['DataAll'].keys(), DicoPuce['DataAll'].keys())
    else :
       print diff(DicInfoCorPuce.keys(), Dico1['DataAll'].keys())
       print diff(Dico1['DataAll'].keys(), DicInfoCorPuce.keys())
    ##on compare pour chaque locus 
    ResParLocus={}
    for Locus in ListeLocus:
        ResParLocus[Locus]=[0]*16
    ResParInd={}
    for Ind in ListeInd :
        ResParInd[Ind]=[0]*16
    for Ind in ListeInd:
        if DicInfoCorPuce!=None :
           IndPuce=DicInfoCorPuce[Ind]
        else :  
           IndPuce=Ind
        for Locus in ListeLocus :
            if Locus not in Dico1['DataAll'][Ind] :
               print "Ind " + Ind + " Locus "+ Locus + " Not present in 1",
               if Locus not in DicoPuce['DataAll'][IndPuce]:
                  print "and in DicoPuce" 
            elif Locus not in DicoPuce['DataAll'][IndPuce] :
               print "Ind " + Ind + " Locus "+ Locus + " Not present in Puce"
            else :
               TmpStat=GetScoreLocus(Dico1['DataAll'][Ind][Locus],DicoPuce['DataAll'][IndPuce][Locus], DicoPuce['InfoSnp'][Locus]["Ref"], DicoPuce['InfoSnp'][Locus]["Alt"][0])
               ResParLocus[Locus][TmpStat]+=1
               ResParInd[Ind][TmpStat]+=1 
    return {'ResParLocus':ResParLocus, 'ResParInd':ResParInd}


def CmpIndWithDico(Dico1, DicoPuce, NomInd):
    ListeLocus=intersect(Dico1['InfoSnp'].keys(), DicoPuce['InfoSnp'].keys())
    DataInd=Dico1['DataAll'][NomInd]
    ListeInd=DicoPuce['DataAll'].keys()
    ResParInd={}
    for Ind in ListeInd :
        ResParInd[Ind]=[0]*16
    for Ind in ListeInd:
        for Locus in ListeLocus :
               TmpStat=GetScoreLocus(DataInd[Locus],DicoPuce['DataAll'][Ind][Locus], DicoPuce['InfoSnp'][Locus]["Ref"])
               ResParInd[Ind][TmpStat]+=1 
    return ResParInd
        

def PrintResInd(FileSortie, ResParInd, Entete, BaliseWiki=False):
    if BaliseWiki :
       sep=" | "
       sepEnteteDebCol=sepEnteteCol=" ^ "
       sepEnteteColFin=" | "
       sepEnteteLigne="^ "
    else :
       sep="\t"
       sepEnteteCol="\t"
       sepEnteteDebCol=sepEnteteColFin=""
       sepEnteteLigne=""
    Chaine=sepEnteteDebCol+Concatener(Entete ,sepEnteteCol)+sepEnteteColFin+"\n"
    for Ind in ResParInd.keys():
        Chaine+=sepEnteteLigne+Ind+sep+Concatener(ResParInd[Ind], sep)+"\n"
    EcrireStatInd=open(FileSortie,"w")
    EcrireStatInd.write(Chaine)
    EcrireStatInd.close()

def CheckDivZero(a,b):
    try :
       c=float(a)/float(b)
    except:
       return -1
    return c
 
def ConcateneStatResumePuce(Resultat) :
    FinalSortie={}
    for Key in Resultat :
        ## Hom et Ref font reference a l'etat de la puce
        NonCouvertPuce= Resultat[Key][0]+ Resultat[Key][2]
        NonCouvertNGS= Resultat[Key][0]+ Resultat[Key][1]
        NonCouvertPar12=Resultat[Key][0]+ Resultat[Key][1]+Resultat[Key][2]
        SumPuceHomRef= Resultat[Key][3] + Resultat[Key][5] + Resultat[Key][12]
        SumPuceHomAlt= Resultat[Key][4] + Resultat[Key][6] + Resultat[Key][13]
        SumFauxHomRef=Resultat[Key][5]
        SumFauxHomAlt=Resultat[Key][6]
        SumPuceHet=Resultat[Key][7]+Resultat[Key][8]+ Resultat[Key][9]+Resultat[Key][10]+Resultat[Key][11]
        SumPuceHetBon=Resultat[Key][7]
        SumNGSHetMauvais=Resultat[Key][8] + Resultat[Key][12] +  Resultat[Key][13] + Resultat[Key][14]
        SumPuceHomNiAltNiRef=Resultat[Key][15]
        Total=sum(Resultat[Key])
        ## % Non couvert puce                   % non Couvert NGS                       % non couvert 2         % Faux Positif hom : Faux en homo et base ref % Faux Negatif : non detecte alors qu'existan %de Het Bon
        FinalSortie[Key]=[CheckDivZero(NonCouvertPuce,Total)*100, CheckDivZero(NonCouvertNGS,Total)*100, CheckDivZero(NonCouvertPar12,Total)*100,  CheckDivZero(SumFauxHomRef,SumPuceHomRef)*100, CheckDivZero(SumFauxHomAlt,SumPuceHomAlt)*100, CheckDivZero(SumPuceHetBon,SumPuceHet)*100,  CheckDivZero(SumPuceHet-SumPuceHetBon,SumPuceHet)*100, CheckDivZero(SumNGSHetMauvais,SumPuceHomRef+SumPuceHomAlt)*100, CheckDivZero(SumPuceHomNiAltNiRef,SumPuceHomRef+SumPuceHomAlt)*100, CheckDivZero(SumPuceHomRef+SumPuceHomAlt+SumPuceHet+NonCouvertPar12+SumPuceHomNiAltNiRef,Total)*100, Total, SumPuceHomAlt, SumPuceHet, SumPuceHomRef-SumFauxHomRef, SumPuceHomAlt-SumFauxHomAlt, SumPuceHetBon, SumPuceHomNiAltNiRef]
    return FinalSortie




def ConcateneStatResumePuceSomme(Resultat) :
    ## Perc Couverture Puce
    NonCouvertPuce=0
    ## Perc Couverture NGS
    NonCouvertNGS=0
    ## Non Couvert par un des deux
    NonCouvertPar12=0
    ## HomozygoteRef 
    SumPuceHomRef=0
    ## HomozygoteRef 
    SumPuceHomAlt=0
    ## Nb Faux Puce Ref Hom
    SumFauxHomRef=0
    ### Hom 
    SumFauxHomAlt=0
    ## Somme des Heterozygote de la puce
    SumPuceHet=0
    ## Somme des Hetereozygote
    SumPuceHetBon=0
    ## Nb Faux Puce Ref Het
    SumNGSHetMauvais=0
    SumPuceHomNiAltNiRef=0
    #Total
    Total=0
    for Key in Resultat :
        ## Hom et Ref font reference a l'etat de la puce
        NonCouvertPuce+= Resultat[Key][0]+ Resultat[Key][2]
        NonCouvertNGS+= Resultat[Key][0]+ Resultat[Key][1]
        NonCouvertPar12+=Resultat[Key][0]+ Resultat[Key][1]+Resultat[Key][2]
        #Nombre de reference sur la puce
        SumPuceHomRef+= Resultat[Key][3] + Resultat[Key][5] + Resultat[Key][12]
        #Nombre d alternatif sur la puce
        SumPuceHomAlt+= Resultat[Key][4] + Resultat[Key][6] + Resultat[Key][13]
        SumFauxHomRef+=Resultat[Key][5]
        SumFauxHomAlt+=Resultat[Key][6]
        SumPuceHet+=Resultat[Key][7]+Resultat[Key][8]+ Resultat[Key][9]+Resultat[Key][10]+Resultat[Key][11]
        SumPuceHetBon+=Resultat[Key][7]
        SumNGSHetMauvais+=Resultat[Key][8] + Resultat[Key][12] +  Resultat[Key][13] + Resultat[Key][14]
        SumPuceHomNiAltNiRef+=Resultat[Key][15]
        Total+=sum(Resultat[Key])
    if Total==0:
       Total=-1
    if SumPuceHomAlt==0:
       SumPuceHomAlt=-2
    if SumPuceHet==0:
       SumPuceHet=-3
    ## % Non couvert puce                   % non Couvert NGS                       % non couvert 2         % Faux Positif hom : Faux en homo et base ref % Faux Negatif : non detecte alors qu'existan %de Het Bon
    return [CheckDivZero(NonCouvertPuce,Total)*100, CheckDivZero(NonCouvertNGS,Total)*100, CheckDivZero(NonCouvertPar12,Total)*100,  CheckDivZero(SumFauxHomRef,SumPuceHomRef)*100, CheckDivZero(SumFauxHomAlt,SumPuceHomAlt)*100, CheckDivZero(SumPuceHetBon,SumPuceHet)*100,  CheckDivZero(SumPuceHet-SumPuceHetBon,SumPuceHet)*100, CheckDivZero(SumNGSHetMauvais,SumPuceHomRef+SumPuceHomAlt)*100, CheckDivZero(SumPuceHomNiAltNiRef,SumPuceHomRef+SumPuceHomAlt)*100, CheckDivZero(SumPuceHomRef+SumPuceHomAlt+SumPuceHet+NonCouvertPar12+SumPuceHomNiAltNiRef,Total)*100, Total, SumPuceHomAlt, SumPuceHet, SumPuceHomRef-SumFauxHomRef, SumPuceHomAlt-SumFauxHomAlt, SumPuceHetBon, SumPuceHomNiAltNiRef]

def GetPosComm(File) :
    read=open(File) 
    FinalListe=[]
    for ligne in read :
        temp=ligne.split()
        FinalListe.append(temp[0]+"_"+temp[1])
    return set(FinalListe)

     

def CompareVCFWithPuce(DataVCF, DicoSnpPuce, DicInfoInd):
    DicoModele={}
    Cmt=0
    for Key in DataVCF:
        Modele=Key
        File=Key
        DicoSnp=DataVCF[Key]
        ResultLocus=CmpLocusPuce(DicoSnp['InfoSnp'], DicoSnpPuce['InfoSnp'])
        ResultatAll=CmpIndLocus(DicoSnp, DicoSnpPuce, DicInfoInd)
        EnteteCmpInfLocus=["Ind","AllNC","Loc1NC","LocPuceNC","RefHomBon","AltHomBon","RefHomFaux","AltHomFaux","DeuxHetBon","DeuxHetFaux","NGSHomRefPuceHet1Bon","NGSHomAltPuceHet1Bon","NGSHomPuceHet2Faux","PuceHomRefNGSHet1Bon","PuceHomAltNGSHet1Bon","PuceHomNGSHet2Faux", "PuceHomAltRefFausse"]
        #PrintResInd(File+".cmppuce.ind.wiki", ResultatAll["ResParInd"], EnteteCmpInfLocus , True)
        PrintResInd(File+".cmppuce.ind.stat", ResultatAll["ResParInd"], EnteteCmpInfLocus)
        EnteteCmpInfLocus=["Locus","AllNC","Loc1NC","LocPuceNC","RefHomBon","AltHomBon","RefHomFaux","AltHomFaux","DeuxHetBon","DeuxHetFaux","NGSHomRefPuceHet1Bon","NGSHomAltPuceHet1Bon","NGSHomPuceHet2Faux","PuceHomRefNGSHet1Bon","PuceHomAltNGSHet1Bon","PuceHomNGSHet2Faux", "PuceHomAltRefFausse"]
        PrintResInd(File+".cmppuce.all.locus.stat", ResultatAll["ResParLocus"], EnteteCmpInfLocus)
        ListeEntetStat=["Ind","% Puce NC", "% NGS NC", " % Non Couvert 1/2", " % FP", "% FN", "% Het Bon", " % FN Het", "% FP", "% Hom Ni Alt Ni Ref","% Somme"]
        ListeEntetStat=["Ind","Perc_Puce_NC", "Perc_NGS_NC", "Perc_Non_Couvert_2", " Perc_FP", "Perc_FN", "Perc_Het_Bon", " Perc_FN_Het", "Perc_FP_Het", "Perc_Hom_Ni_Alt_Ni_Ref", "Perc_Somme", "Total", "NbPuceHomAlt", "NbPuceHet","NbRefBon", "NbAltBon",  "NbHetBon", "NbNiAltNiRef"]
        StatResumeLocus=ConcateneStatResumePuce(ResultatAll["ResParLocus"])
        PrintResInd(File+".cmppuce.locus.resume.stat", StatResumeLocus, ListeEntetStat)
        DicoModele[Modele]=ConcateneStatResumePuceSomme(ResultatAll["ResParInd"])
        PrintResInd(File+".cmppuce.ind.resume.stat", ConcateneStatResumePuce(ResultatAll["ResParInd"]), ListeEntetStat)
        Cmt+=1

def GetDataVCF(ListeFileVCF, BaliseBCF, Entete, SaveInfo, BaliseAll=False) :
    Cmt=0
    Dic={}
    for File in ListeFileVCF :
        print File
        Dic[File]=GetInfoVCF(File,File+".brut.ind.stat" ,File+".brut.locus.stat" , BaliseBCF[Cmt], Entete=Entete, SaveInfo=SaveInfo, BaliseAll=BaliseAll)
        Cmt+=1
    return Dic

def GetEntete(File):
    Lire=open(File)
    Entete=[]
    for Ligne in Lire:
       Entete.append(Ligne.lower().replace("\n",""))
    return Entete

def GetStatutInd(A1,A2): 
    if A1=='.':
      return 2
    if A1==A2:
      return 0
    if A1!=A2:
      return 1

def Getord(Chaine):
    Cmt=0
    for a in Chaine :
        Cmt+=ord(a)
    return Cmt

def CmpDeuxVCF(VCF1, VCF2, File, Nom1, Nom2, BaliseNew):
    ListeLocus=intersect(VCF1['InfoSnp'].keys(), VCF2['InfoSnp'].keys())
    ListeInd = intersect(VCF1['DataAll'].keys(),VCF2['DataAll'].keys())
#ico[ListeInd[CmtInd1]][NomPos]=[A1S,A2S, TmpSplit[CmtInd]]
    CmpVCF=GetMatrix(3,3)
    CmpVCFDiff=GetMatrix(3,3)
    for Ind in ListeInd:
#def NbDifference(A1,B1,A2,B2):
       for Locus in ListeLocus :
           A1=VCF1['DataAll'][Ind][Locus][0]
           A2=VCF1['DataAll'][Ind][Locus][1]
           B1=VCF2['DataAll'][Ind][Locus][0]
           B2=VCF2['DataAll'][Ind][Locus][1]
           StatInd1=GetStatutInd(A1,A2)
           StatInd2=GetStatutInd(B1,B2)
           Diff=Getord(A1)+Getord(A2)-(Getord(B1)+Getord(B2))
           if Diff==0:
              CmpVCF[StatInd1][StatInd2]+=1
           else :
              CmpVCFDiff[StatInd1][StatInd2]+=1
    if BaliseNew==True :
       Ecrire=open(File,'w') 
    else :
       Ecrire=open(File,'a')
    Ecrire.write("FileVCF1 : "+Nom1+"\n")
    Ecrire.write("FileVCF2 : "+Nom2+"\n")
    Ecrire.write("Number Identical locus\n")
    Entete=["HomFile1", "HetFile1","NCFile1"]
    Ecrire.write("\t\tHomFile2\tHetFile2\tNCFile2\n")
    for Cmt in range(0,3):
        Ecrire.write(Entete[Cmt])
        for Cmt2 in range(0,3):
            Ecrire.write("\t"+str(CmpVCF[Cmt][Cmt2]))
        Ecrire.write("\n")
    Ecrire.write("\n")
    Ecrire.write("Number Different locus\n")
    Entete=["HomFile1", "HetFile1","NCFile1"]
    Ecrire.write("\t\tHomFile2\tHetFile2\tNCFile2\n")
    for Cmt in range(0,3):
        Ecrire.write(Entete[Cmt])
        for Cmt2 in range(0,3):
            Ecrire.write("\t"+str(CmpVCFDiff[Cmt][Cmt2]))
        Ecrire.write("\n")
    Ecrire.write("\n")
    Ecrire.close()

def GetFileInfoIndToCompare(FileInfoInd):
    Lire=open(FileInfoInd)
    Dic={}
    for Ligne in Lire :
        SplitL=Ligne.strip().split()
        if len(SplitL)==2 :
           if len(SplitL)>1 and len(SplitL)>1 : 
              Dic[SplitL[0].lower()]=SplitL[1].lower()
    print Dic
    return Dic

## DicoSnpPuce : "/home/jeantristan/EspaceTravail/EspaceTest/AnalyseMultiVCF/DataPuce/Echantillon.vcf"
##
## ListeFile=["bcftools/Snp_Bowtie_20_0_0_2.vcf", "Varscan/Snp_Bowtie_20_0_0.05_0_0_2.vcf"]
def LaunchCommande():
     Commande="exe -s [PuceADN, VCF, opt] -c [VCF1 VCF2...] FileOuptpout -e [ListeEntete opt]  -t [b : bcftools or v:varscan default v * csv number file]\n"
     Fichier50K=GetUnArgument(sys.argv, '-s')
     if len(Fichier50K)>1:
        sys.exit(Commande)
     elif len(Fichier50K)==0 :
        Fichier50K=None
     else :
        Fichier50K=Fichier50K[0]
     ListeFichierVCF=GetUnArgument(sys.argv, '-c')
     if len(ListeFichierVCF) == 0:
        sys.exit(Commande)
     if len(ListeFichierVCF)>1 :
        FichierCompVCF=GetUnArgument(sys.argv,'-oc')
        if len(FichierCompVCF)==0 :
           sys.exit("-oc not found, must be present when there is more one vcf file")
        FichierCompVCF=FichierCompVCF[0]
     ListeTypeFichier=GetUnArgument(sys.argv, '-t')
     if len(ListeTypeFichier)==0:
        BaliseBCF=[False]*len(ListeFichierVCF)
     elif len(ListeTypeFichier)==1 :
        if ListeTypeFichier[0].lower()=='v':
           BaliseBCF=[False]*len(ListeFichierVCF)
        elif (ListeTypeFichier[0].lower())=='b':
           BaliseBCF=[True]*len(ListeFichierVCF)
        else :
           sys.exit("option -t only b or v not "+ListeTypeFichier[0] + "\n"+Commande)
     else :
         if len(ListeTypeFichier)!=len(ListeFichierVCF):
            sys.exit("Number argument of -t option is different of number of csv file\n"+Commande)
         else :
            BaliseBCF=[] 
            for Option in ListeTypeFichier :
               if Option.lower()=='v':
                  BaliseBCF.append(False)
               elif Option.lower()=='b': 
                  BaliseBCF.append(True)
               else :
                  sys.exit("option -t only b or v not "+ListeTypeFichier[0] + "\n"+Commande)
     BaliseAll=GetUnArgument(sys.argv, '-ba')
     if len(BaliseAll)==0 or BaliseAll[0]=='f' :
        BaliseAll=False
     elif BaliseAll[0]=='t' :
        BaliseAll=True
     else :
        sys.exit("option -ba only t or f "+BaliseAll[0] + "\n"+Commande)
     ListeFichierEntet=GetUnArgument(sys.argv, '-e')
     if len(ListeFichierEntet)==0 :
        Entete=None
     else :
        Entete=GetEntete(ListeFichierEntet[0])
     FileInfoInd=GetUnArgument(sys.argv, '-ii')
     DicInfoInd=None
     if FileInfoInd!= None : 
       DicInfoInd=GetFileInfoIndToCompare(FileInfoInd[0])
     SaveInfo=True
     print Entete
     if len(ListeFichierVCF)==1 and Fichier50K==None :
        SaveInfo=False
     DicoVCF=GetDataVCF(ListeFichierVCF,BaliseBCF,Entete, SaveInfo, BaliseAll)
     if Fichier50K!= None :
        DicoSnpPuce=GetInfoVCF(Fichier50K, Fichier50K+".brut.ind.stat",Fichier50K+".brut.locus.stat")
        CompareVCFWithPuce(DicoVCF,  DicoSnpPuce, DicInfoInd)
     Cmt=0
     if len(ListeFichierVCF)>1 :
         for Cmt1 in range(0,len(ListeFichierVCF)):
             if Cmt==0 :
                BaliseNew=True
             else :
                BaliseNew=False
             CmpDeuxVCF(DicoVCF[ListeFichierVCF[Cmt1]], DicoSnpPuce,FichierCompVCF, ListeFichierVCF[Cmt1], "Puce", BaliseNew)
             Cmt+=1
         for Cmt1 in range(0,len(ListeFichierVCF)-1):
             for Cmt2 in range(Cmt1+1, len(ListeFichierVCF)):
                 if Cmt==0:
                    BaliseNew=True
                 else :
                    BaliseNew=False
                 CmpDeuxVCF(DicoVCF[ListeFichierVCF[Cmt1]], DicoVCF[ListeFichierVCF[Cmt2]],FichierCompVCF, ListeFichierVCF[Cmt1], ListeFichierVCF[Cmt2], BaliseNew)
                 Cmt+=1


#
#
def CmpDicoAllModele(Dico1, Dico2, DicoPuce):
    ListeLocus=intersect(intersect(Dico1['InfoSnp'].keys(), DicoPuce['InfoSnp'].keys()),Dico2['InfoSnp'].keys())
    ListeInd=Dico1['DataAll'].keys()
    ##on compare pour chaque locus 
    for Ind in ListeInd:
        for Locus in ListeLocus :
            Score=GetScoreLocus(Dico1['DataAll'][Ind][Locus],DicoPuce['DataAll'][Ind][Locus], Dico2['InfoSnp'][Locus]["Ref"], DicoPuce['InfoSnp'][Locus]["Alt"][0])
            if Score==1 or Score==2:
               print Ind, Locus,Dico1['DataAll'][Ind][Locus][1] ,Dico1['DataAll'][Ind][Locus][2],Dico2['DataAll'][Ind][Locus][1] ,Dico2['DataAll'][Ind][Locus][2], DicoPuce['InfoSnp'][Locus]["Alt"][0]
               raw_input()

if __name__ == '__main__':
   LaunchCommande()
