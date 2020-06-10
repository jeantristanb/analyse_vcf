#!/usr/bin/python
# -*- coding: utf-8 -*-
## Fichier de puce transforme avec 
## Format du fichier 
##CHROM	POS	REF	A1	A2	Panzea
#tfrom re import *
## Analyse seulement un seul VCF
import re
import sys

def intersect(a, b):
     return list(set(a) & set(b))

def OuvrirFichier(Fichier, Type):
    try :
       LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " OuvrirFichier in mode " + Type + "can't OuvrirFichier")
    return LireFich 



def LireFormatPuceV2(Fichier) :
    PosChro=1
    PosPos=2
    PosRef=3
    PosALT1=4
    PosALT2=5
    PosPanzea=6
    LireFich=OuvrirFichier(Fichier, 'r')
    DicChroNonPoly={}
    DicChroPoly={}
    DicChroHet={}
    Stat={}
    cmtSite=0
    for ligne in LireFich:
        if ligne[PosChro]!='#':
            splitligne=ligne.split()
            if splitligne[PosALT1] == splitligne[PosRef] and splitligne[PosALT2] == splitligne[PosRef] :
                DicChroNonPoly[splitligne[PosChro]+"-"+splitligne[PosPos]]=[splitligne[PosRef], splitligne[PosALT1], splitligne[PosALT2]] 
            elif splitligne[PosALT1]!=splitligne[PosALT2] :
                DicChroHet[splitligne[PosChro]+"-"+splitligne[PosPos]]=[splitligne[PosRef], splitligne[PosALT1], splitligne[PosALT2]]
            else :
                DicChroPoly[splitligne[PosChro]+"-"+splitligne[PosPos]]=[splitligne[PosRef], splitligne[PosALT1], splitligne[PosALT2]] 
    DicChroF={}
    DicChroF['Poly']=DicChroPoly
    DicChroF['NonPoly']=DicChroNonPoly
    DicChroF['Het']=DicChroHet
    DicChroF['NomFichier']=Fichier
    LireFich.close()
    return(DicChroF)

### Pour comparer
def ComparerMultiCSV(DicMultiCSV):
    ListeEntete=DicMultiCSV['Entete']
    if len(ListeEntete)==1 :
        return None
    print "Debut de Comparaison des differents fichiers "
    Dico={}
    for Cmt in range(len(ListeEntete)-1):
        Entete1=ListeEntete[Cmt]
        Chro1=DicMultiCSV['Chro'][Entete1]
        Dico[Entete1]={}
        for Cmt2 in range(Cmt+1,len(ListeEntete)):
            CmtBaseIdentique=0
            Entete2=ListeEntete[Cmt2]
            Chro2=DicMultiCSV['Chro'][Entete2]
            Key1=Chro1.keys()
            Key2=Chro2.keys()
            Intersect=intersect(Key1,Key2)
            for PosBase in Intersect :
                if Chro2[PosBase][1] == Chro1[PosBase][1]:
                    CmtBaseIdentique+=1
                else :
                    print  PosBase + " : Polymorphique mais pas meme Alternatif" + Chro1[PosBase][1] + "-" + Chro2[PosBase][1] 
            Dico[Entete1][Entete2]={'NbBaseComm':len(Intersect), 'NbBaseCommAvecMemePoly':CmtBaseIdentique}
    return Dico


def ReadMultiCSV(ListeFichier, ListeEntete) : 
    DicChroF={'Stat':{}, 'NbSite':{}, 'Chro':{}, 'Couv':{}}
    regex=re.compile("ADP=([0-9]+)")
    CmtEntete=0
    ListeBase=['A','T','C','G', 'N']
    for Fichier in ListeFichier :
        Stat={}
        DicChro={}
        for Base in ListeBase:
            Stat[Base]={}
            for Base2 in ListeBase :
                Stat[Base][Base2]=0
        LireFich=OuvrirFichier(Fichier, 'r')
        Entete=ListeEntete[CmtEntete]
        print "Debut lecture de " + Entete + " Fichier " + Fichier
        cmtSite=0
        MeanDP=0
        for ligne in LireFich:
            if ligne[0]!='#':
                splitligne=ligne.split()
                Chro=splitligne[0]
                Pos=splitligne[1]
                PosChro=Chro+"-"+Pos
                Stat[splitligne[3]][splitligne[4]]+=1
                ### Pour la couverture
                DPReg=regex.search(splitligne[7])
                DP=int(DPReg.group(1))
                MeanDP+=DP
                DicChro[PosChro]=(splitligne[3],splitligne[4], splitligne[4],DP)
                cmtSite+=1
        DicChroF['Stat'][Entete]=Stat
        DicChroF['Couv'][Entete]=MeanDP
        DicChroF['NbSite'][Entete]=cmtSite
        DicChroF['Chro'][Entete]=DicChro
        LireFich.close()
        CmtEntete+=1
    DicChroF['Fichier']=ListeFichier
    DicChroF['Entete']=ListeEntete
    return(DicChroF)


def ComparerCSVEt50k(DicMultiCSV, Data50k):
    if Data50k ==None:
       return None
    DicSort={}
    DicSort['NomFichier']=Data50k['NomFichier']
    DicPoly=Data50k['Poly']
    DicNonPoly=Data50k['NonPoly']
    ListePoly=DicPoly.keys()
    ListeNonPoly=DicNonPoly.keys()
    DicSort['NbSnpPuce']=len(ListePoly)+len(ListeNonPoly)
    DicSort['NbSnpPucePoly']=len(ListePoly) 
    DicSort['NbSnpPuceNonPoly']=len(ListeNonPoly)
    ListeEntete=DicMultiCSV['Entete']    
    for Entete in ListeEntete :
        DicSort[Entete]={}
        DicEntete=DicMultiCSV['Chro'][Entete]
        ListePosCSV=DicMultiCSV['Chro'][Entete].keys()
        ListePosCommPoly=intersect(ListePosCSV, ListePoly)
        NbPosCommunPoly=len(ListePosCommPoly)
        NbPosCommunNonPoly=len(intersect(ListePosCSV, ListeNonPoly))
        NbPolyCommun=0
        for PosComm in ListePosCommPoly :
            if DicPoly[PosComm][1]==DicEntete[PosComm][1]:
               NbPolyCommun+=1
        DicSort[Entete]['NbPosCommunePoly']=NbPosCommunPoly
        DicSort[Entete]['NbPosCommuneNonPoly']=NbPosCommunNonPoly
        DicSort[Entete]['NbPosCommunePolyMemeMut']=NbPolyCommun
    return DicSort



def FusionnerDic(DicChroF) :
    DicF={}
    ListeEntete=DicChroF['Entete']
    NbDic=len(DicChroF['Entete'])
    PosRef=NbDic
    ListeLocusI=[]
    for Entete in ListeEntete :
        ListeLocusI.extend(DicChroF['Chro'][Entete].keys())
    print len(ListeLocusI)
    ListeLocus=list(set(ListeLocusI))
    print len(ListeLocus)
    for Locus in ListeLocus:
        DicF[Locus]=["-"]*(NbDic+1)
    for CmtEntete in range(len(ListeEntete)):
        Entete=ListeEntete[CmtEntete] 
        DicE=DicChroF['Chro'][Entete]
        KeyDicE=DicChroF['Chro'][Entete].keys()
        for Locus in KeyDicE:
           DicF[Locus][CmtEntete]=DicE[Locus][1]+'-'+DicE[Locus][2]
           DicF[Locus][PosRef]=DicE[Locus][0]
    return DicF
        

def PrintListe(ListeEntete,SepDeb,Sep,SepFin):
    ChaineF=SepDeb+ListeEntete[0]
    for Entete in ListeEntete[1::] :
         ChaineF+=Sep+str(Entete)
    ChaineF+=SepFin
    return ChaineF


def PrintStatCSVMerge(CSVSortie, TypeSortie, DicComp, DicComp50K):
     Chaine=''
     if TypeSortie=='Wiki':
        SepDebEntete='^  '
        SepEntete='  ^  '
        SepDeb='|  '
        Sep='  |  '
        SepFin='  |\n'
        SepTitreDebut1="====="
        SepTitreFin1="=====\n"
     else :
        SepDebEntete=SepDeb=''
        SepTitreFin1=SepEntete=SepDeb=Sep=SepFin=';'
        SepTitreDebut1=""
     Chaine=SepTitreDebut1+"Descriptif des fichiers"+SepTitreFin1
     ListeEntete=CSVSortie['Entete']
     for CmtEntete in range(len(ListeEntete)) :
         Chaine+=SepDebEntete+ListeEntete[CmtEntete] +Sep+CSVSortie['Fichier'][CmtEntete]+SepFin
     Chaine+=PrintListe(['Nombre de site', 'Couverture Moyenne'],SepDebEntete,SepEntete,SepFin)
     for Entete in CSVSortie['Entete']:
         Chaine+=SepDebEntete+Entete+Sep+str(CSVSortie['NbSite'][Entete])+ Sep+str(CSVSortie['Couv'][Entete]/float(CSVSortie['NbSite'][Entete]))+ SepFin
     Chaine+="\n\n"
     DescriptifBase=CSVSortie['Stat']
     ListeBase=DescriptifBase[ListeEntete[0]].keys()
     EnteteBase=['Ref']
     EnteteBase.extend(ListeEntete)
     Chaine+=PrintListe(EnteteBase,SepDebEntete,SepEntete,SepFin)
     for Base in ListeBase :
         for Base2 in ListeBase :
             Chaine+=SepDebEntete+Base+SepEntete+Base2
             for Entete in ListeEntete:
                  Chaine+=Sep+str(DescriptifBase[Entete][Base][Base2])
             Chaine+=SepFin
     if DicComp!=None:
         Chaine+="\n\n"+SepTitreDebut1+"Comparaison des fichiers"+SepTitreFin1+"\n"
         ListeKeyDicComp=DicComp.keys()
         Chaine+=PrintListe(['Entete 1', 'Entete 2', 'Nombre de position Comm', 'Nombre de position Comm identique'],SepDebEntete,SepEntete,SepFin)
         for KeyDicComp in ListeKeyDicComp:
             for KeyDicComp2 in DicComp[KeyDicComp].keys():
                Chaine+=SepDebEntete+KeyDicComp +" (N="+str(CSVSortie['NbSite'][KeyDicComp]) +")" + SepDebEntete+KeyDicComp2+" (N="+str(CSVSortie['NbSite'][KeyDicComp2]) +")" + Sep + str(DicComp[KeyDicComp][KeyDicComp2]['NbBaseComm']) + Sep + str(DicComp[KeyDicComp][KeyDicComp2]['NbBaseCommAvecMemePoly'])+SepFin

     Chaine+="\n\n"
     if DicComp50K ==None :
         return Chaine
     Chaine+="\n\n"+SepTitreDebut1+"Comparaison avec la puce 50K"+SepTitreFin1+"\n"
     Chaine+=SepDebEntete+"Nom Fichier de la Puce "+Sep+ str(DicComp50K['NomFichier'])+SepFin
     Chaine+=SepDebEntete+"Nombre de SNP sur la puce "+Sep+ str(DicComp50K['NbSnpPuce'])+SepFin
     Chaine+=SepDebEntete+"Nombre de SNP polymorphique non hétérozygote sur la puce "+Sep+ str(DicComp50K['NbSnpPucePoly'])+SepFin
     Chaine+=SepDebEntete+"Nombre de SNP Non polymorphique sur la puce "+Sep+ str(DicComp50K['NbSnpPuceNonPoly'])+SepFin+"\n"
     Chaine+=PrintListe(['', 'Nombre de position Commune polymorphe', 'Nombre de pos com avec même poly', 'Nombre de pos non polymorphe CSV poly'],SepDebEntete,SepEntete,SepFin)
     for Entete in ListeEntete :
         Chaine+=SepDebEntete+Entete+Sep+str(DicComp50K[Entete]['NbPosCommunePoly'])+Sep+str(DicComp50K[Entete]['NbPosCommunePolyMemeMut'])+Sep+str(DicComp50K[Entete]['NbPosCommuneNonPoly']) +SepFin
     Chaine+="\n\n"
     return Chaine



def Test () :
	ListeCSV=['FileVCF/Snp_sup-0_StampyBowtie_20_0.9_0.05.vcf','FileVCF/S_sup-0_Bowtie_20_0.9_0.05.vcf','FileVCF/S_sup-0_Stampy_20_0.9_0.05.vcf']
	ListeEntete=['StampyBowtie', 'Bowtie', 'Stampy']
	CSVSortie=ReadMultiCSV(ListeCSV, ListeEntete)
	DicoPuce=LireFormatPuceV2("Data50k_Ak3-V2.out")
	DicComp=ComparerMultiCSV(CSVSortie)
#DicFusionner=FusionnerDic(CSVSortie)
	DicComp50K=ComparerCSVEt50k(CSVSortie,DicoPuce)
	Chaine=PrintStatCSVMerge(CSVSortie, 'Wiki', DicComp, DicComp50K)
	print Chaine




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
     Commande="exe -s [PuceADN, opt] -c [CSV1 CSV2...] -o FileOuptpout -e [Entete opt]"
     Fichier50K=GetUnArgument(sys.argv, '-s')
     if len(Fichier50K)>1:
        sys.exit(Commande)
     elif len(Fichier50K)==0 :
        Fichier50K=None
     else :
        Fichier50K=Fichier50K[0]
     ListeFichierCSV=GetUnArgument(sys.argv, '-c')
     if len(ListeFichierCSV) == 0:
        sys.exit(Commande)
     ListeOutput=GetUnArgument(sys.argv, '-o')
     if len(ListeOutput)>1:
        sys.exit(Commande)
     elif len(ListeOutput):
         Output=ListeOutput[0]
     else :
         Output=None
     ListeEntete=GetUnArgument(sys.argv, '-e')
     if len(ListeEntete)==0 :
        ListeEntete=ListeFichierCSV
     elif len(ListeEntete)!=len(ListeFichierCSV):
         sys.exit(Commande)
     #ListeCSV=['FileVCF/Snp_sup-0_StampyBowtie_20_0.9_0.05.vcf','FileVCF/S_sup-0_Bowtie_20_0.9_0.05.vcf','FileVCF/S_sup-0_Stampy_20_0.9_0.05.vcf']
     #ListeEntete=['StampyBowtie', 'Bowtie', 'Stampy']
     CSVSortie=ReadMultiCSV(ListeFichierCSV , ListeEntete)
     DicoPuce=LireFormatPuceV2(Fichier50K)
     DicComp=ComparerMultiCSV(CSVSortie)
     #DicFusionner=FusionnerDic(CSVSortie)
     DicComp50K=ComparerCSVEt50k(CSVSortie,DicoPuce)
     Chaine=PrintStatCSVMerge(CSVSortie, 'Wiki', DicComp, DicComp50K)
     if Output==None:
        print Chaine
     else :
        Ecrire=open(Output,'w')
        Ecrire.write(Chaine)
        Ecrire.close()
      

LaunchCommande()


