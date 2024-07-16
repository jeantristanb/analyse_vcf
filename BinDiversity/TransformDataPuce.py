
def Rev(Lettre):
    if Lettre=='A':
       return 'T'
    if Lettre=='C':
       return 'G'
    if Lettre=='T':
       return 'A'
    if Lettre=='G':
       return 'C'




def GetNewGeno(Base ,InfoSnp) :
    if Base[0]=='N' :
       InfoSnp['NC']+=1
       return "./."
    Base1=Base[0]
    Base2=Base[1]
    if InfoSnp['Rev']=='r':
       Base1=Rev(Base1)
       Base2=Rev(Base2)
    if Base1 not in InfoSnp['Base'] :
       InfoSnp['Base'].append(Base1)
       #InfoSnp['Freq'].append(0)
    Pos1=InfoSnp['Base'].index(Base1) 
    if Base2 not in InfoSnp['Base'] :
       InfoSnp['Base'].append(Base2)
       #InfoSnp['Freq'].append(0)
    Pos2=InfoSnp['Base'].index(Base2) 
    if Pos1!=Pos2:
       InfoSnp['Het']+=1
    elif Pos1==0 :
       InfoSnp['WT']+=1
    else : 
       InfoSnp['Hom']+=1
    return str(Pos1)+"/"+str(Pos2)

def GetSnp(FichierSnp, DicoSnp, ListeInd):
    Lire=open(FichierSnp)
    Dico={}
    ListeSnp=set(DicoSnp.keys())
    for Ind in ListeInd :
       Dico[Ind]={}
    for Ligne in Lire:
      TmpSplit=Ligne.split(';')
      NomInd=TmpSplit[3]
      NomPos=TmpSplit[2]
      if NomInd in ListeInd and NomPos in ListeSnp :
          Dico[NomInd][NomPos]=GetNewGeno(TmpSplit[4], DicoSnp[NomPos])
    return Dico
      
def GetSnpV2(FichierSnp, DicoSnp, ListeInd=None):
    Lire=open(FichierSnp)
    ListeNomPos=Lire.readline().replace("\n","").replace("\"","").split()
    Dico={}
    ListeNomSnp=set(DicoSnp)
    for Ligne in Lire:
      LigneSplit=Ligne.replace("\"","").replace("\n","").split()
      NomInd=LigneSplit[0].lower()
      if ListeInd!=None and NomInd in ListeInd :
         Dico[NomInd]={}
         CmtBon=0
         for Cmt in range(1,len(LigneSplit)) :  
             NomPos="\""+ListeNomPos[Cmt-1] +"\""
             if NomPos in ListeNomSnp :
                Dico[NomInd][NomPos]=GetNewGeno(LigneSplit[Cmt], DicoSnp[NomPos])
                CmtBon+=1
          #else :
          #   print NomPos,
    return Dico
      



def GetListeIndividu(FichierInd):
   LireInd=open(FichierInd)
   ListeInd=[]
   for Ligne in LireInd:
      ListeInd.append(Ligne.replace('\n','').lower()) 
   LireInd.close()
   return set(ListeInd)

def GetListeIndividuV2(FichierInd):
   LireInd=open(FichierInd)
   ListeInd={}
   for Ligne in LireInd:
      SplitLigne=Ligne.replace("\"","").replace('\n','').lower().split(';')
      ListeInd[SplitLigne[0]]=SplitLigne[1]
   LireInd.close()
   return ListeInd



def GetReadInfoSnp(FichierInfoSnp):
   Lire=open(FichierInfoSnp)
   Entete=Lire.readline().split()
   PosPass=Entete.index("Filtre")
   PosChro=Entete.index("ChroNew")
   PosPos=Entete.index("PosNew")
   DicoChro={}
   DicoChro2={}
   for Ligne in Lire :
       InfoSnp=Ligne.split()
       if InfoSnp[PosPass]=='TRUE':
         NomChro=InfoSnp[PosChro]
         NomPos=InfoSnp[PosPos]
         DicoChro['"'+InfoSnp[0]+'"']={'Chro':NomChro, 'Pos':NomPos, 'Base':[InfoSnp[5], ], 'Het':0, 'WT':0, 'Hom':0, 'NC':0, 'Rev':InfoSnp[4]}
   Lire.close()
   return DicoChro

def GetInfoInd(FichierInfoInf, ListeIndividu=None):
   Lire=open(FichierInfoInf)
   Dico={}
   DicoReverse={}
   ListeIndividu2=[]
   for Ligne in Lire :
       TmpSplit=Ligne.split(';')
       NomVariete=TmpSplit[2].split('_')[0]
       NomVariete=NomVariete.lower()
       NomInd='"'+TmpSplit[0]+'"'
       if ListeIndividu==None or NomVariete in ListeIndividu:
          Dico[NomInd]=NomVariete
          DicoReverse[NomVariete]=NomInd
          ListeIndividu2.append(NomVariete)
   Lire.close()
   return {'DicoInd':Dico, 'DicoReverse':DicoReverse, 'Ind2':ListeIndividu2}



FichierInfoSnp="/home/dygap/jeantristan/Data/PuceStephane/Info-Ref-V2-V3"
FichierSnp="/home/dygap/jeantristan/Data/PuceStephane/Matrix_Genotyping_Letter_DivZea_50K_HD-SSD-S1P9_Illumina_FW_2013-07-15_.txt"
FichierInfoInd="/home/dygap/jeantristan/AnalyseAll/AnalysePuce/Data/InfoPuceAutreSMH/InfoPuceNew"
FichierListeInd="/home/dygap/jeantristan/AnalyseAll/AnalysePuce/Data/InfoPuceAutreSMH/NomADNEchantillon"

DicoInd=GetListeIndividuV2(FichierInfoInd)
print DicoInd
DicoSnp=GetReadInfoSnp(FichierInfoSnp)
ListeInd=GetListeIndividu(FichierListeInd)

Ecrire=open("AllEchantillon.vcf","w")
Ecrire.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT")


Dico=GetSnpV2(FichierSnp, DicoSnp, ListeInd)

NewNomAccession=[]
NewNomInd=[]

for Accession in Dico.keys() :
   if len(Dico[Accession])!=0:
      NewNomAccession.append(Accession)
      NewNomInd.append(DicoInd[Accession])
   else :
      print DicoInd[Accession]

print len(NewNomAccession), len(NewNomInd), NewNomInd

for Ind in NewNomInd:
    Ecrire.write("\t"+Ind)

Ecrire.write("\n")

for Snp in DicoSnp :
    Ecrire.write(DicoSnp[Snp]['Chro']+"\t"+DicoSnp[Snp]['Pos']+"\t"+Snp + "\t"+DicoSnp[Snp]['Base'][0]+"\t")
    if len(DicoSnp[Snp]['Base'])>1 :
       Ecrire.write(DicoSnp[Snp]['Base'][1])
    if len(DicoSnp[Snp]['Base'])>2 :
       for Base in DicoSnp[Snp]['Base'][2::] :
           print Base
           Ecrire.write(","+Base)
    elif len(DicoSnp[Snp]['Base'])==1 :
       Ecrire.write(".")
    Ecrire.write("\t.\tPASS\tADP=0;WT="+str(DicoSnp[Snp]['WT'])+";HET="+ str(DicoSnp[Snp]['Het']) +";HOM="+str(DicoSnp[Snp]['Hom']) +";NC="+str(DicoSnp[Snp]['NC'])+"\tGT") 
    for Accession in NewNomAccession:
       Ecrire.write("\t"+Dico[Accession][Snp])
    Ecrire.write("\n")

Ecrire.close()
        
    
    
    





